#!/usr/bin/env python3
"""Submit a preparation job, a sample array, and dependent final summaries."""

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

from sample_manifest import check_inputs, discover_samples, stage_sample


REPOSITORY = Path(__file__).resolve().parent
CONFIG_KEYS = (
    "FASTQ_DIR", "FASTQ_GLOB", "OUTPUT_DIR", "SLURM_MAX_SAMPLES",
    "SLURM_MAX_CONCURRENT", "SLURM_ACCOUNT", "SLURM_MAIL_USER", "SLURM_MAIL_TYPE",
    "CONDA_SH", "CONDA_ENV", "PIPELINE_SKIP_CONDA", "SAMPLE_READS", "SAMPLE_SIZE", "SAMPLE_SEED",
    "DECOY_REFERENCE_DIR", "BACTERIA_REFERENCE_DIR", "HOST_REFERENCE_DIR",
    "DECOY_BBSPLIT_REFERENCE", "BACTERIA_BBSPLIT_REFERENCE", "HOST_BBSPLIT_REFERENCE",
    "BBSPLIT_K", "ADAPTER_SEQUENCES", "ADAPTER_REVERSE_COMPLEMENTS", "ADAPTER_MAX_TRIMMING_ROUNDS",
    "CUTADAPT_ERROR_RATE", "CUTADAPT_MIN_OVERLAP", "MIN_READ_LENGTH",
    "THREADS", "GENE_POSITION_BINS", "METAGENE_MIN_FEATURE_READS", "CSV_CONVERSION_SCRIPT",
)


def atomic_json(path, value):
    temporary = path.with_name(path.name + f".{uuid.uuid4().hex}.tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    os.replace(temporary, path)


def config_digest(config):
    return hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()


def positive(config, key, default):
    value = str(config.get(key, default))
    if not re.fullmatch(r"[1-9][0-9]*", value):
        raise ValueError(f"{key} must be a positive integer")
    return int(value)


def read_config(path):
    path = Path(path).expanduser().resolve(strict=True)
    # Source trusted pipeline config as Bash does; read only named settings,
    # without exporting or printing unrelated environment/credential values.
    program = "import json,os,sys; print(json.dumps({k:os.environ[k] for k in sys.argv[1:] if k in os.environ}))"
    result = subprocess.run(["bash", "-c", 'set -ae; source "$1"; shift; exec "$@"',
                             "pipeline-config", str(path), sys.executable, "-c", program, *CONFIG_KEYS],
                            check=True, text=True, capture_output=True)
    config = json.loads(result.stdout)
    for key, default in (("FASTQ_DIR", "."), ("OUTPUT_DIR", ".")):
        config[key] = str(Path(config.get(key, default)).expanduser().resolve())
    for key in ("DECOY_REFERENCE_DIR", "BACTERIA_REFERENCE_DIR", "HOST_REFERENCE_DIR",
                "DECOY_BBSPLIT_REFERENCE", "BACTERIA_BBSPLIT_REFERENCE", "HOST_BBSPLIT_REFERENCE"):
        if config.get(key):
            config[key] = str(Path(config[key]).expanduser().resolve())
    converter = Path(config.get("CSV_CONVERSION_SCRIPT", "bacteria_with_decoys_csvConversion.py"))
    config["CSV_CONVERSION_SCRIPT"] = str(converter if converter.is_absolute() else REPOSITORY / converter)
    config.setdefault("FASTQ_GLOB", "*.fastq*")
    config.setdefault("CONDA_SH", "")
    config.setdefault("CONDA_ENV", "bacteria-RNAseq-pipeline")
    config.setdefault("PIPELINE_SKIP_CONDA", "false")
    config.setdefault("SLURM_MAIL_USER", "")
    config.setdefault("SLURM_MAIL_TYPE", "END,FAIL")
    config["PIPELINE_SKIP_CONDA"] = config["PIPELINE_SKIP_CONDA"].lower()
    if config["PIPELINE_SKIP_CONDA"] not in {"true", "false"}:
        raise ValueError("PIPELINE_SKIP_CONDA must be true or false")
    if config["CONDA_SH"]:
        config["CONDA_SH"] = str(Path(config["CONDA_SH"]).expanduser().resolve())
    if config.get("SAMPLE_READS", "false").lower() not in {"true", "false"}:
        raise ValueError("SAMPLE_READS must be true or false")
    for label in ("DECOY", "BACTERIA"):
        if not config.get(f"{label}_REFERENCE_DIR", config.get(f"{label}_BBSPLIT_REFERENCE", "")):
            raise ValueError(f"Set {label}_REFERENCE_DIR in {path}")
    return config


def validate_job_settings(config):
    if not config.get("SLURM_ACCOUNT", "").strip():
        raise ValueError("Set SLURM_ACCOUNT to your PACE account in the config file")
    if config["PIPELINE_SKIP_CONDA"] != "true":
        if not config["CONDA_SH"] or not Path(config["CONDA_SH"]).is_file():
            raise ValueError("Set CONDA_SH to your conda.sh file in the config file")
        if not config["CONDA_ENV"].strip():
            raise ValueError("Set CONDA_ENV to your conda environment name or path in the config file")
    if config["SLURM_MAIL_USER"] and not config["SLURM_MAIL_TYPE"].strip():
        raise ValueError("Set SLURM_MAIL_TYPE in the config file, or leave SLURM_MAIL_USER empty to disable emails")


def write_config(path, config, **overrides):
    path.parent.mkdir(parents=True, exist_ok=True)
    values = {**config, **{key: str(value) for key, value in overrides.items()}}
    path.write_text("# Snapshot generated for this run; paths are absolute.\n" +
                    "".join(f"{key}={shlex.quote(value)}\n" for key, value in values.items()))


def load_manifest(path):
    manifest = json.loads(Path(path).read_text())
    if manifest.get("version") != 1 or config_digest(manifest["config"]) != manifest["config_digest"]:
        raise ValueError("Unsupported or modified run manifest")
    return manifest


def submit(config_path, sbatch_script):
    config = read_config(config_path)
    validate_job_settings(config)
    maximum = positive(config, "SLURM_MAX_SAMPLES", 40)
    concurrent = positive(config, "SLURM_MAX_CONCURRENT", 8)
    threads = positive(config, "THREADS", 16)
    positive(config, "SAMPLE_SIZE", 5000)
    samples = discover_samples(config["FASTQ_DIR"], config["FASTQ_GLOB"])
    if len(samples) > maximum:
        raise ValueError(f"Found {len(samples)} samples, exceeding SLURM_MAX_SAMPLES={maximum}; increase the setting to process all of them")
    script = Path(sbatch_script).expanduser().resolve(strict=True)
    if not shutil.which("sbatch"):
        raise ValueError("sbatch is not available; submit from your Slurm login node")
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
    root = Path(config["OUTPUT_DIR"]) / "slurm_runs" / run_id
    (root / "logs").mkdir(parents=True)
    manifest_path = root / "run.json"
    manifest = {"version": 1, "run_id": run_id, "run_dir": str(root), "repository": str(REPOSITORY),
                "config": config, "config_digest": config_digest(config), "samples": samples, "jobs": {}}
    atomic_json(manifest_path, manifest)
    write_config(root / "config.env", config)
    print(f"Discovered {len(samples)} samples; at most {concurrent} run concurrently.\nRun directory: {root}", flush=True)
    submitted = []

    def schedule(role, *options):
        command = ["sbatch", "--parsable", f"--chdir={REPOSITORY}", f"--job-name=bacteria_{role}",
                   f"--output={root / 'logs' / (role + '-%A_%a.out')}", f"--cpus-per-task={threads}",
                   f"--account={config['SLURM_ACCOUNT']}"]
        if config["SLURM_MAIL_USER"]:
            command.extend((f"--mail-user={config['SLURM_MAIL_USER']}", f"--mail-type={config['SLURM_MAIL_TYPE']}"))
        else:
            command.append("--mail-type=NONE")
        command.extend(options)
        command.extend((str(script), role, str(manifest_path), str(REPOSITORY)))
        result = subprocess.run(command, cwd=REPOSITORY, check=True, text=True, capture_output=True)
        job_id = result.stdout.strip().splitlines()[-1].split(";", 1)[0]
        if not job_id.isdigit():
            raise ValueError(f"Unexpected sbatch job ID: {result.stdout.strip()}")
        submitted.append(job_id)
        manifest["jobs"][role] = job_id
        atomic_json(manifest_path, manifest)
        return job_id

    try:
        preparation = schedule("prepare")
        array = schedule("sample", f"--array=0-{len(samples) - 1}%{concurrent}",
                         f"--dependency=afterok:{preparation}", "--kill-on-invalid-dep=yes")
        final = schedule("finalize", f"--dependency=afterok:{array}", "--kill-on-invalid-dep=yes",
                         "--cpus-per-task=1", "--mem=8G", "--time=01:00:00")
    except (subprocess.CalledProcessError, ValueError):
        if submitted and shutil.which("scancel"):
            subprocess.run(["scancel", *submitted], check=False)
        manifest["submission_failed"] = True
        atomic_json(manifest_path, manifest)
        raise
    print(f"Preparation job: {preparation}\nSample array: {array}\nFinal CSV job: {final}", flush=True)
    return manifest_path


def job_threads(config):
    value = os.environ.get("SLURM_CPUS_PER_TASK", str(positive(config, "THREADS", 16)))
    if not re.fullmatch(r"[1-9][0-9]*", value):
        raise ValueError("SLURM_CPUS_PER_TASK must be positive")
    return value


def prepare(manifest_path):
    manifest = load_manifest(manifest_path)
    root, config = Path(manifest["run_dir"]), manifest["config"]
    runtime = root / "shared" / "config.env"
    write_config(runtime, config, OUTPUT_DIR=root / "shared", THREADS=job_threads(config))
    subprocess.run(["bash", str(Path(manifest["repository"]) / "map_bacteria_with_decoys.sh"),
                    str(runtime), "--prepare-only"], check=True)
    atomic_json(root / "shared" / "completed.json", {"run_id": manifest["run_id"], "config_digest": manifest["config_digest"]})


def run_sample(manifest_path):
    manifest = load_manifest(manifest_path)
    task = os.environ.get("SLURM_ARRAY_TASK_ID", "")
    if not re.fullmatch(r"[0-9]+", task) or not 0 <= int(task) < len(manifest["samples"]):
        raise ValueError("SLURM_ARRAY_TASK_ID must identify a sample in this run's manifest")
    sample = manifest["samples"][int(task)]
    root, config = Path(manifest["run_dir"]), manifest["config"]
    prepared = json.loads((root / "shared" / "completed.json").read_text())
    if prepared != {"run_id": manifest["run_id"], "config_digest": manifest["config_digest"]}:
        raise ValueError("Shared references do not match this run")
    output = root / "samples" / sample["sample"]
    output.mkdir(parents=True, exist_ok=True)
    marker = output / "completed.json"
    marker.unlink(missing_ok=True)
    (root / "aggregation.completed.json").unlink(missing_ok=True)
    inputs = output / ".inputs"
    stage_sample(sample, inputs, sample_reads=config.get("SAMPLE_READS", "false").lower() == "true",
                 sample_size=positive(config, "SAMPLE_SIZE", 5000), seed=config.get("SAMPLE_SEED", "1"))
    runtime = output / "config.env"
    write_config(runtime, config, FASTQ_DIR=inputs, FASTQ_GLOB="*.fastq*", OUTPUT_DIR=output,
                 SAMPLE_READS="false", THREADS=job_threads(config))
    subprocess.run(["bash", str(Path(manifest["repository"]) / "map_bacteria_with_decoys.sh"), str(runtime),
                    "--prepared-dir", str(root / "shared" / "bbsplit_alignments"), "--skip-summary"], check=True)
    check_inputs(sample)
    from bacteria_with_decoys_csvConversion import build_rows
    from read_breakdown_csv import build_table
    from aggregate_samples import host_enabled
    rows = build_rows(output)
    breakdown = build_table(output, host_enabled=host_enabled(config))
    if len(rows) != 1 or len(breakdown) != 1 or rows[0][0] != sample["sample"] or breakdown[0][0] != sample["sample"]:
        raise ValueError(f"Unexpected sample outputs for {sample['sample']}")
    atomic_json(marker, {"run_id": manifest["run_id"], "sample": sample["sample"],
                         "config_digest": manifest["config_digest"], "reads": sample["reads"]})
    shutil.rmtree(inputs)
    print(f"Completed sample {sample['sample']} (array task {task})", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="role", required=True)
    submission = commands.add_parser("submit")
    submission.add_argument("--config", type=Path, default=Path("config.env"))
    submission.add_argument("--sbatch-script", type=Path, default=REPOSITORY / "run_bacteria_pipeline.sbatch")
    for role in ("prepare", "sample", "finalize"):
        commands.add_parser(role).add_argument("manifest", type=Path)
    args = parser.parse_args()
    try:
        if args.role == "submit":
            submit(args.config, args.sbatch_script)
        elif args.role == "prepare":
            prepare(args.manifest)
        elif args.role == "sample":
            run_sample(args.manifest)
        else:
            from aggregate_samples import aggregate
            load_manifest(args.manifest)
            aggregate(args.manifest)
    except subprocess.CalledProcessError as exc:
        parser.exit(1, f"ERROR: Command failed ({exc.returncode}): {shlex.join(map(str, exc.cmd))}\n{exc.stderr or ''}")
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
