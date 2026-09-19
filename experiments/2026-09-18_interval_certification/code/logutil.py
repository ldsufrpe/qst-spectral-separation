"""Append-only experiment_log.json + environment capture (shared helper, copied
verbatim into each experiment of the 2026-09-18 campaign; SHA-256 logged)."""
import datetime as _dt
import fcntl
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
from pathlib import Path

EXP = Path(__file__).resolve().parent.parent
ROOT = EXP.parent.parent


def now():
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


def sha256_file(p):
    return "sha256:" + hashlib.sha256(Path(p).read_bytes()).hexdigest()


def git_commit():
    try:
        return "git:" + subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "git:unknown"


def capture_env():
    env = EXP / "env"
    env.mkdir(exist_ok=True)
    freeze = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
    (env / "requirements.txt").write_text(freeze)
    cpu = ""
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    mem_kb = 0
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemTotal"):
                mem_kb = int(line.split()[1])
    except OSError:
        pass
    info = dict(os=platform.platform(), kernel=platform.release(),
                python=sys.version.split()[0], cpu=cpu, n_cores=os.cpu_count(),
                ram_gb=round(mem_kb / 2 ** 20, 1), captured_at=now())
    (env / "system_info.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


def peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def load_log(header):
    p = EXP / "experiment_log.json"
    if p.exists():
        return json.loads(p.read_text())
    log = dict(header)
    log.setdefault("runs", [])
    return log


def append_run(header, run):
    with open(EXP / ".experiment_log.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        return _append_run(header, run)


def _append_run(header, run):
    log = load_log(header)
    for key in ("plan_md_hash", "code_commit", "system_info_hash"):
        log[key] = header[key]
    run["run_id"] = f"r{len(log['runs']) + 1:04d}"
    log["runs"].append(run)
    (EXP / "experiment_log.json").write_text(json.dumps(log, indent=2, ensure_ascii=False) + "\n")
    return run["run_id"]


def header(experiment_id, mode, subareas, provenance_note, created_at):
    return dict(
        experiment_id=experiment_id,
        paper_project="2026-09_qst-spectral-separation",
        opportunity_reference="opportunities.md #1",
        mode=mode, subareas_active=subareas,
        target_venue="Linear Algebra and its Applications",
        created_at=created_at,
        provenance_note=provenance_note,
        plan_md_hash=sha256_file(EXP / "plan.md"),
        code_commit=git_commit(),
        system_info_hash=sha256_file(EXP / "env" / "system_info.json"),
        runs=[])
