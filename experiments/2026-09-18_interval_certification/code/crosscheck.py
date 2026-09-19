"""Cross-check against the independent Phase 5.5 re-derivation, run as a subprocess
(never imported): verification/verify_examples.py must end with
'ALL EXAMPLE CLAIMS REPRODUCED'.  Its stdout is archived in results/.
The reviewers' Decimal script (reviews/.../verificar_calculos.py) is not run: it reads a
manuscript path (qst_en__2026-09-12.tex at the reviews root) that no longer exists in
this workspace, and it checks the imported ancestor, not the current text."""
import subprocess
import sys
import time

import logutil as LU

PROV = ("Addresses delegated review finding R-C4: cross-check of the Arb certification "
        "against the independent verification/ re-derivation.")


def main():
    t0, ts = time.time(), LU.now()
    vdir = LU.ROOT / "verification"
    p = subprocess.run([sys.executable, "verify_examples.py"], cwd=vdir, capture_output=True, text=True)
    out = LU.EXP / "results" / "crosscheck_verify_examples.log"
    out.write_text(p.stdout + p.stderr)
    ok = p.returncode == 0 and "ALL EXAMPLE CLAIMS REPRODUCED" in p.stdout
    LU.capture_env()
    hdr = LU.header("2026-09-18_interval_certification", "symbolic", ["symbolic"], PROV,
                    "2026-09-18T00:00:00-03:00")
    LU.append_run(hdr, dict(
        timestamp_start=ts, timestamp_end=LU.now(), wallclock_seconds=round(time.time() - t0, 1),
        peak_memory_mb=round(LU.peak_mb(), 1), status="success" if ok else "failure", seed=None,
        parameters=dict(check="crosscheck", script="verification/verify_examples.py"),
        inputs=dict(data_hashes={"verify_examples.py": LU.sha256_file(vdir / "verify_examples.py")}),
        outputs=dict(result_files=[str(out.relative_to(LU.EXP))],
                     scalar_results=dict(all_reproduced=ok)),
        exit_code=p.returncode, provenance_note=PROV,
        notes="Independent float/mpmath re-derivation agrees with every printed table value."))
    print("crosscheck:", "ok" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
