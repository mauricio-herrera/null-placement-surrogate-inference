from pathlib import Path
import csv
import hashlib

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "FILE_MANIFEST.csv"
SHA = ROOT / "SHA256SUMS.txt"

EXCLUDE_FROM_MANIFEST = {"FILE_MANIFEST.csv", "SHA256SUMS.txt"}
EXCLUDE_FROM_SHA = {"SHA256SUMS.txt"}
EXCLUDE_PARTS = {".git", ".venv", "__pycache__"}


def included(path: Path, exclude_names: set[str]) -> bool:
    rel = path.relative_to(ROOT)
    if path.name in exclude_names:
        return False
    if any(part in EXCLUDE_PARTS for part in rel.parts):
        return False
    if path.suffix == ".pyc":
        return False
    return path.is_file()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

# The manifest deliberately excludes both integrity files, avoiding self-reference.
files = [p for p in sorted(ROOT.rglob("*")) if included(p, EXCLUDE_FROM_MANIFEST)]
with MANIFEST.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["path", "bytes"])
    for p in files:
        w.writerow([p.relative_to(ROOT).as_posix(), p.stat().st_size])

# SHA256SUMS includes FILE_MANIFEST.csv but excludes only itself.
sha_files = [p for p in sorted(ROOT.rglob("*")) if included(p, EXCLUDE_FROM_SHA)]
with SHA.open("w", encoding="utf-8") as f:
    for p in sha_files:
        f.write(f"{sha256(p)}  {p.relative_to(ROOT).as_posix()}\n")

print(f"Wrote {MANIFEST.relative_to(ROOT)} with {len(files)} archived files.")
print(f"Wrote {SHA.relative_to(ROOT)} with {len(sha_files)} SHA256 entries.")
