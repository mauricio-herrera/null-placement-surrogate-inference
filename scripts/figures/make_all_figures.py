from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]

scripts = [
    ROOT / "scripts" / "figures" / "make_figure1.py",
    ROOT / "scripts" / "figures" / "make_figures234.py",
]

for script in scripts:
    print(f"Running {script.name}...")
    subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)

print("Done.")
