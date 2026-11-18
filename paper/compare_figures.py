"""Rasterise regenerated and reference figures and report pixel differences.

    python paper/compare_figures.py results/paper/figures paper/reference/paper_figures [dpi]

Prints, per reference file, the percentage of pixels that differ by more than 10 % of full scale and
the mean absolute difference. Needs pdftoppm (poppler).
"""
import os, subprocess, sys, tempfile, numpy as np
import matplotlib.image as mpimg
NEW, REF = sys.argv[1], sys.argv[2]; DPI = sys.argv[3] if len(sys.argv) > 3 else "60"
rows = []
tmp = tempfile.mkdtemp()
for dp, _, fn in os.walk(REF):
    for f in fn:
        rel = os.path.relpath(os.path.join(dp, f), REF); new = os.path.join(NEW, rel)
        if not os.path.exists(new): rows.append((rel, "MISSING", None, None)); continue
        a, b = os.path.join(tmp, "a"), os.path.join(tmp, "b")
        subprocess.run(["pdftoppm", "-png", "-r", DPI, "-singlefile", os.path.join(dp, f), a], check=True)
        subprocess.run(["pdftoppm", "-png", "-r", DPI, "-singlefile", new, b], check=True)
        A, B = mpimg.imread(a + ".png")[..., :3], mpimg.imread(b + ".png")[..., :3]
        if A.shape != B.shape: rows.append((rel, "SIZE", A.shape, B.shape)); continue
        d = np.abs(A - B).max(axis=2)
        rows.append((rel, "OK", float((d > 0.1).mean() * 100), float(d.mean())))
for r in sorted(rows, key=lambda r: (r[1] != "OK", -(r[2] or 0) if r[1] == "OK" else 0)):
    a, b = (round(r[2], 3), round(r[3], 5)) if r[1] == "OK" else (r[2], r[3])
    print(f"{r[1]:7s} {str(a):>10} {str(b):>10}  {r[0]}")
