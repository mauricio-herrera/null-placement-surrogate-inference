#!/usr/bin/env bash
set -euo pipefail
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error manuscript_final.tex > "build_main_${i}.log"
done
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error Supplementary_Material_Final.tex > "build_supp_${i}.log"
done
if grep -Eq 'undefined|Overfull|Underfull|LaTeX Warning|Package .* Warning' manuscript_final.log Supplementary_Material_Final.log; then
  echo 'Final-pass LaTeX warnings detected:'
  grep -E 'undefined|Overfull|Underfull|LaTeX Warning|Package .* Warning' manuscript_final.log Supplementary_Material_Final.log || true
  exit 2
fi
echo 'Build complete: final-pass logs are clean.'
