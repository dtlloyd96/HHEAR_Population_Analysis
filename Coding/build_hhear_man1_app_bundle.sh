#!/usr/bin/env zsh
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTPUT_DIR="${1:-$ROOT_DIR/HHEAR_Man1_App}"

python3 "$ROOT_DIR/Coding/build_targeted_analyte_mapping.py"
Rscript "$ROOT_DIR/Coding/build_hhear_man1_app_bundle.R" "$OUTPUT_DIR"
