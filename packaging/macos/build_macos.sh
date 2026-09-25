#!/usr/bin/env bash
set -euo pipefail

architecture="${1:-}"
if [[ "$architecture" != "arm64" && "$architecture" != "x86_64" ]]; then
  echo "Usage: build_macos.sh arm64|x86_64" >&2
  exit 2
fi
if [[ "$(uname -s)" != "Darwin" || "$(uname -m)" != "$architecture" ]]; then
  echo "Build on matching macOS hardware; cross-building is not supported." >&2
  exit 2
fi
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$script_dir/../.." && pwd)"
python="$root/.venv/bin/python"
"$python" -m pytest "$root/tests"
"$python" -m compileall -q "$root/src"
"$python" -m PyInstaller --noconfirm \
  --distpath "$root/dist/release-macos" \
  --workpath "$root/build/release-macos" \
  "$script_dir/PichAnalysis.spec"
echo "macOS bundle built for $architecture; runtime and signing validation remain required."
