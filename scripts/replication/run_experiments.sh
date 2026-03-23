#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Edit these paths for your local setup.
IDBC_BINARY="/Users/feb223/projects/coin/intersectionCuts/build-idBC-opt/bin/mibs"
OUTPUT_DIR="$SCRIPT_DIR/results"

# Test-set pairs: NAME PATH NAME PATH ...
TEST_SETS=(
  "iblpSmall" "/Users/feb223/projects/coin/intersectionCuts/Mibs/data/improvingDirectionDatasets/iblpSmall"
  "interSmall" "/Users/feb223/projects/coin/intersectionCuts/Mibs/data/improvingDirectionDatasets/interSmall"
)

DATA_SETS=(
  "iblpSmall" "interSmall"
)

python "$SCRIPT_DIR/make_plots.py" \
  --binariesPath "idBC" "$IDBC_BINARY" \
  --instanceDirs "${TEST_SETS[@]}" \
  --outputDir "$OUTPUT_DIR" \
  --writeParams \
  --testName "idBC" \
  --dataSets "${DATA_SETS[@]}" \

