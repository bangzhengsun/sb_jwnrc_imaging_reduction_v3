#!/usr/bin/env bash
set -eo pipefail
NEXUS_BASE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$NEXUS_BASE/env.sh"
mkdir -p "$NEXUS_BASE/automation/resample_30mas/scratch"
export TMPDIR="$NEXUS_BASE/automation/resample_30mas/scratch"
export PYTHONUNBUFFERED=1
python -u "$NEXUS_BASE/nexus_resample_30mas.py" 2>&1 | tee -a "$NEXUS_BASE/automation/resample_30mas/run.log"
