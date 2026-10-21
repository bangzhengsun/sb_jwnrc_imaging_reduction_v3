#!/bin/bash
set -o pipefail

# Usage check
if [ $# -lt 2 ]; then
  echo "Usage: $0 <filter_name> <round_number>"
  exit 1
fi

filter="$1"
round="$2"
logfile="astrometry_${filter}-Round${round}.log"

# Run Python script and save full output
python SBNRC_05_astrometry.py "$filter" 2>&1 | tee "$logfile"
status=$?
if [ "$status" -ne 0 ]; then exit "$status"; fi

# Search for problematic matches
matches=$(grep -A 5 "Not enough matches" "$logfile")

if [ -n "$matches" ]; then
  echo "⚠️  Some alignment issues found:"
  echo "$matches"
else
  echo "No insufficient-match marker found; numerical astrometry QC remains required."
fi
