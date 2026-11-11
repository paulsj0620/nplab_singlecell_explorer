#!/bin/bash
set -eu
cd -- "$(dirname -- "$0")"
export EXPLORER_SOURCE="${1:-${EXPLORER_SOURCE:-/Users/paul/WorkInProgress/dietary_intervention}}"
export EXPLORER_CACHE="${EXPLORER_CACHE:-$PWD/explorer-data.rds}"
if [[ ! -s "$EXPLORER_CACHE" ]]; then
  Rscript prepare.R "$EXPLORER_SOURCE"
fi
exec Rscript -e 'shiny::runApp(".", host = "127.0.0.1", launch.browser = TRUE)'
