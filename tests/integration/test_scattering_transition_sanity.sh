#!/usr/bin/env bash
set -euo pipefail

exe=${1:?missing executable}
cfg=${2:?missing config}
verify=${3:?missing verifier}

workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

cp "$cfg" "$workdir/input.cfg"

pushd "$workdir" >/dev/null
"$exe" < input.cfg >/dev/null
python3 "$verify" output/scattering_rates.dat
popd >/dev/null
