#!/usr/bin/env bash
set -euo pipefail

exe=${1:?missing executable}
cfg=${2:?missing config}
verify=${3:?missing verifier}

workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

cp "$cfg" "$workdir/input.toml"

pushd "$workdir" >/dev/null
"$exe" < input.toml >/dev/null
python3 "$verify" output/scattering_rates.dat
popd >/dev/null
