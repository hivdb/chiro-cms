#! /bin/bash

set -e

cd $(dirname $0)/..

mkdir -p local/cache-nstrain-aapcnt/
pushd local/cache-nstrain-aapcnt/

rm global.json.gz global-frequencies.json global.json 2> /dev/null || true
curl -sSLR -o global.json.gz "https://nextstrain.org/charon/getDataset?prefix=/ncov/global"
curl -sSLR -o global-frequencies.json "https://nextstrain.org/charon/getDataset?prefix=/ncov/global&type=tip-frequencies"
gzip -d global.json.gz

popd
pipenv run python ./scripts/import_nstrainpcnt.py RdRP S
