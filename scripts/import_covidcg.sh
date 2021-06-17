#! /bin/bash

set -e

cd $(dirname $0)/..

mkdir -p local/cache-covidcg-aapcnt/
pushd local/cache-covidcg-aapcnt/

rm data_package.json 2> /dev/null || true
curl -sSLR -o data_package.json.gz "https://storage.googleapis.com/ve-public/v2.0/data_package.json.gz?nocache=$(shuf -i 1-100000 -n 1)"
gzip -d data_package.json.gz

popd
pipenv run python ./scripts/import_covidcg.py S

