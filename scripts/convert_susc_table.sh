#!/bin/sh

cp ../covid-drdb-reports/susceptibility-data_files/table/*.json ./resources/susc-table/

FOLDER=./resources/susc-table

for f in "$FOLDER"/*.json ; do
    pipenv run python ./scripts/convert_susc_table.py "$f" "$f"
    filename=$(basename $f)
    filename="${filename%.*}"
    echo "$FOLDER/$filename.json"
    pipenv run python ./scripts/json2yaml.py "$FOLDER/$filename.json" "$FOLDER/$filename.yml"
done
