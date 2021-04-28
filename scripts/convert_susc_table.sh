#!/bin/sh

FOLDER=./resources/susc-table
rm -rf "$FOLDER"/*

cp ../covid-drdb-reports/report_tables/*.json "$FOLDER"


for f in "$FOLDER"/*.json ; do
    pipenv run python ./scripts/convert_susc_table.py "$f" "$f"
    filename=$(basename $f)
    filename="${filename%.*}"
    echo "$FOLDER/$filename.json"
    pipenv run python ./scripts/json2yaml.py "$FOLDER/$filename.json" "$FOLDER/$filename.yml"
done
