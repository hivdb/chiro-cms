#!/bin/sh

pipenv run python ./scripts/json2yaml.py ./resources/susc-table/tableS4.json resources/susc-table/tableS4.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/tableS5.json resources/susc-table/tableS5.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/tableS6.json resources/susc-table/tableS6.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/tableS7.json resources/susc-table/tableS7.yml
