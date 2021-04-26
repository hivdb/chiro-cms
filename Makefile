requirements.txt: Pipfile.lock
	@pipenv lock --requirements > requirements.txt

build-docker: requirements.txt
	@docker pull python:3
	@docker build . --no-cache -t hivdb/chiro-cms-builder:latest

push-docker: build-docker
	@docker push hivdb/chiro-cms-builder:latest

pull-docker:
	@docker pull hivdb/chiro-cms-builder:latest

_docker-build:
	@rm -rf build/
	@docker run \
		--mount type=bind,source=$(PWD),target=/app \
		--workdir /app --rm -it \
		hivdb/chiro-cms-builder:latest \
		python build.py

_fast-build:
	@rm -rf build/
	@pipenv run python build.py

local/timestamp/%:
	@mkdir -p $$(dirname $@)
	touch $@

resources/sierra-sars2/outbreak.info/lineages.json: % : local/timestamp/%.$(shell date +"%Y%m%d")
	curl -sSL -o "resources/sierra-sars2/outbreak.info/lineages.json" "https://api.outbreak.info/genomics/lineage?name=*"
	@test "$$(jq .success resources/sierra-sars2/outbreak.info/lineages.json)" = "true" || (echo "success != true" && false)

build: $(shell find . -type f -not -path "./.git*" -a -not -path "*.swp" -a -not -path "*.swo" -a -not -path "*/.DS_Store" -a -not -path "*/.gradle/*" -a -not -path "*/build/*" -a -not -path "*.log" -a -not -path "*/local/*" | sed 's#\([| ]\)#\\\1#g') build.py build_plugins/*.py resources/sierra-sars2/outbreak.info/lineages.json
	@test -e $(shell which pipenv) && make _fast-build || make _docker-build

deploy-dev: build
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/chiro-cms-builder:latest \
		aws s3 sync /app/build s3://cms.hivdb.org/chiro-dev --delete

deploy-dev2: build
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/chiro-cms-builder:latest \
		aws s3 sync /app/build s3://cms.hivdb.org/chiro-dev2 --delete

deploy-prod: build
	@scripts/only-stable.sh make deploy-prod
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/chiro-cms-builder:latest \
		aws s3 sync /app/build s3://cms.hivdb.org/chiro-prod \
		--delete --cache-control max-age=600

deploy-all: deploy-dev deploy-prod

.PHONY: deploy-* _*
