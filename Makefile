build: build.py pages images resources downloads
	@rm -rf build/
	@docker run \
		--mount type=bind,source=$(PWD),target=/app \
		--workdir /app --rm -it \
		hivdb/hivdb-cms-builder:latest \
		python build.py

deploy-dev: build
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/hivdb-cms-builder:latest \
		aws s3 sync /app/build s3://cms.hivdb.org/chiro-dev --delete

deploy-prod: build
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/hivdb-cms-builder:latest \
		aws s3 sync /app/build s3://cms.hivdb.org/chiro-prod \
		--delete --cache-control max-age=600

deploy-all: deploy-dev deploy-prod

.PHONY: deploy-*
