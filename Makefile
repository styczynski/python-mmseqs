ALL_SOURCES := $(shell find ./src -iname *.h -o -iname *.cpp -o -iname *.hpp -o -iname *.c)

lint:
	poetry run black .
	poetry run isort .
	poetry run flakehell lint

install:
	poetry install

update:
	poetry update

test:
	mkdir -p test_reports && poetry run pytest -n 4 --ignore lib --junitxml=test_reports/latest.xml --cov=biosnake --cov-report html && mv htmlcov test_reports/coverage_latest_html

publish:
	poetry build && poetry run s3pypi --bucket pypi.covidgenomics.com --private --region eu-west-1 --dist-path dist

format: lint
	docker run -it -v "$(CURDIR)":/workdir -w /workdir unibeautify/clang-format -style=Google -i $(ALL_SOURCES)

documentation:
	rm -rfd documentation && cd .doc && rm -rfd _build && poetry run make html && mv _build/html ../documentation
