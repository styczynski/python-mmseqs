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
	poetry run pytest -n 4 --ignore lib

publish:
	poetry build && poetry run s3pypi --bucket pypi.covidgenomics.com --private --region eu-west-1 --dist-path dist

format: lint
	docker run -it -v "$(CURDIR)":/workdir -w /workdir unibeautify/clang-format -style=Google -i $(ALL_SOURCES)

documentation:
	rm -rfd documentation && cd .doc && rm -rfd _build && make html && mv _build/html ../documentation
