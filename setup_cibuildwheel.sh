#!/bin/bash

echo "Check if all packages are present..."
sudo apt-get update
sudo apt-get install -y python3 python3-setuptools python3-pkg-resources python3-pip python3-dev libffi-dev build-essential git

echo "Install Poetry..."
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
export PATH="$PATH:$HOME/.poetry/bin"
source $HOME/.poetry/env
echo "Poetry installed."

echo "Export Poetry requirements..."
poetry export -f requirements.txt --output requirements.txt
echo "Install everything..."
pip install -r requirements.txt
