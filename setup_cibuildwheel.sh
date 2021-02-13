#!/bin/bash

echo "Activate EPEL repository"
yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm

echo "Install Python devel"
yum install -y python3-devel.x86_64 || exit 1

echo "Install Poetry..."
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
export PATH="$PATH:$HOME/.poetry/bin"
source $HOME/.poetry/env
echo "Poetry installed."

echo "Export Poetry requirements..."
poetry export -f requirements.txt --output requirements.txt
echo "Install everything..."
pip install -r requirements.txt
