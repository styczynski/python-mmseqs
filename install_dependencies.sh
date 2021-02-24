#!/bin/bash
#
# This is a small script to install all project dependencies on
# Ubuntu-like system.
# @Covid Genomics 2021
#

# Recommended Python version installed by pyenv
PYTHON_VERSION="3.8.3"

# Install libssl and configure repositories
sudo apt update && apt-cache policy libssl1.0-dev
sudo apt-get install -y libssl1.0-dev

# Install standart Python/c++ tools (install standard toolchain)
sudo apt install -y gcc make python3 g++ cmake python3-dev openssl bzip2 \
build-essential zlib1g-dev libbz2-dev \
libreadline-dev libsqlite3-dev wget curl libncurses5-dev libncursesw5-dev \
xz-utils libffi-dev liblzma-dev zlib1g-dev libssl1.1 libssl-dev

# Install and configure pyenv
if pyenv versions ; then
    echo "Pyenv exists. Skipping installation"
else
    echo "Pyenv will be installed"
    curl https://pyenv.run | bash
    export PATH="$HOME/.pyenv/bin:$PATH"
    eval "$(pyenv init -)"
    eval "$(pyenv virtualenv-init -)"
    if pyenv versions ; then
      echo "Pyenv was installed"
    else
      echo "Pyenv installation failed for unknown reasons."
      exit 1
    fi
fi

# Install recommended Python version
if [ "$(pyenv global)" = "$PYTHON_VERSION" ]; then
    echo "Pyenv Python $PYTHON_VERSION is globally configured. Skipping this step"
else
    echo "Installing Python $PYTHON_VERSION via Pyenv"
    pyenv install $PYTHON_VERSION
    pyenv global $PYTHON_VERSION
    if [ "$(pyenv global)" = "$PYTHON_VERSION" ]; then
        echo "Python $PYTHON_VERSION was installed by Pyenv"
    else
        echo "Python $PYTHON_VERSION cannot be installed by Pyenv. For unknown reason installation has failed"
        exit 1
    fi
fi

# Install mmseqs
if mmseqs version ; then
   echo "MMseqs is installed. Skipping installation."
else
    echo "Installing mmseqs"
    wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz
    sudo mv mmseqs-linux-sse41.tar.gz /usr/bin/mmseqs-linux-sse41.tar.gz
    sudo tar xvfz /usr/bin/mmseqs-linux-sse41.tar.gz -C /usr/bin/
    sudo rm /usr/bin/mmseqs-linux-sse41.tar.gz
    export PATH=/usr/bin/mmseqs/bin/:$PATH
    if mmseqs version ; then
        echo "Mmseqs was installed"
    else
        echo "MMseqs installation failed for unknown reasons."
        exit 1
    fi
fi

# Install Poetry
if poetry version ; then
   echo "Poetry is installed. Skipping installation."
else
    echo "Installing Poetry"
    curl -sSL http://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python3 -
    export PATH="$PATH:$HOME/.poetry/bin"
    if mmseqs version ; then
      echo "Poetry was installed"
    else
      echo "Poetry installation failed for unknown reasons."
      exit 1
    fi
fi

# Instal wkhtmltopdf
if wkhtmltopdf --version ; then
   echo "wkhtmltopdf is installed. Skipping installation."
else
    echo "Installing wkhtmltopdf"
    sudo apt-get install -y wkhtmltopdf
    if wkhtmltopdf --version ; then
      echo "wkhtmltopdf was installed"
    else
      echo "wkhtmltopdf installation failed for unknown reasons."
      exit 1
    fi
fi

# Install project
poetry install

# Setup .bashrc
tee -a $HOME/.bashrc <<'EOF'
  # Setup mmseqs
  export PATH=/usr/bin/mmseqs/bin/:$PATH
  # Setup poetry
  export PATH="$PATH:$HOME/.poetry/bin"
  # Setup pyenv
  export PATH="$PATH:$HOME/.pyenv/bin"
  eval "$(pyenv init -)"
  eval "$(pyenv virtualenv-init -)"
  pyenv global 3.8.3
EOF

# Configure current session
source $HOME/.bashrc
poetry shell
