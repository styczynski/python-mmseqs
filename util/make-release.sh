#!/bin/bash -ex
COMMIT="$1"
RELEASE_ID="$2"
RELEASE_MSG="$3"

if [ -z "${GITHUB_TOKEN}" ]; then
	echo "Please set GitHub Token"
  exit 1
fi

function hasCommand() {
	command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

hasCommand github-release
hasCommand echo
hasCommand date
hasCommand wget

# download CI builds
wget -O "biosnake-linux-sse41.tar.gz" "https://biosnake.com/archive/${COMMIT}/biosnake-linux-sse41.tar.gz"
wget -O "biosnake-linux-avx2.tar.gz" "https://biosnake.com/archive/${COMMIT}/biosnake-linux-avx2.tar.gz"
wget -O "biosnake-win64.zip" "https://biosnake.com/archive/${COMMIT}/biosnake-win64.zip"
wget -O "biosnake-osx-sse41.tar.gz" "https://biosnake.com/archive/${COMMIT}/biosnake-osx-sse41.tar.gz"
wget -O "biosnake-osx-avx2.tar.gz" "https://biosnake.com/archive/${COMMIT}/biosnake-osx-avx2.tar.gz"
wget -O "userguide.pdf" "https://biosnake.com/archive/${COMMIT}/userguide.pdf"

# create release tag 
git tag  "${RELEASE_ID}" && git push --tags

# create a formal release
github-release release \
    --user soedinglab \
    --repo biosnake2 \
    --tag "${RELEASE_ID}" \
    --name "Biosnake2 Release $RELEASE_ID" \
    --description "$RELEASE_MSG" \
    --pre-release

# upload AVX2 static binary
github-release upload \
    --user soedinglab \
    --repo biosnake2 \
    --tag "${RELEASE_ID}" \
    --name "Biosnake2-Linux-SSE4_1.tar.gz" \
    --file biosnake-linux-sse41.tar.gz

# upload SSE4.1 static binary
github-release upload \
    --user soedinglab \
    --repo biosnake2 \
    --tag "${RELEASE_ID}" \
    --name "Biosnake2-Linux-AVX2.tar.gz" \
    --file biosnake-linux-avx2.tar.gz

# upload Windows build
github-release upload \
    --user soedinglab \
    --repo biosnake2 \
    --tag "${RELEASE_ID}" \
    --name "Biosnake2-Windows-Unified.zip" \
    --file biosnake-win64.zip

# upload SSE4.1 static binary
github-release upload \
    --user soedinglab \
    --repo biosnake2 \
    --tag "${RELEASE_ID}" \
    --name "Biosnake2-MacOS-SSE4_1.tar.gz" \
    --file biosnake-osx-sse41.tar.gz

# upload AVX2 static binary
github-release upload \
    --user soedinglab \
    --repo biosnake2 \
    --tag "${RELEASE_ID}" \
    --name "Biosnake2-MacOS-AVX2.tar.gz" \
    --file biosnake-osx-avx2.tar.gz

# upload Windows build
github-release upload \
    --user soedinglab \
    --repo biosnake2 \
    --tag "${RELEASE_ID}" \
    --name "Biosnake2-Userguide.pdf" \
    --file userguide.pdf


