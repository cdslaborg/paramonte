#!/usr/bin/env bash

set -o errexit
set -o verbose
set -o pipefail
set -o nounset
set -o errtrace

__file=developer-scripts/travis/install.linux.sh

# Error tracing
# requires `set -o errtrace`
__caf_err_report() {
    error_code=${?}
    echo "Error (code=${error_code}) in ${__file} in function ${1} on line ${2}." >&2
    false
    return ${error_code}
}
# Always provide an error backtrace
trap '__caf_err_report "${FUNCNAME:-.}" ${LINENO}' ERR

if [[ -n "${TRAVIS_TAG}" ]] && ${TRAVIS_SECURE_ENV_VARS} ; then
    brew update > /dev/null
    brew ls --versions gpg2 >/dev/null || brew install gpg2
    brew outdated gpg2 || brew upgrade gpg2
    type -P openssl || brew install openssl
    curl https://izaakbeekman.com/izaak.pubkey.txt | gpg --import
    git tag -v "${TRAVIS_TAG}"
fi
gfortran --version || true
gcc --version || true
g++ --version || true
