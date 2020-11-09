#!/usr/bin/env bash
if ! [[ "${TRAVIS_TAG}" ]] || ! ${TRAVIS_SECURE_ENV_VARS} ; then
    unset encrypted_ef4535c39461_key || true
    unset encrypted_ef4535c39461_iv || true
    rm subkey-328B3A0E-secret.asc{,.enc} || true
else
    echo "Encrypted vars not unset"
fi
if [[ "${OSTYPE}" == [Dd]arwin* ]]; then
    export PATH="${PATH}:${HOME}/Library/Python/2.7/bin"
else
    export PATH="${CACHE}/bin:${PATH}"
    export FC=gfortran-${GCC}
    export CC=gcc-${GCC}
    export CXX=g++-${GCC}
fi
