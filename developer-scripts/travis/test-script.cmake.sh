#!/usr/bin/env bash

set -o errexit
set -o verbose
set -o pipefail
set -o nounset
set -o errtrace

__file=developer-scripts/travis/test-script.cmake.sh

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

echo "Performing Travis-CI script phase for the OpenCoarrays direct cmake build..."

for version in ${GCC}; do
    mkdir "cmake-build-gcc${GCC}"
    export BLD_DIR="cmake-build-gcc${GCC}"
    export FC=gfortran-${version}

    ${FC} --version
    ${CC:-cc} --version || true
    mpif90 --version && mpif90 -show
    mpicc --version && mpicc -show

    # shellcheck disable=SC2153
    for BUILD_TYPE in ${BUILD_TYPES}; do
	# shellcheck disable=SC2015
	[[ -d "${BLD_DIR}" ]] && rm -rf "${BLD_DIR:?}"/* || true
	(
	    cd "${BLD_DIR}"
	    cmake -Wdev \
		  -DCMAKE_INSTALL_PREFIX:PATH="${HOME}/OpenCoarrays" \
		  -DCMAKE_BUILD_TYPE:STRING="${BUILD_TYPE}" \
		  ..
	    make -j 4
	    printf '\nDone compiling OpenCoarrays and tests!\n'
	    CTEST_FLAGS=(--output-on-failure --schedule-random --repeat-until-fail "${NREPEAT:-5}" --timeout "${TEST_TIMEOUT:-200}")
	    printf "\nctest flags: %s\n" "${CTEST_FLAGS[*]}"
	    if [[ "${BUILD_TYPE}" =~ Deb ]]; then
		printf "\nRunning ctest for a debug build...\n\n"
		export NREPEAT=1
	    else
		printf "\nRunning ctest for a non-debug build...\n\n"
	    fi
	    ctest "${CTEST_FLAGS[@]}"
	    make install
	    make uninstall
	)
    done
done
echo "Done."
