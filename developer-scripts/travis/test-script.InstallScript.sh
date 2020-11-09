#!/usr/bin/env bash

set -o errexit
set -o verbose
set -o pipefail
set -o nounset
set -o errtrace

__file=developer-scripts/travis/test-script.InstallScript.sh

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

echo "Performing Travis-CI script phase for the OpenCoarrays installation script..."
export FC=gfortran-${GCC}
export CC=gcc-${GCC}
export CXX=g++-${GCC}
${CC} --version
./install.sh --yes-to-all -i "${HOME}/opencoarrays" -j 4 -f "$(type -P "${FC}")" -c "$(type -P "${CC}")" -C "$(type -P "${CXX}")"
BUILD_LOC=(prerequisites/builds/opencoarrays/*/)
BUILD_LOC_DIR="${BUILD_LOC[${#BUILD_LOC[@]}-1]}"
if [[ -d "${BUILD_LOC_DIR}" ]]; then
    echo "Found opencoarrays build directory created by the install script:"
    echo "   ${BUILD_LOC_DIR}"
    (
	cd "${BUILD_LOC[${#BUILD_LOC[@]}-1]}"
	CTEST_LOC=(../../../installations/cmake/*/bin/ctest)
	INSTALLER_CTEST="${CTEST_LOC[${#CTEST_LOC[@]}-1]}"
	if [[ -x "${INSTALLER_CTEST}" ]] ;then
	    "${INSTALLER_CTEST}" --output-on-failure --schedule-random --repeat-until-fail "${NREPEAT:-5}"
	else
	    ctest --output-on-failure --schedule-random --repeat-until-fail "${NREPEAT:-5}"
	fi
    )
else
    echo "Failed to find install.sh build directory. Contents of prerequisites/builds is:"
    ls prerequisites/builds/*/
    exit 5
fi

echo "Done."
