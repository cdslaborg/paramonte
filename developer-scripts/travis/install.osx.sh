#!/usr/bin/env bash

set -o errexit
set -o verbose
set -o pipefail
set -o nounset
set -o errtrace

__file=developer-scripts/travis/install.osx.sh

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

echo "Performing Travis-CI installation phase on macOS..."

# Update and install via Homebrew on macOS
brew update > /dev/null

brew outdated

{
    mpif90 --version && mpif90 -show
} || echo "No mpif90"
{
    mpicc --version && mpicc -show
} || echo "No mpicc"

type -a cmake || echo "CMake not installed"
cmake --version || true

echo "CC: ${CC:-none defined}"
echo "\$CC --version:"
${CC:-} --version || true

echo "FC: ${FC:-none-defined}"
echo "\$FC --version:"
${FC:-} --version || true

echo "Done."
