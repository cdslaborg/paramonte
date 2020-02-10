#!/usr/bin/env bash
# BASH3 Boilerplate
#
# build.sh
#
#  - Build OpenCoarrays prerequisite packages and their prerequisites
#
# Usage: LOG_LEVEL=7 B3B_USE_CASE=/opt/bash3boilerplate/src/use-case ./my-script.sh -f script_input.txt
#
# More info:
#
#  - https://github.com/kvz/bash3boilerplate
#  - http://kvz.io/blog/2013/02/26/introducing-bash3boilerplate/
#
# Version: 2.0.0
#
# Authors:
#
#  - Kevin van Zonneveld (http://kvz.io)
#  - Izaak Beekman (https://izaakbeekman.com/)
#  - Alexander Rathai (Alexander.Rathai@gmail.com)
#  - Dr. Damian Rouson (http://www.sourceryinstitute.org/) (documentation)
#
# Licensed under MIT
# Copyright (c) 2013 Kevin van Zonneveld (http://kvz.io)

# The invocation of bootstrap.sh below performs the following tasks:
# (1) Import several bash3boilerplate helper functions & default settings.
# (2) Set several variables describing the current file and its usage page.
# (3) Parse the usage information (default usage file name: current file's name with -usage appended).
# (4) Parse the command line using the usage information.

export OPENCOARRAYS_SRC_DIR="${OPENCOARRAYS_SRC_DIR:-${PWD%prerequisites*}}"
export __usage=${OPENCOARRAYS_SRC_DIR}/prerequisites/build.sh-usage
if [[ ! -f "${OPENCOARRAYS_SRC_DIR}/src/libcaf.h" ]]; then
  echo "Please run this script inside the OpenCoarrays source \"prerequisites\" subdirectory"
  echo "or set OPENCOARRAYS_SRC_DIR to the top-level OpenCoarrays source directory path."
  exit 1
fi
export B3B_USE_CASE="${B3B_USE_CASE:-${OPENCOARRAYS_SRC_DIR}/prerequisites/use-case}"
if [[ ! -f "${B3B_USE_CASE:-}/bootstrap.sh" ]]; then
  echo "Please set B3B_USE_CASE to the bash3boilerplate use-case directory path."
  exit 2
fi
# shellcheck source=./use-case/bootstrap.sh
source "${B3B_USE_CASE}/bootstrap.sh" "$@"

# Set up a function to call when receiving an EXIT signal to do some cleanup. Remove if
# not needed. Other signals can be trapped too, like SIGINT and SIGTERM.
function cleanup_before_exit () {
  info "Cleaning up. Done"
}
trap cleanup_before_exit EXIT # The signal is specified here. Could be SIGINT, SIGTERM etc.

### Validation (decide what's required for running your script and error out)
#####################################################################

export __flag_present=1

# shellcheck disable=SC2154
if [[ "${arg_l}" != "${__flag_present}" &&
      "${arg_v}" != "${__flag_present}" &&
      "${arg_h}" != "${__flag_present}" &&
      -z "${arg_D:-${arg_p:-${arg_P:-${arg_U:-${arg_V:-${arg_B}}}}}}" ]]; then
  help "${__base}: Insufficient arguments. Please pass either -B, -D, -h, -l, -L, -p, -P, -U, -v, -V, or a longer equivalent."
fi

# Suppress info and debug messages if -B, -l, -P, -U, -V, -D, or their longer equivalent is present:
 [[ "${arg_l}" == "${__flag_present}" || ! -z "${arg_P:-${arg_U:-${arg_V:-${arg_D:-${arg_B}}}}}" ]] && suppress_info_debug_messages

[ -z "${LOG_LEVEL:-}" ] && emergency "Cannot continue without LOG_LEVEL. "

### Enforce mutual exclusivity of arguments that print single-line output
[ ! -z "${arg_P:-}" ] && [ ! -z "${arg_V:-}" ] && emergency "Only specify one of -P, -U, -V, or their long-form equivalents."
[ ! -z "${arg_P:-}" ] && [ ! -z "${arg_U:-}" ] && emergency "Only specify one of -P, -U, -V, or their long-form equivalents."
[ ! -z "${arg_U:-}" ] && [ ! -z "${arg_V:-}" ] && emergency "Only specify one of -P, -U, -V, or their long-form equivalents."

### Print bootstrapped magic variables to STDERR when LOG_LEVEL
### is at the default value (6) or above.
#####################################################################
# shellcheck disable=SC2154
{
info "__file: ${__file}"
info "__dir: ${__dir}"
info "__base: ${__base}"
info "__os: ${__os}"
info "__usage: ${__usage}"
info "LOG_LEVEL: ${LOG_LEVEL}"

info "-b (--install-branch):   ${arg_b} "
info "-B (--list-branches):    ${arg_B} "
info "-c (--with-c):           ${arg_c} "
info "-C (--with-cxx):         ${arg_C} "
info "-d (--debug):            ${arg_d} "
info "-D (--print-downloader): ${arg_D} "
info "-e (--verbose):          ${arg_e} "
info "-f (--with-fortran):     ${arg_f} "
info "-h (--help):             ${arg_h} "
info "-i (--install-prefix):   ${arg_i} "
info "-I (--install-version):  ${arg_I} "
info "-j (--num-threads):      ${arg_j} "
info "-l (--list-packages):    ${arg_l} "
info "-m (--with-cmake):       ${arg_m} "
info "-M (--with-mpi):         ${arg_M} "
info "-n (--no-color):         ${arg_n} "
info "-o (--only-download):    ${arg_o} "
info "-p (--package):          ${arg_p} "
info "-P (--print-path):       ${arg_P} "
info "-r (--prefix-root):      ${arg_r} "
info "-u (--from-url):         ${arg_u} "
info "-U (--print-url):        ${arg_U} "
info "-v (--version):          ${arg_v} "
info "-V (--print-version):    ${arg_V} "
info "-y (--yes-to-all):       ${arg_y} "
info "-Z (--bootstrap):        ${arg_Z} "
}

if [[ -z "${arg_B}" ]]; then
  # shellcheck source=./build-functions/set_or_print_default_version.sh
  source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/set_or_print_default_version.sh"
  # shellcheck disable=SC2119
  set_or_print_default_version
  export version_to_build="${arg_I:-${arg_b:-${default_version}}}"
fi

# shellcheck source=./build-functions/set_or_print_downloader.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/set_or_print_downloader.sh"
# shellcheck disable=SC2119
set_or_print_downloader $@

# shellcheck source=./build-functions/set_or_print_url.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/set_or_print_url.sh"
set_or_print_url

# shellcheck source=./build-functions/set_or_print_installation_path.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/set_or_print_installation_path.sh"

if [[ -z "${arg_B}" ]]; then
  set_or_print_installation_path
fi

# shellcheck source=./build-functions/download_if_necessary.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/download_if_necessary.sh"
download_if_necessary

# Exit if -o or --only-download was specified when invoking this script
if [[ ${arg_o:-} == "${__flag_present}" ]]; then
   info "No installation to perform: -o or --only-download specified.  Exiting."
   exit 0
fi

# If -Z or --bootstrap was specified, enable bootstrap configure & build
if [[ ${arg_Z:-} != "${__flag_present}" ]]; then
   info "Disabling bootstrap. If the gcc/g++/gfortran build fails, try './install.sh --bootstrap'."
   export bootstrap_configure="--disable-bootstrap"
   export bootstrap_build=""
else
   export bootstrap_configure=""
   export bootstrap_build="bootstrap"
fi

if [[ -z "${arg_B}" ]]; then
  # shellcheck source=./build-functions/unpack_if_necessary.sh
  source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/unpack_if_necessary.sh"
  unpack_if_necessary

  # shellcheck source=./build-functions/set_compilers.sh
  source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/set_compilers.sh"
  set_compilers

  # shellcheck source=./build-functions/build_and_install.sh
  source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/build_and_install.sh"
  build_and_install
fi
