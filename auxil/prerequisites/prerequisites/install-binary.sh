#!/usr/bin/env bash
# BASH3 Boilerplate
#
# install-binary.sh
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


### Start of boilerplate -- do not edit this block #######################
export OPENCOARRAYS_SRC_DIR="${OPENCOARRAYS_SRC_DIR:-${PWD}/../}"
if [[ ! -f "${OPENCOARRAYS_SRC_DIR}/src/libcaf.h" ]]; then
  echo "Please run this script inside the top-level OpenCoarrays source \"prerequisites\" "
  echo "subdirectory or set OPENCOARRAYS_SRC_DIR to the OpenCoarrays source directory path."
  exit 1
fi
export __usage="${OPENCOARRAYS_SRC_DIR}/prerequisites/install-binary.sh-usage"
export B3B_USE_CASE="${B3B_USE_CASE:-${OPENCOARRAYS_SRC_DIR}/prerequisites/use-case}"
if [[ ! -f "${B3B_USE_CASE:-}/bootstrap.sh" ]]; then
  echo "Please set B3B_USE_CASE to the bash3boilerplate use-case directory path."
  exit 2
else
    # shellcheck source=./use-case/bootstrap.sh
    source "${B3B_USE_CASE}/bootstrap.sh" "$@"
fi
### End of boilerplate -- start user edits below #########################


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
if [[ "${arg_l}" != "${__flag_present}" && "${arg_N}" != "${__flag_present}" &&
      "${arg_v}" != "${__flag_present}" && "${arg_h}" != "${__flag_present}" &&
      -z "${arg_D:-${arg_p:-${arg_P:-${arg_U:-${arg_V}}}}}" ]]; then
  help "${__base}: Insufficient arguments. Please pass either -D, -h, -N, -l, -p, -P, -U, -v, -V, or a longer equivalent."
fi

# Suppress info and debug messages if -l, -P, -U, -V, -D, or their longer equivalent is present:
[[ "${arg_l}" == "${__flag_present}" || ! -z "${arg_P:-${arg_U:-${arg_V:-${arg_D}}}}" ]] && suppress_info_debug_messages

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

info "-d (--debug):            ${arg_d} "
info "-D (--print-downloader): ${arg_D} "
info "-e (--verbose):          ${arg_e} "
info "-h (--help):             ${arg_h} "
info "-i (--install-dir):      ${arg_i} "
info "-I (--install-version):  ${arg_I} "
info "-l (--list-packages):    ${arg_l} "
info "-n (--no-color):         ${arg_n} "
info "-p (--package):          ${arg_p}"
info "-P (--print-path):       ${arg_P} "
info "-U (--print-url):        ${arg_U} "
info "-v (--version):          ${arg_v} "
info "-V (--print-version):    ${arg_V} "
}
# shellcheck source=./install-binary-functions/set_or_print_default_version.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/install-binary-functions/set_or_print_default_version.sh"
set_or_print_default_version
export version_to_build="${arg_I:-${default_version}}"

# shellcheck source=./install-binary-functions/set_or_print_downloader.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/install-binary-functions/set_or_print_downloader.sh"
set_or_print_downloader

# shellcheck source=./install-binary-functions/set_or_print_url.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/install-binary-functions/set_or_print_url.sh"
set_or_print_url

# shellcheck source=./build-functions/set_or_print_installation_path.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/set_or_print_installation_path.sh"
set_or_print_installation_path

# shellcheck source=./build-functions/download_if_necessary.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/download_if_necessary.sh"
download_if_necessary

# shellcheck source=./build-functions/unpack_if_necessary.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/unpack_if_necessary.sh"
unpack_if_necessary

# shellcheck source=./install-binary-functions/set_or_print_csv_binary_names.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/install-binary-functions/set_or_print_csv_binary_names.sh"
set_or_print_csv_binary_names

# shellcheck source=./build-functions/set_SUDO_if_needed_to_write_to_directory.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/set_SUDO_if_needed_to_write_to_directory.sh"
set_SUDO_if_needed_to_write_to_directory "${arg_i}"

# shellcheck source=./install-binary-functions/move_binaries_to_install_path.sh
source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/install-binary-functions/move_binaries_to_install_path.sh"
move_binaries_to_install_path
