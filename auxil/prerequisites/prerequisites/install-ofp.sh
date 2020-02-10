#!/usr/bin/env bash
# BASH3 Boilerplate
#
# install-ofp.sh
#
#  - Build the Open Fortran Parser
#
# Usage: ./install-ofp.sh -i /opt
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

export OPENCOARRAYS_SRC_DIR="${OPENCOARRAYS_SRC_DIR:-${PWD}/..}"
if [[ ! -f "${OPENCOARRAYS_SRC_DIR}/src/libcaf.h" ]]; then
  echo "Please run this script inside the OpenCoarrays source \"prerequisites\" subdirectory"
  echo "or set OPENCOARRAYS_SRC_DIR to the top-level OpenCoarrays source directory path."
  exit 1
fi
export __usage="${OPENCOARRAYS_SRC_DIR}/prerequisites/install-ofp.sh-usage"
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

export __flag_present=1

# Verify requirements

[ -z "${LOG_LEVEL:-}" ] && emergency "Cannot continue without LOG_LEVEL. "

# shellcheck disable=SC2154
if [[ "${__os}" != "OSX" ]]; then
   echo "Source tranlsation via OFP is currently supported only on OS X."
   echo "Please submit an issue at https://github.com/sourceryinstitute/opencoarrays/issues."
   emergency "${PWD}/install-ofp.sh: Aborting."
fi

if [[ $(uname) == "Darwin"  ]]; then
  default_ofp_downloader=curl
  args="-LO"
else
  default_ofp_downloader=wget
  args="--no-check-certificate"
fi

# If -D is passed, print the download programs used for OFP and its prerequisites.
# Then exit with normal status.
# shellcheck  disable=SC2154
if [[ "${arg_D}" == "${__flag_present}" ]]; then
  echo "strategoxt-superbundle downloader: $("${OPENCOARRAYS_SRC_DIR}/prerequisites/install-binary.sh" -D strategoxt-superbundle)"
  echo "ofp-sdf default downloader: ${default_ofp_downloader}"
  exit 0
fi

# If -P is passed, print the default installation paths for OFP and its prerequisites.
# Then exit with normal status.
# shellcheck disable=SC2154
install_path="${arg_i}"
strategoxt_superbundle_install_path=$("${OPENCOARRAYS_SRC_DIR}/prerequisites/install-binary.sh" -P strategoxt-superbundle)
# shellcheck disable=SC2154
if [[ "${arg_P}" == "${__flag_present}" ]]; then
  echo "strategoxt-superbundle default installation path: ${strategoxt_superbundle_install_path}"
  echo "ofp default installation path: ${install_path}"
  exit 0
fi

# If -V is passed, print the default versions of OFP and its prerequisites.
# Then exit with normal status.
default_ofp_version=sdf
# shellcheck disable=SC2154
if [[ "${arg_V}" == "${__flag_present}" ]]; then
  echo "strategoxt-superbundle default version: $("${OPENCOARRAYS_SRC_DIR}/prerequisites/install-binary.sh" -V strategoxt-superbundle)"
  echo "ofp default version: ${default_ofp_version}"
  exit 0
fi

# If -U is passed, print the URLs for OFP and its prerequisites.
# Then exit with normal status.
ofp_url_head="https://github.com/sourceryinstitute/opencoarrays/files/213108/"
ofp_url_tail="ofp-sdf.tar.gz"
# shellcheck disable=SC2154
if [[ "${arg_U}" == "${__flag_present}" ]]; then
  echo "strategoxt-superbundle URL: $("${OPENCOARRAYS_SRC_DIR}/prerequisites/install-binary.sh" -U strategoxt-superbundle)"
  echo "ofp URL: ${ofp_url_head}${ofp_url_tail}"
  exit 0
fi

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

info "-d (--debug):            ${arg_d}"
info "-D (--print-downloader): ${arg_D}"
info "-e (--verbose):          ${arg_e}"
info "-h (--help):             ${arg_h}"
info "-i (--install-dir):      ${arg_i}"
info "-I (--install-version):  ${arg_i}"
info "-j (--num-threads):      ${arg_j}"
info "-n (--no-color):         ${arg_n}"
info "-P (--print-path):       ${arg_P}"
info "-U (--print-url):        ${arg_U}"
info "-V (--print-version):    ${arg_V}"
}
# Set OFP installation path to the value of the -i argument if present.
# Otherwise, install OFP in the OpenCoarrays prerequisites/installations directory.
opencoarrays_prerequisites_dir="${OPENCOARRAYS_SRC_DIR}"/prerequisites/
if [[ "${arg_i}" == "${__flag_present}" ]]; then
  install_path="${arg_i}"
else
  install_path="${opencoarrays_prerequisites_dir}"/installations
fi

ofp_prereqs_install_dir="/opt"
# Change present working directory to installation directory
if [[ ! -d "${install_path}" ]]; then
  # shellcheck source=./build-functions/set_SUDO_if_needed_to_write_to_directory.sh
  source "${opencoarrays_prerequisites_dir}/build-functions/set_SUDO_if_needed_to_write_to_directory.sh"
  set_SUDO_if_needed_to_write_to_directory "${install_path}"
  ${SUDO:-} mkdir -p "${install_path}"
fi

# Install OFP prerequisites to /opt (currently the only option)
"${opencoarrays_prerequisites_dir}"/install-binary.sh -p strategoxt-superbundle -i "${strategoxt_superbundle_install_path}"

# Downlaod OFP
pushd "${install_path}"
info "OFP Download command: ${default_ofp_downloader} ${args:-} \"${ofp_url_head}${ofp_url_tail}\""
${default_ofp_downloader} ${args:-} "${ofp_url_head}${ofp_url_tail}" 


# Uncompress OFP
tar xf ofp-sdf.tar.gz
# Return to the original working directory
popd

export SDF2_PATH="${ofp_prereqs_install_dir}"/sdf2-bundle/v2.4/bin
export ST_PATH="${ofp_prereqs_install_dir}"/strategoxt/v0.17/bin
export DYLD_LIBRARY_PATH="${ofp_prereqs_install_dir}"/strategoxt/v0.17/lib:/opt/aterm/v2.5/lib

export OFP_HOME="${install_path}"/ofp-sdf
# shellcheck source=./install-binary-functions/build_parse_table.sh
source "${opencoarrays_prerequisites_dir}"/install-binary-functions/build_parse_table.sh
build_parse_table
