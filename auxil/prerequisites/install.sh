#!/usr/bin/env bash
#
# install.sh
#
# -- This script installs OpenCoarrays and its prerequisites.
#
# OpenCoarrays is distributed under the OSI-approved BSD 3-clause License:
# Copyright (c) 2015-2016, Sourcery, Inc.
# Copyright (c) 2015-2016, Sourcery Institute
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this
#    list of conditions and the following disclaimer in the documentation and/or
#    other materials provided with the distribution.
# 3. Neither the names of the copyright holders nor the names of their contributors
#    may be used to endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Portions of this script derive from BASH3 Boilerplate and are distributed under
# the following license:
#
# The MIT License (MIT)
#
# Copyright (c) 2014 Kevin van Zonneveld
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
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
export OPENCOARRAYS_SRC_DIR="${OPENCOARRAYS_SRC_DIR:-${PWD%/}}"
if [[ ! -f "${OPENCOARRAYS_SRC_DIR}/src/libcaf.h" ]]; then
  echo "Please run this script inside the top-level OpenCoarrays source directory or "
  echo "set OPENCOARRAYS_SRC_DIR to the OpenCoarrays source directory path."
  exit 1
fi
export B3B_USE_CASE="${B3B_USE_CASE:-${OPENCOARRAYS_SRC_DIR}/prerequisites/use-case}"
if [[ ! -f "${B3B_USE_CASE:-}/bootstrap.sh" ]]; then
  echo "Please set B3B_USE_CASE to the bash3boilerplate use-case directory path."
  exit 2
else
    # shellcheck source=./prerequisites/use-case/bootstrap.sh
    source "${B3B_USE_CASE}/bootstrap.sh" "$@"
fi
### End of boilerplate -- start user edits below #########################

# Set expected value of present flags that take no arguments
export __flag_present=1

# Set up a function to call when receiving an EXIT signal to do some cleanup. Remove if
# not needed. Other signals can be trapped too, like SIGINT and SIGTERM.
function cleanup_before_exit () {
  info "Cleaning up. Done"
}
trap cleanup_before_exit EXIT # The signal is specified here. Could be SIGINT, SIGTERM etc.


### Validation (decide what's required for running your script and error out)
#####################################################################

[ -z "${LOG_LEVEL:-}" ] && emergency "Cannot continue without LOG_LEVEL. "

# shellcheck disable=SC2154
if [[ "${arg_v}" == "${__flag_present}" || "${arg_l}" == "${__flag_present}" || ! -z "${arg_P:-${arg_U:-${arg_V:-${arg_D:-${arg_B}}}}}" ]]; then
   print_debug_only=7
   if [ "$(( LOG_LEVEL < print_debug_only ))" -ne 0 ]; then
     debug "Supressing info and debug messages: one of {-l, -v, -P, -U, -V, -D} present."
     suppress_info_debug_messages
   fi
fi

[ ! -z "${arg_D}" ] && [ ! -z "${arg_P:-${arg_U:-${arg_V:-${arg_B}}}}" ] &&
  emergency "Please pass only one of {-B, -D, -p, -P, -U, -V} or a longer equivalent (multiple detected). [exit 101]"

[ ! -z "${arg_P}" ] && [ ! -z "${arg_U:-${arg_V:-${arg_B}}}" ] &&
  emergency "Please pass only one of {-B, -D, -p, -P, -U, -V} or a longer equivalent (multiple detected). [exit 102]"

[ ! -z "${arg_U}" ] && [ ! -z "${arg_V:-${arg_B}}" ] &&
  emergency "Please pass only one of {-B, -D, -p, -P, -U, -V} or a longer equivalent (multiple detected). [exit 103]"

[ ! -z "${arg_V}" ] && [ ! -z "${arg_B}" ] &&
  emergency "Please pass only one of {-B, -D, -p, -P, -U, -V} or a longer equivalent (multiple detected). [exit 104]"

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

info  "-b (--install-branch):   ${arg_b}"
info  "-B (--list-branches):    ${arg_B}"
info  "-c (--with-c):           ${arg_c}"
info  "-C (--with-cxx):         ${arg_C}"
info  "-d (--debug):            ${arg_d}"
info  "-D (--print-downloader): ${arg_D}"
info  "-e (--verbose):          ${arg_e}"
info  "-f (--with-fortran):     ${arg_f}"
info  "-h (--help):             ${arg_h}"
info  "-i (--install-prefix):   ${arg_i}"
info  "-I (--install-version):  ${arg_I}"
info  "-j (--num-threads):      ${arg_j}"
info  "-l (--list-packages):    ${arg_l}"
info  "-m (--with-cmake):       ${arg_m}"
info  "-M (--with-mpi):         ${arg_M}"
info  "-n (--no-color):         ${arg_n}"
info  "-o (--only-download):    ${arg_o}"
info  "-p (--package):          ${arg_p}"
info  "-P (--print-path):       ${arg_P}"
info  "-r (--prefix-root):      ${arg_r}"
info  "-u (--from-url):         ${arg_u}"
info  "-U (--print-url):        ${arg_U}"
info  "-v (--version):          ${arg_v}"
info  "-V (--print-version):    ${arg_V}"
info  "-y (--yes-to-all):       ${arg_y}"
info  "-Z (--bootstrap):        ${arg_Z}"
}
# This file is organized into three sections:
# 1. Command-line argument and environment variable processing.
# 2. Function definitions.
# 3. Main body.
# The script depends on several external programs, including a second script that
# builds prerequisite software.  Building prerequisites requires network access
# unless tar balls of the prerequisites are present.

# TODO:
# 1. Collapse the body of the main conditional branches in the find_or_install function
#    into one new function.
# 2. Verify that script-installed packages meet the minimum version number.
# 3. Add a script_transfer function to collapse the stack_pop x; stack_push z y
#    pattern into one statement
# 4. Consider adding mpich and cmake to the dependency stack before passing them to
#    find_or_install to make the blocks inside find_or_install more uniform.
#    Alternatively, check the dependency stacks for the package before entering the
#    main conditional blocks in find_or_install.
#


# __________ Process command-line arguments and environment variables _____________

this_script="$(basename "$0")"
export this_script

export install_path=${arg_i%/}
export prefix_root="${arg_r:-}"

export num_threads="${arg_j}"
info "num_threads=\"${arg_j}\""

export opencoarrays_src_dir="${OPENCOARRAYS_SRC_DIR}"
info "opencoarrays_src_dir=${OPENCOARRAYS_SRC_DIR}"

export build_path="${opencoarrays_src_dir}"/prerequisites/builds
info "build_path=\"${opencoarrays_src_dir}\"/prerequisites/builds"

export build_script="${opencoarrays_src_dir}"/prerequisites/build.sh
info "build_script=\"${opencoarrays_src_dir}\"/prerequisites/build.sh"

# ___________________ Define functions for use in the Main Body ___________________

# Include stack management functions
#. ./prerequisites/stack.sh
# shellcheck source=./prerequisites/stack.sh
source $opencoarrays_src_dir/prerequisites/stack.sh
stack_new dependency_pkg
stack_new dependency_exe
stack_new dependency_path
stack_new script_installed

# shellcheck source=./prerequisites/install-functions/find_or_install.sh
source $opencoarrays_src_dir/prerequisites/install-functions/find_or_install.sh

# shellcheck source=./prerequisites/install-functions/print_header.sh
source $opencoarrays_src_dir/prerequisites/install-functions/print_header.sh

# shellcheck source=./prerequisites/install-functions/build_opencoarrays.sh
source $opencoarrays_src_dir/prerequisites/install-functions/build_opencoarrays.sh

# shellcheck source=./prerequisites/install-functions/report_results.sh
source $opencoarrays_src_dir/prerequisites/install-functions/report_results.sh

# shellcheck source=./prerequisites/install-functions/install-xcode-clt.sh
source "${opencoarrays_src_dir}/prerequisites/install-functions/install-xcode-clt.sh"

# ___________________ End of function definitions for use in the Main Body __________________


# ________________________________ Start of the Main Body ___________________________________

if [[ "${arg_v}" == "${__flag_present}" || "${arg_V}" == "opencoarrays" ]]; then

  # Print script copyright & version if invoked with -v, -V, or
  # --version argument git attributes handle .VERSION, making it more
  # robust, but fallback version is still manually included. Search
  # for the first version string we encounter and extract it using sed:
  opencoarrays_version=$(sed -n '/[0-9]\{1,\}\(\.[0-9]\{1,\}\)\{1,\}/{s/^\([^.]*\)\([0-9]\{1,\}\(\.[0-9]\{1,\}\)\{1,\}\)\(.*\)/\2/p;q;}' "${opencoarrays_src_dir%/}/.VERSION")
  if [[ "${arg_v}" == "${__flag_present}" ]]; then
    echo "OpenCoarrays ${opencoarrays_version}"
    echo ""
    echo "OpenCoarrays installer"
    echo "Copyright (C) 2015-2016 Sourcery, Inc."
    echo "Copyright (C) 2015-2016 Sourcery Institute"
    echo ""
    echo "OpenCoarrays comes with NO WARRANTY, to the extent permitted by law."
    echo "You may redistribute copies of ${this_script} under the terms of the"
    echo "BSD 3-Clause License.  For more information about these matters, see"
    echo "http://www.sourceryinstitute.org/license.html"
    echo ""
  elif [[ "${arg_V}" == "opencoarrays" ]]; then
    echo "${opencoarrays_version//[[:space:]]/}"
  fi

elif [[ ! -z "${arg_D:-${arg_P:-${arg_U:-${arg_V:-${arg_B}}}}}" ||  "${arg_l}" == "${__flag_present}" ]]; then

  # Delegate to build.sh for the packages it builds
  build_arg=${arg_B:-${arg_D:-${arg_P:-${arg_U:-${arg_V:-${arg_p}}}}}}
  [ ! -z "${arg_B}" ] && build_flag="-B"
  [ ! -z "${arg_D}" ] && build_flag="-D"
  [ ! -z "${arg_P}" ] && build_flag="-P"
  [ ! -z "${arg_U}" ] && build_flag="-U"
  [ ! -z "${arg_V}" ] && build_flag="-V"
  [ "${arg_l}" == "${__flag_present}" ] && build_flag="-l"

  if [[ "${arg_P}" == "opencoarrays" ]]; then

    version="$("${opencoarrays_src_dir}/install.sh" -V opencoarrays)"
    echo "${install_path%/}/opencoarrays/${version}"

  else

    info "Invoking build script with the following command:"
    info "\"${opencoarrays_src_dir}\"/prerequisites/build.sh \"${build_flag}\" \"${build_arg}\""
    "${opencoarrays_src_dir}"/prerequisites/build.sh "${build_flag}" "${build_arg}"

    # Add lines other packages the current script builds
    if [[ "${arg_l}" == "${__flag_present}" ]]; then
      echo "opencoarrays (version $("${opencoarrays_src_dir}/install.sh" -V opencoarrays))"
      echo "ofp (version: ofp-sdf for OS X )"
    fi
  fi

elif [[ "${arg_p:-}" == "opencoarrays" ]]; then

  if [[ "${arg_o}" == "${__flag_present}" ]]; then

    # Download all prerequisites and exit
    info "source ${OPENCOARRAYS_SRC_DIR}/prerequisites/install-functions/download-all-prerequisites.sh"
    source "${OPENCOARRAYS_SRC_DIR}/prerequisites/install-functions/download-all-prerequisites.sh"
    info "download_all_prerequisites"
    download_all_prerequisites

  else

    # Install Xcode command line tools (CLT) if on macOS and if needed
    maybe_install_xcodeCLT
    # Install OpenCoarrays

    cd prerequisites || exit 1
    installation_record=install-opencoarrays.log
    # shellcheck source=./prerequisites/build-functions/set_SUDO_if_needed_to_write_to_directory.sh
    source "${OPENCOARRAYS_SRC_DIR:-}/prerequisites/build-functions/set_SUDO_if_needed_to_write_to_directory.sh"
    version="$("${opencoarrays_src_dir}/install.sh" -V opencoarrays)"

    if [[ -z "${prefix_root:-}" ]]; then
      set_SUDO_if_needed_to_write_to_directory "${install_path}"
    else
      set_SUDO_if_needed_to_write_to_directory "${prefix_root}"
    fi

    if grep -s -q Microsoft /proc/version  ; then
      info "Windows Subsystem for Linux detected.  Invoking windows-install.sh with the following command:"
      info "\"${opencoarrays_src_dir}\"/prerequisites/install-functions/windows-install.sh \"$@\""
      "${opencoarrays_src_dir}"/prerequisites/install-functions/windows-install.sh "$@"
    else
      # Using process substitution "> >(...) -" instead of piping to tee via "2>&1 |" ensures that
      # report_results gets the FC value set in build_opencoarrays
      # Source: http://stackoverflow.com/questions/8221227/bash-variable-losing-its-value-strange
      build_opencoarrays > >( tee ../"${installation_record}" ) -
      report_results 2>&1 | tee -a ../"${installation_record}"
    fi
  fi

elif [[ "${arg_p:-}" == "ofp" ]]; then

  # Install Xcode command line tools (CLT) if on macOS and if needed
  maybe_install_xcodeCLT

  info "Invoking Open Fortran Parser build script with the following command:"
  info "\"${opencoarrays_src_dir}\"/prerequisites/install-ofp.sh"
  "${opencoarrays_src_dir}"/prerequisites/install-ofp.sh

elif [[ ! -z "${arg_p:-}" ]]; then

  # Install Xcode command line tools (CLT) if on macOS and if needed
  maybe_install_xcodeCLT

  info "Invoking build script with the following command:"
  info "\"${opencoarrays_src_dir}\"/prerequisites/build.sh ${*:-}"
  "${opencoarrays_src_dir}"/prerequisites/build.sh "${@:-}"

fi
# ____________________________________ End of Main Body ____________________________________________
