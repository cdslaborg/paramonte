#!/usr/bin/env bash
#
# windows-install.sh
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

[[ -z "${LOG_LEVEL:-}" ]] && emergency "Cannot continue without LOG_LEVEL. "

# shellcheck disable=SC2154
if [[ "${arg_d}" == "${__flag_present}" ]]; then
   print_debug_only=7
   if [[ "$(( LOG_LEVEL < print_debug_only ))" -ne 0 ]]; then
     debug "Supressing info and debug messages: one of {-l, -v, -P, -U, -V, -D} present."
     suppress_info_debug_messages
   fi
fi

# Get linux_distribution name
# shellcheck disable=SC2154
{
info "__file: ${__file}"
info "__dir: ${__dir}"
info "__base: ${__base}"
info "__os: ${__os}"
info "__usage: ${__usage}"
info "LOG_LEVEL: ${LOG_LEVEL}"

info  "-c (--with-c):     [arg] ${arg_c}"
info  "-C (--with-cxx):   [arg] ${arg_C}"
info  "-d (--debug):            ${arg_d}"
info  "-e (--verbose):          ${arg_e}"
info  "-f (--with-fortran):     ${arg_e}"
info  "-i (--install-prefix):   ${arg_i}"
info  "-j (--num-threads):      ${arg_j}"
info  "-m (--with-cmake): [arg] ${arg_m}"
info  "-n (--no-color):         ${arg_n}"
info  "-v (--version):          ${arg_v}"
info  "-V (--version-number):   ${arg_V}"
}

# __________ Process command-line arguments and environment variables _____________

this_script="$(basename "${0}")"
export this_script
debug "this_script=\"${this_script}\""

export install_prefix="${arg_i%/:-${PWD}/prerequisites/installations}"
info "install_prefix=\"${install_prefix}\""

export num_threads="${arg_j}"
info "num_threads=\"${arg_j}\""

opencoarrays_version=$(sed -n '/[0-9]\{1,\}\(\.[0-9]\{1,\}\)\{1,\}/{s/^\([^.]*\)\([0-9]\{1,\}\(\.[0-9]\{1,\}\)\{1,\}\)\(.*\)/\2/p;q;}' "${OPENCOARRAYS_SRC_DIR%/}/.VERSION")

export build_path="${OPENCOARRAYS_SRC_DIR%/}"/prerequisites/builds/opencoarrays/${opencoarrays_version}
info "build_path=\"${build_path}\""

export CMAKE="${arg_m:-cmake}"

verify_this_is_ubuntu()
{
  if [[ ${__os} != "Linux"  ]]; then
    emergency "${__os} not supported: this script is intended for use in Windows Subsystem for Linux "
  fi
  linux_standard_base_i=$(lsb_release -i)
  untrimmed_name=${linux_standard_base_i##*Distributor ID:}
  linux_distribution="${untrimmed_name//[[:space:]]/}"
  info "Linux distribution: ${linux_distribution}"
  if [[ "${linux_distribution:-}" != "Ubuntu" ]]; then
    info "Please run this script inside the Windows Subsystem for Linux (WSL) Ubuntu 16.04 or later,"
    info "which might require joining the Windows Insider Preview program and selecting 'Fast' updates."
    emergency "Error: see above."
  fi
}
verify_this_is_ubuntu

# Ubuntu 16.04 apt-get installs gfortran 5.4.0 or later, which is acceptable for many uses of OpenCoarrays
verify_acceptable_release_number()
{
  linux_standard_base_r=$(lsb_release -r)
  untrimmed_name=${linux_standard_base_r##*Release:}
  release_number="${untrimmed_name//[[:space:]]/}"
  major_release="${release_number%%.*}"
  minor_release="${release_number##*.}"
  info "Release: ${major_release}.${minor_release}"
  if [[ ${major_release} -lt 16  ]]; then
    emergency "Please upgrade to Windows Subsystem for Linux (WSL) Ubuntu 16.04 or later."
  elif [[ ${major_release} -eq 16  ]]; then
    if [[ ${minor_release} -lt "04"  ]]; then
      emergency "Please upgrade to Windows Subsystem for Linux (WSL) Ubuntu 16.04 or later."
    fi
  fi
}
verify_acceptable_release_number

if [[ "${arg_V}" == "${__flag_present}" ]]; then
    # Print just the version number
    info "${opencoarrays_version}"

elif [[ "${arg_v}" == "${__flag_present}" ]]; then

    # Print copyright info and version number
    info "OpenCoarrays ${opencoarrays_version}"
    info ""
    info "OpenCoarrays installer"
    info "Copyright (C) 2015-2017 Sourcery, Inc."
    echo "Copyright (C) 2015-2017 Sourcery Institute"
    info ""
    info "OpenCoarrays comes with NO WARRANTY, to the extent permitted by law."
    info "You may redistribute copies of ${this_script} under the terms of the"
    echo "BSD 3-Clause License.  For more information about these matters, see"
    info "http://www.sourceryinstitute.org/license.html"
    info ""
else

  export FC=${arg_f:-gfortran}
  export CC=${arg_c:-gcc}
  export CXX=${arg_C:-g++}

  # Check for and, if necessary, install OpenCoarrays prerequisites

  if ! type "${CMAKE}" >& /dev/null; then
    sudo apt-get install cmake
  fi
  if ! type "${CXX}" >& /dev/null; then
    sudo apt-get install g++
  fi
  if ! type "${FC}" >& /dev/null; then
    sudo apt-get install gfortran
  fi

  if ! type mpifort >& /dev/null; then
    sudo apt-get install mpich
  fi

  set_SUDO_if_needed_to_write_to_install_dir()
  {
    info "Checking whether the directory ${install_prefix} exists... "
    if [[ -d "${install_prefix}" ]]; then
      info "yes"
      info "Checking whether I have write permissions to ${install_prefix} ... "
      if [[ -w "${install_prefix}" ]]; then
        info "yes"
      else
        info "no"
        SUDO="sudo"
      fi
    else
      info "no"
      info "Checking whether I can create ${install_prefix} ... "
      if mkdir -p "${install_prefix}" >& /dev/null; then
        info "yes."
      else
        info "no."
        SUDO="sudo"
      fi
    fi
  }
  set_SUDO_if_needed_to_write_to_install_dir

  # Install OpenCoarrays

  if [[ -d "${build_path}" ]]; then
    rm -rf "${build_path}"
  fi
  mkdir -p "${build_path}"
  cd "${build_path}" || exit 25
  info "Configuring OpenCoarrays with the following command:"
  info "FC=\"${FC}\" CC=\"${CC}\"  \"${CMAKE}\" \"${OPENCOARRAYS_SRC_DIR}\" -DCMAKE_INSTALL_PREFIX=\"${install_prefix}\""
  FC="${FC}" CC="${CC}" "${CMAKE}" "${OPENCOARRAYS_SRC_DIR}" -DCMAKE_INSTALL_PREFIX="${install_prefix}"
  info "Building OpenCoarrays with the following command:"
  info "make -j ${arg_j}"
  make -j "${arg_j}"
  info "Installing OpenCoarrays with the following command:"
  info "${SUDO:-} make install"
  ${SUDO:-} make install
  if [[ -f "${install_prefix}"/lib/libcaf_mpi.a && -f "${install_prefix}/bin/caf"  && -f "${install_prefix}/bin/cafrun"  ]]; then
    info "OpenCoarrays has been installed in"
    info "${install_prefix}"
  else
    info "Something went wrong. OpenCoarrays is not in the expected location:"
    emergency "${install_prefix}"
  fi
  # See http://stackoverflow.com/questions/31057694/gethostbyname-fail-after-switching-internet-connections/31222970
  loopback_line=$(grep "${NAME}" /etc/hosts)
  if [[ -z "${loopback_line:-}" ]]; then
    info "To ensure the correct functioning of MPI, please add the following line to your /etc/hosts file:"
    info "127.0.0.1 ${NAME}"
  fi
fi
