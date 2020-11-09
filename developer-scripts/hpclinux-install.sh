#!/usr/bin/env bash

# BSD 3-Clause License
#
# -- hpclinux-install.sh
#
#  Install OpenCoarrays inside the HPCLinux distribution from hpclinux.org.
#
#  Usage: cd opencoarrays && ./prerequisites/hpclinux-install.sh"
#
#  Motivation:
#
#    On Fedora-based distributions, the OpenCoarrays installer fails during the stack-based
#    system interrogation (presumably due to problems with the prerequisites/stack.sh script).
#    At least in the case of HPCLinux, system interrogation is unnecessary because HPCLinux is
#    uber-stable and it's  safe bet that all prerequisites need to be installed so this script
#    omits the interrogation and invokes prerequisites/build.sh to build all prerequisites.
#
# Copyright (c) 2016, Sourcery Institute
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if [[ ! -d prerequisites ]]; then
  printf "Please execute this script with your present working directory"
  printf "set to the top level of the OpenCoarrays source tree."
  exit 1
fi

pushd prerequisites

# Build CMake using the system GCC (4.8.1) and prepend its bin subdirectory to the PATH
export cmake_install_path="${PWD}/installations/cmake/cmake"
./build.sh --package cmake --install-prefix "${cmake_install_path}"
if [[ -z "${PATH}" ]]; then
  export PATH="${cmake_install_path}/bin"
else
  export PATH="${cmake_install_path}/bin:${PATH}"
fi

# Build GCC 6.3.0 and prepend its bin subdirectory to the PATH
export gcc_version=6.3.0
export gcc_install_path="${PWD}/installations/gnu/${gcc_version}"
./build.sh --package gcc --install-version "${gcc_version}" --install-prefix "${gcc_install_path}"
export PATH="${gcc_install_path}/bin:${PATH}"
export LD_LIBRARY_PATH="${gcc_install_path}/lib64:${gcc_install_path}/lib:${LD_LIBRARY_PATH}"

# Build MPICH 3.2 and prepend its bin subdirectory to the PATH
export mpich_install_path="${PWD}/installations/mpich"
./build.sh --package mpich --install-prefix "${mpich_install_path}" --num-threads 4
export PATH="${mpich_install_path}/bin:${PATH}"

popd # return to top level of OpenCoarrays source tree

# Build OpenCoarrays
if [[ -d build ]]; then
  printf 'Old build subdirectory found. Ok to delete the "build" subdirectory? (Y/n)'
  read -r delete_build

  if [[ "${delete_build}" == "n" || "${delete_build}" == "no" ]]; then
    printf "n\n"
    printf "Please rename or delete the build subdirectory and restart this script. Aborting. [exit 10]\n"
    exit 10
  else # permission granted to delete build subdirectory
    printf "Y\n"
  fi
  rm -rf build
fi
mkdir build
pushd build

export opencoarrays_install_path="${PWD}/prerequisites/installations/opencoarrays"
FC=gfortran CC=gcc cmake .. -DCMAKE_INSTALL_PREFIX="${opencoarrays_install_path}"
make
ctest --output-on-failure --schedule-random
make install
export PATH="${opencoarrays_install_path}/bin:${PATH}"

popd # return to top level of OpenCoarrays source tree
