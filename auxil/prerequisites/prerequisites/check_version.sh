#!/bin/bash
#
# check_version.sh
#
# -- Verify whether an OpenCoarrays prerequisite meets the required minimum version number.
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

# Interpret the first argument as the package executable name
package="$1"
# Interpret the second argument as the minimimum acceptable package version number
minimum_version="$2"
# Interpret the third argument as indicating the desired verbosity
verbose=$3

this_script=$(basename "$0")

usage()
{
    echo ""
    echo " $this_script - Bash script for verifyin minimum version numbers for OpenCoarrays prerequisites"
    echo ""
    echo " Usage (optional arguments in square brackets): "
    echo "      $this_script [<option>]"
    echo "      $this_script  <package-name> <minimum-version-number> [--verbose]"
    echo ""
    echo " Options:"
    echo "   --help, -h           Show this help message"
    echo "   --list , -l          List the packages whose versions this script can verify"
    echo ""
    echo " Examples:"
    echo ""
    echo "   $this_script cmake 3.4.0"
    echo "   $this_script flex 2.6.0 --verbose"
    echo "   $this_script flex $(./build.sh -V flex )"
    echo "   $this_script --help"
    echo "   $this_script --list"
    echo ""
    echo "[exit 10]"
    exit 10
}

# Print usage information and exit if script is invoked without arguments or with --help or -h as the first argument
if [ $# == 0 ]; then
  usage | less
  exit 20

elif [[ $1 == '--help' || $1 == '-h' ]]; then
  usage | less
  exit 0

elif [[ $1 == '--list' || $1 == '-l' ]]; then
 echo "$this_script currently verifies minimum version numbers for the following OpenCoarrays prerequisites:"
 echo "   cmake, flex, bison, m4"
 exit 40

elif [[ $1 == '-v' || $1 == '-V' || $1 == '--version' ]]; then
  # Print script copyright if invoked with -v, -V, or --version
  # argument Extract version from .VERSION file. Git automatically
  # inserts various things when performing git archive, so ensure we
  # extract just the first version string
  opencoarrays_version=$(sed -n '/[0-9]\{1,\}\(\.[0-9]\{1,\}\)\{1,\}/{s/^\([^.]*\)\([0-9]\{1,\}\(\.[0-9]\{1,\}\)\{1,\}\)\(.*\)/\2/p;q;}' ../.VERSION)
  echo "opencoarrays $opencoarrays_version"
  echo ""
  echo "OpenCoarrays prerequisite version verifier"
  echo "Copyright (C) 2015-2016 Sourcery, Inc."
  echo "Copyright (C) 2015-2016 Sourcery Institute"
  echo ""
  echo "OpenCoarrays comes with NO WARRANTY, to the extent permitted by law."
  echo "You may redistribute copies of $this_script under the terms of the"
  echo "BSD 3-Clause License.  For more information about these matters, see"
  echo "http://www.sourceryinstitute.org/license.html"
  echo ""
fi

package_version_header=$($package --version | head -1)
if [[ "$verbose" == "--verbose" ]]; then
  echo "$package_version_header"
fi

# Extract the text after the final space:
version=${package_version_header##* }
major=${version%%.*}
minor_patch=${version#*.}
minor=${minor_patch%%.*}
patch=${version##*.}
if [[ "$verbose" == "--verbose" ]]; then
  echo "$version = $major . $minor . $patch"
fi

# Extract the text after the final space:
min_major=${minimum_version%%.*}
min_minor_patch=${minimum_version#*.}
min_minor=${min_minor_patch%%.*}
min_patch=${minimum_version##*.}
if [[ "$verbose" == "--verbose" ]]; then
  echo "$minimum_version = $min_major . $min_minor . $min_patch"
fi

if [ "$(( major < min_major ))" -ne 0 ]; then
  if [[ $verbose == "--verbose" ]]; then
    echo "$major < $min_major"
  fi
  exit 10
elif [[ $verbose == "--verbose" ]]; then
  echo "$major >= $min_major"
fi

if [ "$(( minor < min_minor ))" -ne 0 ]; then
  if [[ $verbose == "--verbose" ]]; then
    echo "$minor < $min_minor"
  fi
  exit 20
elif [[ $verbose == "--verbose" ]]; then
  echo "$minor >= $min_minor"
fi

if [ "$(( patch < min_patch ))" -ne 0 ]; then
  if [[ $verbose == "--verbose" ]]; then
    echo "$patch < $min_patch"
  fi
  exit 30
elif [[ $verbose == "--verbose" ]]; then
  echo "$patch >= $min_patch"
fi
exit 0
