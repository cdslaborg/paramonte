#!/usr/bin/env bash
# BASH3 Boilerplate
#
#  bootstrap.sh
#
#  - Exports bash3boilerplate features and variables to the invoking script
#  - Invokes functions containing commands extracted from the bash3boilerplate
#    main.sh as part of a refactoring to facilitate wholesale reuse of main.sh's
#    contents of without modification.
#
# Usage (as invoked in my-script.sh):
#
#   source bootstrap.sh "${@}"
#
# More info:
#
#  - https://github.com/kvz/bash3boilerplate
#  - http://kvz.io/blog/2013/02/26/introducing-bash3boilerplate/
#
# Version: 2.1.0
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

# shellcheck source=./set_environment_and_color.sh
source "${B3B_USE_CASE}"/set_environment_and_color.sh         # turn on errexit, nounset, pipefail, default log level
# shellcheck source=./set_magic_variables.sh
source "${B3B_USE_CASE}"/set_magic_variables.sh "$(caller 0)" # set __dir, __file, __filename, __base, __os
# shellcheck source=./define_functions.sh
source "${B3B_USE_CASE}"/define_functions.sh                  # help/usage function and debug/info output functions
# shellcheck source=./parse_command_line.sh
source "${B3B_USE_CASE}"/parse_command_line.sh "${@:-}"       # parse the command line
# shellcheck source=./set_common_switches.sh
source "${B3B_USE_CASE}"/set_common_switches.sh               # provide defaults for -h, -V, and -d
