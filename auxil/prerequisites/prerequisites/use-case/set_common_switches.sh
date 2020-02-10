# BASH3 Boilerplate
#
#  set_common_switches.sh
#
#  - Sets variables that are useful in conjunction with other bash3boilerplate features
#  - Contains commands extracted from the bash3boilerplate main.sh as part of a refactoring
#    to facilitate wholesale reuse of main.sh's contents of without modification.
#
# Usage (as invoked in bootstrap.sh):
#
#  source set_common_switches.sh
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

### Switches (like -d for debugmode, -h for showing helppage)
#####################################################################

# shellcheck disable=SC2154
{
# debug mode
if [ "${arg_d}" = "1" ]; then
  set -o xtrace
  LOG_LEVEL="7"
  export LOG_LEVEL
fi

# verbose mode
if [ "${arg_e}" = "1" ]; then
  set -o verbose
fi

# help mode
if [ "${arg_h}" = "1" ]; then
  # Help exists with code 1
  help "Help using ${0}"
fi
}
