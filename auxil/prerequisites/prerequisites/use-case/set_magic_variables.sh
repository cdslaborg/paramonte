# BASH3 Boilerplate
#
#  set_magic_variables.sh
#
#  - Sets the variables __dir, __file, __filename, __base, and __os
#  - Defines a function containing commands extracted from the bash3boilerplate
#    main.sh as part of a refactoring to facilitate wholesale reuse of main.sh's
#    contents of without modification.
#
#  Usage (as invoked in bootstrap.sh):
#
#    source set_magic_variables.sh "$(caller 0)"
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

# shellcheck disable=SC2016
[ -z "${1}" ] && echo 'Usage: source set_magic_variables.sh "$(caller 0)"'
# shellcheck disable=SC2034
function set_magic_variables(){
  text_after_final_space="${1##* }"
  __dir="$(cd "$(dirname "${text_after_final_space}")" && pwd)"
  __file="${__dir}/$(basename "${text_after_final_space}")"
  __filename="$(basename "${text_after_final_space}")"
  __base="$(basename "${__file}" .sh)"
  __os="Linux"
  if [[ "${OSTYPE:-}" == "darwin"* ]]; then
    __os="OSX"
  fi
  __usage="${__usage:-${__file}-usage}"
}
set_magic_variables "${@}"
