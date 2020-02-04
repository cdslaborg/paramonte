# BASH3 Boilerplate
#
#  define_functions.sh
#
#  - Defines helper functions containing commands extracted from the
#    bash3boilerplate main.sh as part of a refactoring to facilitate
#    wholesale reuse of main.sh's contents of without modification.
#
# Usage (as invoked in bootstrap.sh):
#
#   source define_functions.sh
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
#
# Licensed under MIT
# Copyright (c) 2013 Kevin van Zonneveld (http://kvz.io)

### Functions
#####################################################################

# shellcheck disable=SC2034
function _fmt ()      {
  local color_debug="\x1b[35m"
  local color_info="\x1b[32m"
  local color_notice="\x1b[34m"
  local color_warning="\x1b[33m"
  local color_error="\x1b[31m"
  local color_critical="\x1b[1;31m"
  local color_alert="\x1b[1;33;41m"
  local color_emergency="\x1b[1;4;5;33;41m"
  local colorvar=color_$1

  local color="${!colorvar:-$color_error}"
  local color_reset="\x1b[0m"
  if [ "${NO_COLOR}" = "true" ] || [[ "${TERM:-}" != "xterm"* ]] || [ -t 1 ]; then
    # Don't use colors on pipes or non-recognized terminals
    color=""; color_reset=""
  fi
  echo -e "$(date -u +"%Y-%m-%d %H:%M:%S UTC") ${color}$(printf "[%9s]" "${1}")${color_reset}";
}

# The block of single-line functions below all print to STDERR,
# leaving STDOUT for piping machine readable information to other
# software. Above each such function is a commented demonstration
# of its usage.  Execution continues after an invocation of each
# function, except the "emergency" function, which causes
# termination with a non-zero exit status.

# shellcheck disable=SC2015
{
# emergency "A \"panic\" condition usually affecting multiple apps/servers/sites. At this level it would usually notify all tech staff on call."
function emergency () {                             echo "$(_fmt emergency) ${*}" 1>&2 || true; exit 1; }
# alert "Should be corrected immediately, therefore notify staff who can fix the problem. An example would be the loss of a primary ISP connection."
function alert ()     { [ "${LOG_LEVEL}" -ge 1 ] && echo "$(_fmt alert) ${*}" 1>&2 || true; }
# critical "Should be corrected immediately, but indicates failure in a primary system, an example is a loss of a backup ISP connection."
function critical ()  { [ "${LOG_LEVEL}" -ge 2 ] && echo "$(_fmt critical) ${*}" 1>&2 || true; }
# error "Non-urgent failures, these should be relayed to developers or admins; each item must be resolved within a given time."
function error ()     { [ "${LOG_LEVEL}" -ge 3 ] && echo "$(_fmt error) ${*}" 1>&2 || true; }
# warning "Warning messages, not an error, but indication that an error will occur if action is not taken, e.g. file system 85% full - each item must be resolved within a given time. This is a debug message"
function warning ()   { [ "${LOG_LEVEL}" -ge 4 ] && echo "$(_fmt warning) ${*}" 1>&2 || true; }
# notice "Events that are unusual but not error conditions - might be summarized in an email to developers or admins to spot potential problems - no immediate action required."
function notice ()    { [ "${LOG_LEVEL}" -ge 5 ] && echo "$(_fmt notice) ${*}" 1>&2 || true; }
# info "Normal operational messages - may be harvested for reporting, measuring throughput, etc. - no action required."
function info ()      { [ "${LOG_LEVEL}" -ge 6 ] && echo "$(_fmt info) ${*}" 1>&2 || true; }
# debug "Info useful to developers for debugging the application, not useful during operations."
function debug ()     { [ "${LOG_LEVEL}" -ge 7 ] && echo "$(_fmt debug) ${*}" 1>&2 || true; }
}
function suppress_debug_messages() { export LOG_LEVEL=6; }
function suppress_info_debug_messages () { export LOG_LEVEL=5; }
function suppress_notice_info_debug_messages () { export LOG_LEVEL=4; }

function help () {
  echo "" 1>&2
  echo " ${*}" 1>&2
  echo "" 1>&2
  # shellcheck disable=SC2154
  cat "${__usage}" 1>&2
  echo "" 1>&2
  exit 1
}
export -f help
export -f emergency
export -f alert
export -f critical
export -f error
export -f warning
export -f notice
export -f info
export -f debug
export suppress_debug_messages
export suppress_info_debug_messages
export suppress_notice_info_debug_messages
