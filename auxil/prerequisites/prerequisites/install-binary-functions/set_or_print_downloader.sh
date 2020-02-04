# If -p, -P, -U, or -V specifies a package, set fetch variable
# If -D specifies a package, print "${fetch}" and exit with normal status
# If -l is present, list all packages and versions and exit with normal status
# shellcheck disable=SC2154
set_or_print_downloader()
{
  # Verify requirements
  if [[ ! -z "${arg_D}" ]]; then
    if  [[ "${arg_N:-}" == "${__flag_present}" || ! -z "${arg_p:-${arg_P:-${arg_U:-${arg_V}}}}" ]]; then
      emergency "Please pass only one of {-D, -N, -p, -P, -U, -V} or a longer equivalent (multiple detected)."
    fi
  fi
  package_name="${arg_D:-${arg_p:-${arg_P:-${arg_U:-${arg_V}}}}}"
  if [[ $(uname) == "Darwin" ]]; then
    wget_or_curl=curl
  else
    wget_or_curl=wget
  fi  
  # This is a bash 3 hack standing in for a bash 4 hash (bash 3 is the lowest common
  # denominator because, for licensing reasons, OS X only has bash 3 by default.)
  # See http://stackoverflow.com/questions/1494178/how-to-define-hash-tables-in-bash
  package_fetch=(
    "strategoxt-superbundle:${wget_or_curl}"
  )
  for package in "${package_fetch[@]}" ; do
     KEY="${package%%:*}"
     VALUE="${package##*:}"
     if [[ "${package_name}" == "${KEY}" ]]; then
       # We recognize the package name so we set the download mechanism:
       fetch=${VALUE}
       # If a printout of the download mechanism was requested, then print it and exit with normal status
       [[ ! -z "${arg_D}" ]] && printf "%s\n" "${fetch}" && exit 0
       break # exit the for loop
     fi
  done
  if [[ -z "${fetch:-}" ]]; then
    emergency "Package ${package_name:-} not recognized.  Use -l or --list-packages to list the allowable names."
  fi
}
