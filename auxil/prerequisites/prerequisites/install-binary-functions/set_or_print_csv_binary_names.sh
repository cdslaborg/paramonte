# If -p, -P, -U, or -V specifies a package, set fetch variable
# If -D specifies a package, print "${fetch}" and exit with normal status
# If -l is present, list all packages and versions and exit with normal status
# shellcheck disable=SC2154
set_or_print_csv_binary_names()
{
  # Verify requirements
  [ "${arg_N:-}" == "${__flag_present}" ] && [ ! -z "${arg_D:-${arg_p:-${arg_P:-${arg_U:-${arg_V}}}}}" ] &&
    emergency "Please pass only one of {-D, -N, -p, -P, -U, -V} or a longer equivalent (multiple detected)."

  # This is a bash 3 hack standing in for a bash 4 hash (bash 3 is the lowest common
  # denominator because, for licensing reasons, OS X only has bash 3 by default.)
  # See http://stackoverflow.com/questions/1494178/how-to-define-hash-tables-in-bash
  package_install_names=(
    "strategoxt-superbundle:aterm,sdf2-bundle,strategoxt"
  )
  for package in "${package_install_names[@]}" ; do
     KEY="${package%%:*}"
     VALUE="${package##*:}"
     if [[ "${package_name}" == "${KEY}" ]]; then
       # We recognize the package name so we set the download mechanism:
       install_names=${VALUE}
       # If a printout of the download mechanism was requested, then print it and exit with normal status
       [[ "${arg_N:-}" == "${__flag_present}" ]] && printf "%s\n" "${install_names}" && exit 0
       break # exit the for loop
     fi
  done
  if [[ -z "${install_names:-}" ]]; then
    emergency "Package ${package_name:-} not recognized.  Use -l or --list-packages to list the allowable names."
  fi
}
