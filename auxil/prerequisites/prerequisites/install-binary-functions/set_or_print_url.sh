# If -p, -D, -P, or -V specifies a package, set package_url
# If -U specifies a package, print the package_url and exit with normal status
# shellcheck disable=SC2154
set_or_print_url()
{
  # Verify requirements
  if [[ ! -z "${arg_U}" ]]; then
    if [[ "${arg_N:-}" == "${__flag_present}" || ! -z "${arg_D:-${arg_p:-${arg_P:-${arg_V}}}}" ]]; then
      emergency "Please pass only one of {-D, -N, -p, -P, -U, -V} or a longer equivalent (multiple detected)."
    fi
  fi
  # Get package name from argument passed with  -p, -D, -P, -V, or -U
  package_name="${arg_D:-${arg_p:-${arg_P:-${arg_U:-${arg_V}}}}}"

  package_url_head=(
    "strategoxt-superbundle;https://github.com/sourceryinstitute/opencoarrays/files/212509/"
  )
  for package in "${package_url_head[@]}" ; do
     KEY="${package%%;*}"
     VALUE="${package##*;}"
     if [[ "${package_name}" == "${KEY}" ]]; then
       # We recognize the package name so we set the URL head:
       url_head="${VALUE}"
       break
     fi
  done

  package_url_tail=(
    "strategoxt-superbundle;strategoxt-superbundle-${version_to_build}.tar.gz"
  )
  for package in "${package_url_tail[@]}" ; do
     KEY="${package%%;*}"
     VALUE="${package##*;}"
     if [[ "${package_name}" == "${KEY}" ]]; then
       # We recognize the package name so we set the URL tail:
       url_tail="${VALUE}"
       break
     fi
  done

  if [[ -z "${url_head:-}" || -z "${url_tail:-}" ]]; then
    emergency "Package ${package_name:-} not recognized.  Use -l or --list-packages to list the allowable names."
  fi

  package_url="${url_head}""${url_tail}"

  # If a printout of the package URL was requested, then print it and exit with normal status
  if [[ ! -z "${arg_U:-}" ]]; then
    printf "%s\n" "${package_url}"
    exit 0
  fi
}
