# Unpack if the unpacked tar ball is not in the present working directory
# shellcheck disable=SC2154
unpack_if_necessary()
{
  if [[ "${fetch}" == "svn" || "${fetch}" == "git" ]]; then
    package_source_directory="${version_to_build}"
  else
    if [[ "${url_tail}" == *tar.gz || "${url_tail}" == *tar.bz2 || "${url_tail}" == *.tar.xz ]]; then
      info "Unpacking ${url_tail}."
      info "pushd ${download_path}"
      pushd "${download_path}"
      info "Unpack command: tar xf ${url_tail}"
      tar xf "${url_tail}"
      info "popd"
      popd
    else
      info "Skipping unpacking because ${url_tail} is not a compressed archive (matching one of *.tar.{gz,bz2,xz}). "
    fi
    # shellcheck disable=SC2034
    package_source_directory="${package_name}-${version_to_build}"
  fi
}
