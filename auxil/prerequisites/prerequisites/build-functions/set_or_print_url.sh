# shellcheck shell=bash
# If -p, -D, -P, or -V specifies a package, set package_url
# If -U specifies a package, print the package_url and exit with normal status
# shellcheck disable=SC2154
set_or_print_url()
{
  # Verify requirements
  [[ -n "${arg_U}" && -n "${arg_D:-${arg_p:-${arg_P:-${arg_V:-${arg_B}}}}}" ]] &&
    emergency "Please pass only one of {-B, -D, -p, -P, -U, -V} or a longer equivalent (multiple detected)."

  # Get package name from argument passed with  -p, -D, -P, -V, or -U
  package_to_build="${arg_p:-${arg_D:-${arg_P:-${arg_U:-${arg_V:-${arg_B}}}}}}"

if [[ -n "${arg_u:-}"  ]]; then
  # User specified a URL from which to download the package
  url_tail="${arg_u##*/}" # set url_tail to text_after_final_slash, greedy expansion needed
  url_head="${arg_u%${url_tail}}" # set url_head text before url_tail
else

  if [[ "${package_to_build}" == 'cmake' ]]; then
    major_minor="${version_to_build%.*}"
  elif [[ "${package_to_build}" == "gcc" ]]; then
    if [[ -z "${arg_b:-${arg_B}}" ]]; then
      gcc_url_head="https://ftpmirror.gnu.org/gcc/gcc-${version_to_build}/"
    else
      gcc_url_head="svn://gcc.gnu.org/svn/gcc/"
    fi
  fi
  package_url_head=(
    "gcc;${gcc_url_head-}"
    "wget;https://ftpmirror.gnu.org/gnu/wget/"
    "m4;https://ftpmirror.gnu.org/gnu/m4/"
    "pkg-config;https://pkgconfig.freedesktop.org/releases/"
    "mpich;https://www.mpich.org/static/downloads/${version_to_build-}/"
    "flex;https://sourceforge.net/projects/flex/files/"
    "make;https://ftpmirror.gnu.org/gnu/make/"
    "bison;https://ftpmirror.gnu.org/gnu/bison/"
    "cmake;https://www.cmake.org/files/v${major_minor:-}/"
    "subversion;https://www.eu.apache.org/dist/subversion/"
  )
  for package in "${package_url_head[@]}" ; do
     KEY="${package%%;*}"
     VALUE="${package##*;}"
     info "KEY=${KEY}  VALUE=${VALUE}"

     if [[ "${package_to_build}" == "${KEY}" ]]; then
       # We recognize the package name so we set the URL head:
       url_head="${VALUE}"
       break
     fi
  done

  # Set differing tails for GCC release downloads versus development branch checkouts
  if [[ "${package_to_build}" == 'gcc' ]]; then
    if [[ "${fetch}" == 'svn' ]]; then
      if [[ "${version_to_build:-}" == "trunk" ]]; then
        gcc_tail="${version_to_build}"
      else
        gcc_tail="branches/${version_to_build:-}"
      fi
    else
      gcc_tail="gcc-${version_to_build}.tar.gz"
    fi
  fi
  package_url_tail=(
    "gcc;${gcc_tail-}"
    "wget;wget-${version_to_build-}.tar.gz"
    "m4;m4-${version_to_build-}.tar.bz2"
    "pkg-config;pkg-config-${version_to_build-}.tar.gz"
    "mpich;mpich-${version_to_build-}.tar.gz"
    "flex;flex-${version_to_build-}.tar.bz2"
    "bison;bison-${version_to_build-}.tar.gz"
    "make;make-${version_to_build-}.tar.bz2"
    "subversion;subversion-${version_to_build-}.tar.gz"
  )
  if [[ $(uname) == "Linux" ]]; then
    package_url_tail+=("cmake;cmake-${version_to_build}-Linux-x86_64.sh")
  else
    package_url_tail+=("cmake;cmake-${version_to_build-}.tar.gz")
  fi
  for package in "${package_url_tail[@]}" ; do
     KEY="${package%%;*}"
     VALUE="${package##*;}"
     if [[ "${package_to_build}" == "${KEY}" ]]; then
       # We recognize the package name so we set the URL tail:
       url_tail="${VALUE}"
       break
     fi
  done


  if [[ -z "${url_head:-}" || -z "${url_tail}" ]]; then
    emergency "Package ${package_name:-} not recognized.  Use --l or --list-packages to list the allowable names."
  fi

fi # end if [[ -n "${arg_u:-}"  ]]; then

  package_url="${url_head}""${url_tail}"

  # If a printout of the package URL was requested, then print it and exit with normal status
  if [[ -n "${arg_U:-}" ]]; then
    printf "%s\n" "${package_url}"
    exit 0
  fi
}
