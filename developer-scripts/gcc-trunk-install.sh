#!/usr/bin/env bash

# Exit on error or use of an unset variable:
set -o errexit
set -o nounset

# Return the highest exit code in a chain of pipes:
set -o pipefail

# Print usage information and exit if this script was invoked with no arugments or with -h or --help
# as the first argument or in a location other than the top-level OpenCoarrays source directory
function usage()
{
  echo "Usage:"
  echo ""
  echo "  cd <opencoarrays-source-directory>"
  echo "  ./developer_scripts/gcc-trunk-install.sh [--patch-file <patch-file-name>] [--packages-root <installation-path>]"
  echo "or"
  echo "  ./developer_scripts/gcc-trunk-install.sh [-p <patch-file-name>] [-r <installation-path>]"
  echo ""
  echo " Square brackets surround optional arguments."
  exit 0
}
[[ "${1:-}" == "-h" || "${1:-}" == "--help"  || ! -f src/libcaf.h ]] && usage

if [[ "${1:-}" == "-r" || "${1:-}" == "--packages-root" ]]; then
  export packages_root="${2}"
  if [[ "${3:-}" == "-r" || "${3:-}" == "--packages-root" ]]; then
    export patch_file="${4}"
  fi
elif [[ "${1:-}" == "-p" || "${1:-}" == "--patch-file" ]]; then
  export patch_file="${2:-}"
  if [[ "${3:-}" == "-r" || "${3:-}" == "--packages-root" ]]; then
    export packages_root="${4}"
  fi
fi
export default_prefix="${HOME}/opt"
export packages_root="${packages_root:-${default_prefix}}"

function set_absolute_path()
{
  : "${1?'set_absolute_path: no argument provided'}"

  arg=${1}
  first_character=$(echo "${arg}" | cut -c1-1)
  if [[ "${first_character}" == "/" ]]; then
    absolute_path="${arg}"
  else
    absolute_path="${PWD%%/}/${arg}"
  fi
}
if [[ ! -z "${patch_file:-}" ]]; then
  set_absolute_path "${patch_file}"
fi

### Define functions
function choose_package_manager()
{
  OS=$(uname)
  case "${OS}" in

    "Darwin" )
      if type brew >& /dev/null; then
        package_manager="brew"
      elif type port >& /dev/null; then
        package_manager="port"
      fi
      ;;

    "Linux" )
      if [[ ! -f /etc/lsb-release ]]; then
        echo "I don't recognize this Linux distribution. Possibly it's too old to have the expected /etc/lsb-release file."
        exit 1
      fi
      . /etc/lsb-release

      case "${DISTRIB_ID:-}" in
        "Ubuntu" )
           package_manager="apt-get"
        ;;
        "Debian" )
           package_manager="dnf"
        ;;
        *)
          echo "I don't recognize the Linux distribution ${DISTRIB_ID:-}"
          exit 1
        ;;
      esac
      ;;

    *)
      echo "I don't recognize the operating system \${OS}"
      exit 1
      ;;
  esac
}
choose_package_manager

### Install any missing prerequisite packages and executable programs ###

### First install the packages for which package management suffices:
function install_if_missing()
{
  # Exit with error message if no first argument:
  : "${1?'\n ERROR: install_if_missing function requires a package name as the first argument.'}"

  package="${1}"          # Interpret the 1st argument as the package name
  executable="${2:-${1}}" # Interpret the 2nd argument as the prerequisite executable program (Default = 1st argument)

  printf "Checking whether ${executable} is in path...  "

  if type ${executable} >& /dev/null; then

    printf "yes.\n"

  else # install package

    printf "no\n"
    echo "sudo ${package_manager} install ${package}"
    sudo "${package_manager}" install ${package}

  fi
}

# Install subversion to check out the GCC trunk if not already in the PATH:
install_if_missing subversion svn

# Install released versions of the GNU compilers if not already in the PATH:
install_if_missing gfortran
install_if_missing g++

# Install build software:
install_if_missing cmake
install_if_missing make
install_if_missing flex

### Install the prerequisites that must be built from source ###

# Download and build the GCC trunk:
echo "Downloading the GCC trunk."
./install.sh --only-download --package gcc --install-branch trunk --yes-to-all

if [[ ! -z "${absolute_path:-}" ]]; then
  # Patch the GCC trunk and rebuild
  echo "Patching the GCC source using ${absolute_path}."
  pushd prerequisites/downloads/trunk
    patch -p0 < "${absolute_path}"
  popd
  export patched="patched"
fi


export GCC_install_prefix=${packages_root}/gnu/trunk
# Build the patched GCC trunk

echo "Rebuilding the patched GCC source."
./install.sh \
  --package gcc \
  --install-branch trunk \
  --yes-to-all \
  --num-threads 4 \
  --disable-bootstrap \
  --install-prefix "${GCC_install_prefix}"

# Verify that GCC installed in the expected path
if ! type "${GCC_install_prefix}"/bin/gfortran >& /dev/null; then
  echo "gfortran is not installed in the expected location ${GCC_install_prefix}."
  exit 1
fi

# Ensure that the just-installed GCC libraries are used when linking
echo "Setting and exporting PATH and LD_LIBRARY_PATH"

function prepend() {

  : "${1?'prepend(): missing prefix'}"
  : "${2?'prepend(): missing environment variable name'}"

  new_prefix="$1"
  export env_var_name="$2"
  eval env_var_val="\$$2"
  if [[ -z $env_var_val ]]; then
    export $env_var_name="$new_prefix"
  else
    export $env_var_name="${new_prefix}":$env_var_val
  fi
}

if type "${GCC_install_prefix}"/bin/gfortran; then
  prepend "${GCC_install_prefix}/bin/" PATH
  prepend "${GCC_install_prefix}/lib64/" LD_LIBRARY_PATH
else
  echo "GCC is not installed in the expected location ${GCC_install_prefix}"
fi

echo "Building MPICH with the ${patched:-} compilers."
export mpich_install_prefix="${packages_root}/mpich/3.2/gnu/trunk"

./install.sh \
  --package mpich \
  --num-threads 4 \
  --yes-to-all \
  --with-fortran "${GCC_install_prefix}/bin/gfortran" \
  --with-c "${GCC_install_prefix}/bin/gcc" \
  --with-cxx "${GCC_install_prefix}/bin/g++" \
  --install-prefix "${mpich_install_prefix}"

# Verify that MPICH installed where expected
if type "${mpich_install_prefix}"/bin/mpifort; then
  prepend "${mpich_install_prefix}/bin/" PATH
else
  echo "MPICH is not installed in the expected location ${mpich_install_prefix}."
  exit 1
fi

# Build OpenCoarrays with the just-built compilers and the just-built MPICH
echo "Building OpenCoarrays."
export opencoarrays_version=$(./install.sh -V opencoarrays)
export opencoarrays_install_prefix="${packages_root}/opencoarrays/${opencoarrays_version}/gnu/trunk"

./install.sh \
  --num-threads 4 \
  --yes-to-all \
  --with-mpi "$mpich_install_prefix" \
  --install-prefix "$opencoarrays_install_prefix"

# Verify that OpenCoarrays installed where expected
if type "${opencoarrays_install_prefix:-}"/bin/caf >& /dev/null; then
  echo ""
  echo "OpenCoarrays is in $opencoarrays_install_prefix"
  echo "To set up your environment for using OpenCoarrays, execute the following command:"
  echo "source ${opencoarrays_install_prefix}/setup.sh"
else
  echo "OpenCoarrays is not in the expected location: ${opencoarrays_install_prefix:-"(no expected location found)"}."
  exit 1
fi
