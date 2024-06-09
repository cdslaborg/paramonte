#!/bin/bash
####################################################################################################################################
####################################################################################################################################
####                                                                                                                            ####
####    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ####
####                                                                                                                            ####
####    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ####
####                                                                                                                            ####
####    This file is part of the ParaMonte library.                                                                             ####
####                                                                                                                            ####
####    LICENSE                                                                                                                 ####
####                                                                                                                            ####
####       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ####
####                                                                                                                            ####
####################################################################################################################################
####################################################################################################################################

#   This Bash script can be called by `install.sh` or any Bash script that needs the following initialization and definitions.
#   This script takes one input argument `caller_name` that is the name of the caller script to set the user notes. default: caller_name=`""`

#IFS=";" # Bash element separator in lists.

paramonte_dir="$( cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd )"
paramonte_dir="$(dirname "${paramonte_dir}")"; export paramonte_dir
paramonte_src_dir="${paramonte_dir}/src"; export paramonte_src_dir
paramonte_auxil_dir="${paramonte_dir}/auxil"; export paramonte_auxil_dir
paramonte_external_dir="${paramonte_dir}/external"; export paramonte_req_dir
paramonte_example_dir="${paramonte_dir}/example"; export paramonte_example_dir
paramonte_benchmark_dir="${paramonte_dir}/benchmark"; export paramonte_benchmark_dir
paramonte_src_fortran_dir="${paramonte_src_dir}/fortran"; export paramonte_src_fortran_dir
paramonte_src_fortran_main_dir="${paramonte_src_fortran_dir}/main"; export paramonte_src_fortran_main_dir
paramonte_src_fortran_test_dir="${paramonte_src_fortran_dir}/test"; export paramonte_src_fortran_test_dir
paramonte_external_doc_dir="${paramonte_external_dir}/paramonted"; export paramonte_external_doc_dir
paramonte_external_doc_out_dir="${paramonte_external_doc_dir}/paramonte"; export paramonte_external_doc_out_dir
paramonte_req_dir="${paramonte_dir}/prerequisites"; export paramonte_req_dir

paramonte_web_github="https://github.com/cdslaborg/paramonte"; export paramonte_req_dir
paramonte_web_github_issues="${paramonte_web_github}/issues/new/choose"; export paramonte_req_dir
build_name="ParaMonte"; export build_name

caller_label=""
if ! [ "${caller_name}" = "" ]; then
    caller_label=" ${caller_name}"
fi

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set up color coding.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

asciiEscVal=27; ESC="$(printf "\\$(printf "%o" "${asciiEscVal}")")";
ColorReset="${ESC}[m"
ColorBold="${ESC}[1m"
Red="${ESC}[31m"
Green="${ESC}[32m"
Yellow="${ESC}[33m"
Blue="${ESC}[34m"
Magenta="${ESC}[35m"
Cyan="${ESC}[36m"
White="${ESC}[37m"
BoldRed="${ESC}[1;31m"
BoldGreen="${ESC}[1;32m"
BoldYellow="${ESC}[1;33m"
BoldBlue="${ESC}[1;34m"
BoldMagenta="${ESC}[1;35m"
BoldCyan="${ESC}[1;36m"
BoldWhite="${ESC}[1;37m"

pmcolor="${BoldCyan}"
pmattn=" ${pmcolor}-- ${build_name}${caller_label} -${ColorReset}"
pmnote="${pmattn} ${BoldYellow}NOTE:${ColorReset}"
pmwarn="${pmattn} ${BoldMagenta}WARNING:${ColorReset}"
pmfatal="${pmattn} ${BoldRed}FATAL ERROR:${ColorReset}"
warning="${BoldMagenta}WARNING${ColorReset}"

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# auxiliary helper functions.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#getLowerCaseChar()
#{
#    if [[ $1 =~ [A-Z] ]];then
#        n=$(printf "%d" "'$1")
#        n=$((n+32))
#        printf \\$(printf "%o" "$n")
#    else
#        printf "%s" "$1"
#    fi
#}
#
#getUpperCaseChar()
#{
#    if [[ $1 =~ [a-z] ]];then
#        n=$(printf "%d" "'$1")
#        n=$((n-32))
#        printf \\$(printf "%o" "$n")
#    else
#        printf "%s" "$1"
#    fi
#}

getLowerCase() {
    echo "$@" | tr '[:upper:]' '[:lower:]'
    #word="$@"
    #for((i=0;i<${#word};i++)); do
    #    ch="${word:$i:1}"
    #    getLowerCaseChar $ch
    #done
}

getUpperCase() {
    echo "$@" | tr '[:lower:]' '[:upper:]'
#    word="$@"
#    for((i=0;i<${#word};i++)); do
#        ch="${word:$i:1}"
#        getUpperCaseChar "$ch"
#    done
}

isnumeric() {
    isNumericValue=true
    word="$@"
    for((i=0;i<${#word};i++)); do
        ch="${word:$i:1}"
        if ! [[ $ch =~ [0-9] ]];then
            isNumericValue=false
            break
        fi
    done
    echo $isNumericValue
}

reportConflict()
{
    usage
    echo >&2 "${pmwarn} conflicting flag values detected."
    echo >&2 "${pmwarn} $1"
    echo >&2 "${pmwarn} The requested build configuration will ignored."
    echo >&2 "${pmwarn} skipping..."
    #exit 1
}

reportBadValue()
{
    usage
    echo >&2 ""
    echo >&2 "${pmfatal} The requested input value $2 specified with "
    echo >&2 "${pmfatal} the input flag $1 is not supported."
    if ! [ -z ${3+x} ]; then
    echo >&2 "${pmfatal} $3"
    fi
    echo >&2 ""
    echo >&2 "${pmfatal} Gracefully exiting."
    echo >&2 ""
    echo >&2 ""
    exit 1
}


# Compare two version strings [$1: version string 1 (v1), $2: version string 2 (v2)]
# Return values:
#   0: v1 == v2
#   1: v1 > v2
#   2: v1 < v2
function compareVersions() {

    # Trivial v1 == v2 test based on string comparison
    [[ "$1" == "$2" ]] && return 0

    # Local variables
    local regex="^(.*)-r([0-9]*)$" va1=() vr1=0 va2=() vr2=0 len i IFS="."

    # Split version strings into arrays, extract trailing revisions
    if [[ "$1" =~ ${regex} ]]; then
        va1=(${BASH_REMATCH[1]})
        [[ -n "${BASH_REMATCH[2]}" ]] && vr1=${BASH_REMATCH[2]}
    else
        va1=($1)
    fi
    if [[ "$2" =~ ${regex} ]]; then
        va2=(${BASH_REMATCH[1]})
        [[ -n "${BASH_REMATCH[2]}" ]] && vr2=${BASH_REMATCH[2]}
    else
        va2=($2)
    fi

    # Bring va1 and va2 to same length by filling empty fields with zeros
    (( ${#va1[@]} > ${#va2[@]} )) && len=${#va1[@]} || len=${#va2[@]}
    for ((i=0; i < len; ++i)); do
        [[ -z "${va1[i]}" ]] && va1[i]="0"
        [[ -z "${va2[i]}" ]] && va2[i]="0"
    done

    # Append revisions, increment length
    va1+=($vr1)
    va2+=($vr2)
    len=$((len+1))

    # *** DEBUG ***
    #echo "TEST: '${va1[@]} (?) ${va2[@]}'"

    # Compare version elements, check if v1 > v2 or v1 < v2
    for ((i=0; i < len; ++i)); do
        if (( 10#${va1[i]} > 10#${va2[i]} )); then
            return 1
        elif (( 10#${va1[i]} < 10#${va2[i]} )); then
            return 2
        fi
    done

    # All elements are equal, thus v1 == v2
    return 0
}

printCopyFailMsg() {
    echo >&2
    echo >&2 "-- ${pmfatal} - Copy action failed. Please resolve the error. Gracefully exiting..."
    echo >&2
    exit 1
}

verify() {
    set +x
    if [ $1 -eq 0 ]; then
        echo >&2 "${pmnote} ${BoldGreen}The ParaMonte $2 appears to have succeeded.${ColorReset}"
    else
        echo >&2
        echo >&2 "${pmfatal} ParaMonte $2 appears to have failed."
        echo >&2 "${pmfatal} If the source of the error cannot be identified,"
        echo >&2 "${pmfatal} consider a fresh installation of ParaMonte's required compilers by calling"
        echo >&2 "${pmfatal} "
        echo >&2 "${pmfatal}     ./install.sh --fresh all"
        echo >&2 "${pmfatal} "
        echo >&2 "${pmfatal} If the error happens during the installation of ParaMonte prerequisites,"
        echo >&2 "${pmfatal} it is possible that the current existing GNU compiler collection installed"
        echo >&2 "${pmfatal} on your system cannot compile the downloaded version of GNU that is required"
        echo >&2 "${pmfatal} for ParaMonte build. In such case, make sure you have a GNU compiler collection"
        echo >&2 "${pmfatal} version ${versionMinGNU} or newer installed on your system, with an updated PATH environmental"
        echo >&2 "${pmfatal} variable, then reinstall ParaMonte."
        echo >&2 "${pmfatal} "
        echo >&2 "${pmfatal} If the error is solely due to the failures of some ParaMonte tests, then"
        echo >&2 "${pmfatal} you may want to skip the testing of the library by specifying \"--test none\" or"
        echo >&2 "${pmfatal} \"--test none\" when calling the ParaMonte installation script."
        echo >&2 "${pmfatal} "
        echo >&2 "${pmfatal}     ./install.sh --test none"
        echo >&2 "${pmfatal} "
        echo >&2 "${pmfatal} To get more help on the usage of the optional install flags, try:"
        echo >&2 "${pmfatal} "
        echo >&2 "${pmfatal}     ./install.sh --help"
        echo >&2 "${pmfatal} "
        echo >&2 "${pmfatal} If all ParaMonte installation attempts fail, please report this issue at"
        echo >&2 "${pmfatal} "
        echo >&2 "${pmfatal}     ${paramonte_web_github_issues}"
        echo >&2 "${pmfatal} "
        echo >&2 "${pmfatal} or by contacting the ParaMonte authors directly (e.g., shahmoradi@utexas.edu)."
        echo >&2
        echo >&2
        echo >&2 "${pmnote} gracefully exiting."
        echo >&2
        exit 1
    fi
}

verifyArgNotKey() {
    # returns 0 if the argument is not a key (beginning with --).
    if [ "${1:0:2}" = "--" ]; then
        echo >&2 "${pmfatal} Key argument (beginning with --) detected as the value of the argument $2."
        echo >&2 "${pmfatal} The specified value: $1"
        echo >&2 "${pmfatal} Gracefully exiting."
        exit 1
    else
        return 0
    fi
}

verifyArgNotEmpty() {
    # returns 0 if the argument is non-empty.
    if [ "$1" = "" ]; then
        echo >&2 "${pmfatal} Empty value for the argument $2 detected."
        echo >&2 "${pmfatal} Gracefully exiting."
        exit 1
    else
        return 0
    fi
}

usage()
{
    echo "$(cat ${paramonte_dir}/install.sh.md)"
    echo "$(cat ${paramonte_dir}/install.config.md)"
}

ulimit -s unlimited || {
    echo >&2 "${pmwarn} Failed to enlarge the stack allocation."
}

if [[ ! -f "${paramonte_src_fortran_main_dir}/pm_kind.F90" ]]; then
  echo >&2
  echo >&2 "${pmfatal} build failed."
  echo >&2 "${pmfatal} Please run this script inside the top-level ParaMonte library root directory."
  echo >&2 "${pmfatal} This is the directory which contains the LICENSE file in the ParaMonte GitHub repository."
  echo >&2
  exit 1
fi

echo "$(cat ${paramonte_dir}/auxil/.paramonte.banner)"

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Define build directories.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if ! [ -d "${paramonte_dir}" ]; then
    echo >&2
    echo >&2 "${pmfatal} The shell variable `paramonte_dir` is unset or is not a directory."
    echo >&2 "${pmfatal} This variable is defined by the install.sh script inside the top-level ParaMonte library root directory."
    echo >&2 "${pmfatal} This is the directory which contains this file in the GitHub repository of ParaMonte."
    echo >&2 "${pmfatal} The ParaMonte build files have been compromised."
    echo >&2 "${pmfatal} Please build a fresh clone of the project from GitHub: ${paramonte_web_github}"
    echo >&2
    exit 1
fi

#paramonte_src_fortran_dir="${paramonte_dir}/src/fortran"; export paramonte_src_fortran_dir
#paramonte_src_fortran_test_dir="${paramonte_dir}/src/fortran/test"; export paramonte_src_fortran_test_dir
ParaMonte_ROOT_BLD_DIR="${paramonte_dir}/bld"; export ParaMonte_ROOT_BLD_DIR
if ! [ -d "${ParaMonte_ROOT_BLD_DIR}" ]; then mkdir -p "${ParaMonte_ROOT_BLD_DIR}"; fi

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set minimum requirements.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

read -r versionMinOpenCoarrays < "${paramonte_dir}"/auxil/version.min.opencoarrays.txt
read -r versionMinOpenMPI < "${paramonte_dir}"/auxil/version.min.openmpi.txt
read -r versionMinMpich < "${paramonte_dir}"/auxil/version.min.mpich.txt
read -r versionMinIntel < "${paramonte_dir}"/auxil/version.min.intel.txt
read -r versionMinGNU < "${paramonte_dir}"/auxil/version.min.gnu.txt

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Determine the operating system.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

uname_os="$(uname -s)"
case "${uname_os}" in
    Darwin*)    os=darwin;;
    CYGWIN*)    os=cygwin;;
    Linux*)     os=linux;;
    MINGW*)     os=mingw;;
    MSYS*)      os=msys;;
    *)          os=""
esac

#if [ "${os}" = "mingw" ]; then
#    uname_os="$(uname -a)"
#    if [[ "${uname_os}" =~ .*[Mm][Ss][Yy][Ss].* ]]; then
#        os="msys"
#    fi
#fi

#if [[ "${os}" =~ .*"os".* ]]; then
if [ "${os}" = "" ]; then
    echo >&2 "      ${pmwarn} Unrecognized system kernel - ${uname_os}"
    echo >&2 "      ${pmwarn} The supported OS kernel names include: linux, darwin, cygwin, mingw"
    echo >&2 "      ${pmwarn} The ParaMonte library build will proceed without guarantee of success."
    echo >&2 "      ${pmwarn} If you believe the error can be resolved or the install script can be improved,"
    echo >&2 "      ${pmwarn} we'd appreciate your reporting of this issue at:"
    echo >&2 "      ${pmwarn} ${paramonte_web_github_issues}"
    echo >&2
    #exit 1
    os="${uname_os}"
fi
os="${os//$'\r'/}" # Remove all carriage returns.
os="${os//$'\n'/}" # Remove all newlines.
export os

# uname_os_FULL="$(uname -a)"
# if [[ "$uname_os_FULL" =~ .*"Microsoft".* && "$uname_os_FULL" =~ .*"Linux".* ]]; then
#     isWSL=true
# else
#     isWSL=false
# fi

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Determine the Make software on the platform and check compatibility.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

makename=make
if  command -v make >/dev/null 2>&1; then
    makeVersion="$(make --version)"
    makeVersionLower="$(getLowerCase "$makeVersion")"
    if [[ "${makeVersionLower}" =~ .*"apple".* ]]; then
        osmake="darwin"
    elif [[ "${makeVersionLower}" =~ .*"cygwin".* ]]; then
        osmake="cygwin"
    elif [[ "${makeVersionLower}" =~ .*"darwin".* ]]; then
        osmake="darwin"
    elif [[ "${makeVersionLower}" =~ .*"linux".* ]]; then
        osmake="linux"
    elif [[ "${makeVersionLower}" =~ .*"mingw".* ]]; then
        osmake="mingw"
    elif [[ "${makeVersionLower}" =~ .*"msys".* ]]; then
        osmake="msys"
    else
        osmake=""
    fi
    echo >&2 "${pmnote} Make software found at: \"$(command -v make)\""
    echo >&2 "${pmnote} Make software version: \"${makeVersion}\""
    if  ! [ "${os}" = "${osmake}" ]; then
        echo >&2 "${pmwarn} Make software version (${osmake}) does is inconsistent with the runtime shell environment (${os})."
        echo >&2 "${pmwarn} ParaMonte library build will proceed with no guarantee of success."
    fi
elif [ "${os}" = "mingw" ] || [ "${os}" = "msys" ]; then
    echo >&2 "${pmwarn} Failed to detect the Make software."
    echo >&2 "${pmwarn} Searching for MINGW/MSYS version of Make binary..."
    if  command -v mingw32-make.exe >/dev/null 2>&1; then
        echo >&2 "${pmnote} mingw32-make.exe Make software found at: \"$(command -v mingw32-make.exe)\""
        makename=mingw32-make.exe
    elif  command -v mingw32-make >/dev/null 2>&1; then
        echo >&2 "${pmnote} mingw32-make Make software found at: \"$(command -v mingw32-make)\""
        makename=mingw32-make
    else
        echo >&2 "${pmwarn} Failed to detect the mingw32-make software. "
        echo >&2 "${pmwarn} Proceeding with the ParaMonte library build with no guarantee of success..."
    fi
fi
export makename

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Determine the architecture.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#arch=$(uname -p) # -p is non-portable. Use -m.
#if [[ "$arch" =~ .*"64".* ]]; then
#    arch="amd64"
#else
#    arch=$(uname -m)
#    if [[ "$arch" =~ .*"64".* ]]; then arch="amd64"; fi
#fi
arch="$(uname -m || uname -p)"
if [[ "${arch}" =~ .*x86_64.* ]]; then
    arch="amd64"
fi
export arch

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set the ParaMonte version (to be used by cmake. UPDATE: CMake read this by itself as of June 2022.)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

paramonte_lang_list=("c" "cpp" "fortran" "matlab" "python" "r")

#paramonte_version_fortran=$(head -n 1 ../../VERSION.md); export paramonte_version_fortran
read -r paramonte_version           < "${paramonte_dir}/VERSION.md"
read -r paramonte_version_c         < "${paramonte_dir}/src/c/VERSION.md"
read -r paramonte_version_cpp       < "${paramonte_dir}/src/cpp/VERSION.md"
read -r paramonte_version_fortran   < "${paramonte_dir}/src/fortran/VERSION.md"
read -r paramonte_version_matlab    < "${paramonte_dir}/src/matlab/VERSION.md"
read -r paramonte_version_python    < "${paramonte_dir}/src/matlab/VERSION.md"

temp="${IFS}"
IFS='.' read -r -a tmparr <<< "${paramonte_version}";           export paramonte_version_major="${tmparr[0]}"
IFS='.' read -r -a tmparr <<< "${paramonte_version_c}";         export paramonte_version_major_c="${tmparr[0]}"
IFS='.' read -r -a tmparr <<< "${paramonte_version_cpp}";       export paramonte_version_major_cpp="${tmparr[0]}"
IFS='.' read -r -a tmparr <<< "${paramonte_version_fortran}";   export paramonte_version_major_fortran="${tmparr[0]}"
IFS='.' read -r -a tmparr <<< "${paramonte_version_matlab}";    export paramonte_version_major_matlab="${tmparr[0]}"
IFS='.' read -r -a tmparr <<< "${paramonte_version_python}";    export paramonte_version_major_python="${tmparr[0]}"
IFS="${temp}"

echo "paramonte_version_major           = ${paramonte_version_major}"
echo "paramonte_version_major_c         = ${paramonte_version_major_c}"
echo "paramonte_version_major_cpp       = ${paramonte_version_major_cpp}"
echo "paramonte_version_major_fortran   = ${paramonte_version_major_fortran}"
echo "paramonte_version_major_matlab    = ${paramonte_version_major_matlab}"
echo "paramonte_version_major_python    = ${paramonte_version_major_python}"

#echo >&2 "${pmnote} ParaMonte version: ${paramonte_version_fortran}" #this is now displayed in the banner below.
#unset fppParaMonteVersion
#unset FPP_PARAMONTE_VERSION_FLAG
#if ! [ -z ${paramonte_version_fortran+x} ]; then
#    export paramonte_version_fortran
#    # NOTE: uncomment the following line, if you wish to insert the ParaMonte version to the source code via Cmake via the compiler preprocessor.
#    # fppParaMonteVersion="${paramonte_version_fortran}"; export fppParaMonteVersion
#    # write the version source include file.
#    if ! [ -f "${ParaMonteVersionFortran_SRC_INC_PATH}" ] || [ "$(grep -c "${paramonte_version_fortran}" "${ParaMonteVersionFortran_SRC_INC_PATH}")" = "0" ]; then
#        echo "! WARNING - DO NOT CHANGE THE CONTENTS OF THIS FILE MANUALLY." > "${ParaMonteVersionFortran_SRC_INC_PATH}"
#        echo "! WARNING - This file is auto-generated by the ParaMonte build scripts." >> "${ParaMonteVersionFortran_SRC_INC_PATH}"
#        echo "self%version = \"${paramonte_version_fortran}\"" >> "${ParaMonteVersionFortran_SRC_INC_PATH}"
#    fi
#fi
##echo "$(cat ./auxil/.paramonte.banner)"
paramonte_version_c="${paramonte_version_c//$'\r'/}" # Remove all carriage returns.
paramonte_version_c="${paramonte_version_c//$'\n'/}" # Remove all newlines.
export paramonte_version_c
paramonte_version_cpp="${paramonte_version_cpp//$'\r'/}" # Remove all carriage returns.
paramonte_version_cpp="${paramonte_version_cpp//$'\n'/}" # Remove all newlines.
export paramonte_version_cpp
paramonte_version_fortran="${paramonte_version_fortran//$'\r'/}" # Remove all carriage returns.
paramonte_version_fortran="${paramonte_version_fortran//$'\n'/}" # Remove all newlines.
export paramonte_version_fortran

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Report build setup.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo >&2
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "                                        ParaMonte library version ${paramonte_version_fortran} build on ${os}                                "
echo >&2 "::::                                                                                                                            ::::"
echo >&2 "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo >&2
echo >&2
echo >&2 "${pmnote}               working directory: $(pwd)"
echo >&2 "${pmnote}        ParaMonte root directory: ${paramonte_dir}"
echo >&2 "${pmnote}   current system's architecture: ${arch}"
echo >&2 "${pmnote} current system's OS kernel name: ${os}"
echo >&2 "${pmnote}                   paramonte_dir: ${paramonte_dir}"
echo >&2 "${pmnote}  paramonte_src_fortran_test_dir: ${paramonte_src_fortran_test_dir}"
echo >&2 "${pmnote}       paramonte_src_fortran_dir: ${paramonte_src_fortran_dir}"
echo >&2 "${pmnote}          ParaMonte_ROOT_BLD_DIR: ${ParaMonte_ROOT_BLD_DIR}"

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Verify CMake installation and version.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

read -r versionMinCMake < "${paramonte_dir}"/auxil/version.min.cmake.txt

versionCurrentCMake="$(cmake --version)"
versionCurrentCMakeArray=($versionCurrentCMake)
versionCurrentCMake="${versionCurrentCMakeArray[2]}"

if [ "${versionCurrentCMake}" = "" ]; then
    cmakeInstallEnabled=true
else
    cmakeInstallEnabled=false
    compareVersions "${versionCurrentCMake}" "${versionMinCMake}"
    if [ "$?" = "2" ]; then
        cmakeInstallEnabled=true
        echo >&2 "${pmwarn} Failed to detect a ParaMonte-compatible installation of cmake."
    else
        echo >&2 "${pmnote} The current cmake installation is ParaMonte compatible."
    fi
fi
#export cmakeInstallEnabled

echo >&2 "${pmnote} cmake path: $(command -v cmake)"
echo >&2 "${pmnote} cmake version current: ${versionCurrentCMake}"
echo >&2 "${pmnote} cmake version required: ${versionMinCMake}"
echo >&2 "${pmnote} cmakeInstallEnabled: ${cmakeInstallEnabled}"

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Setup build functions.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

getPathFC() {
    # Return the compiler absolute path given its name or path.
    local fc="$1"
    if [ "${fc}" = "" ] || [ "${fc}" = "none" ]; then
        local pathFC=""
    else
        if [ "${os}" = "mingw" ] || [ "${os}" = "msys" ] || [ "${os}" = "cygwin" ]; then
            extensions=(".exe" "")
        else
            extensions=("")
        fi
        for extension in "${extensions[@]}"; do
            #echo >&2 "${pmnote} Checking compiler executable with extension: \"${extension}\""
            if [ -f "${fc}${extension}" ]; then # check if the specified compiler can be found in the environment.
                local pathFC="$(ls --indicator-style=none ${fc}${extension}*)"
                break
            else
                local pathFC="$(command -v "${fc}${extension}")"
                if [ -f "${pathFC}" ]; then
                    break
                else
                    local pathFC="$(ls --indicator-style=none $(command -v "${fc}${extension}")*)"
                    if [ -f "${pathFC}" ]; then
                        break
                    fi
                fi
            fi
        done
    fi
    echo "${pathFC}"
}

getCSID() {
    # Return the compiler suite name (ID) given its path.
    local csvs="csid"
    local pathFC="$1"
    if [ "${pathFC}" = "" ]; then
        local csid="csid"
    else
        if [[ "${pathFC}" =~ .*"ifort".* || "${pathFC}" =~ .*"ifx".* ]]; then
            local csid="intel"
        elif [[ "${pathFC}" =~ .*"gfortran".* ]]; then
            local csid="gnu"
        fi
    fi
    echo "${csid}"
}

getCSVS() {
    # Return the compiler suite major version given its path and csid.
    local csvs="csvs"
    local pathFC="$1"
    local csid="$2"
    if [ -f "${pathFC}" ]; then
        if [ "${csid}" = "intel" ] || [ "${csid}" = "gnu" ]; then
            tempsrc="$(mktemp).F90"
            tempdir="$(mktemp -d)"
            if ! [ -d "${tempdir}" ]; then
                tempsrc="getCompilerVersion.F90"
                tempdir="."
            fi
            versionExtractionFailed=true
            #echo >&2 "${pmnote} Changing directory to: ${tempdir}"
            cd "${tempdir}" && cp "${paramonte_dir}/auxil/getCompilerVersion.F90" "${tempsrc}"
            if "${pathFC}" "${tempsrc}" -o "${tempsrc}".exe; then
                chmod +x "${tempsrc}".exe && "${tempsrc}".exe "${tempsrc}".tmp && {
                    local csvs=$(head -n 1 "${tempsrc}".tmp)
                    versionExtractionFailed=false
                    rm *.tmp *.exe *.o >/dev/null 2>&1 # || continue
                }
            fi
            if [ "${versionExtractionFailed}" = "true" ]; then
                echo >&2 "${pmwarn} Failed to fetch the compiler suite name and version."
                echo >&2 "${pmwarn} Proceeding with no guarantee of build success..."
            fi
            #echo >&2 "${pmnote} Changing directory to: ${paramonte_dir}"
            cd "${paramonte_dir}"
        fi
    fi
    echo "${csvs}"
}
