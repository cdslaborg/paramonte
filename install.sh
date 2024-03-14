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

####    See the file install.sh.usage in the same folder for usage guidelines of this Batch script.
####    
####    NOTE: Do not change the contents of this file unless you know what the consequences are.
####    This is the Bash script file that builds objects, shared libraries,
####    as well as the test and example binaries of the ParaMonte library in
####    POSIX-like (non-Windows) environments (such as Unix Shell, Git Bash).
####    Upon invocation of this file from a Bash command-line interface,
####    this script will parse the user-provided flags and their values
####    to build the ParaMonte library.
####    
####    To redirect the output to an external file (e.g., install.sh.out), try:
####    
####        install.sh >install.sh.out 2>&1
####    
####    to redirect output to the external file install.sh.out and run the installation in background, try:
####    
####        install.sh >install.sh.out 2>&1 &; jobs; disown
####        jobs; disown

# The following lines must appear in the specified order before anything else in the script.
#paramonte_dir="${paramonte_dir:-${PWD%/}}"
FILE_DIR="$( cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd )"
caller_name="$(basename "${BASH_SOURCE[0]}")"
source ./auxil/install.init.sh
#workingDir="$(pwd)"

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Reset the script arguments.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if [ 0 -lt 1 ]; then # just to allow toggling in notepad++.

    unset bdir
    export FOR_COARRAY_NUM_IMAGES=3
    ddir="${paramonte_dir}/bin"
    flag_ddir="-Dddir=${ddir}"

    unset list_build
    unset list_checking
    unset list_fc
    unset list_lang
    unset list_lib
    unset list_mem
    unset list_par

    unset flag_bench
    unset flag_benchpp
    unset flag_blas
    unset flag_codecov
    unset flag_cfi
    unset flag_deps
    unset flag_exam
    unset flag_exampp
    unset flag_fpp
    unset flag_fresh
    unset flag_j
    unset flag_lapack
    unset flag_matlabdir
    unset flag_me
    unset flag_mod
    unset flag_nproc
    unset flag_perfprof
    unset flag_pdt
    unset flag_purity
    unset flag_test

    unset flag_ski
    unset flag_iki
    unset flag_lki
    unset flag_cki
    unset flag_rki

    while [ "$1" != "" ]; do
        case "$1" in

            #### list args

            --build )       shift
                            verifyArgNotKey "$1" --build
                            verifyArgNotEmpty "$1" --build
                            list_build="$1"
                            ;;

            --checking )    shift
                            verifyArgNotKey "$1" --checking
                            verifyArgNotEmpty "$1" --checking
                            list_checking="$1"
                            ;;

            --fc )          shift
                            verifyArgNotKey "$1" --fc
                            verifyArgNotEmpty "$1" --fc
                            list_fc="$1"
                            ;;

            --lang )        shift
                            verifyArgNotKey "$1" --lang
                            verifyArgNotEmpty "$1" --lang
                            list_lang="$1"
                            ;;

            --lib )         shift
                            verifyArgNotKey "$1" --lib
                            verifyArgNotEmpty "$1" --lib
                            list_lib="$1"
                            ;;

            --mem )         shift
                            verifyArgNotKey "$1" --mem
                            verifyArgNotEmpty "$1" --mem
                            list_mem="$1"
                            ;;

            --par )         shift
                            verifyArgNotKey "$1" --par
                            verifyArgNotEmpty "$1" --par
                            list_par="$1"
                            ;;

            #### flag args

            --bench )       shift
                            verifyArgNotKey "$1" --bench
                            verifyArgNotEmpty "$1" --bench
                            flag_bench="-Dbench=$1"
                            ;;

            --benchpp )     shift
                            verifyArgNotKey "$1" --benchpp
                            verifyArgNotEmpty "$1" --benchpp
                            flag_benchpp="-Dbenchpp=$1"
                            ;;

            --blas )        shift
                            verifyArgNotKey "$1" --blas
                            verifyArgNotEmpty "$1" --blas
                            flag_blas="-Dblas=$1"
                            ;;

            --codecov )     shift
                            verifyArgNotKey "$1" --codecov
                            verifyArgNotEmpty "$1" --codecov
                            flag_codecov="-Dcodecov=$1"
                            ;;

            --deps )        shift
                            verifyArgNotKey "$1" --deps
                            verifyArgNotEmpty "$1" --deps
                            flag_deps="-Ddeps=$1"
                            ;;

            --exam )        shift
                            verifyArgNotKey "$1" --exam
                            verifyArgNotEmpty "$1" --exam
                            flag_exam="-Dexam=$1"
                            ;;

            --exampp )      shift
                            verifyArgNotKey "$1" --exampp
                            verifyArgNotEmpty "$1" --exampp
                            flag_exampp="-Dexampp=$1"
                            ;;

            --cfi )         shift
                            verifyArgNotKey "$1" --cfi
                            verifyArgNotEmpty "$1" --cfi
                            flag_cfi="-Dcfi=$1"
                            ;;

            --fpp )         shift
                            verifyArgNotKey "$1" --fpp
                            verifyArgNotEmpty "$1" --fpp
                            flag_fpp="-Dfpp=$1"
                            ;;

            --fresh )       shift
                            verifyArgNotKey "$1" --fresh
                            verifyArgNotEmpty "$1" --fresh
                            flag_fresh="-Dfresh=$1"
                            ;;

            --lapack )      shift
                            verifyArgNotKey "$1" --lapack
                            verifyArgNotEmpty "$1" --lapack
                            flag_lapack="-Dlapack=$1"
                            ;;

            --matlabdir )   shift
                            verifyArgNotKey "$1" --matlabdir
                            verifyArgNotEmpty "$1" --matlabdir
                            flag_matlabdir="-Dmatlabdir=\"$1\""
                            ;;

            --me )          shift
                            verifyArgNotKey "$1" --me
                            verifyArgNotEmpty "$1" --me
                            flag_me="-Dme=\"$1\""
                            ;;

            --mod )         shift
                            verifyArgNotKey "$1" --mod
                            verifyArgNotEmpty "$1" --mod
                            flag_mod="-Dmod=$1"
                            ;;

            --nproc )       shift
                            verifyArgNotKey "$1" --nproc
                            verifyArgNotEmpty "$1" --nproc
                            isNumericValue="$(isnumeric ${nproc})"
                            if ! [ "${isNumericValue}" = "true" ]; then
                                reportBadValue "--nproc" $nproc "The spoecified number of processors must be a positive integer."
                            fi
                            FOR_COARRAY_NUM_IMAGES="$1"
                            flag_nproc="-Dnproc=$1"
                            ;;

            --perfprof )    shift
                            verifyArgNotKey "$1" --perfprof
                            verifyArgNotEmpty "$1" --perfprof
                            flag_perfprof="-Dperfprof=$1"
                            ;;

            --pdt )         shift
                            verifyArgNotKey "$1" --pdt
                            verifyArgNotEmpty "$1" --pdt
                            flag_pdt="-Dpdt=$1"
                            ;;

            --purity )      shift
                            verifyArgNotKey "$1" --purity
                            verifyArgNotEmpty "$1" --purity
                            flag_purity="-Dpurity=$1"
                            ;;

            --test )        shift
                            verifyArgNotKey "$1" --test
                            verifyArgNotEmpty "$1" --test
                            flag_test="-Dtest=$1"
                            ;;

            #### flag args: type kind

            --ski )         shift
                            verifyArgNotKey "$1" --ski
                            verifyArgNotEmpty "$1" --ski
                            flag_ski="-Dski=$1"
                            ;;

            --iki )         shift
                            verifyArgNotKey "$1" --iki
                            verifyArgNotEmpty "$1" --iki
                            flag_iki="-Diki=$1"
                            ;;

            --lki )         shift
                            verifyArgNotKey "$1" --lki
                            verifyArgNotEmpty "$1" --lki
                            flag_lki="-Dlki=$1"
                            ;;

            --cki )         shift
                            verifyArgNotKey "$1" --cki
                            verifyArgNotEmpty "$1" --cki
                            flag_cki="-Dcki=$1"
                            ;;

            --rki )         shift
                            verifyArgNotKey "$1" --rki
                            verifyArgNotEmpty "$1" --rki
                            flag_rki="-Drki=$1"
                            ;;

            #### other args

            --bdir )        shift
                            verifyArgNotKey "$1" --bdir
                            verifyArgNotEmpty "$1" --bdir
                            bdir="$1"
                            ;;

            --ddir )        shift
                            verifyArgNotKey "$1" --ddir
                            verifyArgNotEmpty "$1" --ddir
                            ddir="$1"
                            flag_ddir="-Dddir=${ddir}"
                            ;;

            --help )        usage
                            echo >&2 ""
                            echo >&2 ""
                            exit
                            ;;

            -j )            shift
                            verifyArgNotKey "$1" -j
                            verifyArgNotEmpty "$1" -j
                            flag_j="-j $1"
                            ;;

            * )             usage
                            echo >&2 ""
                            echo >&2 "-- ParaMonte - FATAL: The specified flag $1 does not exist."
                            echo >&2 "-- ParaMonte - FATAL: The specified flag $1 does not exist."
                            echo >&2 ""
                            echo >&2 "-- ParaMonte - gracefully exiting."
                            echo >&2 ""
                            echo >&2 ""
                            exit 1
        esac
        shift
    done

fi

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set the default values for the input command line arguments.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if [ 0 -lt 1 ]; then # just to allow toggling in notepad++.

    #### list args

    if  [ -z ${list_build+x} ]; then
        list_build="release"
    else
        list_build="$(getLowerCase "$list_build")"
    fi

    if  [ -z ${list_checking+x} ]; then
        list_checking="nocheck"
    else
        list_checking="$(getLowerCase "$list_checking")"
    fi

    if [ -z ${list_fc+x} ]; then
        #   On Windows OS, particularly with MinGW, CMake fails if the specified compiler name or path does not have the
        #   the file extension ".exe". Given that this extension is unlikely to change in the future, and that it is not used
        #   on Unix systems, try suffixing the extension to the compiler file path. If it fails, use the default specified path.
        compilers=("ifort" "ifx" "gfortran-13" "gfortran-12" "gfortran-11" "gfortran-10" "gfortran")
        if [ "${os}" = "mingw" ] || [ "${os}" = "msys" ] || [ "${os}" = "cygwin" ]; then
            extensions=(".exe" "")
        else
            extensions=("")
        fi
        for compiler in "${compilers[@]}"; do
            for extension in "${extensions[@]}"; do
                echo >&2 "${pmnote} Checking existence of compiler \"${compiler}\" executable with extension \"${extension}\""
                # check if the specified compiler can be found in the environment.
                if command -v "${compiler}${extension}" >/dev/null 2>&1; then
                    list_fc="$(command -v "${compiler}${extension}")"
                    break 2
                fi
            done
        done
        if  [ "${list_fc}" = "" ]; then
            echo >&2 "${pmwarn} Failed to detect any compatible Fortran compiler in the environment."
            echo >&2 "${pmwarn} You can manually specify the Fortran compiler or its path via the install script flag \"--fc\"."
            echo >&2 "${pmwarn} The build will proceed with no guarantee of success."
            list_fc="default"
        else
            echo >&2 "${pmnote} The identified Fortran compiler path is: fc=\"${list_fc}\""
        fi
    fi

    if  [ -z ${list_lang+x} ]; then
        list_lang="fortran"
    else
        list_lang="$(getLowerCase "$list_lang")"
        list_lang="${list_lang/c++/cpp}"
    fi

    if  [ -z ${list_lib+x} ]; then
        list_lib="shared"
    else
        # Replace `dynamic` with `shared`.
        list_lib=${list_lib/dynamic/shared}
        list_lib="$(getLowerCase "$list_lib")"
    fi

    if  [ -z ${list_mem+x} ]; then
        list_mem="heap"
    else
        list_mem="$(getLowerCase "$list_mem")"
    fi

    if  [ -z ${list_par+x} ]; then
        list_par="serial"
    else
        # Replace `none` with `serial`.
        list_par=${list_par/none/serial}
        list_par="$(getLowerCase "$list_par")"
    fi

    if  [ -z ${flag_j+x} ]; then
        flag_j="-j"
    fi

fi

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set CMake default flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if [ 0 -lt 1 ]; then # just to allow toggling in notepad++.

    if  [[ "${flag_fresh}" =~ .*"prereq".* || "${flag_fresh}" =~ .*"all".* ]]; then
        if  [ -d "${paramonte_req_dir}" ]; then
            echo >&2 "${pmnote} Removing the old prerequisites of the ParaMonte library build at: paramonte_req_dir=\"${paramonte_req_dir}\""
            rm -rf "${paramonte_req_dir}"
        fi
    fi

    # Set the CMake build generator.

    flag_g="-G"
    unset cmakeBuildGenerator
    if [ "${os}" = "mingw" ]; then
        cmakeBuildGenerator="MinGW Makefiles"
    elif [ "${os}" = "msys" ]; then
        cmakeBuildGenerator="MSYS Makefiles"
    else
        cmakeBuildGenerator="Unix Makefiles"
    fi

fi

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Configure and build the ParaMonte library.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

for fc in ${list_fc//;/$'\n'}; do

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Set up the CMake fc flag.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    unset flag_fc
    if [ -f "${fc}" ]; then
        fcpath="${fc}"
        flag_fc="-Dfc=${fcpath}"
    elif ! [[ "${fc}" =~ [Dd][Ee][Ff][Aa][Uu][Ll][Tt] ]]; then
        fcpath="$(getPathFC "${fc}")"
        if [ -f "${fcpath}" ]; then
            flag_fc="-Dfc=${fcpath}"
        else
            echo >&2 "${pmfatal} Failed to detect the full path for the specified compiler: fcpath=\"${fcpath}\""
            reportBadValue "--fc" "${fc}"
        fi
        echo >&2 "${pmnote} Fortran compiler path: fcpath=\"${fcpath}\""
    fi

    # Get the compiler ID and version to be used in the build path.

    csid=$(getCSID "${fcpath}")
    csvs=$(getCSVS "${fcpath}" "${csid}")
    echo >&2 "${pmnote} compiler suite: ${csid}"
    echo >&2 "${pmnote} compiler version: ${csvs}"

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Build the ParaMonte library with the specified compiler.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    for lang in $list_lang; do

        flag_lang="-Dlang=${lang}"

        #lang_is_dynamic="false"
        #if [ "${lang}" = "julia" ] || [ "${lang}" = "matlab" ] || [ "${lang}" = "mathematica" ] || [ "${lang}" = "python" ] || [ "${lang}" = "r" ]; then
        #    lang_is_dynamic="true"
        #fi

        ## Set source preprocessing output on by default for Fortran.
        #if [ "${lang}" = "fortran" ] && [ "${flag_fpp}" = "" ]; then
        #    flag_fpp="-Dfpp=generic"
        #fi

        for build in $list_build; do

            flag_build="-Dbuild=${build}"

            for lib in $list_lib; do

                flag_lib="-Dlib=${lib}"

                for mem in $list_mem; do

                    flag_mem="-Dmem=${mem}"

                    for par in $list_par; do

                        flag_par="-Dpar=${par}"

                        for checking in $list_checking; do

                            flag_checking="-Dchecking=${checking}"

                            #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                            # Set the ParaMonte CMake build directory.
                            #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                            if [ 0 -lt 1 ]; then # just to allow toggling in notepad++.

                                # First, determine the MPI library name to be used in the ParaMonte library name and build directory.

                                if [ "${par}" = "mpi" ]; then
                                    if [ -z ${me+x} ]; then
                                        me=mpiexec
                                    fi
                                    if ! [[ -f "${me}" ]]; then
                                        me="$(command -v ${me})"
                                        if ! [[ -f "${me}" ]]; then
                                            echo >&2 "${pmwarn} The specified mpiexec path does not appear to be a valid path to an mpiexec executable: me=\"${me}\""
                                        fi
                                    fi
                                    parname=mpi
                                    mpiVersionInfo="$(${me} --version)"
                                    if [[ "${mpiVersionInfo}" =~ .*"Intel".* ]]; then
                                        parname="impi"
                                    elif [[ "${mpiVersionInfo}" =~ .*[oO][pP][eE][nN][rR][tT][eE].* ]] || [[ "${mpiVersionInfo}" =~ .*[oO][pP][eE][nN]-?[mM][pP][iI].* ]]; then
                                        parname="openmpi"
                                    elif [[ "${mpiVersionInfo}" =~ .*[mM][pP][iI][cC][hH].* ]]; then
                                        parname="mpich"
                                    else # look for mpichversion
                                        mpichversion_path="$(dirname ${me})"/mpichversion
                                        if [ -f "${mpichversion_path}" ] && [[ "$(mpichversion_path)" =~ .*[mM][pP][iI][cC][hH].* ]]; then
                                            parname="mpich"
                                        fi
                                    fi
                                    if [ "${parname}" = "mpi" ]; then
                                        echo >&2 "${pmwarn} The MPI library vendor could not be identified."
                                        echo >&2 "${pmwarn} The MPI library behavior does not match the Intel, MPICH, or OpenMPI libraries."
                                        echo >&2 "${pmwarn} The ParaMonte library name will be suffixed with the generic \"${parname}\" label."
                                    fi
                                elif [ "${par}" = "omp" ]; then
                                    parname="openmp"
                                else
                                    parname="${par}"
                                fi

                                if [ -z ${bdir+x} ]; then
                                    paramonte_bld_dir="${paramonte_dir}/bld/${os}/${arch}/${csid}/${csvs}/${build}/${lib}/${mem}/${parname}/${lang}/${checking}"
                                    if [[ "${flag_perfprof}" =~ .*"all".* ]]; then
                                        paramonte_bld_dir="${paramonte_bld_dir}/perfprof"
                                    fi
                                    if [[ "${flag_codecov}" =~ .*"all".* ]]; then
                                        paramonte_bld_dir="${paramonte_bld_dir}/codecov"
                                    fi
                                    echo >&2 "${pmnote} The ParaMonte library build directory: paramonte_bld_dir=\"${paramonte_bld_dir}\""
                                else
                                    echo >&2 "${pmnote} User-specified library build directory detected: bdir=\"${bdir}\""
                                    paramonte_bld_dir="${bdir}"
                                fi

                                # Make the build directory if needed.

                                if ! [ -d "${paramonte_bld_dir}" ]; then
                                    echo >&2 "${pmnote} Generating the ParaMonte build directory..."
                                    mkdir -p "${paramonte_bld_dir}/"
                                fi

                            fi

                            #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                            #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                            # Configure and build the library via CMake.
                            #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                            #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                            if [ 0 -lt 1 ]; then # just to allow toggling in notepad++.

                                echo >&2 "${pmnote} All generated build files will be stored at: ${paramonte_bld_dir}"
                                echo >&2 "${pmnote} Changing directory to: ${paramonte_bld_dir}"
                                echo >&2 ""
                                echo >&2 "########################################################################################%%"
                                echo >&2 ""
                                echo >&2 "${pmnote} Invoking CMake as:"
                                echo >&2 ""

                                set -x
                                (cd "${paramonte_bld_dir}" && \
                                cmake \
                                "${paramonte_dir}" \
                                ${flag_g} "${cmakeBuildGenerator}" \
                                "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON" \
                                "${flag_ddir}" \
                                "${flag_build}" \
                                "${flag_checking}" \
                                "${flag_lang}" \
                                "${flag_lib}" \
                                "${flag_mem}" \
                                "${flag_par}" \
                                ${flag_fc} \
                                ${flag_bench} \
                                ${flag_benchpp} \
                                ${flag_blas} \
                                ${flag_codecov} \
                                ${flag_cfi} \
                                ${flag_deps} \
                                ${flag_exam} \
                                ${flag_exampp} \
                                ${flag_fpp} \
                                ${flag_fresh} \
                                ${flag_lapack} \
                                ${flag_matlabdir} \
                                ${flag_me} \
                                ${flag_mod} \
                                ${flag_nproc} \
                                ${flag_perfprof} \
                                ${flag_pdt} \
                                ${flag_purity} \
                                ${flag_test} \
                                ${flag_ski} \
                                ${flag_iki} \
                                ${flag_lki} \
                                ${flag_cki} \
                                ${flag_rki} \
                                )
                                verify $? "configuration with cmake"

                                echo >&2 ""
                                echo >&2 "########################################################################################%%"
                                echo >&2 ""

                               #(cd "${paramonte_bld_dir}" && $makename ${flag_j})
                                (cd "${paramonte_bld_dir}" && cmake --build "${paramonte_bld_dir}" ${flag_j})
                                verify $? "build with make"

                               #(cd "${paramonte_bld_dir}" && $makename install)
                                (cd "${paramonte_bld_dir}" && cmake --build "${paramonte_bld_dir}" --target install ${flag_j})
                                verify $? "installation"

                               #(cd "${paramonte_bld_dir}" && $makename deploy)
                                (cd "${paramonte_bld_dir}" && cmake --build "${paramonte_bld_dir}" --target deploy ${flag_j})
                                verify $? "deployment"

                               #(cd "${paramonte_bld_dir}" && $makename test && echo)
                                (cd "${paramonte_bld_dir}" && cmake --build "${paramonte_bld_dir}" --target test)
                                verify $? "testing"

                               #(cd "${paramonte_bld_dir}" && $makename example)
                                (cd "${paramonte_bld_dir}" && cmake --build "${paramonte_bld_dir}" --target example)
                                verify $? "examples build and run"

                               #(cd "${paramonte_bld_dir}" && $makename benchmark)
                                (cd "${paramonte_bld_dir}" && cmake --build "${paramonte_bld_dir}" --target benchmark)
                                verify $? "benchmarks build and run"

                            fi

                            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                            #### CODECOVE: Generate *.gcda *.gcno codecov files for the Fortran test source files
                            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                            if [[ "${flag_codecov}" =~ .*"all".* ]]; then

                                # \todo
                                # \warning
                                # This version extraction relies on "version " appearing before the lcov version number.
                                # This must be adjusted to any potential future changes in the LCOV version string.

                                if  command -v lcov >/dev/null 2>&1; then
                                    lcovVersion=$(lcov -v | grep -Po '(?<=version )[^;]+')
                                else
                                    unset lcovVersion
                                fi

                                unset htmlSubDir
                                if [ "${par}" = "mpi" ]; then
                                    parallelismText="MPI Parallel"
                                    htmlSubDir="mpi"
                                elif [[ "${par}" =~ .*"caf".* ]]; then
                                    parallelismText="Coarray Parallel"
                                    htmlSubDir="caf"
                                else
                                    parallelismText="Serial"
                                    htmlSubDir="serial"
                                fi
                                htmlDir="${paramonte_external_codecov_fortran_dir}/${htmlSubDir}"
                                #htmlDir="${paramonte_external_codecov_fortran_dir}/${paramonte_version_fortran}/${PMLIB_BASE_NAME}"
                                htmlTitleCodeCov="ParaMonte ${paramonte_version_fortran} :: ${parallelismText} Fortran - Code Coverage Report"

                                if [[ ${csid} == [gG][nN][uU] ]]; then

                                    if command -v gcov >/dev/null 2>&1; then

                                        gcovPath="$(command -v gcov)"
                                        echo >&2 "${pmnote} GNU gcov detected at: ${gcovPath}"
                                        echo >&2 "${pmnote} Invoking gcov to generate coverage report..."

                                        gcovFortranDataDir=$(find "${paramonte_bld_dir}" -name pm_blas*.o)
                                        gcovFortranDataDir=$(dirname "${gcovFortranDataDir}")

                                        if [ -d "${gcovFortranDataDir}" ]; then

                                            gcovFortranDir="${paramonte_bld_dir}"/gcov
                                            if [ -d "${gcovFortranDir}" ]; then
                                                echo >&2 "${pmnote} Removing the old existing gcov files directory: ${gcovFortranDir}"
                                                rm -rf "${gcovFortranDir}"
                                            fi
                                            echo >&2 "${pmnote} Generating the gcov files directory: ${gcovFortranDir}"
                                            mkdir -p "${gcovFortranDir}"
                                            cd "${gcovFortranDir}"

                                            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                                            #### GCOV Fortran source: Generate *.gcda *.gcno codecov files for the Fortran source files
                                            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                                            for srcFileName in "${paramonte_src_fortran_dir}"/*.F90; do
                                                if ! [[ "${srcFileName}" =~ .*"pm_array.F90".* ]]; then

                                                    srcFileNameBase=$(basename -- "${srcFileName}")
                                                    #if ! [[ "${srcFileName}" =~ .*".inc.F90".* ]]; then
                                                    #objFilePath=$(find "${gcovFortranDataDir}" -name ${srcFileNameBase}.o)
                                                    # The following assumes that cmake names the object files with full source file name (including file extension).
                                                    objFilePath="${gcovFortranDataDir}/${srcFileNameBase}.o"
                                                    echo >&2 Now running: gcov "${paramonte_src_fortran_dir}/${srcFileName}" -o "${objFilePath}"

                                                    gcov "${paramonte_src_fortran_dir}/${srcFileName}" -o "${objFilePath}" || {
                                                        echo >&2 "${pmfatal} Fatal Error: Code Coverage analysis via GNU gcov tool failed."
                                                        exit 1
                                                    }
                                                    #fi

                                                fi
                                            done

                                            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                                            #### GCOV TEST SOURCE: Generate *.gcda *.gcno codecov files for the Fortran test source files
                                            #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                                            #### first attempt to infer the cmake object files' directory

                                            gcovFortranTestDataDir=$(find "${paramonte_bld_dir}" -name main.F90.o)
                                            gcovFortranTestDataDir=$(dirname "${gcovFortranTestDataDir}")

                                            if [ -d "${gcovFortranTestDataDir}" ]; then
                                                gcovFortranTestDir="${paramonte_bld_dir}"/gcov
                                                #rm -rf "${gcovFortranTestDir}" This would also the new Fortran coverage files.
                                                if ! [ -d "${gcovFortranTestDir}" ]; then
                                                    mkdir -p "${gcovFortranTestDir}"
                                                fi
                                                cd "${gcovFortranTestDir}"
                                                for srcFileName in "${paramonte_src_fortran_test_dir}"/*.F90; do
                                                    if ! [[ "${srcFileName}" =~ .*"main.F90".* ]]; then

                                                        srcFileNameBase=$(basename -- "${srcFileName}")
                                                        objFilePath="${gcovFortranTestDataDir}/${srcFileNameBase}.o"
                                                        echo >&2 gcov "${paramonte_src_fortran_test_dir}/${srcFileName}" -o "${objFilePath}"

                                                        gcov "${paramonte_src_fortran_test_dir}/${srcFileName}" -o "${objFilePath}" \
                                                        || {
                                                            echo >&2 "${pmwarn} The ParaMonte Code Coverage analysis of the test source files via GNU gcov tool failed."
                                                            echo >&2 "${pmwarn} The ParaMonte test source file: ${paramonte_src_fortran_test_dir}/${srcFileName}"
                                                            echo >&2 "${pmwarn} Skipping..."
                                                        }
                                                        #fi
                                                    fi
                                                done
                                            else
                                                echo >&2 "${pmwarn} The directory for the *.gcda *.gcno codecov data files does not exist."
                                                echo >&2 "${pmwarn} The expected directory path: ${gcovFortranTestDataDir}"
                                                echo >&2 "${pmwarn} Skipping code coverage report generation for the test files..."
                                            fi

                                            #::::::::::::::::::::::::::::::::::::::::
                                            # LCOV: generate lcov summary report file
                                            #::::::::::::::::::::::::::::::::::::::::

                                            if command -v lcov >/dev/null 2>&1; then

                                                #### generate the Fortran code coverage report file

                                                lcovFortranDir="${paramonte_bld_dir}"/lcov
                                                rm -rf "${lcovFortranDir}"
                                                mkdir -p "${lcovFortranDir}"
                                                cd "${lcovFortranDir}"

                                                lcovOutputFortranFilePath="${lcovFortranDir}/paramonte.fortran.coverage.info"

                                                unset branchCoverageFlag
                                                # Add the following flag to lcov to enable branch coverage:
                                                # branchCoverageFlag="--rc lcov_branch_coverage=1"
                                                # "${branchCoverageFlag}"

                                                lcov --capture \
                                                --directory "${gcovFortranDataDir}" \
                                                --output-file "${lcovOutputFortranFilePath}" \
                                                || {
                                                    echo >&2 "${pmwarn} Code Coverage report generation via LCOV tool failed."
                                                    #exit 1
                                                }

                                                #### generate the fortran code coverage report file

                                                lcovOutputCombinedFilePath="${lcovFortranDir}/paramonte.combined.coverage.info"

                                                unset lcovOutputTestFilePath
                                                if ls "${gcovFortranTestDir}"/*.gcov 1> /dev/null 2>&1; then
                                                    echo >&2 "${pmnote} generating the code coverage report file for the ParaMonte test files..."
                                                    lcovOutputTestFilePath="${lcovFortranDir}/paramonte.test.coverage.info"
                                                    lcov --capture \
                                                    --directory "${gcovFortranTestDataDir}" \
                                                    --output-file "${lcovOutputTestFilePath}" \
                                                    && {
                                                        echo >&2 "${pmnote} Combining all LCOV code coverage report files as a single final report file..."
                                                        lcov --add-tracefile "${lcovOutputFortranFilePath}" -a "${lcovOutputTestFilePath}" -o "${lcovOutputCombinedFilePath}"
                                                    } || {
                                                        echo >&2 "${pmwarn} Code Coverage report generation for the ParaMonte test source files via lcov tool failed."
                                                        echo >&2 "${pmwarn} Skipping..."
                                                    }
                                                else
                                                    echo >&2 "${pmwarn} Failed to detect the *.gcda *.gcno codecov data files for the ParaMonte test source files."
                                                    echo >&2 "${pmwarn} The expected directory path for the files: ${gcovFortranTestDir}"
                                                    echo >&2 "${pmwarn} The coverage report for the ParaMonte test source file will not be included."
                                                    echo >&2 "${pmwarn} Skipping the lcov code coverage report generation for the test files..."
                                                fi

                                                #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                                                # HTML: convert the lcov summary file to the final html report files
                                                #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                                                if command -v genhtml >/dev/null 2>&1; then

                                                    if [ -d "${htmlDir}" ]; then
                                                        rm -rf "${htmlDir}"
                                                    else
                                                        mkdir -p "${htmlDir}"
                                                    fi

                                                    if ! [ -f "${lcovOutputCombinedFilePath}" ]; then
                                                        cp "${lcovOutputFortranFilePath}" "${lcovOutputCombinedFilePath}" || {
                                                            echo >&2 "${pmwarn} Copy action failed:"
                                                            echo >&2 "${pmwarn} from: ${lcovOutputFortranFilePath}"
                                                            echo >&2 "${pmwarn}   to: ${lcovOutputCombinedFilePath}"
                                                        }
                                                    fi

                                                    genhtml \
                                                    "${lcovOutputCombinedFilePath}" \
                                                    --output-directory "${htmlDir}" \
                                                    --legend \
                                                    --title "${htmlTitleCodeCov}" \
                                                    && {

                                                        echo >&2 "${pmnote} The code coverage build files are stored at: ${paramonte_bld_dir}"
                                                        echo >&2 "${pmnote} The code coverage report files are stored at: ${htmlDir}"

                                                        # postprocess the html files

                                                        pmlinkopen='<a href="https:\/\/www.cdslab.org\/paramonte\/" target="_blank">'
                                                        pmlinklogo='<img alt="The ParaMonte Documentation Website" src="https:\/\/cdslaborg.github.io\/paramonted/fortran\/html\/logo.png"\/>'
                                                        pmlinkclose='<\/a>'

                                                        original='<tr><td class="title">LCOV - code coverage report<\/td><\/tr>'
                                                        modified='<tr><td class="title">'
                                                        modified+="${pmlinkopen}"
                                                        modified+="${pmlinklogo}"
                                                        modified+="${pmlinkclose}"
                                                        modified+='<\/td><\/tr>'

                                                        footer='<tr><td class="versionInfo">'
                                                        footer+='<a href="https:\/\/www.cdslab.org\/paramonte"><b>ParaMonte: Plain Powerful Parallel Monte Carlo Library<\/b><\/a>&nbsp;<br>'
                                                        footer+='<a href="https:\/\/www.cdslab.org" target="_blank"><b>The Computational Data Science Lab<\/b><\/a><br>'
                                                        footer+="&copy; Copyright 2012 - $(date +%Y)"
                                                        footer+='<\/td><\/tr>'

                                                        shopt -s globstar
                                                        for htmlFilePath in "${htmlDir}"/**/*.html; do # Whitespace-safe and recursive
                                                            sed -i "s/${original}/${modified}/g" "${htmlFilePath}"
                                                            sed -i "/<tr><td class=\"versionInfo\">/c\\${footer}" "${htmlFilePath}"
                                                        done

                                                        #scfile="${paramonte_dir}/auxil/sc.html"
                                                        #if [ -f "${scfile}" ]; then
                                                        #    echo >&2
                                                        #    echo >&2 "${pmnote} processing the sc file contents..."
                                                        #    echo >&2
                                                        #    sccontents=`cat "${paramonte_dir}/auxil/sc.html"`
                                                        #    sed -e '/<\/body>/r${scfile}' "${htmlFilePath}"
                                                        #else
                                                        #    echo >&2
                                                        #    echo >&2 "${pmnote} ${warning} ${paramonte_dir}/auxil/sc.html is missing in your clone."
                                                        #    echo >&2 "${pmnote} ${warning} This is not critical, unless you are a ParaMonte developer and"
                                                        #    echo >&2 "${pmnote} ${warning} aim to publicly release this code coverage report. To obtain a "
                                                        #    echo >&2 "${pmnote} ${warning} copy of the file, contact the ParaMonte lead developer at"
                                                        #    echo >&2 "${pmnote} ${warning} "
                                                        #    echo >&2 "${pmnote} ${warning} shahmoradi@utexas.edu"
                                                        #    echo >&2
                                                        #fi

                                                    } || {

                                                        echo >&2 "${pmwarn} Code Coverage report generation via genhtml failed."

                                                    }
                                                    # "${branchCoverageFlag}" \
                                                    #--title "<a href=\"https://github.com/cdslaborg/paramonte\" target=\"_blank\">ParaMonte Fortran</a> code coverage report" \

                                                    ## generate test files code coverage
                                                    #
                                                    #gcovFortranTestDataDir=$(find "${paramonte_bld_obj_dir}" -name test_pm_*.o)
                                                    #gcovFortranTestDataDir=$(dirname "${gcovFortranTestDataDir}")
                                                    #
                                                    #lcovFortranTestDir="${paramonte_bld_dir}"/test/lcov
                                                    #mkdir -p "${lcovFortranTestDir}"
                                                    #cd "${lcovFortranTestDir}"
                                                    #
                                                    #lcov --capture --directory "${gcovFortranTestDataDir}" --output-file ./paramonte.coverage.info
                                                    #
                                                    #genhtml paramonte.coverage.info --output-directory "${lcovFortranTestDir}/html"

                                                else
                                                    echo >&2 "${pmwarn} Failed to find the GENHTML test coverage summarizer."
                                                    echo >&2 "${pmwarn} The genhtml program is required to generate the coverage report."
                                                    echo >&2 "${pmwarn} If you believe genhtml is already installed on your system,"
                                                    echo >&2 "${pmwarn} please make sure the path its directory is added to the"
                                                    echo >&2 "${pmwarn} PATH environmental variable of your terminal."
                                                    echo >&2 "${pmwarn} Once added, rerun the ParaMonte code coverage."
                                                fi

                                            else
                                                echo >&2 "${pmwarn} Failed to find the LCOV test coverage summarizer."
                                                echo >&2 "${pmwarn} The lcov program is required to generate the coverage report."
                                                echo >&2 "${pmwarn} If you believe lcov is already installed on your system,"
                                                echo >&2 "${pmwarn} please make sure the path its directory is added to the"
                                                echo >&2 "${pmwarn} PATH environmental variable of your terminal."
                                                echo >&2 "${pmwarn} Once added, rerun the ParaMonte code coverage."
                                            fi

                                            cd "${paramonte_dir}"

                                        else
                                            echo >&2 "${pmfatal} Failed to find the ParaMonte library objects directory."
                                            echo >&2 "${pmfatal} "
                                            echo >&2 "${pmfatal} gracefully exiting The ParaMonte build script."
                                            exit 1
                                        fi

                                    else
                                        echo >&2 "${pmfatal} Fatal Error: Failed to find the GNU gcov test coverage program."
                                        echo >&2 "${pmfatal} The gcov program is required to generate the coverage report."
                                        echo >&2 "${pmfatal} If you believe gcov is already installed on your system,"
                                        echo >&2 "${pmfatal} please make sure the path its directory is added to the"
                                        echo >&2 "${pmfatal} PATH environmental variable of your terminal."
                                        echo >&2 "${pmfatal} Once added, rerun the ParaMonte code coverage."
                                        echo >&2 "${pmfatal} "
                                        echo >&2 "${pmfatal} Gracefully exiting The ParaMonte build script."
                                        exit 1
                                    fi

                                else

                                    echo >&2 "${pmfatal} Code coverage with compilers other than GNU gfortran is currently unsupported."
                                    exit 1
                                fi

                            fi # codecov

                        done

                    done

                done

            done

        done

    done

done

echo >&2 ""
echo >&2 "${pmnote} All build files for all requested build configurations are stored at: \"${paramonte_dir}/bld/\""
if [[ ! "${flag_codecov}" =~ .*"true".* ]]; then
    echo >&2 "${pmnote} The installed binary files for all requested build configurations are ready to use at: \"${ddir}/\""
fi
echo >&2 ""
echo >&2 "${pmnote} ${BoldGreen}mission accomplished.${ColorReset}"
echo >&2 ""
