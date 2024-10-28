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
    unset flag_matlabroot
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

    ####
    #### Flag to skip timely production steps in development and documentation modes.
    #### if `--dev` is specified  by the user, the macro `dev_enabled` will be set to `1`,
    #### preventing the package copy to the final binary installation directory.
    ####
    #### The variable `ntry` is to bypass the need for duplicate build with CMake for development and testing times.
    #### The duplicate build with CMake is required to ensure the generation of FPP source files in the output package.
    #### This need for duplicate builds is an issue within the current CMake build scripts of the ParaMonte library that
    #### must be resolved in the future with a better solution.
    ####

    flag_dev="-Ddev_enabled=0"
    ntry=2

    while [ "$1" != "" ]; do
        case "$1" in

            #### list args

            --build )       shift
                            verifyArgNotKey "$1" --build
                            #verifyArgNotEmpty "$1" --build
                            list_build="$1"
                            ;;

            --checking )    shift
                            verifyArgNotKey "$1" --checking
                            #verifyArgNotEmpty "$1" --checking
                            list_checking="$1"
                            ;;

            --fc )          shift
                            verifyArgNotKey "$1" --fc
                            #verifyArgNotEmpty "$1" --fc
                            list_fc="$1"
                            ;;

            --lang )        shift
                            verifyArgNotKey "$1" --lang
                            #verifyArgNotEmpty "$1" --lang
                            list_lang="$1"
                            ;;

            --lib )         shift
                            verifyArgNotKey "$1" --lib
                            #verifyArgNotEmpty "$1" --lib
                            list_lib="$1"
                            ;;

            --mem )         shift
                            verifyArgNotKey "$1" --mem
                            #verifyArgNotEmpty "$1" --mem
                            list_mem="$1"
                            ;;

            --par )         shift
                            verifyArgNotKey "$1" --par
                            #verifyArgNotEmpty "$1" --par
                            list_par="$1"
                            ;;

            #### flag args

            --bench )       shift
                            verifyArgNotKey "$1" --bench
                            #verifyArgNotEmpty "$1" --bench
                            flag_bench="-Dbench=$1"
                            ;;

            --benchpp )     shift
                            verifyArgNotKey "$1" --benchpp
                            #verifyArgNotEmpty "$1" --benchpp
                            flag_benchpp="-Dbenchpp=$1"
                            ;;

            --blas )        shift
                            verifyArgNotKey "$1" --blas
                            #verifyArgNotEmpty "$1" --blas
                            flag_blas="-Dblas=$1"
                            ;;

            --codecov )     shift
                            verifyArgNotKey "$1" --codecov
                            #verifyArgNotEmpty "$1" --codecov
                            flag_codecov="-Dcodecov=$1"
                            ;;

            --deps )        shift
                            verifyArgNotKey "$1" --deps
                            #verifyArgNotEmpty "$1" --deps
                            flag_deps="-Ddeps=$1"
                            ;;

            --exam )        shift
                            verifyArgNotKey "$1" --exam
                            #verifyArgNotEmpty "$1" --exam
                            flag_exam="-Dexam=$1"
                            ;;

            --exampp )      shift
                            verifyArgNotKey "$1" --exampp
                            #verifyArgNotEmpty "$1" --exampp
                            flag_exampp="-Dexampp=$1"
                            ;;

            --cfi )         shift
                            verifyArgNotKey "$1" --cfi
                            #verifyArgNotEmpty "$1" --cfi
                            flag_cfi="-Dcfi=$1"
                            ;;

            --fpp )         shift
                            verifyArgNotKey "$1" --fpp
                            #verifyArgNotEmpty "$1" --fpp
                            flag_fpp="-Dfpp=$1"
                            ;;

            --fresh )       shift
                            verifyArgNotKey "$1" --fresh
                            #verifyArgNotEmpty "$1" --fresh
                            flag_fresh="-Dfresh=$1"
                            ;;

            --lapack )      shift
                            verifyArgNotKey "$1" --lapack
                            #verifyArgNotEmpty "$1" --lapack
                            flag_lapack="-Dlapack=$1"
                            ;;

            --matlabroot )  shift
                            verifyArgNotKey "$1" --matlabroot
                            #verifyArgNotEmpty "$1" --matlabroot
                            flag_matlabroot="-Dmatlabroot=\"$1\""
                            ;;

            --me )          shift
                            verifyArgNotKey "$1" --me
                            #verifyArgNotEmpty "$1" --me
                            flag_me="-Dme=\"$1\""
                            ;;

            --mod )         shift
                            verifyArgNotKey "$1" --mod
                            #verifyArgNotEmpty "$1" --mod
                            flag_mod="-Dmod=$1"
                            ;;

            --nproc )       shift
                            verifyArgNotKey "$1" --nproc
                            #verifyArgNotEmpty "$1" --nproc
                            isNumericValue="$(isnumeric ${nproc})"
                            if ! [ "${isNumericValue}" = "true" ]; then
                                reportBadValue "--nproc" $nproc "The spoecified number of processors must be a positive integer."
                            fi
                            FOR_COARRAY_NUM_IMAGES="$1"
                            flag_nproc="-Dnproc=$1"
                            ;;

            --perfprof )    shift
                            verifyArgNotKey "$1" --perfprof
                            #verifyArgNotEmpty "$1" --perfprof
                            flag_perfprof="-Dperfprof=$1"
                            ;;

            --pdt )         shift
                            verifyArgNotKey "$1" --pdt
                            #verifyArgNotEmpty "$1" --pdt
                            flag_pdt="-Dpdt=$1"
                            ;;

            --purity )      shift
                            verifyArgNotKey "$1" --purity
                            #verifyArgNotEmpty "$1" --purity
                            flag_purity="-Dpurity=$1"
                            ;;

            --test )        shift
                            verifyArgNotKey "$1" --test
                            #verifyArgNotEmpty "$1" --test
                            flag_test="-Dtest=$1"
                            ;;

            #### flag args: type kind

            --ski )         shift
                            verifyArgNotKey "$1" --ski
                            #verifyArgNotEmpty "$1" --ski
                            flag_ski="-Dski=$1"
                            ;;

            --iki )         shift
                            verifyArgNotKey "$1" --iki
                            #verifyArgNotEmpty "$1" --iki
                            flag_iki="-Diki=$1"
                            ;;

            --lki )         shift
                            verifyArgNotKey "$1" --lki
                            #verifyArgNotEmpty "$1" --lki
                            flag_lki="-Dlki=$1"
                            ;;

            --cki )         shift
                            verifyArgNotKey "$1" --cki
                            #verifyArgNotEmpty "$1" --cki
                            flag_cki="-Dcki=$1"
                            ;;

            --rki )         shift
                            verifyArgNotKey "$1" --rki
                            #verifyArgNotEmpty "$1" --rki
                            flag_rki="-Drki=$1"
                            ;;

            #### other args

            --bdir )        shift
                            verifyArgNotKey "$1" --bdir
                            #verifyArgNotEmpty "$1" --bdir
                            bdir="$1"
                            ;;

            --ddir )        shift
                            verifyArgNotKey "$1" --ddir
                            #verifyArgNotEmpty "$1" --ddir
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
                            #verifyArgNotEmpty "$1" -j
                            flag_j="-j $1"
                            ;;

            --dev )         flag_dev="-Ddev_enabled=1"
                            ntry=1
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

    ####
    #### list args
    ####

    if  [ "${list_build}" = "" ]; then
        list_build="release"
    else
        list_build="$(getLowerCase "$list_build")"
    fi

    if  [ "${list_checking}" = "" ]; then
        list_checking="nocheck"
    else
        list_checking="$(getLowerCase "$list_checking")"
    fi

    ####
    #### Set the Fortran compiler.
    ####

    if [ "${list_fc}" = "" ]; then
        #   On Windows OS, particularly with MinGW, CMake fails if the specified compiler name or path does not have the
        #   the file extension ".exe". Given that this extension is unlikely to change in the future, and that it is not used
        #   on Unix systems, try suffixing the extension to the compiler file path. If it fails, use the default specified path.
        #fccompilers=("ifort" "ifx" "gfortran" "gfortran-10" "gfortran-11" "gfortran-12" "gfortran-13" "gfortran-14" "gfortran-15" "gfortran-16" "gfortran-17" "gfortran-18" "gfortran-19" "gfortran-20")
        fccompilers=("ifort" "ifx" "gfortran-20" "gfortran-19" "gfortran-18" "gfortran-17" "gfortran-16" "gfortran-15" "gfortran-14" "gfortran-13" "gfortran-12" "gfortran-11" "gfortran-10" "gfortran")
        if [ "${os}" = "mingw" ] || [ "${os}" = "msys" ] || [ "${os}" = "cygwin" ]; then
            extensions=(".exe" "")
        else
            extensions=("")
        fi
        for fccompiler in "${fccompilers[@]}"; do
            for extension in "${extensions[@]}"; do
                echo >&2 "${pmnote} Checking for the existence of the Fortran compiler \"${fccompiler}\" executable with extension \"${extension}\""
                # check if the specified compiler can be found in the environment.
                if  command -v "${fccompiler}${extension}" >/dev/null 2>&1; then
                    list_fc="$(command -v "${fccompiler}${extension}")"
                    break 2
                fi
            done
        done

        if  [ "${list_fc}" = "" ]; then
            echo >&2 "${pmwarn} Failed to detect any compatible Fortran compiler in the environment."
            echo >&2 "${pmwarn} You can manually specify the Fortran compiler or its path via the install.sh script flag \"--fc\"."
            echo >&2 "${pmwarn} The build will proceed with no guarantee of success."
            list_fc="default"
        else
            echo >&2 "${pmnote} The identified Fortran compiler path is: fc=\"${list_fc}\""
        fi
    fi

    ####
    #### Set the library language.
    ####

    if  [ "${list_lang}" = "" ]; then
        list_lang="fortran"
    else
        list_lang="$(getLowerCase "$list_lang")"
        list_lang="${list_lang/c++/cpp}"
    fi

    if  [ "${list_lib}" = "" ]; then
        list_lib="shared"
    else
        # Replace `dynamic` with `shared`.
        list_lib=${list_lib/dynamic/shared}
        list_lib="$(getLowerCase "$list_lib")"
    fi

    if  [ "${list_mem}" = "" ]; then
        list_mem="heap"
    else
        list_mem="$(getLowerCase "$list_mem")"
    fi

    if  [ "${list_par}" = "" ]; then
        list_par="serial"
    else
        # Replace `none` with `serial`.
        list_par=${list_par/none/serial}
        list_par="$(getLowerCase "$list_par")"
    fi

    if  [ "${flag_j}" = "" ]; then
        flag_j="-j"
    fi

    #### Optional arguments

    if  [ "${flag_exampp}" = "" ] && ! [ "${flag_exam}" = "" ] ; then
        flag_exampp="${flag_exam/-Dexam=/-Dexampp=}"
    fi

    if  [ "${flag_benchpp}" = "" ] && ! [ "${flag_bench}" = "" ] ; then
        flag_benchpp="${flag_bench/-Dbench=/-Dbenchpp=}"
    fi

fi

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set CMake default flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if [ 0 -lt 1 ]; then # just to allow toggling in notepad++.

    if  [[ "${flag_fresh}" =~ .*"prereq".* ]]; then
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

if [ 0 -lt 1 ]; then # just to allow further control.
for fc in ${list_fc//;/$'\n'}; do # replace `;` with newline character.

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

    for lang in ${list_lang//;/$'\n'}; do # replace `;` with newline character.

        flag_lang="-Dlang=${lang}"

        #lang_is_dynamic="false"
        #if [ "${lang}" = "julia" ] || [ "${lang}" = "matlab" ] || [ "${lang}" = "mathematica" ] || [ "${lang}" = "python" ] || [ "${lang}" = "r" ]; then
        #    lang_is_dynamic="true"
        #fi

        ## Set source preprocessing output on by default for Fortran.
        #if [ "${lang}" = "fortran" ] && [ "${flag_fpp}" = "" ]; then
        #    flag_fpp="-Dfpp=generic"
        #fi

        for build in ${list_build//;/$'\n'}; do # replace `;` with newline character.

            flag_build="-Dbuild=${build}"

            for lib in ${list_lib//;/$'\n'}; do # replace `;` with newline character.

                flag_lib="-Dlib=${lib}"

                for mem in ${list_mem//;/$'\n'}; do # replace `;` with newline character.

                    flag_mem="-Dmem=${mem}"

                    for par in ${list_par//;/$'\n'}; do # replace `;` with newline character.

                        flag_par="-Dpar=${par}"

                        for checking in ${list_checking//;/$'\n'}; do # replace `;` with newline character.

                            flag_checking="-Dchecking=${checking}"

                            #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                            # Set the ParaMonte CMake build directory.
                            #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                            if [ 0 -lt 1 ]; then # just to allow toggling in notepad++.

                                # First, determine the MPI library name to be used in the ParaMonte library name and build directory.

                                if [ "${par}" = "mpi" ]; then
                                   #if [ -z ${me+x} ]; then
                                    if [ "${me}" = "" ]; then
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

                               #if [ -z ${bdir+x} ]; then
                                if [ "${bdir}" = "" ]; then
                                    paramonte_bld_dir="${paramonte_dir}/bld/${os}/${arch}/${csid}/${csvs}/${build}/${lib}/${mem}/${parname}/${checking}/${lang}"
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

                                if  [[ "${flag_fresh}" =~ .*"all".* ]]; then
                                    if  [ -d "${paramonte_bld_dir}" ]; then
                                        echo >&2 "${pmnote} Removing the old prerequisites of the ParaMonte library build at: paramonte_bld_dir=\"${paramonte_bld_dir}\""
                                        rm -rf "${paramonte_bld_dir}"
                                    fi
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

                                #### The following loop temporarily bypasses an existing bug where the first fresh installation
                                #### does not copy the FPP source files to the deployment and installation directories.

                                for ((n=0; n<$ntry; n++)); do

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
                                    "${flag_dev}" \
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
                                    ${flag_matlabroot} \
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

                                done

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
                                htmlDir="${paramonte_external_doc_out_dir}/codecov/fortran/${paramonte_version_major_fortran}/${htmlSubDir}"
                                htmlTitleCodeCov="ParaMonte ${paramonte_version_fortran} :: ${parallelismText} Fortran - Code Coverage Report"

                                if [[ ${csid} == [gG][nN][uU] ]] && command -v gcov >/dev/null 2>&1 && command -v lcov >/dev/null 2>&1; then

                                    #   gcov should be run with the current directory the same as that when you invoked the compiler.
                                    #   Otherwise it will not be able to locate the source files.
                                    #   gcov produces files called mangledname.gcov in the current directory.
                                    #   These contain the coverage information of the source file they correspond to.
                                    #   One .gcov file is produced for each source (or header) file containing code, which was compiled to produce the data files.
                                    #   The mangledname part of the output file name is usually simply the source file name, but can be something more complicated
                                    #   if the `-l` or `-p` options are given. Refer to those options for details.
                                    #
                                    #   -o directory|file
                                    #   --object-directory directory
                                    #   --object-file file
                                    #
                                    #       Specify either the directory containing the gcov data files, or the object path name.
                                    #       The .gcno, and .gcda data files are searched for using this option.
                                    #       If a directory is specified, the data files are in that directory and named after the input file name, without its extension.
                                    #       If a file is specified here, the data files are named after that file, without its extension.

                                    gcovPath="$(command -v gcov)"
                                    echo >&2 "${pmnote} GNU gcov detected at: ${gcovPath}"
                                    echo >&2 "${pmnote} Invoking gcov to generate coverage report..."

                                    paramonte_bld_lcov_dir="${paramonte_bld_dir}"/lcov
                                    rm -rf "${paramonte_bld_lcov_dir}"
                                    mkdir -p "${paramonte_bld_lcov_dir}"

                                    lcovCombinedTraceFileName="all.tracefile.info"
                                    lcovCombinedTraceFilePath="${paramonte_bld_lcov_dir}/${lcovCombinedTraceFileName}"

                                    #unset branchCoverageFlag
                                    # Add the following flag to lcov to enable branch coverage:
                                    # branchCoverageFlag="--rc lcov_branch_coverage=1"
                                    # "${branchCoverageFlag}"

                                    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                                    #### GCOV main source: Generate *.gcda *.gcno codecov files for the Fortran main source files
                                    #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                                    if [ 0 -lt 1 ]; then # just to allow toggling in notepad++. Disable conditional to prevent example coverage report.

                                        colections=("main" "test")

                                        for collection in ${colections[@]}; do

                                            if [ "${collection}" = "main" ]; then
                                                # pm_blas is a unique file name. That is why we choose and use it.
                                                sobjFileDir=$(find "${paramonte_bld_dir}/obj" -name pm_blas*.o)
                                                gcovFileDir=$(find "${paramonte_bld_dir}/obj" -name pm_blas*.gcno)
                                            elif [ "${collection}" = "test" ]; then
                                                # main is a unique file name. That is why we choose and use it.
                                                sobjFileDir=$(find "${paramonte_bld_dir}/test/obj" -name main*.o)
                                                gcovFileDir=$(find "${paramonte_bld_dir}/test/obj" -name main*.gcno)
                                            else
                                                echo >&2 "${pmfatal} Unrecognized coverage report collection."
                                                echo >&2 "${pmfatal} This is an internal ParaMonte build script error."
                                                echo >&2 "${pmfatal} Please report this problem to the ParaMonte library developers on GitHub issues."
                                                exit 1
                                            fi

                                            sobjFileDir=$(dirname "${sobjFileDir}")
                                            gcovFileDir=$(dirname "${gcovFileDir}")
                                            echo >&2 "${pmnote} sobjFileDir=${sobjFileDir}"
                                            echo >&2 "${pmnote} gcovFileDir=${gcovFileDir}"

                                            if [ -d "${sobjFileDir}" ] && [ -d "${gcovFileDir}" ]; then

                                                cd "${gcovFileDir}"
                                                for objFilePath in "${sobjFileDir}"/*.o; do
                                                    srcFileName="${objFilePath##*/}" # extract the basename.
                                                    srcFileName="${srcFileName%.o}" # extract the filename (without extension).
                                                    srcFilePath="${paramonte_src_fortran_main_dir}/${srcFileName}"
                                                    echo >&2 "${pmnote} running gcov  ${srcFilePath}  -o  ${gcovFileDir}/${srcFileName}"
                                                                                gcov "${srcFilePath}" -o "${gcovFileDir}/${srcFileName}.null" || {
                                                                                    echo >&2 "${pmwarn} The ParaMonte Code Coverage analysis of the ${collection} object files via GNU gcov tool failed."
                                                                                    echo >&2 "${pmwarn} The ParaMonte ${collection} source file: ${srcFilePath}"
                                                                                    echo >&2 "${pmwarn} The ParaMonte ${collection} object file: ${objFilePath}"
                                                                                }
                                                done

                                                cd "${paramonte_bld_lcov_dir}"
                                                lcovCurrentTraceFilePath="${paramonte_bld_lcov_dir}/lcov.${collection}.tracefile.info"
                                                echo >&2 "${pmnote} Invoking LCOV report file for the ParaMonte Fortran ${collection} files: ${lcovCurrentTraceFilePath}"
                                                lcov --capture --directory "${gcovFileDir}" --output-file "${lcovCurrentTraceFilePath}" && {
                                                    if [ -f "${lcovCombinedTraceFilePath}" ]; then
                                                        echo >&2 "${pmnote} Merging lcov report file with the main report file: ${lcovCombinedTraceFilePath}"
                                                        mv "${lcovCombinedTraceFilePath}" "${lcovCombinedTraceFilePath}".temp
                                                        lcov --add-tracefile "${lcovCombinedTraceFilePath}".temp -a "${lcovCurrentTraceFilePath}" -o "${lcovCombinedTraceFilePath}"
                                                    else
                                                        cp -arf "${lcovCurrentTraceFilePath}" "${lcovCombinedTraceFilePath}"
                                                    fi
                                                } || {
                                                    echo >&2 "${pmwarn} Code Coverage report generation for the ParaMonte ${collection} source files via lcov tool failed."
                                                    #exit 1
                                                }

                                            else

                                                echo >&2 "${pmwarn} Failed to detect the *.o, *.gcda, and *.gcno data files for the ParaMonte ${collection} source files."
                                                echo >&2 "${pmwarn} The expected directory path for the object files: ${sobjFileDir}"
                                                echo >&2 "${pmwarn} The expected directory path for the GCOV   files: ${gcovFileDir}"
                                                echo >&2 "${pmwarn} The coverage report for the ParaMonte ${collection} source file will not be included."
                                                echo >&2 "${pmwarn} Skipping the GCOV code coverage report generation for the ${collection} files..."

                                            fi

                                        done

                                    fi

                                    ####
                                    #### Add the ParaMonte example gcov report files.
                                    ####

                                    if [ 1 -lt 0 ]; then # just to allow toggling in notepad++. Disable conditional to prevent example coverage report.

                                        collection="example"
                                        for modpath in "${paramonte_bld_dir}"/pkg/example/pm_arrayCenter*/; do
                                            for exppath in "${modpath}"/*/; do

                                                modname=$(basename "${modpath}")
                                                expname=$(basename "${exppath}")

                                                #### First attempt to infer the cmake object files' directory.

                                                sobjFileDir=$(find "${exppath}" -name main*.o)
                                                sobjFileDir=$(dirname "${sobjFileDir}")

                                                gcovFileDir=$(find "${exppath}" -name main*.gcno)
                                                gcovFileDir=$(dirname "${gcovFileDir}")

                                                if [ -d "${sobjFileDir}" ] && [ -d "${gcovFileDir}" ]; then

                                                    cd "${gcovFileDir}"
                                                    for objFilePath in "${sobjFileDir}"/*.o; do
                                                        srcFilePath="${paramonte_example_dir}/${lang}/${modname}/${expname}/main.F90"
                                                        echo >&2 "${pmnote} running gcov  ${srcFilePath}  --relative-only --source-prefix  ${paramonte_bld_dir}/pkg/  -o  ${gcovFileDir}/main.F90"
                                                                                    gcov "${srcFilePath}" --relative-only --source-prefix "${paramonte_bld_dir}/pkg/" -o "${gcovFileDir}/main.F90.null" || {
                                                                                        echo >&2 "${pmwarn} The ParaMonte Code Coverage analysis of the test object files via GNU gcov tool failed."
                                                                                        echo >&2 "${pmwarn} The ParaMonte test source file: ${srcFilePath}"
                                                                                        echo >&2 "${pmwarn} The ParaMonte test object file: ${objFilePath}"
                                                                                    }
                                                    done

                                                    cd "${paramonte_bld_lcov_dir}"
                                                    lcovCurrentTraceFilePath="${paramonte_bld_lcov_dir}/lcov.${collection}.${modname}.${expname}.tracefile.info"
                                                    echo >&2 "${pmnote} Invoking LCOV report file for the ParaMonte Fortran ${collection} files: ${lcovCurrentTraceFilePath}"
                                                    lcov --capture --directory "${gcovFileDir}" --output-file "${lcovCurrentTraceFilePath}" && {
                                                        if [ -f "${lcovCombinedTraceFilePath}" ]; then
                                                            echo >&2 "${pmnote} Merging lcov report file with the main report file: ${lcovCombinedTraceFilePath}"
                                                            mv "${lcovCombinedTraceFilePath}" "${lcovCombinedTraceFilePath}".temp
                                                            lcov --add-tracefile "${lcovCombinedTraceFilePath}".temp -a "${lcovCurrentTraceFilePath}" -o "${lcovCombinedTraceFilePath}"
                                                        else
                                                            cp -arf "${lcovCurrentTraceFilePath}" "${lcovCombinedTraceFilePath}"
                                                        fi
                                                    } || {
                                                        echo >&2 "${pmwarn} Code Coverage report generation failed for the ParaMonte ${collection}${modname}.${expname} source file: ${exppath}"
                                                        #exit 1
                                                    }

                                                else

                                                    echo >&2 "${pmwarn} Failed to detect the *.o, *.gcda, and *.gcno data files for the ParaMonte test source files."
                                                    echo >&2 "${pmwarn} The expected directory path for the object files: ${sobjFileDir}"
                                                    echo >&2 "${pmwarn} The expected directory path for the GCOV   files: ${gcovFileDir}"
                                                    echo >&2 "${pmwarn} The coverage report for the ParaMonte test source file will not be included."
                                                    echo >&2 "${pmwarn} Skipping the GCOV code coverage report generation for the test files..."

                                                fi

                                            done
                                        done

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

                                        #if ! [ -f "${lcovCombinedTraceFilePath}" ]; then
                                        #    cp "${lcovOutputFortranFilePath}" "${lcovCombinedTraceFilePath}" || {
                                        #        echo >&2 "${pmwarn} Copy action failed:"
                                        #        echo >&2 "${pmwarn} from: ${lcovOutputFortranFilePath}"
                                        #        echo >&2 "${pmwarn}   to: ${lcovCombinedTraceFilePath}"
                                        #    }
                                        #fi

                                        genhtml "${lcovCombinedTraceFilePath}" --output-directory "${htmlDir}" --legend --title "${htmlTitleCodeCov}" && {

                                            echo >&2 "${pmnote} The code coverage build files are stored at: ${paramonte_bld_dir}"
                                            echo >&2 "${pmnote} The code coverage report files are stored at: ${htmlDir}"

                                            # postprocess the html files.

                                            # The variable `original` contains what LCOV software generates.
                                            # The variable `modified` contains what what we intend to replace the LCOV header with.

                                            original="<tr><td class=\"title\">LCOV - code coverage report<\/td><\/tr>"
                                            modified="<tr><td class=\"title\">"
                                            modified+="<a href=\"https:\/\/www.cdslab.org\/paramonte\/fortran\/${paramonte_version_major_fortran}\" target=\"_blank\">"
                                            modified+="<img alt=\"https:\/\/www.cdslab.org\/paramonte\/fortran\/${paramonte_version_major_fortran}\" src=\"https:\/\/raw.githubusercontent.com\/cdslaborg\/paramonte\/16e8fc2abd6a4263c9f6ee25f7a9f45435443688\/img\/banner.png\"\/>"
                                            modified+="<\/a>"
                                            modified+="<\/td><\/tr>"

                                            footer="<tr><td class=\"versionInfo\">"
                                            footer+="<a href=\"https:\/\/www.cdslab.org\/paramonte\"><b>ParaMonte: Parallel Monte Carlo and Machine Learning Library<\/b><\/a>&nbsp;<br>"
                                            footer+="<a href=\"https:\/\/www.cdslab.org\" target=\"_blank\"><b>The Computational Data Science Lab<\/b><\/a><br>"
                                            footer+="&copy; Copyright 2012 - $(date +%Y)"
                                            footer+="<\/td><\/tr>"

                                            shopt -s globstar
                                            for htmlFilePath in "${htmlDir}"/**/*.html; do # Whitespace-safe and recursive
                                                #sed -i "s/${paramonte_bld_dir//\//\\/}\/pkg\///g" "${htmlFilePath}"
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
                                        #sobjFortranTestDir=$(find "${paramonte_bld_obj_dir}" -name test_pm_*.o)
                                        #sobjFortranTestDir=$(dirname "${sobjFortranTestDir}")
                                        #
                                        #lcovFortranTestDir="${paramonte_bld_dir}"/test/lcov
                                        #mkdir -p "${lcovFortranTestDir}"
                                        #cd "${lcovFortranTestDir}"
                                        #
                                        #lcov --capture --directory "${sobjFortranTestDir}" --output-file ./paramonte.coverage.info
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

                                    cd "${paramonte_dir}"

                                else

                                    echo >&2 "${pmfatal} Note: Code coverage with compilers other than GNU gfortran is currently unsupported."
                                    echo >&2 "${pmfatal} Fatal Error: Failed to find the GNU compilers, GNU GCOV and/or LCOV coverage software."
                                    echo >&2 "${pmfatal} The GCOV and LCOV software are required to generate the coverage report."
                                    echo >&2 "${pmfatal} If you believe GCOV/LCOV are already installed on your system,"
                                    echo >&2 "${pmfatal} please make sure the path to their binary directories are"
                                    echo >&2 "${pmfatal} added to the PATH environmental variable of your terminal."
                                    echo >&2 "${pmfatal} Once added, rerun the ParaMonte code coverage."
                                    echo >&2 "${pmfatal} "
                                    echo >&2 "${pmfatal} Gracefully exiting The ParaMonte build script."
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
fi

echo >&2 ""
echo >&2 "${pmnote} All build files for all requested build configurations are stored at: \"${paramonte_dir}/bld/\""
if [[ ! "${flag_codecov}" =~ .*"true".* ]]; then
    echo >&2 "${pmnote} The installed binary files for all requested build configurations are ready to use at: \"${ddir}/\""
fi

####
#### Compress all binary folders.
####

if [ -d "${ddir}" ]; then
    if  command -v "tar" >/dev/null 2>&1; then
        echo >&2
        echo >&2 "${pmnote} Compressing all subdirectories in the directory: ${ddir}"
        echo >&2
        cd "${ddir}"
        for subdir in ./*; do
            if [ -d "${subdir}" ]; then
                compressionEnabled=false
                if [[ "${subdir}" =~ .*"_c".* ]]; then compressionEnabled=true; fi
                if [[ "${subdir}" =~ .*"_cpp".* ]]; then compressionEnabled=true; fi
                if [[ "${subdir}" =~ .*"_fortran".* ]]; then compressionEnabled=true; fi
                if [[ "${subdir}" =~ .*"_matlab".* ]]; then compressionEnabled=true; fi
                if [[ "${subdir}" =~ .*"_python".* ]]; then compressionEnabled=true; fi
                if [ "${compressionEnabled}" = "true" ]; then
                    tarfile="${subdir}.tar.gz"
                    #cd "${subdir}"
                    if [ -f "${tarfile}" ]; then
                        echo >&2 "${pmwarn} Compressed subdirectory already exists: ${tarfile}"
                        echo >&2 "${pmwarn} Overwriting the existing archive file..."
                    fi
                    echo >&2 "${pmnote} Compressing subdirectory: ${subdir}"
                    tar -cvzf ${tarfile} --exclude="${subdir}/setup.sh" "${subdir}" && {
                        echo >&2 "${pmnote} Subdirectory compressed: ${tarfile}"
                    }|| {
                        echo >&2
                        echo >&2 "${pmfatal} Compression failed for subdirectory: ${subdir}"
                        echo >&2 "${pmfatal} Gracefully exiting."
                        echo >&2
                        exit 1
                    }
                    #cd ..
                fi
            else
                echo >&2 "${pmnote} Non-directory object detected: ${subdir}"
            fi
        done
    else
        echo >&2
        echo >&2 "${pmwarn} The tar archive maker cannot be found on your system."
        echo >&2 "${pmwarn} Skipping the archive file generation from the output binary files..."
        echo >&2
    fi
elif ! [ "${flag_deps}" = "" ]; then
    echo >&2
    echo >&2 "${pmwarn} The requested input target directory ${ddir} specified with the input flag --dir does not exist."
    echo >&2
fi

echo >&2 ""
echo >&2 "${pmnote} ${BoldGreen}mission accomplished.${ColorReset}"
echo >&2 ""
