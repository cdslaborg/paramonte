!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define a quoting macro for better illustration in the output file.
!> \cond excluded
#ifdef  __GFORTRAN__
# define QSTART(X) "&
# define QEND(X) &X"
#else   /* default quoting macro */
# define GET_QUOTED(X) #X
# define QSTART(X) &
# define QEND(X) GET_QUOTED(X)
#endif
!> \endcond excluded

!>  \brief
!>  This module contains procedures and data that provide general information about the ParaMonte library, its interfaces, and its build.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday 00:01 AM, January 1, 2018, Institute for Computational Engineering and Sciences, University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_paramonte

    use pm_kind, only: SK, IK, LK
    use pm_str, only: NLC, UNDEFINED
    use iso_fortran_env, only: compiler_options, compiler_version

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar constant of type `character` of default kind \SK,
    !>  containing the web portal address for the ParaMonte library repository.<br>
    !>
    !>  \see
    !>  [paramonte](@ref pm_paramonte::paramonte)<br>
    !>  [PARAMONTE_WEB_DOC](@ref pm_paramonte::PARAMONTE_WEB_DOC)<br>
    !>  [PARAMONTE_WEB_ISSUES](@ref pm_paramonte::PARAMONTE_WEB_ISSUES)<br>
    !>  [PARAMONTE_WEB_REPOSITORY](@ref pm_paramonte::PARAMONTE_WEB_REPOSITORY)<br>
    !>  [PARAMONTE_COMPILER_OPTIONS](@ref pm_paramonte::PARAMONTE_COMPILER_OPTIONS)<br>
    !>  [PARAMONTE_COMPILER_VERSION](@ref pm_paramonte::PARAMONTE_COMPILER_VERSION)<br>
    !>  [paramonte_type](@ref pm_paramonte::paramonte_type)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>
    !>  \final{PARAMONTE_WEB_REPOSITORY}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: PARAMONTE_WEB_REPOSITORY = SK_"https://github.com/cdslaborg/paramonte"

    !>  \brief
    !>  The scalar constant of type `character` of default kind \SK,
    !>  containing the web portal address for reporting the ParaMonte library issues.<br>
    !>
    !>  \see
    !>  [paramonte](@ref pm_paramonte::paramonte)<br>
    !>  [PARAMONTE_WEB_DOC](@ref pm_paramonte::PARAMONTE_WEB_DOC)<br>
    !>  [PARAMONTE_WEB_ISSUES](@ref pm_paramonte::PARAMONTE_WEB_ISSUES)<br>
    !>  [PARAMONTE_WEB_REPOSITORY](@ref pm_paramonte::PARAMONTE_WEB_REPOSITORY)<br>
    !>  [PARAMONTE_COMPILER_OPTIONS](@ref pm_paramonte::PARAMONTE_COMPILER_OPTIONS)<br>
    !>  [PARAMONTE_COMPILER_VERSION](@ref pm_paramonte::PARAMONTE_COMPILER_VERSION)<br>
    !>  [paramonte_type](@ref pm_paramonte::paramonte_type)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>
    !>  \final{PARAMONTE_WEB_ISSUES}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: PARAMONTE_WEB_ISSUES = SK_"https://github.com/cdslaborg/paramonte/issues"

    !>  \brief
    !>  The scalar constant of type `character` of default kind \SK,
    !>  containing the web portal address for reporting the ParaMonte library issues.<br>
    !>
    !>  \see
    !>  [paramonte](@ref pm_paramonte::paramonte)<br>
    !>  [PARAMONTE_WEB_DOC](@ref pm_paramonte::PARAMONTE_WEB_DOC)<br>
    !>  [PARAMONTE_WEB_ISSUES](@ref pm_paramonte::PARAMONTE_WEB_ISSUES)<br>
    !>  [PARAMONTE_WEB_REPOSITORY](@ref pm_paramonte::PARAMONTE_WEB_REPOSITORY)<br>
    !>  [PARAMONTE_COMPILER_OPTIONS](@ref pm_paramonte::PARAMONTE_COMPILER_OPTIONS)<br>
    !>  [PARAMONTE_COMPILER_VERSION](@ref pm_paramonte::PARAMONTE_COMPILER_VERSION)<br>
    !>  [paramonte_type](@ref pm_paramonte::paramonte_type)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>
    !>  \final{PARAMONTE_WEB_DOC}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: PARAMONTE_WEB_DOC = SK_"https://www.cdslab.org/paramonte/"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar constant of type `character` of default kind \SK,
    !>  containing the version of the compiler with which the ParaMonte library has been built.<br>
    !>
    !>  \see
    !>  [paramonte](@ref pm_paramonte::paramonte)<br>
    !>  [PARAMONTE_WEB_DOC](@ref pm_paramonte::PARAMONTE_WEB_DOC)<br>
    !>  [PARAMONTE_WEB_ISSUES](@ref pm_paramonte::PARAMONTE_WEB_ISSUES)<br>
    !>  [PARAMONTE_WEB_REPOSITORY](@ref pm_paramonte::PARAMONTE_WEB_REPOSITORY)<br>
    !>  [PARAMONTE_COMPILER_OPTIONS](@ref pm_paramonte::PARAMONTE_COMPILER_OPTIONS)<br>
    !>  [PARAMONTE_COMPILER_VERSION](@ref pm_paramonte::PARAMONTE_COMPILER_VERSION)<br>
    !>  [paramonte_type](@ref pm_paramonte::paramonte_type)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>
    !>  \final{PARAMONTE_COMPILER_VERSION}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: PARAMONTE_COMPILER_VERSION = compiler_version()

    !>  \brief
    !>  The scalar constant of type `character` of default kind \SK,
    !>  containing the compiler options used to build the ParaMonte library.<br>
    !>
    !>  \see
    !>  [paramonte](@ref pm_paramonte::paramonte)<br>
    !>  [PARAMONTE_WEB_DOC](@ref pm_paramonte::PARAMONTE_WEB_DOC)<br>
    !>  [PARAMONTE_WEB_ISSUES](@ref pm_paramonte::PARAMONTE_WEB_ISSUES)<br>
    !>  [PARAMONTE_WEB_REPOSITORY](@ref pm_paramonte::PARAMONTE_WEB_REPOSITORY)<br>
    !>  [PARAMONTE_COMPILER_OPTIONS](@ref pm_paramonte::PARAMONTE_COMPILER_OPTIONS)<br>
    !>  [PARAMONTE_COMPILER_VERSION](@ref pm_paramonte::PARAMONTE_COMPILER_VERSION)<br>
    !>  [paramonte_type](@ref pm_paramonte::paramonte_type)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>
    !>  \final{PARAMONTE_COMPILER_OPTIONS}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: PARAMONTE_COMPILER_OPTIONS = compiler_options()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This a derived type devoid of any components or methods whose instantiated objects are used within the ParaMonte library
    !>  to signify the use of the ParaMonte style (vs. alternative approaches) in performing various actions within the library.<br>
    !>
    !>  \see
    !>  [paramonte](@ref pm_paramonte::paramonte)<br>
    !>  [PARAMONTE_WEB_DOC](@ref pm_paramonte::PARAMONTE_WEB_DOC)<br>
    !>  [PARAMONTE_WEB_ISSUES](@ref pm_paramonte::PARAMONTE_WEB_ISSUES)<br>
    !>  [PARAMONTE_WEB_REPOSITORY](@ref pm_paramonte::PARAMONTE_WEB_REPOSITORY)<br>
    !>  [PARAMONTE_COMPILER_OPTIONS](@ref pm_paramonte::PARAMONTE_COMPILER_OPTIONS)<br>
    !>  [PARAMONTE_COMPILER_VERSION](@ref pm_paramonte::PARAMONTE_COMPILER_VERSION)<br>
    !>  [paramonte_type](@ref pm_paramonte::paramonte_type)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>
    !>  \final{paramonte_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: paramonte_type
    end type

    !>  \brief
    !>  The scalar module constant of type [paramonte_type](@ref pm_paramonte::paramonte_type) that can be used within the ParaMonte library
    !>  to signify the use of the ParaMonte style (vs. alternative approaches) in performing various actions within the library.<br>
    !>
    !>  \see
    !>  [paramonte](@ref pm_paramonte::paramonte)<br>
    !>  [PARAMONTE_WEB_DOC](@ref pm_paramonte::PARAMONTE_WEB_DOC)<br>
    !>  [PARAMONTE_WEB_ISSUES](@ref pm_paramonte::PARAMONTE_WEB_ISSUES)<br>
    !>  [PARAMONTE_WEB_REPOSITORY](@ref pm_paramonte::PARAMONTE_WEB_REPOSITORY)<br>
    !>  [PARAMONTE_COMPILER_OPTIONS](@ref pm_paramonte::PARAMONTE_COMPILER_OPTIONS)<br>
    !>  [PARAMONTE_COMPILER_VERSION](@ref pm_paramonte::PARAMONTE_COMPILER_VERSION)<br>
    !>  [paramonte_type](@ref pm_paramonte::paramonte_type)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>
    !>  \final{paramonte}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(paramonte_type), parameter :: paramonte = paramonte_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: paramonte
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! The ParaMonte programming language interface.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \cond excluded
#if     C_ENABLED

#define PROGLANG c
#define PROGNAME"C"

#elif   COBOL_ENABLED

#define PROGLANG cobol
#define PROGNAME"COBOL"

#elif   CPP_ENABLED

#define PROGLANG cpp
#define PROGNAME"C++"

#elif   CSHARP_ENABLED

#define PROGLANG csharp
#define PROGNAME"C#"

#elif   GO_ENABLED

#define PROGLANG go
#define PROGNAME"Go"

#elif   FORTRAN_ENABLED

#define PROGLANG fortran
#define PROGNAME"Fortran"

#elif   JAVA_ENABLED

#define PROGLANG java
#define PROGNAME"Java"

#elif   JAVASCRIPT_ENABLED

#define PROGLANG javascript
#define PROGNAME"JavaScript"

#elif   JULIA_ENABLED

#define PROGLANG julia
#define PROGNAME"Julia"

#elif   MATLAB_ENABLED

#define PROGLANG matlab
#define PROGNAME"MATLAB"

#elif   MATTHEMATICA_ENABLED

#define PROGLANG mathematica
#define PROGNAME"Wolfram Mathematica"

#elif   PYTHON_ENABLED

#define PROGLANG python
#define PROGNAME"Python"

#elif   R_ENABLED

#define PROGLANG r
#define PROGNAME"R"

#elif   RUST_ENABLED

#define PROGLANG rust
#define PROGNAME"Rust"

#elif   SAS_ENABLED

#define PROGLANG sas
#define PROGNAME"SAS"

#elif   SWIFT_ENABLED

#define PROGLANG swift
#define PROGNAME"Swift"

#else
! Assume default Fortran.
#define PROGLANG fortran
#define PROGNAME"Fortran"
!#error  "Unrecognized interface."
#endif
!>  \endcond excluded

    !>  \brief
    !>  The scalar constant of type `character` of default kind \SK containing the full generic name
    !>  of the programming language environment for which the ParaMonte library has been configured.
    !>
    !>  \see
    !>  [envis](@ref pm_paramonte::envis)<br>
    !>  [envis_type](@ref pm_paramonte::envis_type)<br>
    !>  [envname](@ref pm_paramonte::envname)<br>
    !>
    !>  \final{filext_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: envname = PROGNAME//" programming language"

    !>  \brief
    !>  This is the derived type for generating objects containing scalar components of type `logical` named after
    !>  different programming language environments that are currently or will be supported for access to the ParaMonte library.
    !>
    !>  \details
    !>  This derived type is of minimal usage outside the ParaMonte library internal routines.<br>
    !>  If needed, the constant object instances of this derived type [envis](@ref pm_paramonte::envis) can be used.<br>
    !>
    !>  \see
    !>  [envis](@ref pm_paramonte::envis)<br>
    !>  [envis_type](@ref pm_paramonte::envis_type)<br>
    !>  [envname](@ref pm_paramonte::envname)<br>
    !>
    !>  \final{envis_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: envis_type
        logical(LK) :: c            = .false._LK
        logical(LK) :: cobol        = .false._LK
        logical(LK) :: cpp          = .false._LK
        logical(LK) :: csharp       = .false._LK
        logical(LK) :: go           = .false._LK
        logical(LK) :: fortran      = .false._LK
        logical(LK) :: java         = .false._LK
        logical(LK) :: javascript   = .false._LK
        logical(LK) :: julia        = .false._LK
        logical(LK) :: matlab       = .false._LK
        logical(LK) :: mathematica  = .false._LK
        logical(LK) :: python       = .false._LK
        logical(LK) :: r            = .false._LK
        logical(LK) :: rust         = .false._LK
        logical(LK) :: sas          = .false._LK
        logical(LK) :: swift        = .false._LK
    end type

    !>  \brief
    !>  The scalar constant object of type [envis_type](@ref pm_paramonte::envis_type) whose components
    !>  are all `.false.` except the environment for which the ParaMonte library is currently configured.
    !>
    !>  \see
    !>  [envis](@ref pm_paramonte::envis)<br>
    !>  [envis_type](@ref pm_paramonte::envis_type)<br>
    !>  [envname](@ref pm_paramonte::envname)<br>
    !>
    !>  \final{envis}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(envis_type), parameter :: envis = envis_type(PROGLANG = .true._LK)

!>  \cond excluded
#undef PROGLANG
#undef PROGNAME
!>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `public` scalar `character` constant of default kind \SK containing the ParaMonte library build date.
    !>
    !>  \details
    !>  The generated date is,
    !>  <ol>
    !>      <li>    the value of the predefined compiler preprocessor macro `__TIMESTAMP__` if the compiler suite is Intel or GNU.<br>
    !>      <li>    the value of CMake build generator `TIMESTAMP` with format `"%a %b %d %H:%M:%S %Y"` if the build is configured with CMake.<br>
    !>  </ol>
    !>
    !>  \interface{PARAMONTE_BUILD_DATE}
    !>  \code{.F90}
    !>
    !>      use pm_paramonte, only: PARAMONTE_BUILD_DATE
    !>
    !>      print *, PARAMONTE_BUILD_DATE
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [PARAMONTE_SPLASH](@ref pm_paramonte::PARAMONTE_SPLASH)<br>
    !>  [PARAMONTE_VERSION](@ref pm_paramonte::PARAMONTE_VERSION)<br>
    !>  [PARAMONTE_BUILD_DATE](@ref pm_paramonte::PARAMONTE_BUILD_DATE)<br>
    !>  [getParaMonteSplash](@ref pm_paramonte::getParaMonteSplash)<br>
    !>
    !>  \example{PARAMONTE_BUILD_DATE}
    !>  \include{lineno} example/pm_paramonte/PARAMONTE_BUILD_DATE/main.F90
    !>  \compilef{PARAMONTE_BUILD_DATE}
    !>  \output{PARAMONTE_BUILD_DATE}
    !>  \include{lineno} example/pm_paramonte/PARAMONTE_BUILD_DATE/main.out.F90
    !>
    !>  \final{PARAMONTE_BUILD_DATE}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if INTEL_ENABLED || GNU_ENABLED
    character(*, SK), parameter :: PARAMONTE_BUILD_DATE = __TIMESTAMP__
    !>  \cond excluded
#elif defined PARAMONTE_BUILD_TIMESTAMP
    character(*, SK), parameter :: PARAMONTE_BUILD_DATE = PARAMONTE_BUILD_TIMESTAMP
#else
    character(*, SK), parameter :: PARAMONTE_BUILD_DATE = UNDEFINED
#endif
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `public` scalar `character` constant of default kind \SK containing the ParaMonte library version.
    !>
    !>  \details
    !>  The generated version is the value reported in the
    !>  [VERSION.md](https://raw.githubusercontent.com/cdslaborg/paramonte/main/src/fortran/VERSION.md)
    !>  in the root directory of the ParaMonte library git repository.
    !>
    !>  \interface{PARAMONTE_VERSION}
    !>  \code{.F90}
    !>
    !>      use pm_paramonte, only: PARAMONTE_VERSION
    !>
    !>      print *, PARAMONTE_VERSION
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [PARAMONTE_SPLASH](@ref pm_paramonte::PARAMONTE_SPLASH)<br>
    !>  [PARAMONTE_VERSION](@ref pm_paramonte::PARAMONTE_VERSION)<br>
    !>  [PARAMONTE_BUILD_DATE](@ref pm_paramonte::PARAMONTE_BUILD_DATE)<br>
    !>  [getParaMonteSplash](@ref pm_paramonte::getParaMonteSplash)<br>
    !>
    !>  \example{PARAMONTE_VERSION}
    !>  \include{lineno} example/pm_paramonte/PARAMONTE_VERSION/main.F90
    !>  \compilef{PARAMONTE_VERSION}
    !>  \output{PARAMONTE_VERSION}
    !>  \include{lineno} example/pm_paramonte/PARAMONTE_VERSION/main.out.F90
    !>
    !>  \final{PARAMONTE_VERSION}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
#if defined __PARAMONTE_VERSION__
    character(*, SK), parameter :: PARAMONTE_VERSION = __PARAMONTE_VERSION__
    !>  \cond excluded
#else
    ! WARNING: The following ParaMonte library version tag will be taken from the declaration
    ! WARNING: that is generated by the preprocessing build script of the ParaMonte library.
    ! WARNING: This is superior to the above method of using the compiler Preprocessor
    ! WARNING: since CMAKE triggers a complete lengthy rebuild of the library when the
    ! WARNING: only the library version has changed.
#include "pm_paramonte@version.inc.F90"
#endif
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `public` scalar `character` constant of default kind \SK containing the ParaMonte splash information.
    !>
    !>  \details
    !>  The splash information includes the following items,
    !>  <ol>
    !>      <li>    The library name.
    !>      <li>    The library version.
    !>      <li>    The library build date.
    !>      <li>    The developing organization name.
    !>      <li>    The developing organization affiliations.
    !>      <li>    The project PI contact information.
    !>      <li>    The project GitHub web address.
    !>  </ol>
    !>
    !>  \interface{PARAMONTE_SPLASH}
    !>  \code{.F90}
    !>
    !>      use pm_paramonte, only: PARAMONTE_SPLASH
    !>
    !>      print *, PARAMONTE_SPLASH
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [PARAMONTE_SPLASH](@ref pm_paramonte::PARAMONTE_SPLASH)<br>
    !>  [PARAMONTE_VERSION](@ref pm_paramonte::PARAMONTE_VERSION)<br>
    !>  [PARAMONTE_BUILD_DATE](@ref pm_paramonte::PARAMONTE_BUILD_DATE)<br>
    !>  [getParaMonteSplash](@ref pm_paramonte::getParaMonteSplash)<br>
    !>
    !>  \example{PARAMONTE_SPLASH}
    !>  \include{lineno} example/pm_paramonte/PARAMONTE_SPLASH/main.F90
    !>  \compilef{PARAMONTE_SPLASH}
    !>  \output{PARAMONTE_SPLASH}
    !>  \include{lineno} example/pm_paramonte/PARAMONTE_SPLASH/main.out.F90
    !>
    !>  \final{PARAMONTE_SPLASH}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    character(*, SK), parameter :: PARAMONTE_SPLASH =   NLC//NLC// &
                                                        SK_"ParaMonte"//NLC// &
                                                        SK_"Parallel Monte Carlo and"//NLC// &
                                                        SK_"Machine Learning Library"//NLC// &
                                                        NLC// &
                                                        SK_"Version "//PARAMONTE_VERSION//NLC// &
                                                        NLC// &
                                                        SK_"Build: "//PARAMONTE_BUILD_DATE//NLC// &
                                                        NLC// &
                                                        SK_"Developed by"//NLC// &
                                                        NLC// &
                                                        SK_"The computational Data Science Lab"//NLC// &
                                                        NLC// &
                                                        SK_"at"//NLC// &
                                                        NLC// &
                                                        SK_"Department of Physics"//NLC// &
                                                        SK_"Data Science Program, College of Science"//NLC// &
                                                        SK_"The University of Texas at Arlington"//NLC// &
                                                        NLC// &
                                                        SK_"Multiscale Modeling Group"//NLC// &
                                                        SK_"Center for Computational Oncology"//NLC// &
                                                        SK_"Oden Institute for Computational Engineering and Sciences"//NLC// &
                                                        SK_"Department of Aerospace Engineering and Engineering Mechanics"//NLC// &
                                                        SK_"Department of Neurology, Dell-Seton Medical School"//NLC// &
                                                        SK_"Department of Biomedical Engineering"//NLC// &
                                                        SK_"The University of Texas at Austin"//NLC// &
                                                        NLC// &
                                                        SK_"For questions and further information, please contact the PI:"//NLC// &
                                                        NLC// &
                                                        SK_"Amir Shahmoradi"//NLC// &
                                                        !NLC// &
                                                        SK_"shahmoradi@utexas.edu"//NLC// &
                                                        !SK_"amir.shahmoradi@uta.edu"//NLC// &
                                                        !SK_"ashahmoradi@gmail.com"//NLC// &
                                                        !SK_"amir@physics.utexas.edu"//NLC// &
                                                        !SK_"amir@austin.utexas.edu"//NLC// &
                                                        !SK_"amir@ph.utexas.edu"//NLC// &
                                                        !NLC// &
                                                        !NLC// &
                                                        !SK_"https://www.cdslab.org/pm"//NLC// &
                                                        SK_"cdslab.org/pm"//NLC// &
                                                        NLC//NLC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an `allocatable` scalar `character` of default kind \SK
    !>  containing the ParaMonte library splash information centered with the requested line width, margins, and fillings.
    !>
    !>  \details
    !>  See the documentation of [PARAMONTE_SPLASH](@ref pm_paramonte::PARAMONTE_SPLASH) for further splash information details.<br>
    !>  The output splash begins and ends with a newline character returned by the intrinsic `new_line("a")` for the requested `character` kind.<br>
    !>
    !>  \param[in]      width   :   The input scalar `integer` of default kind \IK, containing the width of the splash **including** the lengths of the left and right margins.<br>
    !>                              (**optional**, default = `132`.)
    !>  \param[in]      lwsize  :   The input scalar `integer` of default kind \IK representing the width of the left wrapper margin of each line of the output `splash`.<br>
    !>                              (**optional**, default = `4`.)
    !>  \param[in]      twsize  :   The input scalar `integer` of default kind \IK representing the width of the top wrapper margin of the output `splash`.<br>
    !>                              (**optional**, default = `2`.)
    !>  \param[in]      rwsize  :   The input scalar `integer` of default kind \IK representing the width of the right wrapper margin of each line of the output `splash`.<br>
    !>                              (**optional**, default = `4`.)
    !>  \param[in]      bwsize  :   The input scalar `integer` of default kind \IK representing the width of the bottom wrapper margin of the output `splash`.<br>
    !>                              (**optional**, default = `2`.)
    !>  \param[in]      fill    :   The input scalar of the same type and kind as the output `slpash`, of length type parameter `1`,
    !>                              containing the value to fill the new elements (if any, excluding the margins) surrounding the splash text in each line of the output splash.<br>
    !>                              (**optional**, default = `" "`, the whitespace character.)
    !>  \param[in]      lwfill  :   The input scalar of the same type and kind as the output `slpash`, of length type parameter `1`,
    !>                              containing the value to fill the left wrapper margin (if any) of each line of the output `splash`.<br>
    !>                              (**optional**, default = `"%"`.)
    !>  \param[in]      twfill  :   The input scalar of the same type and kind as the output `slpash`, of length type parameter `1`,
    !>                              containing the value to fill the top wrapper margin (if any) of the output `splash`.<br>
    !>                              (**optional**, default = `"%"`.)
    !>  \param[in]      rwfill  :   The input scalar of the same type and kind as the output `slpash`, of length type parameter `1`,
    !>                              containing the value to fill the right wrapper margin (if any) of each line of the output `splash`.<br>
    !>                              (**optional**, default = `"%"`.)
    !>  \param[in]      bwfill  :   The input scalar of the same type and kind as the output `slpash`, of length type parameter `1`,
    !>                              containing the value to fill the bottom wrapper margin (if any) of the output `splash`.<br>
    !>                              (**optional**, default = `"%"`.)
    !>
    !>  \return
    !>  `splash`                :   The output `allocatable` scalar of type `character` of default kind \SK, containing the ParaMonte splash.<br>
    !>                              Lines within the `splash` are separated via the newline character returned by the intrinsic `new_line(SK_"a")`.<br>
    !>
    !>  \interface{getParaMonteSplash}
    !>  \code{.F90}
    !>
    !>      use pm_paramonte, only: getParaMonteSplash
    !>      character(:, SK), allocatable :: splash
    !>
    !>      splash = getParaMonteSplash(width = width, lwsize = lwsize, twsize = twsize, rwsize = rwsize, bwsize = bwsize, fill = fill, lwfill = lwfill, twfill = twfill, rwfill = rwfill, bwfill = bwfill)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [PARAMONTE_SPLASH](@ref pm_paramonte::PARAMONTE_SPLASH)<br>
    !>  [PARAMONTE_VERSION](@ref pm_paramonte::PARAMONTE_VERSION)<br>
    !>  [PARAMONTE_BUILD_DATE](@ref pm_paramonte::PARAMONTE_BUILD_DATE)<br>
    !>  [getParaMonteSplash](@ref pm_paramonte::getParaMonteSplash)<br>
    !>
    !>  \example{getParaMonteSplash}
    !>  \include{lineno} example/pm_paramonte/getParaMonteSplash/main.F90
    !>  \compilef{getParaMonteSplash}
    !>  \output{getParaMonteSplash}
    !>  \include{lineno} example/pm_paramonte/getParaMonteSplash/main.out.F90
    !>
    !>  \final{getParaMonteSplash}
    !>
    !>  \author
    !>  \AmirShahmoradi, March 22, 2012, 00:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    PURE function getParaMonteSplash(width, lwsize, twsize, rwsize, bwsize, fill, lwfill, twfill, rwfill, bwfill) result(splash)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getParaMonteSplash
#endif
        use pm_option, only: getOption
        use pm_arrayFind, only: setLoc
        use pm_arrayResize, only: setResized
        use pm_arrayCenter, only: setCentered
        integer(IK)     , intent(in)    , optional  :: width, lwsize, twsize, rwsize, bwsize
        character(1, SK), intent(in)    , optional  :: fill , lwfill, twfill, rwfill, bwfill
        character(:, SK), allocatable               :: splash
        integer(IK)     , allocatable               :: locNLC(:)
        character(1, SK)                            :: fill_def , lwfill_def, twfill_def, rwfill_def, bwfill_def
        integer(IK)                                 :: splashLen, iline, ibeg, iend, sbeg, lenNLC, linewidth, lenLocNLC
        integer(IK)                                 :: lwsize_def, twsize_def, rwsize_def, bwsize_def
        lwsize_def = getOption(4_IK, lwsize)
        twsize_def = getOption(2_IK, twsize)
        rwsize_def = getOption(4_IK, rwsize)
        bwsize_def = getOption(2_IK, bwsize)
        lwfill_def = getOption(SK_"%", lwfill)
        twfill_def = getOption(SK_"%", twfill)
        rwfill_def = getOption(SK_"%", rwfill)
        bwfill_def = getOption(SK_"%", bwfill)
        fill_def = getOption(SK_" ", fill)
        lenNLC = len(NLC, IK)
        call setResized(locNLC, 47_IK)
        call setLoc(locNLC, lenLocNLC, PARAMONTE_SPLASH, NLC, blindness = lenNLC)
        ! Ensure the width is minimally set to the maximum line length (~61) of the splash message.
        !print *, getOption(132_IK, width)
        linewidth = max(maxval(locNLC(2 : lenLocNLC) - locNLC(1 : lenLocNLC - 1)) - lenNLC, getOption(132_IK, width))! + lwsize_def + rwsize_def
        splashLen = lenNLC + (twsize_def + lenLocNLC + bwsize_def) * (linewidth + lenNLC)
        allocate(character(splashLen, SK) :: splash)
        splash(1 : lenNLC) = NLC
        ibeg = lenNLC + 1_IK
        do iline = 1, twsize_def
            iend = ibeg + linewidth - 1_IK
            splash(ibeg : iend) = repeat(twfill_def, linewidth)
            ibeg = iend + lenNLC + 1
            splash(iend + 1 : ibeg - 1) = NLC
        end do
        if (lenLocNLC > 0_IK) then
            sbeg = 1_IK ! splash line beginning.
            do iline = 1, lenLocNLC
                iend = ibeg + linewidth - 1_IK
                call setCentered( splash(ibeg : iend) & ! LCOV_EXCL_LINE
                                , PARAMONTE_SPLASH(sbeg : locNLC(iline) - 1) & ! LCOV_EXCL_LINE
                                , lmsize = lwsize_def & ! LCOV_EXCL_LINE
                                , rmsize = rwsize_def & ! LCOV_EXCL_LINE
                                , lmfill = lwfill_def & ! LCOV_EXCL_LINE
                                , rmfill = rwfill_def & ! LCOV_EXCL_LINE
                                , fill = fill_def & ! LCOV_EXCL_LINE
                                )
                sbeg = locNLC(iline) + lenNLC
                ibeg = iend + lenNLC + 1_IK
                splash(iend + 1 : ibeg - 1) = NLC
            end do
        else
            call setCentered( splash & ! LCOV_EXCL_LINE
                            , PARAMONTE_SPLASH & ! LCOV_EXCL_LINE
                            , lmsize = lwsize_def & ! LCOV_EXCL_LINE
                            , rmsize = rwsize_def & ! LCOV_EXCL_LINE
                            , lmfill = lwfill_def & ! LCOV_EXCL_LINE
                            , rmfill = rwfill_def & ! LCOV_EXCL_LINE
                            , fill = fill_def & ! LCOV_EXCL_LINE
                            ) ! LCOV_EXCL_LINE
        end if
        do iline = 1, bwsize_def
            iend = ibeg + linewidth - 1_IK
            splash(ibeg : iend) = repeat(bwfill_def, linewidth)
            ibeg = iend + lenNLC + 1
            splash(iend + 1 : ibeg - 1) = NLC
        end do
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_paramonte ! LCOV_EXCL_LINE