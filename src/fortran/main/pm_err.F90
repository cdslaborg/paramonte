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

!>  \brief
!>  This module contains classes and procedures for reporting and handling errors.<br>
!>
!>  \note
!>  There is a subtle difference between the impure [setAsserted](@ref pm_err::setAsserted)
!>  and the impure [setAborted](@ref pm_err::setAborted) generic interfaces.<br>
!>  <ol>
!>      <li>    The [setAsserted](@ref pm_err::setAsserted) procedures are frequently used internally
!>              within the various components of the ParaMonte library to assert a runtime condition.<br>
!>              Many of the assertions are related to the testing of the library.<br>
!>              Frequently however, the library procedures contain runtime checks that call
!>              [setAsserted](@ref pm_err::setAsserted) when the library is built with the preprocessing flag `CHECK_ENABLED=1`.<br>
!>              The failure of an assertion within [setAsserted](@ref pm_err::setAsserted) usually implies a more severe low-level error in the library.<br>
!>              The default behavior of [setAsserted](@ref pm_err::setAsserted) in such cases is to call `error stop` to halt the program unless it is
!>              specifically overridden by the user (via the input argument `renabled`).<br>
!>      <li>    By contrast, [setAborted](@ref pm_err::setAborted) should be usually called
!>              when an assertion is known to have failed **before** calling [setAborted](@ref pm_err::setAborted).<br>
!>              Furthermore, the default behavior of the program upon encountering an error depends on the compile-time
!>              constant [SOFT_EXIT_ENABLED](@ref pm_err::SOFT_EXIT_ENABLED) whose value depends on the build type and target
!>              programming language environment of the library (C, C++, Fortran, Java, Julia, MATLAB, Python, R, ...).<br>
!>              This behavior is critical for graceful handling of runtime errors in dynamic languages where a
!>              call to `error stop` shuts down the entire dynamic session causing data loss and wasting user time.<br>
!>              Additionally, [setAborted](@ref pm_err::setAborted) allows flexible formatting of
!>              the output message in both the stdout and an output unit specified by the user.<br>
!>  </ol>
!>
!>  [test_pm_err](@ref test_pm_err)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_err

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter         :: MODULE_NAME = "@pm_err"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#if CODECOV_ENABLED || BASIC_TEST_ENABLED || SAMPLER_TEST_ENABLED || ((MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED) && !CAF_ENABLED && !MPI_ENABLED)
#if (MATLAB_ENABLED || PYTHON_ENABLED || R_ENABLED) && !(CAF_ENABLED || MPI_ENABLED)
    !>  \cond excluded
    logical(LK)     , parameter         :: SOFT_EXIT_ENABLED = .true._LK
    !>  \endcond excluded
#else
    !>  \brief
    !>  The scalar `logical` of default kind \LK.<br>
    !>  If `.true.`, the [setAborted](@ref pm_err::setAborted) generic interface avoid calling the Fortran `error stop` statement upon
    !>  detecting the occurrence of an error.<br>
    !>  Instead it return the program control to the caller procedure for further actions.<br>
    !>  This flag is important for testing the parallel implementations of the routines of the ParaMonte library.<br>
    !>
    !>  \final{SOFT_EXIT_ENABLED}
    !>
    !>  \author
    !>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    logical(LK)     , parameter         :: SOFT_EXIT_ENABLED = .false._LK
#endif

    !  \brief
    !  The scalar `public` module variable of type `logical` of default kind \LK.<br>
    !  If `.true.`, the [setAborted](@ref pm_err::setAborted) subroutines avoid calling the Fortran `error stop` statement upon
    !  detecting the occurrence of an error. Instead they return the program control to the caller procedure for further actions.<br>
    !  This flag is important for testing the parallel implementations of the routines of the ParaMonte library.<br>
    !logical(LK) :: mv_isTestingMode = .false._LK

    !>  \brief
    !>  The scalar integer of default kind \IK containing the null value initially
    !>  assigned to the `stat` component of the [err_type](@ref pm_err::err_type).<br>
    !>  It is used to detect and control the printing of the compiler error `stat` number if any exists.<br>
    !>
    !>  \final{STATNULL}
    !>
    !>  \author
    !>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    integer(IK)                         :: STATNULL = -huge(0_IK)

   !logical(LK)     , parameter :: ERR_HANDLING_REQUESTED = .false._LK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating objects to gracefully and
    !>  verbosely handle runtime unexpected behavior in the ParaMonte library.<br>
    !>
    !>  \details
    !>  Always initialize objects of type [err_type](@ref pm_err::err_type)
    !>  if you prefer to avoid the manual allocation of the `msg` component of the object.<br>
    !>  The pre-allocated component can be readily passed to various Fortran routines for catching error messages.<br>
    !>
    !>  \param[in]  occurred    :   The input scalar `logical` of default kind \LK that is
    !>                              used to initialize the corresponding component of the output object.<br>
    !>                              (**optional**, default = `.false.`)
    !>  \param[in]  stat        :   The input scalar `integer` of default kind \IK that contains the error flag or status.<br>
    !>                              It is used to to initialize the corresponding component of the output object.<br>
    !>                              (**optional**, default = `-huge(0_IK)`)
    !>  \param[in]  msg         :   The input `allocatable` scalar `character` of default kind \SK containing the error message.<br>
    !>                              It is used to to initialize the corresponding component of the output object.<br>
    !>                              (**optional**, default = `repeat(" ", 255)`)
    !>
    !>  \return
    !>  `err`                   :   The output scalar object of type [err_type](@ref pm_err::err_type).<br>
    !>
    !>  \interface{err_type}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: err_type
    !>      type(err_type) :: err
    !>
    !>      err = err_type()
    !>      print *, err%occurred
    !>      print *, err%stat
    !>      print *, err%msg
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>
    !>  \example{err_type}
    !>  \include{lineno} example/pm_err/err_type/main.F90
    !>  \compilef{err_type}
    !>  \output{err_type}
    !>  \include{lineno} example/pm_err/err_type/main.out.F90
    !>
    !>  \final{err_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    type :: err_type
        logical(LK)                     :: iswarned = .false._LK    !<  The scalar `logical` of default kind \LK that is `.true.` if a runtime warning occurs.<br>
        logical(LK)                     :: occurred = .false._LK    !<  The scalar `logical` of default kind \LK that is `.true.` if a runtime error occurs.<br>
        integer(IK)                     :: stat     = -huge(0_IK)   !<  The scalar `integer` of default kind \IK to contain the error flag or status
                                                                    !<  code returned by the compiler or program upon encountering an error.<br>
        character(:, SK), allocatable   :: msg                      !<  The `allocatable` scalar `character` of default kind \SK containing the error message.<br>
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface err_type
    pure module function err_typer(occurred, stat, msg) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: err_typer
#endif
        use pm_kind, only: SK, IK
        logical(LK)         , intent(in)    , optional  :: occurred
        integer(IK)         , intent(in)    , optional  :: stat
        character(*, SK)    , intent(in)    , optional  :: msg
        type(err_type)                                  :: err
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !>  This is the derived type for generating objects to gracefully and
!    !>  verbosely handle runtime unexpected behavior in the ParaMonte library.<br>
!    !>
!    !>  \see
!    !>  [setNoted](@ref pm_err::setNoted)<br>
!    !>  [setWarned](@ref pm_err::setWarned)<br>
!    !>  [setAborted](@ref pm_err::setAborted)<br>
!    !>  [setAsserted](@ref pm_err::setAsserted)<br>
!    !>
!    !>  \example{message_type}
!    !>  \include{lineno} example/pm_err/err_type/main.F90
!    !>  \compilef{message_type}
!    !>  \output{message_type}
!    !>  \include{lineno} example/pm_err/err_type/main.out.F90
!    !>
!    !>  \final{message_type}
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
!    type :: message_type
!        logical(LK)                     :: occurred = .false._LK        !<  The scalar `logical` of default kind \LK that is `.true.` if a runtime error occurs.<br>
!        integer(IK)                     :: stat     = -huge(0_IK)       !<  The scalar `integer` of default kind \IK to contain the error flag or status
!                                                                        !<  code returned by the compiler or program upon encountering an error.<br>
!        character(:, SK), allocatable   :: msg                          !<  The scalar `character` of default kind \SK containing the error message.<br>
!       !logical(LK)                     , private   :: opened = .false._LK  !<  \private    The scalar `logical` of default kind \LK that is internally set to `.false.` once the object is constructed.<br>
!       !                                                                    !!              It is later set to `.true.` by all other methods of the class to allow closing of the file in the `final` subroutine.<br>
!       !                                                                    !!              This flag is needed to avoid unintentional closing of the opened file by the finalizer routine of the derived type.<br>
!        integer(IK)                     , private   :: tmsize = 0_IK        !<  \private    The scalar `integer` of default kind \IK representing the number of records to skip on top of the current record.<br>
!        integer(IK)                     , private   :: count = 1_IK         !<  \private    The scalar `integer` of default kind \IK representing the number of times to repeat the current record.<br>
!        integer(IK)                     , private   :: bmsize = 0_IK        !<  \private    The scalar `integer` of default kind \IK representing the number of records to skip on the bottom of the current record.<br>
!        integer(IK)                                 :: unit = output_unit   !<  \public     The scalar `integer` of default kind \IK representing the unit number of the file to which the record must be written.<br>
!        logical(LK)                     , private   :: sticky = .false._LK  !<  \public     The scalar `logical` of default kind \LK. If `.true.` the relevant properties of the current object can be later changed by its methods.<br>
!        logical(LK)                     , private   :: generic = .true._LK  !<  \public     The scalar `logical` of default kind \LK that is `.true.` only if the generic formatting of routine is used for display (unspecified by the user).<br>
!        character(3, SK)                , private   :: advance = SK_"YES"   !<  \private    The scalar `character` of default kind \SK representing the default IO `advance` mode for reading/writing actions.<br>
!        character(:, SK), allocatable   , private   :: format               !<  \private    The scalar `character` of default kind \SK representing the default IO format with which the records should be displayed.<br>
!        character(:, SK), allocatable   , private   :: delim                !<  \private    The scalar `character` of default kind \SK representing the default delimiter to be used for wrapping strings before displaying them.<br>
!        character(:, SK), allocatable   , private   :: file                 !<  \private    The scalar `character` of default kind \SK representing the file path to which the record must be written.<br>
!    !> \cond excluded
!    contains
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if     SK5_ENABLED
!        procedure, nopass                       ::          setNoted_D0_SK5
!        generic                                 :: note =>  setNoted_D0_SK5
!#endif
!#if     SK4_ENABLED
!        procedure, nopass                       ::          setNoted_D0_SK4
!        generic                                 :: note =>  setNoted_D0_SK4
!#endif
!#if     SK3_ENABLED
!        procedure, nopass                       ::          setNoted_D0_SK3
!        generic                                 :: note =>  setNoted_D0_SK3
!#endif
!#if     SK2_ENABLED
!        procedure, nopass                       ::          setNoted_D0_SK2
!        generic                                 :: note =>  setNoted_D0_SK2
!#endif
!#if     SK1_ENABLED
!        procedure, nopass                       ::          setNoted_D0_SK1
!        generic                                 :: note =>  setNoted_D0_SK1
!#endif
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if     SK5_ENABLED
!        procedure, nopass                       ::          setWarned_D0_SK5
!        generic                                 :: warn =>  setWarned_D0_SK5
!#endif
!#if     SK4_ENABLED
!        procedure, nopass                       ::          setWarned_D0_SK4
!        generic                                 :: warn =>  setWarned_D0_SK4
!#endif
!#if     SK3_ENABLED
!        procedure, nopass                       ::          setWarned_D0_SK3
!        generic                                 :: warn =>  setWarned_D0_SK3
!#endif
!#if     SK2_ENABLED
!        procedure, nopass                       ::          setWarned_D0_SK2
!        generic                                 :: warn =>  setWarned_D0_SK2
!#endif
!#if     SK1_ENABLED
!        procedure, nopass                       ::          setWarned_D0_SK1
!        generic                                 :: warn =>  setWarned_D0_SK1
!#endif
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if     SK5_ENABLED
!        procedure, nopass                       ::          stop_D0_SK5
!        generic                                 :: stop =>  stop_D0_SK5
!#endif
!#if     SK4_ENABLED
!        procedure, nopass                       ::          stop_D0_SK4
!        generic                                 :: stop =>  stop_D0_SK4
!#endif
!#if     SK3_ENABLED
!        procedure, nopass                       ::          stop_D0_SK3
!        generic                                 :: stop =>  stop_D0_SK3
!#endif
!#if     SK2_ENABLED
!        procedure, nopass                       ::          stop_D0_SK2
!        generic                                 :: stop =>  stop_D0_SK2
!#endif
!#if     SK1_ENABLED
!        procedure, nopass                       ::          stop_D0_SK1
!        generic                                 :: stop =>  stop_D0_SK1
!#endif
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a decorated string resulting from the concatenation of [getFile(FILE)](@ref pm_err::getFile)
    !>  and [getLine(LINE)](@ref pm_err::getLine) where `FILE` and `LINE` are replaced by the specified input source file and line of interest.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface are meant to facilitate tracing of origin of an error in error message when it occurs.<br>
    !>  Although the task performed by this generic interface is rather trivial, the goal is to offer consistent style over the entire library
    !>  for traceback information collection.<br>
    !>
    !>  \param[in]  file    :   The input scalar `character` of processor default kind `kind("")`,
    !>                          containing the path of the file within which an error or an interesting event has occurred.<br>
    !>                          The path of a source file is frequently returned by the compiler-defined macro `__FILE__`.<br>
    !>                          At least two compilers (Intel and GNU) define the macro `__FILE__` at compile time.<br>
    !>  \param[in]  line    :   The input scalar `integer` of processor default kind `kind(0)`,
    !>                          containing the line at which an error or an interesting event has occurred.<br>
    !>                          The number of a line within a source file is frequently returned by the compiler-defined macro `__LINE__`.<br>
    !>                          At least two compilers (Intel and GNU) define the macro `__LINE__` at compile time.<br>
    !>
    !>  \return
    !>  `str`               :   The output scalar `character` of default kind \SK containing a decoration of the input `file` and `line`.<br>
    !>
    !>  \interface{getFine}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK
    !>      use pm_err, only: getFine
    !>      character(:, SK), allocatable :: str
    !>
    !>      str = getFine(file, line)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \devnote
    !>  The name `Fine` results from merging `File` with `Line`.<br>
    !>
    !>  \see
    !>  [getFile](@ref pm_err::getFile)<br>
    !>  [getFine](@ref pm_err::getFine)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_err/getFine/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_err/getFine/main.out.F90
    !>
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final{getFine}
    !>
    !>  \todo
    !>  \pmed Adding an example usage to this interface can be helpful.<br>
    !>
    !>  \author
    !>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface getFine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getFine_SK5(file, line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFine_SK5
#endif
        use pm_kind, only: SK
        character(*)    , intent(in)            :: file
        integer         , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

#if SK4_ENABLED
    pure module function getFine_SK4(file, line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFine_SK4
#endif
        use pm_kind, only: SK
        character(*)    , intent(in)            :: file
        integer         , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

#if SK3_ENABLED
    pure module function getFine_SK3(file, line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFine_SK3
#endif
        use pm_kind, only: SK
        character(*)    , intent(in)            :: file
        integer         , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

#if SK2_ENABLED
    pure module function getFine_SK2(file, line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFine_SK2
#endif
        use pm_kind, only: SK
        character(*)    , intent(in)            :: file
        integer         , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

#if SK1_ENABLED
    pure module function getFine_SK1(file, line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFine_SK1
#endif
        use pm_kind, only: SK
        character(*)    , intent(in)            :: file
        integer         , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a decorated string as `"@file(FILE)"` where `FILE` is replaced by the input source file of interest.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface are meant to facilitate tracing of origin of an error in error message when it occurs.<br>
    !>  Although the task performed by this generic interface is rather trivial, the goal is to offer consistent style over the entire library
    !>  for traceback information collection.<br>
    !>
    !>  \param[in]  file    :   The input scalar `character` of kind \SKALL, containing the path of the file within which an error or an interesting event has occurred.<br>
    !>                          The path of a source file is frequently returned by the compiler-defined macro `__FILE__`.<br>
    !>                          At least two compilers (Intel and GNU) define the macro `__FILE__` at compile time.<br>
    !>
    !>  \return
    !>  `str`               :   The output scalar `character` of default kind \SK containing a decoration of the input `file`.<br>
    !>
    !>  \interface{getFile}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK
    !>      use pm_err, only: getFile
    !>      character(:, SK), allocatable :: str
    !>
    !>      str = getFile(file)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [getFile](@ref pm_err::getFile)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_err/getFile/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_err/getFile/main.out.F90
    !>
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final{getFile}
    !>
    !>  \todo
    !>  \pmed Adding an example usage to this interface can be helpful.<br>
    !>
    !>  \author
    !>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface getFile

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getFile_SK5(file) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFile_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)            :: file
        character(len(file) + 7, SK)            :: str
    end function
#endif

#if SK4_ENABLED
    pure module function getFile_SK4(file) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFile_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)            :: file
        character(len(file) + 7, SK)            :: str
    end function
#endif

#if SK3_ENABLED
    pure module function getFile_SK3(file) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFile_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)            :: file
        character(len(file) + 7, SK)            :: str
    end function
#endif

#if SK2_ENABLED
    pure module function getFile_SK2(file) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFile_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)            :: file
        character(len(file) + 7, SK)            :: str
    end function
#endif

#if SK1_ENABLED
    pure module function getFile_SK1(file) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFile_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)            :: file
        character(len(file) + 7, SK)            :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a decorated string as `"@line(LINE)"` where `LINE` is replaced by the input source line of interest.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface are meant to facilitate tracing of origin of an error in error message when it occurs.<br>
    !>  Although the task performed by this generic interface is rather trivial, the goal is to offer consistent style over the entire library
    !>  for traceback information collection.<br>
    !>
    !>  \param[in]  line    :   The input scalar `integer` of kind \IKALL, containing the line at which an error or an interesting event has occurred.<br>
    !>                          The number of a line within a source file is frequently returned by the compiler-defined macro `__LINE__`.<br>
    !>                          At least two compilers (Intel and GNU) define the macro `__LINE__` at compile time.<br>
    !>
    !>  \return
    !>  `str`               :   The output scalar `character` of default kind \SK containing a decoration of the input `line`.<br>
    !>
    !>  \interface{getLine}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK
    !>      use pm_err, only: getLine
    !>      character(:, SK), allocatable :: str
    !>
    !>      str = getLine(line)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [getFile](@ref pm_err::getFile)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_err/getLine/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_err/getLine/main.out.F90
    !>
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final{getLine}
    !>
    !>  \todo
    !>  \pmed Adding an example usage to this interface can be helpful.<br>
    !>
    !>  \author
    !>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface getLine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function getLine_IK5(line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLine_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)    , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

#if IK4_ENABLED
    pure module function getLine_IK4(line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLine_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)    , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

#if IK3_ENABLED
    pure module function getLine_IK3(line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLine_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)    , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

#if IK2_ENABLED
    pure module function getLine_IK2(line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLine_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)    , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

#if IK1_ENABLED
    pure module function getLine_IK1(line) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLine_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)    , intent(in)            :: line
        character(:, SK), allocatable           :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Verify the input `assertion` holds and if it does not, print the (optional) input message on `stdout` and halt
    !>  the program via `error stop` or (optionally) return the program control to the caller routine, if requested.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface are meant to be used primarily within the ParaMonte library
    !>  routines to verify the runtime conditions that must hold when the library procedures are called by the end-user.<br>
    !>
    !>  \param[in]  assertion   :   The input scalar `logical` of kind \LKALL. If `.true.` it means that the desired condition holds.<br>
    !>                              Otherwise, the program halts by calling `error stop` unless `renabled = .true.`.<br>
    !>  \param[in]  msg         :   The input scalar `character` of default kind \SK containing a descriptive message of the condition
    !>                              that must hold that leads to `assertion = .true.`.<br>
    !>                              If present, the message if printed out on the output unit `output_unit`
    !>                              before halting the program (or returning the control to the caller).<br>
    !>                              (**optional**, default = `""`)
    !>  \param[in]  renabled    :   The input scalar `logical` of default kind \LK. If `.true.`, the program-control will be
    !>                              returned to the calling routine instead of terminating the program via `error stop`.<br>
    !>                              This capability is particularly needed during the testing of the various components
    !>                              of the ParaMonte library, when fatal error must be handled gracefully.<br>
    !>                              (**optional**, default = `.false.`)
    !>
    !>  \interface{setAsserted}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: setAsserted
    !>
    !>      call setAsserted(assertion, msg = msg, renabled = renabled)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final{setAsserted}
    !>
    !>  \todo
    !>  \pmed Adding an example usage to this interface can be helpful.<br>
    !>
    !>  \author
    !>  \AmirShahmoradi, 3:43 AM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface setAsserted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setAsserted_LK5(assertion, msg, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsserted_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)    , intent(in)            :: assertion
        character(*, SK), intent(in), optional  :: msg
        logical(LK)     , intent(in), optional  :: renabled
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setAsserted_LK4(assertion, msg, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsserted_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)    , intent(in)            :: assertion
        character(*, SK), intent(in), optional  :: msg
        logical(LK)     , intent(in), optional  :: renabled
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setAsserted_LK3(assertion, msg, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsserted_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)    , intent(in)            :: assertion
        character(*, SK), intent(in), optional  :: msg
        logical(LK)     , intent(in), optional  :: renabled
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setAsserted_LK2(assertion, msg, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsserted_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)    , intent(in)            :: assertion
        character(*, SK), intent(in), optional  :: msg
        logical(LK)     , intent(in), optional  :: renabled
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setAsserted_LK1(assertion, msg, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsserted_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)    , intent(in)            :: assertion
        character(*, SK), intent(in), optional  :: msg
        logical(LK)     , intent(in), optional  :: renabled
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the base derived type for constructing subclasses that contain the specifications of the generic message interfaces.<br>
    !>
    !>  \details
    !>  The generic message interfaces of this module include,
    !>  <ol>
    !>      <li>    [setMarked](@ref pm_err::setMarked)<br>
    !>      <li>    [setNoted](@ref pm_err::setNoted)<br>
    !>      <li>    [setWarned](@ref pm_err::setWarned)<br>
    !>      <li>    [setAborted](@ref pm_err::setAborted)<br>
    !>  </ol>
    !>  The subclasses of this derived contain a dynamic method `show()` that acts as a convenience wrapper around the corresponding static generic interfaces mentioned above.<br>
    !>  The main use of these classes is in situations where arbitrary user-specified messages must be delivered to an output unit **repeatedly with a uniform identical style**.<br>
    !>  In such cases, the settings for the message can be specified once via the default constructor of this derived type and used subsequently via its wrapper method [show](@ref pm_err::show).<br>
    !>  This base derived type minimally defines the following type components listed below as input arguments to the default constructor of the type.<br>
    !>
    !>  \param[in]  prefix      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**, see the corresponding default in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  indent      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**, see the corresponding default in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  break       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**, see the corresponding default in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  newline     :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**, see the corresponding default in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  width       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**, see the corresponding default in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  maxwidth    :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**, see the corresponding default in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  tmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines preceding the message to be displayed.<br>
    !>                              (**optional**, see the corresponding default in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  bmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines succeeding the message to be displayed.<br>
    !>                              (**optional**, see the corresponding default in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  unit        :   The input scalar `integer` of default kind \IK containing the output file unit.<br>
    !>                              (**optional**, default = `output_unit` taken from the intrinsic module `iso_fortran_env`.)
    !>  \param[in]  sticky      :   The input scalar of type `logical` of default kind \LK.<br>
    !>                              If `.true.`, the global properties (i.e., components) of the object of class [message_type](@ref pm_err::message_type)
    !>                              can be overridden by the corresponding arguments passed to the dynamic methods of the object.<br>
    !>                              (**optional**, default = `.false.`)
    !>
    !>  \return
    !>  `message`               :   An object of class [message_type](@ref pm_err::message_type) containing the attributes and wrapper method for outputting a message.<br>
    !>
    !>  \interface{message_type}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: message_type, mark_type, note_type, warn_type, stop_type
    !>      class(message_type), allocatable :: message
    !>
    !>      message = message_type()
    !>      message = mark_type()
    !>      message = note_type()
    !>      message = warn_type()
    !>      message = stop_type()
    !>
    !>  \endcode
    !>
    !>  \note
    !>  This base derived type is meant to be primary provide the scaffolding required to define the subclasses,
    !>  <ol>
    !>      <li>    [mark_type](@ref pm_err::mark_type)<br>
    !>      <li>    [note_type](@ref pm_err::note_type)<br>
    !>      <li>    [warn_type](@ref pm_err::warn_type)<br>
    !>      <li>    [stop_type](@ref pm_err::stop_type)<br>
    !>  </ol>
    !>
    !>  \see
    !>  [show](@ref pm_err::show)<br>
    !>  [mark_type](@ref pm_err::mark_type)<br>
    !>  [note_type](@ref pm_err::note_type)<br>
    !>  [warn_type](@ref pm_err::warn_type)<br>
    !>  [stop_type](@ref pm_err::stop_type)<br>
    !>  [message_type](@ref pm_err::message_type)<br>
    !>  [display_type](@ref pm_io::display_type)<br>
    !>  [getStrWrapped](@ref pm_str::getStrWrapped)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>
    !>  \test
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final{message_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: message_type
        character(:, SK), allocatable   , public    :: prefix                       !<  \public  See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
        character(:, SK), allocatable   , public    :: indent                       !<  \public  See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
        character(:, SK), allocatable   , public    :: break                        !<  \public  See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
        character(:, SK), allocatable   , public    :: newline                      !<  \public  See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
        integer(IK)     , allocatable   , public    :: width                        !<  \public  See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
        integer(IK)     , allocatable   , public    :: maxwidth                     !<  \public  See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
        integer(IK)     , allocatable   , public    :: tmsize                       !<  \public  See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
        integer(IK)     , allocatable   , public    :: bmsize                       !<  \public  See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
        integer(IK)     , allocatable   , public    :: unit                         !<  \public  See the corresponding definition in the documentation of [message_type](@ref pm_err::message_type).<br>
        logical(LK)                     , public    :: sticky = .false._LK          !<  \public  See the corresponding argument in the documentation of [message_type](@ref pm_err::message_type).<br>
       !logical(LK)                     , private   :: uninit = .true._LK           !<  \private The scalar logical that is `.true.` if and only if the type instance has not been initialized with the constructor.<br>
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for constructing objects that contain the specifications of the generic interface [setMarked](@ref pm_err::setMarked)
    !>  along with a dynamic method [show](@ref pm_err_show) that acts as a convenience wrapper around the generic interface [setMarked](@ref pm_err::setMarked).<br>
    !>
    !>  \details
    !>  See the documentation of the parent type [message_type](@ref pm_err::message_type) for more details and inherited derived type components.<br>
    !>
    !>  \return
    !>  `mark`                  :   An object of type [mark_type](@ref pm_err::mark_type) containing the attributes and wrapper method for outputting a remark.<br>
    !>
    !>  \interface{mark_type}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: mark_type
    !>      type(mark_type) :: spec
    !>
    !>      spec = mark_type( prefix = prefix       &
    !>                      , indent = indent       &
    !>                      , break = break         &
    !>                      , newline = newline     &
    !>                      , width = width         &
    !>                      , maxwidth = maxwidth   &
    !>                      , tmsize = tmsize       &
    !>                      , bmsize = bmsize       &
    !>                      , unit = unit           &
    !>                      , sticky = sticky       &
    !>                      )
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  The allocation status of the components is used to signify the `optional` argument presence within the methods.<br>
    !>
    !>  \see
    !>  [mark_type](@ref pm_err::mark_type)<br>
    !>  [note_type](@ref pm_err::note_type)<br>
    !>  [warn_type](@ref pm_err::warn_type)<br>
    !>  [stop_type](@ref pm_err::stop_type)<br>
    !>  [display_type](@ref pm_io::display_type)<br>
    !>  [message_type](@ref pm_err::message_type)<br>
    !>  [getStrWrapped](@ref pm_str::getStrWrapped)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>
    !>  \example{mark_type}
    !>  \include{lineno} example/pm_err/mark_type/main.F90
    !>  \compilef{mark_type}
    !>  \output{mark_type}
    !>  \include{lineno} example/pm_err/mark_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11-13}
    !>  \desc
    !>  \gfortran cannot properly construct the `allocatable` scalar non-`character` components
    !>  of objects of type [warn_type](@ref pm_err::warn_type) using the default constructor.<br>
    !>  For example, when `unit` is set via the default constructor, the program behaves as if the `unit`
    !>  component of the object is allocated but unset, yielding a `segmentation fault` error.<br>
    !>  \remedy
    !>  For now, the custom constructor bypasses \gfortran bug.<br>
    !>
    !>  \final{mark_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(message_type) :: mark_type
    contains
        !>  \cond excluded
#if     SK5_ENABLED
        procedure, pass ::          setMarkedMethod_D0_SK5
        generic         :: show =>  setMarkedMethod_D0_SK5
#endif
#if     SK4_ENABLED
        procedure, pass ::          setMarkedMethod_D0_SK4
        generic         :: show =>  setMarkedMethod_D0_SK4
#endif
#if     SK3_ENABLED
        procedure, pass ::          setMarkedMethod_D0_SK3
        generic         :: show =>  setMarkedMethod_D0_SK3
#endif
#if     SK2_ENABLED
        procedure, pass ::          setMarkedMethod_D0_SK2
        generic         :: show =>  setMarkedMethod_D0_SK2
#endif
#if     SK1_ENABLED
        procedure, pass ::          setMarkedMethod_D0_SK1
        generic         :: show =>  setMarkedMethod_D0_SK1
#endif
        !>  \endcond excluded
    end type

    !>  \cond excluded
    interface mark_type
        module procedure :: mark_typer
    end interface
    !>  \endcond excluded

    interface
    !>  \brief
    !>  Generate and return an object of type [mark_type](@ref pm_err::mark_type) with the user-specified input attributes.
    !>
    !>  \details
    !>  This generic interface serves as the custom constructor for objects of type [mark_type](@ref pm_err::mark_type).<br>
    !>  See the documentation of [mark_type](@ref pm_err::mark_type) for the documentation of the input arguments.<br>
    !>  This custom constructor exists because of a \gfortran bug in the implementation of the default constructor.<br>
    !>  See the bug description below.<br>
    !>
    !>  \see
    !>  [mark_type](@ref pm_err::mark_type)<br>
    !>
    !>  \test
    !>  [test_pm_err](@ref test_pm_err)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11-13}
    !>  \desc
    !>  \gfortran cannot properly construct the `allocatable` scalar non-`character` components
    !>  of objects of type [mark_type](@ref pm_err::mark_type) using the default constructor.<br>
    !>  For example, when `unit` is set via the default constructor, the program behaves as if the `unit`
    !>  component of the object is allocated but unset, yielding a `segmentation fault` error.<br>
    !>  \remedy
    !>  For now, the custom constructor bypasses \gfortran bug.<br>
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{mark_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    pure module function mark_typer(prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: mark_typer
#endif
        character(*, SK), intent(in), optional  :: prefix
        character(*, SK), intent(in), optional  :: indent
        character(*, SK), intent(in), optional  :: break
        character(*, SK), intent(in), optional  :: newline
        integer(IK)     , intent(in), optional  :: width
        integer(IK)     , intent(in), optional  :: maxwidth
        integer(IK)     , intent(in), optional  :: tmsize
        integer(IK)     , intent(in), optional  :: bmsize
        integer(IK)     , intent(in), optional  :: unit
        logical(LK)     , intent(in), optional  :: sticky
        type(mark_type)                         :: self
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for constructing objects that contain the specifications of the generic interface [setNoted](@ref pm_err::setNoted)
    !>  along with a dynamic method [show](@ref pm_err_show) that acts as a convenience wrapper around the generic interface [setNoted](@ref pm_err::setNoted).<br>
    !>
    !>  \details
    !>  See the documentation of the parent type [message_type](@ref pm_err::message_type) for more details and inherited derived type components.<br>
    !>
    !>  \return
    !>  `note`                  :   An object of type [note_type](@ref pm_err::note_type) containing the attributes and wrapper method for outputting a notification message.<br>
    !>
    !>  \interface{note_type}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: note_type
    !>      type(note_type) :: spec
    !>
    !>      spec = note_type( prefix = prefix       &
    !>                      , indent = indent       &
    !>                      , break = break         &
    !>                      , newline = newline     &
    !>                      , width = width         &
    !>                      , maxwidth = maxwidth   &
    !>                      , tmsize = tmsize       &
    !>                      , bmsize = bmsize       &
    !>                      , unit = unit           &
    !>                      , sticky = sticky       &
    !>                      )
    !>
    !>  \endcode
    !>
    !>  \devnote
    !>  The allocation status of the components is used to signify the `optional` argument presence within the methods.<br>
    !>
    !>  \see
    !>  [mark_type](@ref pm_err::mark_type)<br>
    !>  [note_type](@ref pm_err::note_type)<br>
    !>  [warn_type](@ref pm_err::warn_type)<br>
    !>  [stop_type](@ref pm_err::stop_type)<br>
    !>  [display_type](@ref pm_io::display_type)<br>
    !>  [message_type](@ref pm_err::message_type)<br>
    !>  [getStrWrapped](@ref pm_str::getStrWrapped)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>
    !>  \example{note_type}
    !>  \include{lineno} example/pm_err/note_type/main.F90
    !>  \compilef{note_type}
    !>  \output{note_type}
    !>  \include{lineno} example/pm_err/note_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11-13}
    !>  \desc
    !>  \gfortran cannot properly construct the `allocatable` scalar non-`character` components
    !>  of objects of type [warn_type](@ref pm_err::warn_type) using the default constructor.<br>
    !>  For example, when `unit` is set via the default constructor, the program behaves as if the `unit`
    !>  component of the object is allocated but unset, yielding a `segmentation fault` error.<br>
    !>  \remedy
    !>  For now, the custom constructor bypasses \gfortran bug.<br>
    !>
    !>  \final{note_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(message_type) :: note_type
    contains
        !>  \cond excluded
#if     SK5_ENABLED
        procedure, pass ::          setNotedMethod_D0_SK5
        generic         :: show =>  setNotedMethod_D0_SK5
#endif
#if     SK4_ENABLED
        procedure, pass ::          setNotedMethod_D0_SK4
        generic         :: show =>  setNotedMethod_D0_SK4
#endif
#if     SK3_ENABLED
        procedure, pass ::          setNotedMethod_D0_SK3
        generic         :: show =>  setNotedMethod_D0_SK3
#endif
#if     SK2_ENABLED
        procedure, pass ::          setNotedMethod_D0_SK2
        generic         :: show =>  setNotedMethod_D0_SK2
#endif
#if     SK1_ENABLED
        procedure, pass ::          setNotedMethod_D0_SK1
        generic         :: show =>  setNotedMethod_D0_SK1
#endif
        !>  \endcond excluded
    end type

    !>  \cond excluded
    interface note_type
        module procedure :: note_typer
    end interface
    !>  \endcond excluded

    interface
    !>  \brief
    !>  Generate and return an object of type [note_type](@ref pm_err::note_type) with the user-specified input attributes.
    !>
    !>  \details
    !>  This generic interface serves as the custom constructor for objects of type [note_type](@ref pm_err::note_type).<br>
    !>  See the documentation of [note_type](@ref pm_err::note_type) for the documentation of the input arguments.<br>
    !>  This custom constructor exists because of a \gfortran bug in the implementation of the default constructor.<br>
    !>  See the bug description below.<br>
    !>
    !>  \see
    !>  [note_type](@ref pm_err::note_type)<br>
    !>
    !>  \test
    !>  [test_pm_err](@ref test_pm_err)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11-13}
    !>  \desc
    !>  \gfortran cannot properly construct the `allocatable` scalar non-`character` components
    !>  of objects of type [note_type](@ref pm_err::note_type) using the default constructor.<br>
    !>  For example, when `unit` is set via the default constructor, the program behaves as if the `unit`
    !>  component of the object is allocated but unset, yielding a `segmentation fault` error.<br>
    !>  \remedy
    !>  For now, the custom constructor bypasses \gfortran bug.<br>
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{note_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    pure module function note_typer(prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: note_typer
#endif
        character(*, SK), intent(in), optional  :: prefix
        character(*, SK), intent(in), optional  :: indent
        character(*, SK), intent(in), optional  :: break
        character(*, SK), intent(in), optional  :: newline
        integer(IK)     , intent(in), optional  :: width
        integer(IK)     , intent(in), optional  :: maxwidth
        integer(IK)     , intent(in), optional  :: tmsize
        integer(IK)     , intent(in), optional  :: bmsize
        integer(IK)     , intent(in), optional  :: unit
        logical(LK)     , intent(in), optional  :: sticky
        type(note_type)                         :: self
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for constructing objects that contain the specifications of the generic interface [setWarned](@ref pm_err::setWarned)
    !>  along with a dynamic method [show](@ref pm_err_show) that acts as a convenience wrapper around the generic interface [setWarned](@ref pm_err::setWarned).<br>
    !>
    !>  \details
    !>  See the documentation of the parent type [message_type](@ref pm_err::message_type) for more details and inherited derived type components.<br>
    !>
    !>  \return
    !>  `warn`                  :   An object of type [warn_type](@ref pm_err::warn_type) containing the attributes and wrapper method for outputting a warning message.<br>
    !>
    !>  \interface{warn_type}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: warn_type
    !>      type(warn_type) :: spec
    !>
    !>      spec = warn_type( prefix = prefix       &
    !>                      , indent = indent       &
    !>                      , break = break         &
    !>                      , newline = newline     &
    !>                      , width = width         &
    !>                      , maxwidth = maxwidth   &
    !>                      , tmsize = tmsize       &
    !>                      , bmsize = bmsize       &
    !>                      , unit = unit           &
    !>                      , sticky = sticky       &
    !>                      )
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [mark_type](@ref pm_err::mark_type)<br>
    !>  [note_type](@ref pm_err::note_type)<br>
    !>  [warn_type](@ref pm_err::warn_type)<br>
    !>  [stop_type](@ref pm_err::stop_type)<br>
    !>  [display_type](@ref pm_io::display_type)<br>
    !>  [message_type](@ref pm_err::message_type)<br>
    !>  [getStrWrapped](@ref pm_str::getStrWrapped)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>
    !>  \example{warn_type}
    !>  \include{lineno} example/pm_err/warn_type/main.F90
    !>  \compilef{warn_type}
    !>  \output{warn_type}
    !>  \include{lineno} example/pm_err/warn_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11-13}
    !>  \desc
    !>  \gfortran cannot properly construct the `allocatable` scalar non-`character` components
    !>  of objects of type [warn_type](@ref pm_err::warn_type) using the default constructor.<br>
    !>  For example, when `unit` is set via the default constructor, the program behaves as if the `unit`
    !>  component of the object is allocated but unset, yielding a `segmentation fault` error.<br>
    !>  \remedy
    !>  For now, the custom constructor bypasses \gfortran bug.<br>
    !>
    !>  \final{warn_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(message_type) :: warn_type
    contains
        !>  \cond excluded
#if     SK5_ENABLED
        procedure, pass ::          setWarnedMethod_D0_SK5
        generic         :: show =>  setWarnedMethod_D0_SK5
#endif
#if     SK4_ENABLED
        procedure, pass ::          setWarnedMethod_D0_SK4
        generic         :: show =>  setWarnedMethod_D0_SK4
#endif
#if     SK3_ENABLED
        procedure, pass ::          setWarnedMethod_D0_SK3
        generic         :: show =>  setWarnedMethod_D0_SK3
#endif
#if     SK2_ENABLED
        procedure, pass ::          setWarnedMethod_D0_SK2
        generic         :: show =>  setWarnedMethod_D0_SK2
#endif
#if     SK1_ENABLED
        procedure, pass ::          setWarnedMethod_D0_SK1
        generic         :: show =>  setWarnedMethod_D0_SK1
#endif
        !>  \endcond excluded
    end type

    !>  \cond excluded
    interface warn_type
        module procedure :: warn_typer
    end interface
    !>  \endcond excluded

    interface
    !>  \brief
    !>  Generate and return an object of type [warn_type](@ref pm_err::warn_type) with the user-specified input attributes.
    !>
    !>  \details
    !>  This generic interface serves as the custom constructor for objects of type [warn_type](@ref pm_err::warn_type).<br>
    !>  See the documentation of [warn_type](@ref pm_err::warn_type) for the documentation of the input arguments.<br>
    !>  This custom constructor exists because of a \gfortran bug in the implementation of the default constructor.<br>
    !>  See the bug description below.<br>
    !>
    !>  \see
    !>  [warn_type](@ref pm_err::warn_type)<br>
    !>
    !>  \test
    !>  [test_pm_err](@ref test_pm_err)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11-13}
    !>  \desc
    !>  \gfortran cannot properly construct the `allocatable` scalar non-`character` components
    !>  of objects of type [warn_type](@ref pm_err::warn_type) using the default constructor.<br>
    !>  For example, when `unit` is set via the default constructor, the program behaves as if the `unit`
    !>  component of the object is allocated but unset, yielding a `segmentation fault` error.<br>
    !>  \remedy
    !>  For now, the custom constructor bypasses \gfortran bug.<br>
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{warn_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    pure module function warn_typer(prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: warn_typer
#endif
        character(*, SK), intent(in), optional  :: prefix
        character(*, SK), intent(in), optional  :: indent
        character(*, SK), intent(in), optional  :: break
        character(*, SK), intent(in), optional  :: newline
        integer(IK)     , intent(in), optional  :: width
        integer(IK)     , intent(in), optional  :: maxwidth
        integer(IK)     , intent(in), optional  :: tmsize
        integer(IK)     , intent(in), optional  :: bmsize
        integer(IK)     , intent(in), optional  :: unit
        logical(LK)     , intent(in), optional  :: sticky
        type(warn_type)                         :: self
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for constructing objects that contain the specifications of the generic interface [setAborted](@ref pm_err::setAborted)
    !>  along with a dynamic method [show](@ref pm_err_show) that acts as a convenience wrapper around the generic interface [setAborted](@ref pm_err::setAborted).<br>
    !>
    !>  \details
    !>  See the documentation of the parent type [message_type](@ref pm_err::message_type) for more details and inherited derived type components.<br>
    !>
    !>  \param[in]  help        :   The input scalar of type `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              For details, see the description of the corresponding argument of the [setAborted](@ref pm_err::setAborted) dynamic method of the derived type.<br>
    !>                              (**optional**. The default value is set as described in the documentation of the corresponding argument of [setAborted](@ref pm_err::setAborted).)
    !>  \param[in]  stat        :   The input scalar `integer` of default kind \IK.<br>
    !>                              For details, see the description of the corresponding argument of the [setAborted](@ref pm_err::setAborted) dynamic method of the derived type.<br>
    !>                              (**optional**. The default value is set as described in the documentation of the corresponding argument of [setAborted](@ref pm_err::setAborted).)
    !>  \param[in]  renabled    :   The input scalar `logical` of default kind \LK.<br>
    !>                              For details, see the description of the corresponding argument of the [setAborted](@ref pm_err::setAborted) dynamic method of the derived type.<br>
    !>                              (**optional**. The default value is set as described in the documentation of the corresponding argument of [setAborted](@ref pm_err::setAborted).)
    !>
    !>  \return
    !>  `stop`                  :   An object of type [stop_type](@ref pm_err::stop_type) containing the attributes and wrapper method for aborting a program with an abortion message.<br>
    !>
    !>  \interface{stop_type}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: stop_type
    !>      type(stop_type) :: stop
    !>
    !>      stop = stop_type( prefix = prefix       &
    !>                      , indent = indent       &
    !>                      , break = break         &
    !>                      , newline = newline     &
    !>                      , width = width         &
    !>                      , maxwidth = maxwidth   &
    !>                      , tmsize = tmsize       &
    !>                      , bmsize = bmsize       &
    !>                      , unit = unit           &
    !>                      , sticky = sticky       &
    !>                      , stat = stat           &
    !>                      , help = help           &
    !>                      , renabled = renabled   &
    !>                      )
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [mark_type](@ref pm_err::mark_type)<br>
    !>  [note_type](@ref pm_err::note_type)<br>
    !>  [warn_type](@ref pm_err::warn_type)<br>
    !>  [stop_type](@ref pm_err::stop_type)<br>
    !>  [display_type](@ref pm_io::display_type)<br>
    !>  [message_type](@ref pm_err::message_type)<br>
    !>  [getStrWrapped](@ref pm_str::getStrWrapped)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>
    !>  \example{stop_type}
    !>  \include{lineno} example/pm_err/stop_type/main.F90
    !>  \compilef{stop_type}
    !>  \output{stop_type}
    !>  \include{lineno} example/pm_err/stop_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11-13}
    !>  \desc
    !>  \gfortran cannot properly construct the `allocatable` scalar non-`character` components
    !>  of objects of type [stop_type](@ref pm_err::stop_type) using the default constructor.<br>
    !>  For example, when `unit` is set via the default constructor, the program behaves as if the `unit`
    !>  component of the object is allocated but unset, yielding a `segmentation fault` error.<br>
    !>  \remedy
    !>  For now, the custom constructor bypasses \gfortran bug.<br>
    !>
    !>  \final{stop_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(message_type) :: stop_type
        character(:, SK)    , allocatable, public   :: help     !<  \public See the corresponding definition in the documentation of the default constructor of the derived type [stop_type](@ref pm_err::stop_type).<br>
        integer(IK)         , allocatable, public   :: stat     !<  \public See the corresponding definition in the documentation of the default constructor of the derived type [stop_type](@ref pm_err::stop_type).<br>
        logical(LK)         , allocatable, public   :: renabled !<  \public See the corresponding definition in the documentation of the default constructor of the derived type [stop_type](@ref pm_err::stop_type).<br>
    contains
        !>  \cond excluded
#if     SK5_ENABLED
        procedure, pass ::          setAbortedMethod_D0_SK5
        generic         :: show =>  setAbortedMethod_D0_SK5
#endif
#if     SK4_ENABLED
        procedure, pass ::          setAbortedMethod_D0_SK4
        generic         :: show =>  setAbortedMethod_D0_SK4
#endif
#if     SK3_ENABLED
        procedure, pass ::          setAbortedMethod_D0_SK3
        generic         :: show =>  setAbortedMethod_D0_SK3
#endif
#if     SK2_ENABLED
        procedure, pass ::          setAbortedMethod_D0_SK2
        generic         :: show =>  setAbortedMethod_D0_SK2
#endif
#if     SK1_ENABLED
        procedure, pass ::          setAbortedMethod_D0_SK1
        generic         :: show =>  setAbortedMethod_D0_SK1
#endif
        !>  \endcond excluded
    end type

    !>  \cond excluded
    interface stop_type
        module procedure :: stop_typer
    end interface
    !>  \endcond excluded

    interface
    !>  \brief
    !>  Generate and return an object of type [stop_type](@ref pm_err::stop_type) with the user-specified input attributes.
    !>
    !>  \details
    !>  This generic interface serves as the custom constructor for objects of type [stop_type](@ref pm_err::stop_type).<br>
    !>  See the documentation of [stop_type](@ref pm_err::stop_type) for the documentation of the input arguments.<br>
    !>  This custom constructor exists because of a \gfortran bug in the implementation of the default constructor.<br>
    !>  See the bug description below.<br>
    !>
    !>  \see
    !>  [stop_type](@ref pm_err::stop_type)<br>
    !>
    !>  \test
    !>  [test_pm_err](@ref test_pm_err)<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11-13}
    !>  \desc
    !>  \gfortran cannot properly construct the `allocatable` scalar non-`character` components
    !>  of objects of type [stop_type](@ref pm_err::stop_type) using the default constructor.<br>
    !>  For example, when `unit` is set via the default constructor, the program behaves as if the `unit`
    !>  component of the object is allocated but unset, yielding a `segmentation fault` error.<br>
    !>  \remedy
    !>  For now, the custom constructor bypasses \gfortran bug.<br>
    !>
    !>  \todo
    !>  \plow A test should be implemented for arrays of size that can be represented *only* by an \IKD integer.
    !>
    !>  \final{stop_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    pure module function stop_typer(prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky, help, stat, renabled) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: stop_typer
#endif
        character(*, SK), intent(in), optional  :: prefix
        character(*, SK), intent(in), optional  :: indent
        character(*, SK), intent(in), optional  :: break
        character(*, SK), intent(in), optional  :: newline
        integer(IK)     , intent(in), optional  :: width
        integer(IK)     , intent(in), optional  :: maxwidth
        integer(IK)     , intent(in), optional  :: tmsize
        integer(IK)     , intent(in), optional  :: bmsize
        integer(IK)     , intent(in), optional  :: unit
        logical(LK)     , intent(in), optional  :: sticky
        character(*, SK), intent(in), optional  :: help
        integer(IK)     , intent(in), optional  :: stat
        logical(LK)     , intent(in), optional  :: renabled
        type(stop_type)                         :: self
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_err_show
    !>  Write the input string `msg` in the format of a notification/warning/abortion message to the output.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface are dynamic member methods of objects
    !>  whose class is dictated by the implicitly-passed `class` argument `self` as described below.<br>
    !>  to obtain a wrapped text within the specified `width` and `maxwidth` where the input `prefix` is suffixed with ` - NOTE: `
    !>  and then prepended to the beginning of each wrapped line and then printed on the display.<br>
    !>
    !>  See the documentations of,
    !>  <ol>
    !>      <li>    [mark_type](@ref pm_err::mark_type)<br>
    !>      <li>    [note_type](@ref pm_err::note_type)<br>
    !>      <li>    [warn_type](@ref pm_err::warn_type)<br>
    !>      <li>    [stop_type](@ref pm_err::stop_type)<br>
    !>  </ol>
    !>  for example usage.<br>
    !>
    !>  \param[inout]   self        :   The input/output **implicitly-passed** scalar that can be of,<br>
    !>                                  <ol>
    !>                                      <li>    class [mark_type](@ref pm_err::mark_type)<br>
    !>                                      <li>    class [note_type](@ref pm_err::note_type)<br>
    !>                                      <li>    class [warn_type](@ref pm_err::warn_type)<br>
    !>                                      <li>    class [stop_type](@ref pm_err::stop_type)<br>
    !>                                  </ol>
    !>  \param[in]      msg         :   The input scalar `character` of kind \SKALL representing the text to be wrapped and displayed.<br>
    !>  \param[in]      prefix      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                                  (**optional**)
    !>  \param[in]      indent      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                                  (**optional**)
    !>  \param[in]      break       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                                  (**optional**)
    !>  \param[in]      newline     :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                                  (**optional**)
    !>  \param[in]      width       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                                  (**optional**)
    !>  \param[in]      maxwidth    :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                                  (**optional**)
    !>  \param[in]      tmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines preceding the text in the display.<br>
    !>                                  (**optional**, default = `1_IK`)
    !>  \param[in]      bmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines succeeding the text in the display.<br>
    !>                                  (**optional**, default = `0_IK`)
    !>  \param[in]      unit        :   The input scalar `integer` of default kind \IK containing the output file unit.<br>
    !>                                  If the input `unit` is different from the standard output (`output_unit`), all messages will also be printed on stdout recursively.<br>
    !>                                  (**optional**, default = `output_unit` taken from the intrinsic module `iso_fortran_env`.)
    !>  \param[in]      sticky      :   The input scalar of type `logical` of default kind \LK.<br>
    !>                                  If `.true.`, the global properties (i.e., components) of the object of class [message_type](@ref pm_err::message_type)
    !>                                  can be overridden by the corresponding arguments passed to the dynamic methods of the object.<br>
    !>                                  (**optional**, default = `.false.`)
    !>  \param[in]      stat        :   The input scalar `integer` of default kind \IK containing the error stat code, frequently returned by the compiler.<br>
    !>                                  (**optional**. If missing, no error code will be suffixed to the output error message. It can be present only if `self` is of class [stop_type](@ref pm_err::stop_type).)
    !>  \param[in]      help        :   The input scalar `character` of the same kind as `msg` containing a descriptive help message to display separately after the fatal error message.<br>
    !>                                  (**optional**, default = `""`. It can be present only if `self` is of class [stop_type](@ref pm_err::stop_type).)
    !>  \param[in]      renabled    :   The input scalar `logical` of default kind \LK that if `.true.`, the program-control will be returned to the calling routine instead of terminating the program via `error stop`.<br>
    !>                                  This capability is particularly needed during the testing of the various components of the ParaMonte library, when fatal error must be handled gracefully.<br>
    !>                                  (**optional**, default = [SOFT_EXIT_ENABLED](@ref pm_err::SOFT_EXIT_ENABLED). It can be present only if `self` is of class [stop_type](@ref pm_err::stop_type).)
    !>
    !>  \interface{show}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: show
    !>
    !>      call self%show(msg, prefix = prefix, indent = indent, break = break, newline = newline, width = width, maxwidth = maxwidth, tmsize = tmsize, bmsize = bmsize, unit = unit)
    !>      call self%show(msg, prefix = prefix, indent = indent, break = break, newline = newline, width = width, maxwidth = maxwidth, tmsize = tmsize, bmsize = bmsize, unit = unit, stat = stat, help = help, renabled = renabled) ! only if self is of class `stop_type`.
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final{show}
    !>
    !>  \author
    !>  \AmirShahmoradi, 9:49 PM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface show

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setMarkedMethod_D0_SK5(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedMethod_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        class(mark_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setMarkedMethod_D0_SK4(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedMethod_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        class(mark_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setMarkedMethod_D0_SK3(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedMethod_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        class(mark_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setMarkedMethod_D0_SK2(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedMethod_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        class(mark_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setMarkedMethod_D0_SK1(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedMethod_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        class(mark_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setNotedMethod_D0_SK5(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedMethod_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        class(note_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setNotedMethod_D0_SK4(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedMethod_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        class(note_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setNotedMethod_D0_SK3(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedMethod_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        class(note_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setNotedMethod_D0_SK2(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedMethod_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        class(note_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setNotedMethod_D0_SK1(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedMethod_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        class(note_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setWarnedMethod_D0_SK5(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedMethod_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        class(warn_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setWarnedMethod_D0_SK4(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedMethod_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        class(warn_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setWarnedMethod_D0_SK3(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedMethod_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        class(warn_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setWarnedMethod_D0_SK2(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedMethod_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        class(warn_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setWarnedMethod_D0_SK1(self, msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedMethod_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        class(warn_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setAbortedMethod_D0_SK5(self, msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedMethod_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        class(stop_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setAbortedMethod_D0_SK4(self, msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedMethod_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        class(stop_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setAbortedMethod_D0_SK3(self, msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedMethod_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        class(stop_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setAbortedMethod_D0_SK2(self, msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedMethod_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        class(stop_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setAbortedMethod_D0_SK1(self, msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, sticky, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedMethod_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        class(stop_type)            , intent(inout)                 :: self
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
        logical(LK)                 , intent(in)    , optional      :: sticky
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Write the input string in the format of a custom notification, warning, or other types of messages to the output.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface call [getStrWrapped](@ref pm_str::getStrWrapped) with the input string
    !>  to obtain a wrapped text within the specified `width` and `maxwidth` where the input `prefix` is prepended
    !>  to the beginning of **each** wrapped line and then printed on the display.<br>
    !>
    !>  \param[in]  msg         :   The input scalar `character` of kind \SKALL representing the text to be wrapped and displayed.<br>
    !>  \param[in]  prefix      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**, default = `""`)
    !>  \param[in]  indent      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  break       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  newline     :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  width       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  maxwidth    :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  tmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines preceding the text in the display.<br>
    !>                              (**optional**, default = `1_IK`)
    !>  \param[in]  bmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines succeeding the text in the display.<br>
    !>                              (**optional**, default = `0_IK`)
    !>  \param[in]  unit        :   The input scalar `integer` of default kind \IK containing the output file unit.<br>
    !>                              (**optional**, default = `output_unit` taken from the intrinsic module `iso_fortran_env`.)
    !>
    !>  \interface{setMarked}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: setMarked
    !>
    !>      call setMarked(msg, prefix = prefix, indent = indent, break = break, newline = newline, width = width, maxwidth = maxwidth, tmsize = tmsize, bmsize = bmsize, unit = unit)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setMarked](@ref pm_err::setMarked)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_err/setMarked/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_err/setMarked/main.out.F90
    !>
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, 9:49 PM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface setMarked

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setMarkedStatic_D0_SK5(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedStatic_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setMarkedStatic_D0_SK4(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedStatic_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setMarkedStatic_D0_SK3(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedStatic_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setMarkedStatic_D0_SK2(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedStatic_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setMarkedStatic_D0_SK1(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMarkedStatic_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Write the input string in the format of a notification to the output.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface call [getStrWrapped](@ref pm_str::getStrWrapped) with the input string
    !>  to obtain a wrapped text within the specified `width` and `maxwidth` where the input `prefix` is suffixed with ` - NOTE: `
    !>  and then prepended to the beginning of each wrapped line and then printed on the display.<br>
    !>
    !>  \param[in]  msg         :   The input scalar `character` of kind \SKALL representing the text to be wrapped and displayed.<br>
    !>  \param[in]  prefix      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**, default = `""`)
    !>  \param[in]  indent      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  break       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  newline     :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  width       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  maxwidth    :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**. The default is set by [getStrWrapped](@ref pm_str::getStrWrapped).)
    !>  \param[in]  tmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines preceding the text in the display.<br>
    !>                              (**optional**, default = `1_IK`)
    !>  \param[in]  bmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines succeeding the text in the display.<br>
    !>                              (**optional**, default = `0_IK`)
    !>  \param[in]  unit        :   The input scalar `integer` of default kind \IK containing the output file unit.<br>
    !>                              (**optional**, default = `output_unit` taken from the intrinsic module `iso_fortran_env`.)
    !>
    !>  \interface{setNoted}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: setNoted
    !>
    !>      call setNoted(msg, prefix = prefix, indent = indent, break = break, newline = newline, width = width, maxwidth = maxwidth, tmsize = tmsize, bmsize = bmsize, unit = unit)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_err/setNoted/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_err/setNoted/main.out.F90
    !>
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, 9:49 PM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface setNoted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setNotedStatic_D0_SK5(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedStatic_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setNotedStatic_D0_SK4(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedStatic_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setNotedStatic_D0_SK3(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedStatic_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setNotedStatic_D0_SK2(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedStatic_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setNotedStatic_D0_SK1(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNotedStatic_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Write the input string in the format of a warning to the output.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface call [getStrWrapped](@ref pm_str::getStrWrapped) with the input string
    !>  to obtain a wrapped text within the specified `width` and `maxwidth` where the input `prefix` is suffixed with ` - WARNING: `
    !>  and then prepended to the beginning of each wrapped line and then printed on the display.<br>
    !>
    !>  \param[in]  msg         :   The input scalar `character` of kind \SKALL representing the text to be wrapped and displayed.<br>
    !>  \param[in]  prefix      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  indent      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  break       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  newline     :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  width       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  maxwidth    :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  tmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines preceding the text in the display.<br>
    !>                              (**optional**, default = `1_IK`)
    !>  \param[in]  bmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines succeeding the text in the display.<br>
    !>                              (**optional**, default = `0_IK`)
    !>  \param[in]  unit        :   The input scalar `integer` of default kind \IK containing the output file unit.<br>
    !>                              (**optional**, default = `output_unit` taken from the intrinsic module `iso_fortran_env`.)
    !>
    !>  \interface{setWarned}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: setWarned
    !>
    !>      call setWarned(msg, prefix = prefix, indent = indent, break = break, newline = newline, width = width, maxwidth = maxwidth, tmsize = tmsize, bmsize = bmsize, unit = unit)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_err/setWarned/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_err/setWarned/main.out.F90
    !>
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, 9:49 PM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface setWarned

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setWarnedStatic_D0_SK5(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedStatic_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setWarnedStatic_D0_SK4(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedStatic_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setWarnedStatic_D0_SK3(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedStatic_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setWarnedStatic_D0_SK2(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedStatic_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setWarnedStatic_D0_SK1(msg, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnedStatic_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Write the input string in the format of a fatal error to the output, then call `error stop` or return the control to the caller if requested.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface call [getStrWrapped](@ref pm_str::getStrWrapped) with the input string
    !>  to obtain a wrapped text within the specified `width` and `maxwidth` where the input `prefix` is suffixed with ` - FATAL: `
    !>  and then prepended to the beginning of each wrapped line and then printed on the display.<br>
    !>
    !>  \param[in]  msg         :   The input scalar `character` of kind \SKALL representing the text to be wrapped and displayed.<br>
    !>  \param[in]  help        :   The input scalar `character` of the same kind as `msg` containing a descriptive help message to display separately after the fatal error message.<br>
    !>                              (**optional**, default = `""`)
    !>  \param[in]  prefix      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  indent      :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  break       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  newline     :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  width       :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  maxwidth    :   See the corresponding definition in the documentation of [getStrWrapped](@ref pm_str::getStrWrapped).<br>
    !>                              (**optional**)
    !>  \param[in]  tmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines preceding the text in the display.<br>
    !>                              (**optional**, default = `1_IK`)
    !>  \param[in]  bmsize      :   The input scalar of `integer` of default kind \IK representing the number of empty lines succeeding the text in the display.<br>
    !>                              (**optional**, default = `0_IK`)
    !>  \param[in]  unit        :   The input scalar `integer` of default kind \IK containing the output file unit.<br>
    !>                              If the input `unit` is different from the standard output (`output_unit`), all messages will also be printed on stdout recursively.<br>
    !>                              (**optional**, default = `output_unit` taken from the intrinsic module `iso_fortran_env`.)
    !>  \param[in]  stat        :   The input scalar `integer` of default kind \IK containing the error stat code, frequently returned by the compiler.<br>
    !>                              (**optional**. If missing, no error code will be suffixed to the output error message.)
    !>  \param[in]  renabled    :   The input scalar `logical` of default kind \LK.<br>
    !>                              If `.true.`, the program-control will be returned to the calling
    !>                              routine instead of terminating the program via `error stop`.<br>
    !>                              This capability is particularly needed during the testing of the various components
    !>                              of the ParaMonte library, when fatal error must be handled gracefully.<br>
    !>                              (**optional**, default = [SOFT_EXIT_ENABLED](@ref pm_err::SOFT_EXIT_ENABLED))
    !>
    !>  \interface{setAborted}
    !>  \code{.F90}
    !>
    !>      use pm_err, only: setAborted
    !>
    !>      call setAborted(msg, help = help, prefix = prefix, indent = indent, break = break, newline = newline, width = width, maxwidth = maxwidth, tmsize = tmsize, bmsize = bmsize, unit = unit, stat = stat, renabled = renabled)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setNoted](@ref pm_err::setNoted)<br>
    !>  [setWarned](@ref pm_err::setWarned)<br>
    !>  [setAborted](@ref pm_err::setAborted)<br>
    !>  [setAsserted](@ref pm_err::setAsserted)<br>
    !>  [getLine](@ref pm_err::getLine)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_err/setAborted/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_err/setAborted/main.out.F90
    !>
    !>  [test_pm_err](@ref test_pm_err)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, 9:49 PM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
    interface setAborted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setAbortedStatic_D0_SK5(msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedStatic_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setAbortedStatic_D0_SK4(msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedStatic_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setAbortedStatic_D0_SK3(msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedStatic_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setAbortedStatic_D0_SK2(msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedStatic_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setAbortedStatic_D0_SK1(msg, help, prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit, stat, renabled)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAbortedStatic_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: msg
        character(*,SKG)            , intent(in)    , optional      :: help, prefix, indent, break, newline
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        integer(IK)                 , intent(in)    , optional      :: tmsize, bmsize, unit, stat
        logical(LK)                 , intent(in)    , optional      :: renabled
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_err ! LCOV_EXCL_LINE