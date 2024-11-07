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
!>  This module contains a simple unit-testing framework for the Fortran libraries, including the ParaMonte library.
!>
!>  \details
!>  Tests within the ParaMonte testing framework are based on modules.<br>
!>  Each module `pm_*` has a corresponding test module `test_pm_*` that implements tests of all objects and routines of the target module.<br>
!>  Each test module must contain a public `subroutine` named `setTest()` with the following interface,
!>  \code{.F90}
!>      subroutine setTest()
!>      end subroutine
!>  \endcode
!>  within which all test are implemented as desired.<br>
!>  These module subroutines are subsequently called in the ParaMonte `main.F90` file.<br>
!>  See the documentation of [test_type](@ref pm_test::test_type) for further details on the testing approach and assertion verification.<br>
!>
!>  \see
!>  [pm_err](@ref pm_err)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_test

    use pm_io, only: display_type
    use pm_kind, only: SK, IK, LK, RKG => RK
    use, intrinsic  :: iso_fortran_env, only: output_unit
    use pm_strANSI, only: bright, fcyan, fgreen, fyellow, fred, underlined, reset
    use pm_parallelism, only: image_type
    use pm_arrayResize, only: setResized
    use pm_container, only: css_type
    use pm_timer, only: timer_type
    use pm_val2str, only: getStr
    use pm_err, only: err_type

    implicit none

    private
    public :: SK, IK, LK, test_type, setSummary

    !>  \cond excluded
    character(*, SK)        , parameter     :: NLC = new_line(SK_"a")
    character(*, SK)        , parameter     :: MODULE_NAME = SK_"pm_test"
   !integer(IK)             , parameter     :: TEST_COUNTER_LEN = 8_IK
    integer(IK)             , parameter     :: TIME_FIELD_LEN = 9_IK
    integer(IK)                             :: mv_testCounterOld
    integer(IK)                             :: mv_testCounter
    integer(IK)                             :: mv_nfail
    integer(IK)                             :: mv_npass
    type(css_type)          , allocatable   :: mv_failedTestFuncName(:)
    character(:, SK)        , allocatable   :: mc_passedString
    character(:, SK)        , allocatable   :: mc_failedString
    logical(LK)                             :: mv_uninit = .true._LK
    class(timer_type)       , allocatable   :: mv_timer !<  \public This timer is exclusively used to time tests at the global level.
    type(image_type)                        :: mv_image
    !>  \endcond excluded

    !>  \brief
    !>  This is the module `private` derived type for constructing objects that contain the input and output directory paths for tests.
    !>
    !>  \details
    !>  This derived type is meant to be used internally within the parent module.<br>
    !>
    !>  \see
    !>  [dir_type](@ref pm_test::dir_type)<br>
    !>  [file_type](@ref pm_test::file_type)<br>
    !>  [test_type](@ref pm_test::test_type)<br>
    !>
    !>  \final{dir_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: dir_type
        character(:, SK)    , allocatable   :: inp !<  \public The scalar `allocatable` of type `character` of default kind \SK containing the path to the input  directory for the current test.
        character(:, SK)    , allocatable   :: out !<  \public The scalar `allocatable` of type `character` of default kind \SK containing the path to the output directory for the current test.
    end type


    !>  \brief
    !>  This is the module `private` derived type for constructing objects that contain the path to a file and its unit.
    !>
    !>  \details
    !>  This derived type is meant to be used internally within the parent module.<br>
    !>
    !>  \see
    !>  [dir_type](@ref pm_test::dir_type)<br>
    !>  [file_type](@ref pm_test::file_type)<br>
    !>  [test_type](@ref pm_test::test_type)<br>
    !>
    !>  \final{file_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: file_type
        integer(IK)                         :: unit
        character(:, SK)    , allocatable   :: path
    end type

    !type                                    :: strSplit_type
    !    integer(IK)         , allocatable   :: sindex(:,:)
    !    character(:, SK)    , allocatable   :: val
    !end type

    !type                                    :: subject_type
    !    type(strSplit_type)                 :: name             !<  \public An object containing information about the name and its parts (separated by underscore) of the current test subject (a procedure, type, ...).
    !end type

    !type                                    :: testFunc_type
    !    character(:, SK)    , allocatable   :: id               !<  \public The current number (ID) of the test function for a given object to be tested.
    !    character(:, SK)    , allocatable   :: name             !<  \public The name of the test function for the currently given object that is to be tested.
    !    class(timer_type)   , allocatable   :: timer            !<  \public This timer is exclusively used to time test cases within a given test internally.
    !    integer(IK)                         :: counter          !<  \public The counter for the test cases within a given test function for a given test subject.
    !end type

    !type                                    :: current_type
    !    integer(IK)                         :: counter          !<  \public The `protected` scalar object of type `integer` of default kind \IK, containing count of individual assertions made within a specified test function or test module.<br>
    !                                                            !!          The procedure name is automatically extracted from the input procedure name to the [run](@ref pm_test::setTestFunc) method of the object of type [test_type](@ref pm_test::test_type).<br>
    !end type

    type                                    :: scope_type
        class(timer_type)   , allocatable   :: timer            !<  \public This timer is exclusively used to time tests at the target level/scope (that is, to time all tests within a given module, or within a given test function).<br>
        character(:, SK)    , allocatable   :: name             !<  \public The `public` scalar `allocatable` of type `character` of default kind \SK, containing the name of the current scope of testing (e.g., host name of the procedure being tested, or the procedure name).<br>
    end type

    type, extends(scope_type)               :: func_type
        integer(IK)                         :: counter = 0_IK   !<  \private    The `private` scalar object of type `integer` of default kind \IK, containing count of individual assertions made within a
                                                                !!              specified test function passed, via a single call, to the `run()` method of an object of type [test_type](@ref pm_test::testt_type).<br>
        character(:, SK)    , allocatable   :: id               !<  \public     The `public` scalar `allocatable` of type `character` of default kind \SK, containing the ID of the current test for a given procedure to be tested.<br>
                                                                !!              The presence of this component is historical when individual tests for a single target procedure were implemented in separate test functions.<br>
                                                                !!              In such a case, the test functions of a given procedure can be suffixed with a unique ID as in `test_*_ID()` where `ID` is replaced with the test function number.<br>
                                                                !!              Test function IDs are not necessary anymore as of ParaMonvte V.2 because individual test functions can now implement as many tests as desired.<br>
                                                                !!              If a test function does not have an ID, the `id` component of the parent object of type [test_type](@ref pm_test::test_type) is set to the default value of `1`.<br>
    end type

    !>  \brief
    !>  This is the derived type [test_type](@ref pm_test::test_type) for generating objects
    !>  that facilitate testing of a series of procedures or concepts within a unified testing framework.
    !>
    !>  \details
    !>  This procedure is a constructor of [test_type](@ref pm_test::test_type).<br>
    !>  See also the documentation of [test_typer](@ref pm_test::test_typer), the constructor of [test_type](@ref pm_test::test_type).<br>
    !>
    !>  \param[in]  host        :   The input scalar of type `character` of default kind \SK containing the name of the host of the procedure or code that is being tested.<br>
    !>                              Frequently, in modern Fortran, this is the name of the module containing the procedure being tested.<br>
    !>                              This name is used only for information display purposes.<br>
    !>                              (**optional**, default = `SK_""`, that is, the target procedures to be tested are assumed to be `external` and not in a module or in any specific host scope.)
    !>  \param[in]  inp         :   The input scalar of type `character` of arbitrary length type parameter, containing the input directory
    !>                              for the current testing suite (i.e., all tests performed with the output object of this constructor).<br>
    !>                              The specified value for this argument filles the corresponding `inp` component of the `dir` component of the output `test` object.<br>
    !>                              (**optional**, default = `SK_"./input"`)
    !>  \param[in]  out         :   The input scalar of type `character` of arbitrary length type parameter, containing the output directory
    !>                              for the current testing suite (i.e., all tests performed with the output object of this constructor).<br>
    !>                              The specified value for this argument filles the corresponding `out` component of the `out` component of the output `test` object.<br>
    !>                              (**optional**, default = `SK_"./output"`)
    !>  \param[in]  traceable   :   The input scalar of type `logical` of default kind \LK, that can be used to control the display of debugging information for failed tests.<br>
    !>                              Specifying this argument directly sets the corresponding `traceable` component of the output `test` object.<br>
    !>                              By convention, setting `traceable` to `.true.` is used within the test implementations to allow displaying information about assertion failures should any occur.<br>
    !>                              The ability to display descriptive information about the source of a test failure and its origins if often desired at the time of debugging.<br>
    !>                              Therefore, this component is set to `.true.` when runtime checks are enabled (i.e., FFP macro `CHECK_ENABLED=1` is set, frequently corresponding to the debug mode)
    !>                              and `.false.` otherwise (frequently corresponding to the release mode).<br>
    !>                              This optional input argument overrides the default behavior.<br>
    !>                              **Note that if `traceable` is set to `.true._LK`, the program will `error stop` by encountering the first test failure**.<br>
    !>                              (**optional**, default = `.true._LK` if the library is built with preprocessor macro `CHECK_ENABLED=1`, otherwise `.false._LK`)
    !>
    !>  \return
    !>  `test`                  :   The output scalar object of type [test_type](@ref pm_test::test_type) containing the specifics of the runtime system shell.<br>
    !>
    !>  \interface{test_typer}
    !>  \code{.F90}
    !>
    !>      use pm_test, only: test_type
    !>      type(test_type) :: test
    !>
    !>      test = test_type(host = host, inp = inp, out = out, traceable = traceable)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure must be called by all images/processes in parallel mode as it contains a call to `image_type()%%sync`.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [test_type](@ref pm_test::test_type)<br>
    !>
    !>  \final{test_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: test_type
        integer(IK)                         , private   :: counter = 0_IK           !<  \private    The `private` scalar object of type `integer` of default kind \IK, containing count of individual assertions made by the `run()` method of an object of type [test_type](@ref pm_test::test_type).<br>
        logical(LK)                         , private   :: asserted = .true._LK     !<  \private    The `private` scalar of type `logical` of default kind \LK that is `.true.` <b>if and only if</b> all tests within the current test suite performed by an object of [test_type](@ref pm_test::test_type) have passed successfully. It is used for test summary.<br>
        logical(LK)                         , public    :: traceable                !<  \public     The `protected` scalar object of type `logical` of default kind \LK, that can be used to control the display of debugging information for failed tests. See also the corresponding `optional` argument of the [constructor](@ref pm_test::test_type).<br>
        type(dir_type)                      , public    :: dir                      !<  \public     The `public` scalar object of type [dir_type](@ref pm_test::dir_type) containing the input and output directories for the test data.<br>
        type(func_type)                     , public    :: func                     !<  \public     The scalar object of type [func_type](@ref pm_test::func_type) containing the name of the current procedure being tested and the timer starting at the beginning of the current test (corresponding to the most recent call to the `run()` method of the object of type [test_type](@ref pm_test::test_type).<br>
        type(scope_type)                    , public    :: host                     !<  \public     The scalar object of type [scope_type](@ref pm_test::scope_type) containing the name of the current host (scoping unit) being tested and the timer starting at the beginning of the creating the object of type [test_type](@ref pm_test::test_type).<br>
        type(display_type)                  , public    :: disp                     !<  \public     The `public` scalar of type [display_type](@ref pm_io::display_type) containing tools for displaying information in the desired output, for example, debugging information when a test fails.<br>
        type(image_type)                    , public    :: image                    !<  \public     The `public` scalar of type [image_type](@ref pm_parallelism::image_type) containing information about the parallel processes in parallel mode.<br>
        type(file_type)                     , public    :: file                     !<  \public     The `public` scalar of type [file_type](@ref pm_test::file_type) that serves as convenience to store the unit number and file path of any working file during the testing.<br>
       !type(err_type)                      , private   :: err                      !<  \private    The `private` scalar of type [err_type](@ref pm_err::err_type) containing information about any errors hat may occur during the testing.<br>
    contains
       !procedure                           , pass      :: openFile
        procedure                           , pass      :: summarize => setTestSummary
        procedure                           , pass      :: assert => setTestAsserted
        procedure                           , pass      :: run => setTestFunc
    end type

    interface test_type
        module procedure :: test_typer
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Initialize the global module variables for testing.<br>
    !>
    !>  \details
    !>  This subroutine is called internally by the procedures of this module and methods of objects of type [test_type](@ref pm_test::test_type)
    !>  at the beginning of the first test, before any testing is done to ensure information for
    !>  the final global summary by [setSummary](@ref pm_test::setSummary) is properly collected.
    !>
    !>  \interface{setInitial}
    !>  \code{.F90}
    !>
    !>      call setInitial()
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [test_type](@ref pm_test::test_type)<br>
    !>  [setSummary](@ref pm_test::setSummary)<br>
    !>
    !>  \final{setInitial}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    subroutine setInitial()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInitial
#endif
        ! set the mode of error handling to testing mode
        !mv_isTestingMode = .true._LK
        mv_uninit = .false._LK
        mv_image = image_type()
        mv_timer = timer_type()
        mc_passedString = bright//fgreen//SK_"passed"//reset
        mc_failedString = bright//fred  //SK_"FAILED"//reset
        mv_testCounterOld = 0_IK
        mv_testCounter = 0_IK
        mv_npass = 0_IK
        mv_nfail = 0_IK
        ! preallocate the names of failed tests.
        call setResized(mv_failedTestFuncName, mv_nfail + 1_IK)
        write(output_unit, "(*(g0,:,' '))")
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Summarize the collection of all tests performed on all modules (or scoping units).<br>
    !>
    !>  \details
    !>  Calling this procedure at the end of a series of tests performed by a collection of objects of type [test_type](@ref pm_test::test_type)
    !>  will lead to display of relevant summary of all passed and failed tests performed and further details about the failed tests and test timing.<br>
    !>
    !>  \interface{setSummary}
    !>  \code{.F90}
    !>
    !>      call setSummary()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure must be called by all parallel images/processes/threads.<br>
    !>
    !>  \see
    !>  [test_type](@ref pm_test::test_type)<br>
    !>  [setInitial](@ref pm_test::setInitial)<br>
    !>
    !>  \final{setSummary}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    subroutine setSummary()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSummary
#endif
        use pm_io, only: TABEQV
        use pm_err, only: setAborted
        use pm_mathNumSys, only: getCountDigit
        real(RKG) :: percentageTestPassed
        character(:, SK), allocatable :: msg, color
        character(*), parameter :: format = "(*(g0,:,' '))"
        integer(IK) :: ntotal, ifail, ndigit
        if (mv_uninit) call setInitial()
        ntotal = mv_npass + mv_nfail
        percentageTestPassed = 0._RKG
        if (ntotal /= 0_IK) percentageTestPassed = 100_IK * mv_npass / real(ntotal, RKG)
        if (mv_image%is%first) then
            msg = bright//fyellow//getStr(ntotal)//SK_" tests performed. "//reset
            if (nint(percentageTestPassed, IK) == 0_IK) then
                color = fred
            else
                color = fgreen
            end if
            msg = msg//bright//color//getStr(percentageTestPassed,SK_"(f0.2)")//SK_"% of "//getStr(ntotal)//SK_" tests passed. "//reset
            if (0_IK < mv_nfail) then
                color = fred
            else
                color = fgreen
            end if
            msg = msg//bright//color//getStr(mv_nfail, SK_"(g0)")//SK_" tests failed. "//reset
            msg = msg//bright//fyellow//SK_"The total elapsed wall-clock time: "//getStr(mv_timer%time(since = mv_timer%start), SK_"(f0.6)")//SK_" seconds."//reset
            write(output_unit, format) NLC//msg//NLC
            ndigit = getCountDigit(mv_nfail)
            if (mv_nfail > 0_IK) then
                write(output_unit, format) NLC//bright//fred//SK_"The following tests FAILED:"//reset//NLC
                do ifail = 1, mv_nfail
                    write(output_unit, format) bright//fred//TABEQV//adjustr(getStr(ifail, length = ndigit))//SK_") "//mv_failedTestFuncName(ifail)%val//reset
                end do
                write(output_unit, format) NLC//bright//fred//SK_"Errors occurred while running the ParaMonte tests."
                write(output_unit, format) SK_"Please report this issue at: "//bright//fcyan//underlined//SK_"https://github.com/cdslaborg/paramonte/issues"//reset//NLC//NLC
#if             !CHECK_ENABLED
                write(output_unit, format) NLC//bright//fred//SK_"To get information about the cause of failure, compile the library with FPP runtime checks enabled and rerun the tests."//reset//NLC
#endif
            end if
        end if
        !if (mv_nfail > 0_IK) !error stop !call setAborted()
        call mv_image%finalize()
        mv_uninit = .true._LK
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an object of type [test_type](@ref pm_test::test_type).
    !>
    !>  \details
    !>  This procedure is a constructor of [test_type](@ref pm_test::test_type).<br>
    !>  See also the documentation of [test_type](@ref pm_test::test_type).<br>
    !>
    !>  \param[in]  host        :   The input scalar of type `character` of default kind \SK containing the name of the host of the procedure or code that is being tested.<br>
    !>                              Frequently, in modern Fortran, this is the name of the module containing the procedure being tested.<br>
    !>                              This name is used only for information display purposes.<br>
    !>                              (**optional**, default = `SK_"@unknown_scope"`, that is, the target procedures to be tested are assumed to be `external` and not in a module or in any specific host scope.)
    !>  \param[in]  inp         :   The input scalar of type `character` of arbitrary length type parameter, containing the input directory
    !>                              for the current testing suite (i.e., all tests performed with the output object of this constructor).<br>
    !>                              The specified value for this argument filles the corresponding `inp` component of the `dir` component of the output `test` object.<br>
    !>                              (**optional**, default = `SK_"./input"`)
    !>  \param[in]  out         :   The input scalar of type `character` of arbitrary length type parameter, containing the output directory
    !>                              for the current testing suite (i.e., all tests performed with the output object of this constructor).<br>
    !>                              The specified value for this argument filles the corresponding `out` component of the `out` component of the output `test` object.<br>
    !>                              (**optional**, default = `SK_"./output"`)
    !>  \param[in]  traceable   :   The input scalar of type `logical` of default kind \LK, that can be used to control the display of debugging information for failed tests.<br>
    !>                              Specifying this argument directly sets the corresponding `traceable` component of the output `test` object.<br>
    !>                              By convention, setting `traceable` to `.true.` is used within the test implementations to allow displaying information about assertion failures should any occur.<br>
    !>                              The ability to display descriptive information about the source of a test failure and its origins if often desired at the time of debugging.<br>
    !>                              Therefore, this component is set to `.true.` when runtime checks are enabled (i.e., FFP macro `CHECK_ENABLED=1` is set, frequently corresponding to the debug mode)
    !>                              and `.false.` otherwise (frequently corresponding to the release mode).<br>
    !>                              This optional input argument overrides the default behavior.<br>
    !>                              **Note that if `traceable` is set to `.true._LK`, the program will `error stop` by encountering the first test failure**.<br>
    !>                              (**optional**, default = `.true._LK` if the library is built with preprocessor macro `CHECK_ENABLED=1`, otherwise `.false._LK`)
    !>
    !>  \return
    !>  `test`                  :   The output scalar object of type [test_type](@ref pm_test::test_type) containing the specifics of the runtime system shell.<br>
    !>
    !>  \interface{test_typer}
    !>  \code{.F90}
    !>
    !>      use pm_test, only: test_type
    !>      type(test_type) :: test
    !>
    !>      test = test_type(host = host, inp = inp, out = out, traceable = traceable)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure must be called by all images/processes in parallel mode as it contains a call to [image_type()%sync()](@ref pm_parallelism::image_type)`.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [test_type](@ref pm_test::test_type)<br>
    !>
    !>  \final{test_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    function test_typer(host, inp, out, traceable) result(test)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: test_typer
#endif
        use pm_err, only: err_type, setAborted
        use pm_sysPath, only: isFailedMakeDir, isDir
        character(*, SK), intent(in), optional :: host, inp, out
        logical(LK), intent(in), optional :: traceable
        type(test_type) :: test
        if (mv_uninit) call setInitial()
        test%image = mv_image
        test%disp = display_type()
        test%func%timer = timer_type()
        test%host%timer = timer_type()
        if (present(traceable)) then
            test%traceable = traceable
        else
#if         CHECK_ENABLED
            test%traceable = .true._LK
#else
            test%traceable = .false._LK
#endif
        end if
        if (present(host)) then
            test%host%name = trim(adjustl(host))
        else
            test%host%name = SK_"@unknown_scope"
        end if
        if (present(inp)) then
            test%dir%inp = inp
        else
            test%dir%inp = SK_"./input"
        end if
        if (present(out)) then
            test%dir%out = out
        else
            test%dir%out = SK_"./output"
        end if
        ! mkdir the output directory if it does not exists.
        if (test%image%is%first .and. .not. isDir(test%dir%out)) then
            if (isFailedMakeDir(test%dir%out)) error stop MODULE_NAME//SK_"@test_typer(): Failed to generate the test module output directory: "//test%dir%out
        end if
        call test%image%sync()
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Run the input test function and verify the assertion returned by the test function.
    !>
    !>  \details
    !>  This procedure is a method of [test_type](@ref pm_test::test_type).<br>
    !>  See also the documentation of [test_type](@ref pm_test::test_type).<br>
    !>
    !>  \param[inout]   self            :   The input/output object of class [test_type](@ref pm_test::test_type) that is passed implicitly.
    !>  \param          getAssertion    :   The `external` user-specified test function that takes no input arguments and returns a scalar of type `logical` of default kind \LK.<br>
    !>                                      The returned value must be `.true.` <b>if and only if</b> the test assertion is verified within the specified test function.<br>
    !>                                      The following illustrates the generic interface of `getAssertion()`,
    !>                                      \code{.F90}
    !>                                          function getAssertion() result(assertion)
    !>                                              use pm_kind, only: LK
    !>                                              logical(LK) :: assertion
    !>                                          end function
    !>                                      \endcode
    !>  \param[in]      name            :   The input scalar of type `character` of default kind \SK containing the name of the test function as is,
    !>                                      that is, the actual name of target of the dummy argument `getAssertion()`.<br>
    !>                                      The test function name is assumed to follow the template `test_*()` or `test_*_ID()` where,<br>
    !>                                      <ol>
    !>                                          <li>    `*` is taken as the name of the target procedure that is being tested.<br>
    !>                                                  This name is used for displaying testing details on the specified output.<br>
    !>                                          <li>    `ID` is the optional integer in test function name representing the ID of the current test function for the target procedure being tested.<br>
    !>                                                  The default value for `ID` is `1`, if missing.<br>
    !>                                      </ol>
    !>                                      If the specified function name is not prefixed with `test_`, the whole name (excluding any `ID` suffix) is taken as the name of the procedure being tested.<br>
    !>                                      (**optional**, default = `SK_"@test_unknown"`.)
    !>
    !>  \interface{setTestFunc}
    !>  \code{.F90}
    !>
    !>      use pm_test, only: test_type
    !>      type(test_type) :: test
    !>
    !>      test = test_type()
    !>      call test%run(getAssertion, name = name)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [test_type](@ref pm_test::test_type)<br>
    !>
    !>  \final{setTestFunc}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    subroutine setTestFunc(self, getAssertion, name)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setTestFunc
#endif
        use pm_strASCII, only: isStrInteger
        procedure(logical(LK)) :: getAssertion
        class(test_type), intent(inout) :: self
        character(*, SK), intent(in), optional :: name
        logical(LK) :: assertion
        if (present(name)) then
            self%func%name = trim(adjustl(name))
            if (4_IK < len(self%func%name, IK)) then
                if (self%func%name(1:5) == SK_"test_") self%func%name = self%func%name(6:)
            end if
            block
                use pm_arrayFind, only: setLoc
                integer(IK), allocatable :: loc(:)
                integer(IK) :: nloc
                call setResized(loc, 2_IK)
                call setLoc(loc, nloc, self%func%name, SK_"_", blindness = 1_IK)
                if (nloc == 0_IK) then
                    self%func%id = SK_"#1"
                elseif (isStrInteger(self%func%name(loc(nloc) + 1 :))) then
                    self%func%id = SK_"#"//self%func%name(loc(nloc) + 1 :)
                    self%func%name = self%func%name(1 : loc(nloc) - 1)
                else
                    self%func%id = SK_"#1"
                    self%func%name = name
                end if
            end block
        else
            self%func%id = SK_"#1"
            self%func%name = SK_"unknown"
        end if
        self%counter = 0_IK
       !self%func%name = SK_"Test_"//self%host%name(2:)//SK_"@"//trim(adjustl(name))
        self%func%name = self%host%name//SK_"@"//trim(adjustl(self%func%name))
        self%func%timer = timer_type()
        assertion = getAssertion()
        self%asserted = self%asserted .and. assertion
        call self%assert(assertion)
        self%func%counter = 0_IK
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Test the validity of the input `assertion` and if it does not hold,<br>
    !>  <ol>
    !>      <li>    in `debug` and `testing` modes, call `error stop` after displaying the assertion description `desc`.<br>
    !>      <li>    in `release` mode, report the test failure and continue with the rest of the tests.<br>
    !>  </ol>
    !>
    !>  \details
    !>  This procedure is a method of the class [test_type](@ref pm_test::test_type).<br>
    !>
    !>  \param[inout]   self        :   The input/output object of class [test_type](@ref pm_test::test_type) that is passed implicitly.
    !>  \param[in]      assertion   :   The input scalar of type `logical` of default kind \LK containing the assertion to be verified.
    !>  \param[in]      desc        :   The input scalar of type `character` of default kind \SK of arbitrary length type parameter containing
    !>                                  a description of assertion for display if the assertion fails and `self%traceable` is `.true.`.<br>
    !>                                  (**optional**, default = `"unavailable"`)
    !>  \param[in]      line        :   The input scalar of type `integer` of default kind \IK,
    !>                                  containing the line number for the source of the assertion within the test.<br>
    !>                                  Supplying this information helps trace the source of the test failure should it occur.<br>
    !>                                  (**optional**)
    !>
    !>  \interface{setTestAsserted}
    !>  \code{.F90}
    !>
    !>      use pm_test, only: test_type
    !>      type(test_type) :: test
    !>
    !>      test = test_type()
    !>      call test%run(assertion = assertion, desc = desc, line = line)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [test_type](@ref pm_test::test_type)<br>
    !>
    !>  \final{setTestAsserted}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 20, 2015, 1:11 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    subroutine setTestAsserted(self, assertion, desc, line)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setTestAsserted
#endif
        use pm_err, only: setAborted
        use pm_arrayInit, only: getCoreHalo
        logical(LK), intent(in) :: assertion
        class(test_type), intent(inout) :: self
        integer(IK), intent(in), optional :: line
        character(*, SK), intent(in), optional :: desc
        character(:, SK), allocatable :: statusMsg
        character(:, SK), allocatable :: desc_def
        character(:, SK), allocatable :: traceback
        character(:, SK), allocatable :: assertionID
        character(:, SK), allocatable :: errmsg
        integer(IK) :: iid
        self%func%timer%delta = self%func%timer%time(since = self%func%timer%start) - self%func%timer%clock
        self%func%timer%clock = self%func%timer%time(since = self%func%timer%start)
        self%func%counter = self%func%counter + 1_IK
        if (present(desc)) then
            desc_def = desc
        else
            desc_def = SK_"unavailable"
        end if
        if (present(line)) then
            traceback = SK_"traceback test source line: "//getStr(line)
        else
            traceback = SK_""
        end if
        assertionID = self%func%name//self%func%id//SK_"-"//getStr(self%func%counter)//SK_" "
        if (assertion) then
            mv_npass = mv_npass + 1_IK
            statusMsg = mc_passedString
        else
            mv_nfail = mv_nfail + 1_IK
            statusMsg = mc_failedString
            if (size(mv_failedTestFuncName, 1, IK) < mv_nfail) call setResized(mv_failedTestFuncName)
            mv_failedTestFuncName(mv_nfail)%val = assertionID
        end if
        mv_testCounter = mv_testCounter + 1_IK
        do iid = 1, self%image%count
            if (iid == self%image%id) then
                write(self%disp%unit, "(*(A))") & ! LCOV_EXCL_LINE
                SK_"["//adjustr(getStr(mv_testCounter))//SK_"] testing "// & ! LCOV_EXCL_LINE
                getCoreHalo(79_IK, assertionID, SK_".", 0_IK)//SK_" "//statusMsg// & ! LCOV_EXCL_LINE
                SK_" in "//adjustr(getStr(self%func%timer%delta, SK_"(f0.4)", TIME_FIELD_LEN))// & ! LCOV_EXCL_LINE
                SK_" out of "//adjustr(getStr(self%func%timer%clock, SK_"(f0.4)", TIME_FIELD_LEN))// & ! LCOV_EXCL_LINE
                SK_" seconds on image "//getStr(self%image%id)
            end if
#if         CAF_ENABLED || MPI_ENABLED
            block
                integer :: stat
                call execute_command_line(" ", cmdstat = stat)
                flush(output_unit)
                call self%image%sync()
            end block
#endif
        end do
        call self%func%timer%wait(0.0001_RKG)
        if (self%traceable .and. .not. assertion) then
            errmsg = SK_"The test assertion is FALSE."//NLC//SK_"The assertion description: "//trim(adjustl(desc_def))//NLC//traceback ! LCOV_EXCL_LINE
            call setAborted(msg = errmsg, prefix = SK_"ParaMonteTest"//self%func%name) ! LCOV_EXCL_LINE
            error stop ! LCOV_EXCL_LINE
        end if
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Summarize the suite of tests performed by the parent object of type [test_type](@ref pm_test::test_type) of this method.
    !>
    !>  \details
    !>  This procedure is a method of [test_type](@ref pm_test::test_type).<br>
    !>  See also the documentation of [test_type](@ref pm_test::test_type).<br>
    !>  This procedure can be effectively considered as the destructor of an object of type [test_type](@ref pm_test::test_type).<br>
    !>
    !>  \param[inout]   self    :   The input/output scalar object of class [test_type](@ref pm_test::test_type) that is to be destructed.<br>
    !>
    !>  \warning
    !>  This method must be called by all images/processes/threads as it contains a global synchronization upon exit.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [test_type](@ref pm_test::test_type)<br>
    !>
    !>  \final{setTestSummary}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    subroutine setTestSummary(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setTestSummary
#endif
        use pm_arrayInit, only: getCoreHalo
        use pm_mathNumSys, only: getCountDigit
        class(test_type), intent(inout) :: self
        character(:, SK), allocatable :: statusMsg
        if (self%asserted) then
            statusMsg = mc_passedString
        else
            statusMsg = mc_failedString
        end if
        self%host%timer%delta = self%host%timer%time(since = self%host%timer%start)
        write(self%disp%unit, "(*(g0,:,' '))") bright//fyellow// & ! LCOV_EXCL_LINE
        SK_"["//adjustr(getStr(mv_testCounter - mv_testCounterOld, length = getCountDigit(mv_testCounter)))//SK_"] testing "// & ! LCOV_EXCL_LINE
        getCoreHalo(79_IK, self%host%name//SK_" ", SK_".", 0_IK)//SK_" "//statusMsg// & ! LCOV_EXCL_LINE
        SK_" in "//adjustr(getStr(self%host%timer%delta, SK_"(f0.4)", TIME_FIELD_LEN))// & ! LCOV_EXCL_LINE
        SK_" out of "//adjustr(getStr(mv_timer%time(since = mv_timer%start), SK_"(f0.4)", TIME_FIELD_LEN))// & ! LCOV_EXCL_LINE
        SK_" seconds on image "//getStr(self%image%id)//reset
        mv_testCounterOld = mv_testCounter
        call self%image%sync()
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if 0
    !>  \cond excluded
    function openFile(self, path, label, status, position) result(file)
        class(test_type), intent(in) :: self
        character(*, SK), intent(in), optional :: path, label, status, position
        character(:, SK), allocatable :: prefix_def, status_def, position_def
        type(file_type) :: file
        if (present(position)) then
            position_def = position
        else
            position_def = SK_"asis"
        end if
        if (present(status)) then
            status_def = status
        else
            status_def = SK_"unknown"
        end if
        if (present(path)) then
            file%path = path
        else
            if (present(label)) then
                prefix_def = SK_"@"//label
            else
                prefix_def = SK_""
            end if
            file%path = self%dir%out//SK_"/"//self%func%name//prefix_def//SK_"@"//getStr(self%image%id)//SK_".txt"
        end if
        !if (getStrLower(position_def)=="append") status_def = SK_"old"
#if     INTEL_ENABLED && WINDOWS_ENABLED
#define SHARED, shared
#else
#define SHARED
#endif
        open(file = file%path, newunit = file%unit, status = status_def, position = position_def SHARED)
#undef  SHARED
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    ! Return the negative natural logarithm of MVN distribution evaluated at the input vector `Point` of length `ndim`.
    function getLogFuncMVN(ndim,Point) result(logFunc) BIND_C
        implicit none
        integer(IK) , intent(in)        :: ndim
        real(RKG)   , intent(in)        :: Point(ndim)
        real(RKG)                       :: logFunc
#if     CFI_ENABLED
        value                           :: ndim
#endif
        real(RKG)   , parameter         :: LOG_INVERSE_SQRT_TWO_PI = log(1._RKG / sqrt(2._RKG * acos(-1._RKG)))

        !block
        !    use pm_sysShell, only: sleep
        !    use pm_err, only: err_type
        !    type(err_type) :: err
        !    call sleep(seconds=5000.e-6_RKG,err=err)
        !end block

        !block
        !    real(RKG), allocatable :: unifrnd(:,:)
        !    allocate(unifrnd(200,20))
        !    call random_number(unifrnd)
        !    logFunc = sum(unifrnd) - 0.5_RKG * sum(Point**2) - sum(unifrnd)
        !    deallocate(unifrnd)
        !end block

        logFunc = ndim * LOG_INVERSE_SQRT_TWO_PI - 0.5_RKG * sum(Point**2_IK)

        !block
        !    integer(IK), parameter :: nmode = 2_IK
        !    real(RKG):: LogAmplitude(nmode), mean(nmode), invCov(nmode), logSqrtDetInvCovMat(nmode)
        !    LogAmplitude           = [1._RKG, 1._RKG]
        !    mean                   = [0._RKG, 7._RKG]
        !    invCov              = [1._RKG,1._RKG]
        !    logSqrtDetInvCovMat    = [1._RKG,1._RKG]
        !    logFunc = getLogProbGausMix ( nmode = 2_IK &
        !                                , nd = 1_IK &
        !                                , np = 1_IK &
        !                                , LogAmplitude = LogAmplitude &
        !                                , mean = mean &
        !                                , invCov = invCov &
        !                                , logSqrtDetInvCovMat = logSqrtDetInvCovMat &
        !                                , Point = Point(1) &
        !                                )
        !end block

    end function getLogFuncMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the negative natural logarithm of MVN distribution evaluated at the input vector `Point` of length `ndim`.
    function getLogFuncBanana2D(ndim,Point) result(logFunc) BIND_C
        use pm_distMultiNorm, only: getMultiNormLogPDFNF
        use pm_distMultiNorm, only: getMultiNormLogPDF
        use pm_mathLogSumExp, only: getLogSumExp
        implicit none
        integer(IK) , intent(in)        :: ndim
        real(RKG)   , intent(in)        :: Point(ndim)
        real(RKG)                       :: logFunc
#if     CFI_ENABLED
        value                           :: ndim
#endif
        integer(IK) , parameter         :: NPAR = 2_IK ! sum(Banana,gaussian) normalization factor
        real(RKG)   , parameter         :: normfac = 0.3_RKG ! 0.6_RKG ! sum(Banana,gaussian) normalization factor
        real(RKG)   , parameter         :: lognormfac = log(normfac) ! sum(Banana,gaussian) normalization factor
        real(RKG)   , parameter         :: a = 0.7_RKG, b = 1.5_RKG ! parameters  of Banana function
        real(RKG)   , parameter         :: MeanB(NPAR) = [ -5.0_RKG , 0._RKG ] ! mean vector of Banana function
        real(RKG)   , parameter         :: MeanG(NPAR) = [  3.5_RKG , 0._RKG ] ! mean vector of Gaussian function
        real(RKG)   , parameter         :: CovMatB(NPAR,NPAR) = reshape([0.25_RKG,0._RKG,0._RKG,0.81_RKG],shape=shape(CovMatB)) ! Covariance matrix of Banana function
        real(RKG)   , parameter         :: CovMatG(NPAR,NPAR) = reshape([0.15_RKG,0._RKG,0._RKG,0.15_RKG],shape=shape(CovMatB)) ! Covariance matrix of Gaussian function
        real(RKG)   , parameter         :: InvCovMatB(NPAR,NPAR) = reshape([4._RKG,0._RKG,0._RKG,1.23456790123457_RKG],shape=shape(InvCovMatB)) ! Inverse Covariance matrix of Banana function
        real(RKG)   , parameter         :: InvCovMatG(NPAR,NPAR) = reshape([6.66666666666667_RKG,0._RKG,0._RKG,6.66666666666667_RKG],shape=shape(InvCovMatG)) ! Inverse Covariance matrix of Gaussian function
        real(RKG)   , parameter         :: logSqrtDetInvCovB = log(sqrt(4.93827160493827_RKG)) ! Determinant of the Inverse Covariance matrix of Banana function
        real(RKG)   , parameter         :: logSqrtDetInvCovG = log(sqrt(44.4444444444445_RKG)) ! Determinant of the Inverse Covariance matrix of Gaussian function
        real(RKG)                       :: PointSkewed(NPAR) ! transformed parameters that transform the Gaussian to the Banana function
        real(RKG)                       :: LogProb(2) ! logProbMVN, logProbBanana

        PointSkewed(1) = -Point(1)
        PointSkewed(2) = +Point(2)

        ! Gaussian function

        LogProb(1) = lognormfac + getMultiNormLogPDF(PointSkewed, MeanG, InvCovMatG, getMultiNormLogPDFNF(ndim, logSqrtDetInvCovG)) ! logProbMVN

        ! Do variable transformations for the Skewed-Gaussian (banana) function.

        PointSkewed(2)  = a * PointSkewed(2)
        PointSkewed(1)  = PointSkewed(1)/a - b*(PointSkewed(2)**2 + a**2)

        ! Banana function

        LogProb(2) = lognormfac + getMultiNormLogPDF(PointSkewed, MeanB, InvCovMatB, getMultiNormLogPDFNF(ndim, logSqrtDetInvCovB)) ! logProbBanana

        !MeanBnew(2) = a * MeanBnew(2)
        !MeanBnew(1) = MeanBnew(1)/a - b*(MeanBnew(2)**2 + a**2)

        !logFunc = LogProb(2)
        !logFunc = lognormfac + LogProb(1)
        logFunc = getLogSumExp(LogProb, maxval(LogProb))

    end function getLogFuncBanana2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the 2D EggBox function in the cubic unit domain [0,1].
    !>  \remark
    !>  The log(integral) of this function in the unit domain is `235.88`.
    function getLogFuncEggBox2D(ndim,Point) result(logFunc) BIND_C
        integer(IK) , intent(in)        :: ndim
        real(RKG)   , intent(in)        :: Point(ndim)
        real(RKG)                       :: logFunc
        real(RKG)   , parameter         :: PI = acos(-1._RKG)
#if     CFI_ENABLED
        value :: ndim
#endif
        !logFunc = (2._RKG + cos(0.5_RKG*Point(1)) * cos(0.5_RKG*Point(2)) )**5_IK
        logFunc = (2._RKG + cos(5_IK * PI * Point(1) - 2.5_RKG * PI) * cos(5_IK * PI * Point(2) - 2.5_RKG * PI)) ** 5.0
    end function getLogFuncEggBox2D
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#endif

end module pm_test ! LCOV_EXCL_LINE
