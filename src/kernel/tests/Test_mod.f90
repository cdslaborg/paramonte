!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This is the driver module of the unit testing framework of the ParaMonte kernel routines.
!>  \author Amir Shahmoradi

! LCOV_EXCL_START
module Test_mod

    use, intrinsic  :: iso_fortran_env, only: output_unit
    use JaggedArray_mod, only: CharVec_type
    use Parallelism_mod, only: Image_type
    use Timer_mod, only: Timer_type
    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type, mv_isTestingMode

    implicit none

    character(*)    , parameter     :: MODULE_NAME = "Test_mod"

    real(RK)                        :: mv_timeTotal
    integer(IK)                     :: mv_testCounterOld
    integer(IK)                     :: mv_testCounter
    integer(IK)                     :: mv_nfailMax
    integer(IK)                     :: mv_nfail
    integer(IK)                     :: mv_npass
    type(CharVec_type), allocatable :: mv_FailedTestFuncName(:)

    character(:), allocatable       :: mc_passedString
    character(:), allocatable       :: mc_failedString

    type(Image_type)                :: mv_Image

    type :: Test_type
#if defined DEBUG_ENABLED || defined TESTING_ENABLD || defined CODECOV_ENABLD
        logical                     :: isVerboseMode = .true.
#else
        logical                     :: isVerboseMode = .false.
#endif
        integer(IK)                 :: outputUnit
        character(:), allocatable   :: moduleName
        character(:), allocatable   :: outDir
        character(:), allocatable   :: inDir
        type(Timer_type)            :: Timer
        type(Image_type)            :: Image
        type(Err_type)              :: Err
    contains
        procedure, pass :: run => runTest
        procedure, pass :: finalize => finalizeTest
    end type Test_type

    interface Test_type
        module procedure constructTest
    end interface Test_type

    abstract interface
        function test_func_proc() result(assertion)
            implicit none
            logical :: assertion
        end function test_func_proc
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Initialize the test for the first time. This function must be called
    !> at the beginning of the main test file, before any testing is done.
    subroutine setup()

        use Decoration_mod, only: style

        implicit none

        ! set the mode of error handling to testing mode

        mv_isTestingMode = .true.

        ! set up image counts

        call mv_Image%query()

        mc_passedString = style("passed", "bright", "green")
        mc_failedString = style("FAILED", "bright", "red")

        mv_timeTotal = 0._RK
        mv_testCounterOld = 0_IK
        mv_testCounter = 0_IK
        mv_nfailMax = 100_IK
        mv_npass = 0_IK
        mv_nfail = 0_IK

        ! preallocate the names of failed tests

        if (allocated(mv_FailedTestFuncName)) deallocate(mv_FailedTestFuncName)
        allocate(mv_FailedTestFuncName(mv_nfailMax))

        write(output_unit, "(*(g0,:,' '))")

    end subroutine setup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine finalize()

        use Constants_mod, only: RK, IK
        use Decoration_mod, only: style
        use String_mod, only: num2str
        use Err_mod, only: abort
        implicit none
        real(RK)                    :: percentageTestPassed
        integer(IK)                 :: ntotal
        character(:), allocatable   :: color
        character(:), allocatable   :: msg
        integer(IK)                 :: i

        ntotal = mv_npass + mv_nfail

        percentageTestPassed = 0._RK
        if (ntotal/=0_IK) percentageTestPassed = 100_IK * mv_npass / real(ntotal, kind=RK)

        if (mv_Image%isFirst) then

            msg = style( num2str(ntotal)//" tests performed. ", "bright", "yellow")

            color = "green"; if (nint(percentageTestPassed,kind=IK)==0_IK) color = "red"
            msg = msg // style( num2str(percentageTestPassed,"(f0.2)")//"% of "//num2str(ntotal)//" tests passed. ", "bright", color)

            color = "green"; if (mv_nfail>0_IK) color = "red"
            msg = msg // style( num2str(mv_nfail,"(g0)")//" tests failed. ", "bright", color)

            msg = msg // style( "The total elapsed wall-clock time: " // num2str(mv_timeTotal,"(f0.6)")//" seconds.", "bright", "yellow" )

            write(output_unit, "(*(g0,:,' '))")
            write(output_unit, "(*(g0,:,' '))") msg
            write(output_unit, "(*(g0,:,' '))")

            if (mv_nfail>0_IK) then

                write(output_unit, "(*(g0,:,' '))")
                write(output_unit, "(*(g0,:,' '))") style("The following tests FAILED:", "bright", "red")
                write(output_unit, "(*(g0,:,' '))")
                do i = 1, mv_nfail
                    write(output_unit, "(*(g0,:,' '))") style("    "//num2str(i,"(g0)")//" - "//mv_FailedTestFuncName(i)%record, "bright", "red")
                end do
                write(output_unit, "(*(g0,:,' '))")
                write(output_unit, "(*(g0,:,' '))") style("Errors occurred while running the ParaMonte tests. Please report this issue at:", "bright", "red")
                write(output_unit, "(*(g0,:,' '))")
                write(output_unit, "(*(g0,:,' '))") style("    https://github.com/cdslaborg/paramonte/issues", "bright", "cyan")
#if ! defined DEBUG_ENABLED
                write(output_unit, "(*(g0,:,' '))")
                write(output_unit, "(*(g0,:,' '))") style("To get more information about the failure, you may also compile and run the test in debug mode.", "bright", "red")
#endif
                write(output_unit, "(*(g0,:,' '))")

            end if

        end if

        if (mv_nfail>0_IK) then
            call abort()
        else
            call mv_Image%finalize()
        end if

    end subroutine finalize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructTest(moduleName, inDir, outDir) result(Test)

        use Decoration_mod, only: style
        use Constants_mod, only: NLC
        use String_mod, only: num2str
        use Path_mod, only: Path_type
        use Err_mod, only: abort

        implicit none

        character(*), intent(in)            :: moduleName
        character(*), intent(in), optional  :: inDir, outDir
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME//"@constructTest()"
        type(Test_type)                     :: Test
        type(Path_type)                     :: Path

        ! set up image counts

        Test%Image          = mv_Image
        Test%outputUnit     = output_unit
        Test%moduleName     = moduleName

        ! set up the timer

        Test%Timer = Timer_type(Test%Err)
        if (Test%Err%occurred) then
            Test%Err%msg = PROCEDURE_NAME // ": Error occurred while setting up the test timer."//NLC// Test%Err%msg
            call abort( Err = Test%Err, prefix = "ParaMonteUnitTesting", newline = NLC, outputUnit = Test%outputUnit )
            return
        end if

        ! set up the tests input path

        if (present(inDir)) then
            Test%inDir = inDir
        else
            Test%inDir = "./input"
        end if
        Path = Path_type(Test%inDir)
        if (Path%Err%occurred) then
            call abort( Err = Path%Err &
                      , newline = "\n" &
                      , outputUnit = Test%outputUnit &
                      )
        else
            Test%inDir = Path%modified
        end if

        ! set up the tests output path

        if (present(outDir)) then
            Test%outDir = outDir
        else
            Test%outDir = "./output"
        end if
        Path = Path_type(Test%outDir)
        if (Path%Err%occurred) then
            call abort( Err = Path%Err &
                     !, prefix = "FATAL: " &
                      , newline = "\n" &
                      , outputUnit = Test%outputUnit &
                      )
        else
            Test%outDir = Path%modified
            if (Test%Image%isFirst) then
                Path%Err = Path%mkdir(Test%outDir)
                if (Path%Err%occurred) then
                    call abort  ( Err = Path%Err &
                                !, prefix = "FATAL: " &
                                , newline = "\n" &
                                , outputUnit = Test%outputUnit &
                                )
                end if
            end if
        end if

        call mv_Image%sync()

    end function constructTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure must be called by all MPI/CAF processors.
    subroutine runTest(Test, test_func, test_func_name)

        use String_mod, only: String_type, num2str, isInteger, replaceStr, padString
        use Err_mod, only: abort
        implicit none
        class(Test_type), intent(inout) :: Test
        procedure(test_func_proc)       :: test_func
        character(*), intent(in)        :: test_func_name
        character(:), allocatable       :: counterStr, funcName
        type(String_type)               :: Message
        logical                         :: assertion
        integer                         :: imageID

        funcName = test_func_name(6:) ! remove "test_"
        Message%Parts = Message%split(string = funcName, delim = "_", npart = Message%nPart)
        if ( Message%nPart==1_IK .or. .not. isInteger(Message%Parts(Message%nPart)%record) ) then
            counterStr = "#1"
        else
            funcName = funcName(1:len(funcName)-len("_"//Message%Parts(Message%nPart)%record))
            counterStr = "#"//Message%Parts(Message%nPart)%record
        !else
        !    Test%Err%occurred = .true.
        !    Test%Err%msg = "Invalid function name passed to runTest() method of the Test object: " // test_func_name
        !    call abort  ( Err = Test%Err &
        !                , newline = "\n" &
        !                , outputUnit = Test%outputUnit &
        !                )
        end if

        call Test%Timer%toc()
        assertion = test_func()
        call Test%Timer%toc()

        if (assertion) then

            mv_npass = mv_npass + 1
            Message%value = mc_passedString

        else

            Message%value = mc_failedString
            mv_nfail = mv_nfail + 1

            ! set up the allocation for the names of failed tests

            if (mv_nfail==mv_nfailMax) then
                block
                    type(CharVec_type), allocatable :: mv_FailedTestFuncNameTemp(:)
                    mv_nfailMax = 2 * mv_nfailMax
                    allocate(mv_FailedTestFuncNameTemp(mv_nfailMax))
                    mv_FailedTestFuncNameTemp(1:mv_nfail) = mv_FailedTestFuncName
                    call move_alloc(from=mv_FailedTestFuncNameTemp, to=mv_FailedTestFuncName)
                end block
            end if

            mv_FailedTestFuncName(mv_nfail)%record = "Test_"//Test%moduleName(2:)//"@"//test_func_name

        end if

        mv_testCounter = mv_testCounter + 1

        do imageID = 1, Test%Image%count

            if (imageID==Test%Image%id) then
                write(Test%outputUnit,"(*(g0,:,' '))") "["//adjustr(num2str(mv_testCounter,minLen=4_IK))//"] testing" &
                                                     , padString( Test%moduleName//"@"//funcName//" "//counterStr//" ","." , 79_IK ) &
                                                     , Message%value, "in", adjustr(num2str(Test%Timer%Time%delta,"(f0.4)",8_IK)), "seconds on image "//num2str(Test%Image%id)

            end if

#if defined CAF_ENABLED || MPI_ENABLED
            call execute_command_line(" ", cmdstat = Test%Err%stat)
            flush(output_unit)
            call Test%Image%sync()
#endif

        end do

        block
            use System_mod, only: sleep
            use Err_mod, only: Err_type
            type(Err_type) :: Err
            call sleep(0.02_RK, Err)
        end block

    end subroutine runTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine finalizeTest(Test)
        use Decoration_mod, only: style
        use String_mod, only: num2str, padString
        implicit none
        class(Test_type), intent(inout) :: Test
        if (Test%Image%isFirst) then
            Test%Err%msg = "["//adjustr(num2str(mv_testCounter-mv_testCounterOld,minLen=4_IK))//"] testing " &
                         //padString(Test%moduleName//" ",".",79)//" isdone in " &
                         //adjustr(num2str(Test%Timer%Time%total,"(f0.4)",8_IK))//" seconds on image "//num2str(Test%Image%id)
            write(Test%outputUnit, "(*(g0,:,' '))") style(Test%Err%msg, "bright", "yellow")
        end if
        mv_testCounterOld = mv_testCounter
        if (allocated(Test%moduleName)) deallocate(Test%moduleName)
        mv_timeTotal = mv_timeTotal + Test%Timer%Time%total
        call mv_Image%sync()
    end subroutine finalizeTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the density function value of the uncorrelated Multivariate Normal distribution for the given input
    !> location specified by the vector `Point` of length `ndim`.
#if defined CFI_ENABLED
    function getLogFuncMVN(ndim,Point) result(logFunc) bind(C)
#else
    function getLogFuncMVN(ndim,Point) result(logFunc)
#endif
        ! This function returns the probability density function of the standard multivariate normal distribution of ndim dimensions.
        !use Statistics_mod, only: getLogProbMixMVN
        use Constants_mod, only : IK, RK
        implicit none

        !! Standard MultiVariate Normal (SMVN) specifications: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
        !real(RK), parameter :: LOG_SMVN_COEF = NDIM*log(1._RK/sqrt(2._RK*acos(-1._RK))) ! log(1/sqrt(2*Pi)^ndim)

#if defined CFI_ENABLED
        integer(IK), intent(in), value  :: ndim
#else
        integer(IK), intent(in)         :: ndim
#endif
        real(RK), intent(in)            :: Point(ndim)
        real(RK)                        :: logFunc

        !block
        !    use System_mod, only: sleep
        !    use Err_mod, only: Err_type
        !    type(Err_type) :: Err
        !    call sleep(seconds=5000.e-6_RK,Err=Err)
        !end block

        !block
        !    real(RK), allocatable :: unifrnd(:,:)
        !    allocate(unifrnd(200,20))
        !    call random_number(unifrnd)
        !    logFunc = sum(unifrnd) - 0.5_RK * sum(Point**2) - sum(unifrnd)
        !    deallocate(unifrnd)
        !end block

        logFunc = -sum(Point**2)

       !block
       !    integer(IK), parameter :: nmode = 2_IK
       !    real(RK) :: LogAmplitude(nmode), MeanVec(nmode), InvCovMat(nmode), LogSqrtDetInvCovMat(nmode)
       !    LogAmplitude        = [1._RK, 1._RK]
       !    MeanVec             = [0._RK, 7._RK]
       !    InvCovMat           = [1._RK,1._RK]
       !    LogSqrtDetInvCovMat = [1._RK,1._RK]
       !    logFunc = getLogProbGausMix ( nmode = 2_IK &
       !                                , nd = 1_IK &
       !                                , np = 1_IK &
       !                                , LogAmplitude = LogAmplitude &
       !                                , MeanVec = MeanVec &
       !                                , InvCovMat = InvCovMat &
       !                                , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
       !                                , Point = Point(1) &
       !                                )
       !end block
    end function getLogFuncMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_mod ! LCOV_EXCL_LINE
! LCOV_EXCL_STOP