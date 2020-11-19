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
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Test_mod

    use, intrinsic  :: iso_fortran_env, only: output_unit
    use JaggedArray_mod, only: CharVec_type
    use Timer_mod, only: Timer_type
    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*)    , parameter     :: MODULE_NAME = "Test_mod"

    real(RK)                        :: mv_timeTotal
    logical                         :: mv_imageIsFirst
    integer(IK)                     :: mv_imageCount
    integer(IK)                     :: mv_imageID
    integer(IK)                     :: mv_nfailMax = 1000_IK
    integer(IK)                     :: mv_nfail = 0_IK
    integer(IK)                     :: mv_npass = 0_IK
    type(CharVec_type), allocatable :: mv_FailedTestFuncName(:)

    character(:), allocatable       :: mc_passedString
    character(:), allocatable       :: mc_failedString

    type                            :: Image_type
        integer(IK)                 :: id, count
        logical                     :: isFirst = .false.
        logical                     :: isMaster = .false.
        logical                     :: isNotFirst = .false.
        logical                     :: isNotMaster = .false.
        character(:), allocatable   :: name
    end type Image_type

    type :: Test_type
#if defined DBG_ENABLED
        logical                     :: isDebugMode = .true.
#else
        logical                     :: isDebugMode = .false.
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

        ! set up image counts

#if defined CAF_ENABLED
        mv_imageID = this_image()
        mv_imageCount = num_images()
#elif defined MPI_ENABLED
        block
            use mpi
            integer(IK) :: ierrMPI
            logical     :: isInitialized
            call mpi_initialized( isInitialized, ierrMPI )
            if (.not. isInitialized) call mpi_init(ierrMPI)
            call mpi_comm_rank(mpi_comm_world, mv_imageID, ierrMPI)
            call mpi_comm_size(mpi_comm_world, mv_imageCount, ierrMPI)
            mv_imageID = mv_imageID + 1_IK ! make the ranks consistent with Fortran coarray indexing conventions
        end block
#else
        mv_imageID = 1_IK
        mv_imageCount = 1_IK
#endif
        mv_imageIsFirst = mv_imageID == 1_IK

        ! preallocate the names of failed tests

        if (allocated(mv_FailedTestFuncName)) deallocate(mv_FailedTestFuncName)
        allocate(mv_FailedTestFuncName(mv_nfailMax))

        mc_passedString = style("passed", "bright", "green")
        mc_failedString = style("failed", "bright", "red")

        mv_timeTotal = 0._RK

        write(output_unit, "(*(g0,:,' '))")

    end subroutine setup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine finalize()

        use Constants_mod, only: RK, IK
        use Decoration_mod, only: style
        use String_mod, only: num2str
        implicit none
        integer(IK)                 :: percentageTestPassed
        real(RK)                    :: ntotal
        character(:), allocatable   :: color
        character(:), allocatable   :: msg
        integer(IK)                 :: i

        ntotal = mv_npass + mv_nfail

        percentageTestPassed = nint(100_IK * mv_npass / ntotal)

        if (mv_imageIsFirst) then

            color = "green"; if (percentageTestPassed==0_IK) color = "red"
            msg = style( num2str(percentageTestPassed,"(g0)")//"% tests passed. ", "bright", color)

            color = "green"; if (mv_nfail>0_IK) color = "red"
            msg = msg // style( num2str(mv_nfail,"(g0)")//" tests failed. ", "bright", color)

            msg = msg // style( "Total wall clock time: " // num2str(mv_timeTotal,"(f0.6)")//" seconds.", "bright", "yellow" )

            write(output_unit, "(*(g0,:,' '))")
            write(output_unit, "(*(g0,:,' '))") msg

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
                write(output_unit, "(*(g0,:,' '))")

                error stop

            end if

        end if

#if defined MPI_ENABLED
        call finalizeMPI()
#elif defined CAF_ENABLED
        sync all
#endif

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

        ! set up image counts

        Test%Image%id       = mv_imageID
        Test%Image%count    = mv_imageCount
        Test%Image%isFirst  = mv_imageIsFirst
        Test%Image%name     = "@process(" // num2str(Test%Image%id) // ")"
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
            Test%inDir = "./input/"
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
            Test%outDir = "./output/"
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
            Path%Err = Path%mkdir(Test%outDir)
            if (Path%Err%occurred) then
                call abort  ( Err = Path%Err &
                            !, prefix = "FATAL: " &
                            , newline = "\n" &
                            , outputUnit = Test%outputUnit &
                            )
            end if
        end if

        ! announce the module being tested

        if (Test%Image%isFirst) then
            Test%Err%msg = "Testing "//Test%moduleName//" ..."
            write(Test%outputUnit, "(*(g0,:,' '))") style(Test%Err%msg, "bright", "yellow")
#if defined CAF_ENABLED
            ! sync images
            sync images(*)
        else
            sync images(1)
#endif
        end if

#if defined MPI_ENABLED
        block
            use mpi
            integer(IK) :: ierrMPI
            call mpi_barrier(mpi_comm_world,ierrMPI)
        end block
#endif

    end function constructTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine runTest(Test, test_func, funcName)

        use String_mod, only: String_type, num2str
        use Err_mod, only: abort
        implicit none
        class(Test_type), intent(inout) :: Test
        procedure(test_func_proc)       :: test_func
        character(*), intent(in)        :: funcName
        character(:), allocatable       :: counterStr
        type(String_type)               :: Message
        logical                         :: assertion

        Message%Parts = Message%splitStr(string=funcName(6:), delimiter="_", nPart=Message%nPart)
        if (Message%nPart==2_IK) then
            counterStr = "#"//Message%Parts(2)%record
        elseif (Message%nPart==1_IK) then
            counterStr = "#1"
        else
            Test%Err%occurred = .true.
            Test%Err%msg = "Invalid function name passed to runTest() method of the Test object: " // funcName
            call abort  ( Err = Test%Err &
                        , newline = "\n" &
                        , outputUnit = Test%outputUnit &
                        )
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
            mv_FailedTestFuncName(mv_nfail)%record = "Test_"//Test%moduleName(2:)//"@"//funcName
        end if

        write(Test%outputUnit,"(*(g0,:,' '))") "testing", Test%moduleName//"@"//Message%Parts(1)%record, counterStr &
                                             , "...", Message%value, "in", num2str(Test%Timer%Time%delta,"(f0.6)"), "seconds"

    end subroutine runTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine finalizeTest(Test)
        use Decoration_mod, only: style
        use String_mod, only: num2str
        implicit none
        class(Test_type), intent(inout) :: Test
        if (Test%Image%isFirst) then
            Test%Err%msg = "Testing "//Test%moduleName//" completed in "//num2str(Test%Timer%Time%total,"(f0.6)")//" seconds."
            write(Test%outputUnit, "(*(g0,:,' '))") style(Test%Err%msg, "bright", "yellow")
        end if
        if (allocated(Test%moduleName)) deallocate(Test%moduleName)
        mv_timeTotal = mv_timeTotal + Test%Timer%Time%total
        call sync()
    end subroutine finalizeTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine sync()
        implicit none
#if defined MPI_ENABLED
        block
            use mpi
            integer(IK) :: ierrMPI
            logical     :: isFinalized
            call mpi_finalized( isFinalized, ierrMPI )
            if (.not. isFinalized) then
                call mpi_barrier(mpi_comm_world,ierrMPI)
            end if
        end block
#elif defined CAF_ENABLED
        sync all
#endif
    end subroutine sync

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined MPI_ENABLED
    subroutine finalizeMPI()
        use mpi
        integer(IK) :: ierrMPI
        logical     :: isFinalized
        call mpi_finalized( isFinalized, ierrMPI )
        if (.not. isFinalized) then
            call mpi_barrier(mpi_comm_world,ierrMPI)
            call mpi_finalize(ierrMPI)
        end if
    end subroutine finalizeMPI
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_mod
