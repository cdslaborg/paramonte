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

    use Constants_mod, only: IK
    use Decoration_mod, only: writeDecoratedText, write
    implicit none

    type                            :: Image_type
        integer(IK)                 :: id, count
        logical                     :: isFirst = .false.
        logical                     :: isMaster = .false.
        logical                     :: isNotFirst = .false.
        logical                     :: isNotMaster = .false.
        character(:), allocatable   :: name
    end type Image_type

    type :: Test_type
        logical                     :: assertion = .false.
#if defined DBG_ENABLED
        logical                     :: isDebugMode = .true.
#else
        logical                     :: isDebugMode = .false.
#endif
        integer(IK)                 :: outputUnit
        character(:), allocatable   :: moduleName
        character(:), allocatable   :: testName
        character(:), allocatable   :: outDir
        character(:), allocatable   :: inDir
        type(Image_type)            :: Image
    contains
        procedure, pass :: testing
        procedure, pass :: skipping
        procedure, pass :: verify
        procedure, pass :: checkForErr
        procedure, pass :: finalize
    end type Test_type

    interface Test_type
        module procedure constructTest
    end interface Test_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructTest(moduleName) result(Test)
        use, intrinsic  :: iso_fortran_env, only: output_unit
        use String_mod, only: num2str
        use Path_mod, only: Path_type
        use Err_mod, only: abort
        implicit none
        character(*), intent(in)    :: moduleName
        type(Test_type)             :: Test
        type(Path_type)             :: Path
#if defined CAF_ENABLED
        Test%Image%id        = this_image()
        Test%Image%count     = num_images()
#elif defined MPI_ENABLED
        block
            use mpi
            integer(IK) :: ierrMPI
            logical     :: isInitialized
            call mpi_initialized( isInitialized, ierrMPI )
            if (.not. isInitialized) call mpi_init(ierrMPI)
            call mpi_comm_rank(mpi_comm_world, Test%Image%id, ierrMPI)
            call mpi_comm_size(mpi_comm_world, Test%Image%count, ierrMPI)
            Test%Image%id = Test%Image%id + 1_IK ! make the ranks consistent with Fortran coarray indexing conventions
        end block
#else
        Test%Image%id        = 1_IK
        Test%Image%count     = 1_IK
#endif
        Test%Image%name     = "@process(" // num2str(Test%Image%id) // ")"
        Test%Image%isFirst  = Test%Image%id == 1_IK
        Test%assertion      = .false.
        Test%outputUnit     = output_unit
        Test%moduleName     = moduleName

        Path = Path_type("./input/")
        if (Path%Err%occurred) then
            call abort( Err = Path%Err &
                     !, prefix = "FATAL: " &
                      , newline = "\n" &
                      , outputUnit = Test%outputUnit &
                      )
        else
            Test%inDir = Path%modified
        end if

        Path = Path_type("./out/")
        if (Path%Err%occurred) then
            call abort( Err = Path%Err &
                     !, prefix = "FATAL: " &
                      , newline = "\n" &
                      , outputUnit = Test%outputUnit &
                      )
        else
            Test%outDir = Path%modified
        end if

        if (Test%Image%isFirst) then
            call writeDecoratedText( text = "\nTesting module " // moduleName // "\n" &
                                   , marginTop = 2, marginBot = 1, outputUnit = Test%outputUnit, newLine = "\n" )
            Path%Err = Path%mkdir(Test%outDir)
            if (Path%Err%occurred) then
                call abort  ( Err = Path%Err &
                            !, prefix = "FATAL: " &
                            , newline = "\n" &
                            , outputUnit = Test%outputUnit &
                            )
            end if
#if defined CAF_ENABLED
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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine testing(Test,testName)
        implicit none
        class(Test_type), intent(inout) :: Test
        character(*), intent(in) :: testName
        Test%testName = testName
        if (Test%Image%isFirst) call write(Test%outputUnit,2,0,1, "Testing " // Test%testName // " ..." )
    end subroutine testing

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine verify(Test,message)
        implicit none
        class(Test_type), intent(inout)     :: Test
        character(*), intent(in), optional  :: message
        if (Test%assertion) then
            Test%assertion = .false.
#if defined DBG_ENABLED
            call write(Test%outputUnit,0,0,1, "Successful " // Test%Image%name // "." )
#endif
        else
            call write(Test%outputUnit,1,1,1, "Fatally Failed " // Test%Image%name // "." )
            if (present(message)) call write(Test%outputUnit,1,0,1, message )
            call write(Test%outputUnit)
            error stop
        end if
    end subroutine verify

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForErr(Test,Err)
        use Err_mod, only: Err_type, abort
        use String_mod, only: num2str
        implicit none
        class(Test_type), intent(inout) :: Test
        type(Err_type), intent(inout)   :: Err
        if (Err%occurred) then
            call abort( Err = Err &
                      , newline = "\n" &
                      , outputUnit = Test%outputUnit &
                      )
        end if
    end subroutine checkForErr

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine skipping(Test)
        implicit none
        class(Test_type), intent(inout) :: Test
        if (Test%Image%isFirst) then
            call write(Test%outputUnit,0,0,1, "TEST WARNING: No result-verification test available for " // Test%testName // ". skipping... " )
        end if
        Test%assertion = .false.
    end subroutine skipping

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine finalize(Test)
        class(Test_type), intent(inout) :: Test
        if (allocated(Test%moduleName)) deallocate(Test%moduleName)
        if (allocated(Test%testName)) deallocate(Test%testName)
#if defined MPI_ENABLED
        block
            use mpi
            integer(IK) :: ierrMPI
            logical     :: isFinalized
            call mpi_finalized( isFinalized, ierrMPI )
            if (.not. isFinalized) then
                call mpi_barrier(mpi_comm_world,ierrMPI)
                !call mpi_finalize(ierrMPI)
            end if
        end block
#elif defined CAF_ENABLED
        sync all
#endif
    end subroutine finalize

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_mod
