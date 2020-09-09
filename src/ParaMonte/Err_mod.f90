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

module Err_mod

    character(*), parameter :: MODULE_NAME = "@Err_mod"

    logical     , parameter :: ERR_HANDLING_REQUESTED = .false.

#if (defined MATLAB_ENABLED || defined PYTHON_ENABLED) && !defined CAF_ENABLED && !defined MPI_ENABLED
    logical     , parameter :: SOFT_EXIT_ENABLED = .true.
#else
    logical     , parameter :: SOFT_EXIT_ENABLED = .false.
#endif

    type :: Err_type
        logical                     :: occurred = .false.
        integer                     :: stat     = -huge(0)
        integer                     :: statNull = -huge(0)
        character(:), allocatable   :: msg
    end type Err_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine abort(Err,prefix,newline,outputUnit,returnEnabled)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: abort
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Decoration_mod, only: write
        use Constants_mod, only: NLC
        implicit none
        type(Err_type), intent(in)          :: Err
        character(*), intent(in), optional  :: prefix, newline
        integer     , intent(in), optional  :: outputUnit
        logical     , intent(in), optional  :: returnEnabled

        logical                             :: returnEnabledDefault
        character(:), allocatable           :: pfx, msg, nlstr
        character(63)                       :: dummyChar1, imageChar !, dummyChar2

        if (present(returnEnabled)) then
            returnEnabledDefault = returnEnabled
        else
            returnEnabledDefault = SOFT_EXIT_ENABLED
        end if

#if defined CAF_ENABLED
        write(imageChar ,"(g0)") this_image()
#elif defined MPI_ENABLED
        block
            use mpi
            integer :: imageID, ierrMPI
            call mpi_comm_rank(mpi_comm_world, imageID, ierrMPI)
            write(imageChar ,"(g0)") imageID + 1
        end block
#else
        imageChar =  "1"
#endif

        if (present(newline)) then
            nlstr = newline
        else
            nlstr = NLC
        end if

        if (Err%stat==Err%statNull) then    ! it is a null error code, ignore it and do not report the error code
            msg = Err%msg
        else
            write(dummyChar1,"(g0)") Err%stat
           !write(dummyChar2,"(g0)") Err%statNull
           !msg =   Err%msg // nlstr // "Error Code: " // trim(adjustl(dummyChar1)) // ". Null Error Code: " // trim(adjustl(dummyChar2)) // "."
            msg =   Err%msg // nlstr // "Error Code: " // trim(adjustl(dummyChar1)) // "."
        end if

        if (present(prefix)) then
            call informUser(msg,prefix//" - FATAL: ",nlstr,outputUnit)
            pfx = prefix
        else
            call informUser(msg," - ",nlstr,outputUnit)
            pfx = ""
        end if

        if (present(outputUnit)) then
            if (outputUnit/=output_unit) then
                call write(outputUnit,1,0,1, pfx // " - Please Correct the error(s) and rerun the simulation." )
                call write(outputUnit,0,0,1, pfx // " - If the cause of the error cannot be diagnosed, please report it at:" )
                call write(outputUnit,0,0,1, pfx // " -" )
                call write(outputUnit,0,0,1, pfx // " -     https://github.com/cdslaborg/paramonte/issues" )
                call write(outputUnit,0,0,1, pfx // " -" )
                call write(outputUnit,0,2,1, pfx // " - Gracefully Exiting on image " // trim(adjustl(imageChar)) // "." )
            end if
        end if

        ! notify the user on screen too

        call write(output_unit,1,0,1, pfx // " - FATAL: Runtime error occurred." )
        call write(output_unit,0,0,1, pfx // " - FATAL: For more information, see the output '*_report.txt' file (if generated)." )
        call write(output_unit,0,2,1, pfx // " - FATAL: Gracefully Exiting on image " // trim(adjustl(imageChar)) // "." )

        flush(output_unit) ! call execute_command_line(" ")
        flush(outputUnit)

        ! wait for one second:
        block
            use Constants_mod, only: RK
            use, intrinsic  :: iso_fortran_env, only: int64
            integer(int64)  :: countOld, countNew, countMax
            real(RK)        :: countRate
            call system_clock( count=countOld, count_rate=countRate, count_max=countMax )
            if (countOld/=-huge(0_int64) .and. countRate/=0._RK .and. countMax==0_int64) then
                loopWait: do
                    call system_clock( count=countNew )
                    if (countNew==countMax) then
                        if (returnEnabledDefault) return
                        error stop
                    elseif ( real(countNew-countOld,kind=RK) / countRate >= 2._RK ) then
                        exit loopWait
                    end if
                    cycle
                end do loopWait
            end if
        end block

#if defined MPI_ENABLED
        block
            use mpi
            integer :: ierrMPI, errcode
            errcode = 1; call mpi_abort(mpi_comm_world, errcode, ierrMPI)
        end block
#endif
        if (returnEnabledDefault) return
        error stop

    end subroutine abort

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine warn(msg,prefix,newline,outputUnit,marginTop,marginBot)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: warn
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in)            :: msg
        character(*), intent(in), optional  :: prefix, newline
        integer(IK) , intent(in), optional  :: outputUnit,marginTop,marginBot
        if (present(prefix)) then
            call informUser ( msg           = msg                       &
                            , prefix        = prefix // " - WARNING: "  &
                            , newline       = newline                   &
                            , outputUnit    = outputUnit                &
                            , marginTop     = marginTop                 &
                            , marginBot     = marginBot                 &
                            )
        else
            call informUser ( msg           = msg                       &
                            , prefix        = " - WARNING: "            &
                            , newline       = newline                   &
                            , outputUnit    = outputUnit                &
                            , marginTop     = marginTop                 &
                            , marginBot     = marginBot                 &
                            )
        end if
    end subroutine warn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine note(msg,prefix,newline,outputUnit,marginTop,marginBot)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: note
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in)            :: msg
        character(*), intent(in), optional  :: prefix, newline
        integer(IK) , intent(in), optional  :: outputUnit,marginTop,marginBot
        if (present(prefix)) then
            call informUser ( msg           = msg                   &
                            , prefix        = prefix // " - NOTE: " &
                            , newline       = newline               &
                            , outputUnit    = outputUnit            &
                            , marginTop     = marginTop             &
                            , marginBot     = marginBot             &
                            )
        else
            call informUser ( msg           = msg                   &
                            , prefix        = " - NOTE: "           &
                            , newline       = newline               &
                            , outputUnit    = outputUnit            &
                            , marginTop     = marginTop             &
                            , marginBot     = marginBot             &
                            )
        end if
    end subroutine note

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine informUser(msg,prefix,newline,outputUnit,wrapSplit,wrapWidth,marginTop,marginBot)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: informUser
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Decoration_mod, only: write, getListOfLines, wrapText
        use JaggedArray_mod, only: CharVec_type
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in)            :: msg
        character(*), intent(in), optional  :: prefix, newline, wrapSplit
        integer(IK) , intent(in), optional  :: outputUnit, wrapWidth, marginTop, marginBot

        integer(IK)                         :: stdout, sizeList, sizeListJustified
        integer(IK)                         :: irecord, ijustified, width
        integer(IK)                         :: padTop, padBot, padTopCurrent, padBotCurrent
        character(:), allocatable           :: pfx, split
        type(CharVec_type), allocatable     :: List(:), ListJustified(:)

        if (present(outputUnit)) then
            stdout = outputUnit
        else
            stdout = output_unit
        end if
        if (present(prefix)) then
            pfx = prefix
        else
            pfx = ""
        end if
        if (present(wrapSplit)) then
            split = wrapSplit
        else
            split = " "
        end if
        if (present(wrapWidth)) then
            width = wrapWidth
        else
            width = 100_IK
        end if
        if (present(marginTop)) then
            padTop = marginTop
        else
            padTop = 1_IK
        end if
        if (present(marginBot)) then
            padBot = marginBot
        else
            padBot = 1_IK
        end if

        List = getListOfLines(string=msg,delimiter=newline)
        sizeList = size(List)
        do irecord = 1, sizeList
            ListJustified = wrapText( string    = List(irecord)%record  &
                                    , width     = width                 &
                                    , split     = split                 &
                                    , pad       = " "                   &
                                    )
            sizeListJustified = size(ListJustified)
            do ijustified = 1, sizeListJustified
                padTopCurrent = 0_IK
                padBotCurrent = 0_IK
                if (irecord==1 .and. ijustified==1_IK) padTopCurrent = padTop ! the very first line
                if (irecord==sizeList .and. ijustified==sizeListJustified) padBotCurrent = padBot ! the very last line
                call write(stdout, padTopCurrent, padBotCurrent, 1_IK, pfx // ListJustified(ijustified)%record )
            end do
        end do
        if (.not.present(marginBot)) call write(stdout)

    end subroutine informUser

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Err_mod