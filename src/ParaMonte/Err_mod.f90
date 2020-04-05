!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of the ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module Err_mod

    character(*), parameter :: MODULE_NAME = "@Err_mod"

    logical, parameter :: ERR_HANDLING_REQUESTED = .false.

    type :: Err_type
        logical                     :: occurred = .false.
        integer                     :: stat     = -huge(0)
        integer                     :: statNull = -huge(0)
        character(:), allocatable   :: msg
    end type Err_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine abort(Err,prefix,newline,outputUnit)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: abort
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: NLC
        use Decoration_mod, only: write
        implicit none
        type(Err_type), intent(in)          :: Err
        character(*), intent(in), optional  :: prefix, newline
        integer     , intent(in), optional  :: outputUnit

        character(:), allocatable           :: pfx, msg, nlstr
        character(63)                       :: dummyChar1, imageChar !, dummyChar2

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

        if (Err%stat==Err%statNull) then    ! it's a null error code, ignore it and do not report the error code
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
                call write(outputUnit,1,0,1, pfx // " - For further help, contact Amir Shahmoradi via:" )
                call write(outputUnit,0,0,1, pfx // " - a.shahmoradi@gmail.com" )
                call write(outputUnit,0,0,1, pfx // " - shahmoradi@utexas.edu" )
                call write(outputUnit,0,0,1, pfx // " - cdslab.org/ParaMonte/" )
                call write(outputUnit,1,2,1, pfx // " - Gracefully Exiting on image " // trim(adjustl(imageChar)) // "." )

            end if
        end if

        if (outputUnit/=output_unit) then
            ! notify the user on screen too
            call write(output_unit,1,0,1, pfx // " - FATAL: Runtime error occurred." )
            call write(output_unit,0,0,1, pfx // " - FATAL: For more information please see the report file." )
            call write(output_unit,0,2,1, pfx // " - FATAL: Gracefully Exiting on image " // trim(adjustl(imageChar)) // "." )
        end if

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
        error stop

    end subroutine abort

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine warn(msg,prefix,newline,outputUnit)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: warn
#endif
        implicit none
        character(*), intent(in)            :: msg
        character(*), intent(in), optional  :: prefix, newline
        integer     , intent(in), optional  :: outputUnit
        if (present(prefix)) then
            call informUser(msg,prefix // " - WARNING: ",newline,outputUnit)
        else
            call informUser(msg," - WARNING: ",newline,outputUnit)
        end if
    end subroutine warn

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Err_mod