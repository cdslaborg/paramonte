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
!>  This file contains procedure implementations of [pm_except](@ref pm_except).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%
#if     getFine_ENABLED
        !%%%%%%%%%%%%%%

        str = getFile(file)//getLine(line)

        !%%%%%%%%%%%%%%
#elif   getFile_ENABLED
        !%%%%%%%%%%%%%%

        str = SKC_'@file('//file//SKC_")"

        !%%%%%%%%%%%%%%
#elif   getLine_ENABLED
        !%%%%%%%%%%%%%%

        str = repeat(SK_" ", range(0_IKC) + 3) ! sign. 2 is essential. extra 1 is cautionary.
        write(str, "(I0)") line
        str = '@line('//trim(str)//")"

        !%%%%%%%%%%%%%%%%%%
#elif   setAsserted_ENABLED
        !%%%%%%%%%%%%%%%%%%

        character(1, SK), parameter :: BEL = achar(7, SK)
        character(*, SK), parameter :: NLC = new_line(SK_"a")
        if (.not. assertion) then
            if (present(msg)) then
#if             MEXPRINT_ENABLED
                call mexPrintf(msg//repeat(BEL, 3)//NLC)
#else
                write(output_unit,"(A)") msg//repeat(BEL, 3)
#endif
            end if
            if (present(renabled)) then
                if (renabled) return
            end if
            error stop "Assertion failed."
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMarked_ENABLED && Static_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_str, only: getStrWrapped

        character(:,SKC), allocatable   :: def_prefix, remark
        character(1,SKC), parameter     :: LF = new_line(SKC_"a") ! char(10, SKC)
        integer(IK)                     :: def_unit, def_tmsize, def_bmsize

        ! Set the prefix.

        if (present(prefix)) then
            def_prefix = prefix
        else
            def_prefix = SKC_" - REMARK: "
        end if

        ! Set the default.

        def_unit = int(output_unit, IK)
        if (present(unit)) def_unit = unit

        def_tmsize = 1_IK
        if (present(tmsize)) def_tmsize = tmsize

        def_bmsize = 0_IK
        if (present(bmsize)) def_bmsize = bmsize

        ! Wrap the message and write the text to the output

        remark = repeat(LF, def_tmsize)//getStrWrapped(msg, prefix = def_prefix, indent = indent, break = break, newline = newline, linefeed = LF, width = width, maxwidth = maxwidth)//repeat(LF, def_bmsize)
#if     MEXPRINT_ENABLED
        if (def_unit == output_unit) then
            call mexPrintf(remark//new_line(SKC_"a"))
            return
        end if
#endif
        write(def_unit, "(a)") remark

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setNoted_ENABLED && Static_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(*,SKC), parameter :: REMARK = SKC_" - NOTE: "
        if (present(prefix)) then
            call setMarked(msg, prefix//REMARK, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
        else
            call setMarked(msg, REMARK, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setWarned_ENABLED && Static_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(*,SKC), parameter :: REMARK = SKC_" - WARNING: "
        if (present(prefix)) then
            call setMarked(msg, prefix//REMARK, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
        else
            call setMarked(msg, REMARK, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setAborted_ENABLED && Static_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(1,SKC), parameter :: LF = new_line(SKC_"a") ! char(10, SKC)
        character(*,SKC), parameter :: REMARK = SKC_" - FATAL: "
        character(:,SKC), allocatable :: def_prefix, def_msg
        character(31,SKC) :: val2str

        if (present(prefix)) then
            def_prefix = prefix//REMARK
        else
            def_prefix = REMARK
        end if

        ! Set the error code.

        if (present(stat)) then
            write(val2str, "(2(g0))") LF//SKC_"ERROR CODE: ", stat
            def_msg = msg//trim(adjustl(val2str))
        else
            def_msg = msg
        end if

        ! Set the processor ID.

        block
            use pm_parallelism, only: getImageID
            write(val2str, "(g0)") getImageID()
        end block

        ! Report the final troubleshooting info.

        def_msg = def_msg//LF//SKC_"Please correct the error(s) and rerun the program."//LF
        if (present(help)) def_msg = def_msg//help//LF
        def_msg = def_msg//SKC_"Gracefully exiting on image/process "//trim(adjustl(val2str))//SKC_"."//LF//LF

        call setMarked(msg, def_prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, unit)

        ! Report to stdout and flush the output.

        block
            use iso_fortran_env, only: output_unit
            if (present(unit)) then
                if (unit /= output_unit) call setMarked(msg//repeat(achar(7), 3), def_prefix, indent, break, newline, width, maxwidth, tmsize, bmsize, int(output_unit, IK)) ! Set off the alarm via BEL character.
                flush(unit)
            end if
            flush(output_unit) ! call execute_command_line(" ")
        end block

        ! LCOV_EXCL_START

        ! Return or halt the program.

        block
            logical(LK) :: def_renabled
            def_renabled = SOFT_EXIT_ENABLED
            if (present(renabled)) def_renabled = renabled
            if (def_renabled) return
        end block

        ! Wait for one second before aborting the program.

        block
            use pm_kind, only: IKD
            integer(IKD) :: countOld, countNew, countMax, countRate
            call system_clock(countOld, countRate, countMax)
            if (countOld /= -huge(0_IKD) .and. countRate /= 0_IKD .and. countMax /= 0_IKD) then
                loopWait: do
                    call system_clock(countNew)
                    if (real(abs(countNew - countOld)) / real(countRate) >= 1.) exit loopWait
                end do loopWait
            end if
        end block

        ! abort.

#if     MPI_ENABLED
        block
            use mpi !mpi_f08, only: mpi_abort, mpi_comm_world
            integer :: ierrMPI
            call mpi_abort(mpi_comm_world, 1, ierrMPI)
        end block
#else
        error stop
#endif
        ! LCOV_EXCL_STOP

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   Method_ENABLED && (setMarked_ENABLED || setNoted_ENABLED || setWarned_ENABLED || setAborted_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setNoted_ENABLED
#define SET_MARKED setNoted
#elif   setMarked_ENABLED
#define SET_MARKED setMarked
#elif   setWarned_ENABLED
#define SET_MARKED setWarned
#elif   setAborted_ENABLED
#define SET_MARKED setAborted
        integer(IK), allocatable :: def_stat
        logical(LK), allocatable :: def_renabled
        character(:,SKC), allocatable :: def_help
#else
#error  "Unrecognized interface."
#endif
        integer(IK), allocatable :: def_width, def_maxwidth, def_tmsize, def_bmsize, def_unit
        character(:,SKC), allocatable :: def_prefix, def_indent, def_break, def_newline
        if (present(sticky)) self%sticky = sticky
        if (self%sticky) then
            if (present(prefix  ))self%prefix   = prefix    ;
            if (present(indent  ))self%indent   = indent    ;
            if (present(break   ))self%break    = break     ;
            if (present(newline ))self%newline  = newline   ;
            if (present(width   ))self%width    = width     ;
            if (present(maxwidth))self%maxwidth = maxwidth  ;
            if (present(tmsize  ))self%tmsize   = tmsize    ;
            if (present(bmsize  ))self%bmsize   = bmsize    ;
            if (present(unit    ))self%unit     = unit      ;
#if         setAborted_ENABLED
            if (present(help    ))self%help     = help      ;
            if (present(stat    ))self%stat     = stat      ;
            if (present(renabled))self%renabled = renabled  ;
#elif       !(setMarked_ENABLED || setNoted_ENABLED || setWarned_ENABLED)
#error      "Unrecognized interface."
#endif
        end if
        if (present(prefix  )) then; def_prefix   = prefix    ; elseif (allocated(self%prefix     )) then; def_prefix   = self%prefix    ; end if
        if (present(indent  )) then; def_indent   = indent    ; elseif (allocated(self%indent     )) then; def_indent   = self%indent    ; end if
        if (present(break   )) then; def_break    = break     ; elseif (allocated(self%break      )) then; def_break    = self%break     ; end if
        if (present(newline )) then; def_newline  = newline   ; elseif (allocated(self%newline    )) then; def_newline  = self%newline   ; end if
        if (present(width   )) then; def_width    = width     ; elseif (allocated(self%width      )) then; def_width    = self%width     ; end if
        if (present(maxwidth)) then; def_maxwidth = maxwidth  ; elseif (allocated(self%maxwidth   )) then; def_maxwidth = self%maxwidth  ; end if
        if (present(tmsize  )) then; def_tmsize   = tmsize    ; elseif (allocated(self%tmsize     )) then; def_tmsize   = self%tmsize    ; end if
        if (present(bmsize  )) then; def_bmsize   = bmsize    ; elseif (allocated(self%bmsize     )) then; def_bmsize   = self%bmsize    ; end if
        if (present(unit    )) then; def_unit     = unit      ; elseif (allocated(self%unit       )) then; def_unit     = self%unit      ; end if
#if     setAborted_ENABLED
        if (present(help    )) then; def_help     = help      ; elseif (allocated(self%help       )) then; def_help     = self%help      ; end if
        if (present(stat    )) then; def_stat     = stat      ; elseif (allocated(self%stat       )) then; def_stat     = self%stat      ; end if
        if (present(renabled)) then; def_renabled = renabled  ; elseif (allocated(self%renabled   )) then; def_renabled = self%renabled  ; end if
#elif   !(setMarked_ENABLED || setNoted_ENABLED || setWarned_ENABLED)
#error  "Unrecognized interface."
#endif
        call SET_MARKED( msg &
                        , prefix    = def_prefix &
                        , indent    = def_indent &
                        , break     = def_break &
                        , newline   = def_newline &
                        , width     = def_width &
                        , maxwidth  = def_maxwidth &
                        , tmsize    = def_tmsize &
                        , bmsize    = def_bmsize &
                        , unit      = def_unit &
#if                     setAborted_ENABLED
                        , stat      = def_stat &
                        , help      = def_help &
                        , renabled  = def_renabled &
#elif                   !(setMarked_ENABLED || setNoted_ENABLED || setWarned_ENABLED)
#error                  "Unrecognized interface."
#endif
                        )
        if (allocated(def_prefix   )) deallocate(def_prefix     )
        if (allocated(def_indent   )) deallocate(def_indent     )
        if (allocated(def_break    )) deallocate(def_break      )
        if (allocated(def_newline  )) deallocate(def_newline    )
        if (allocated(def_width    )) deallocate(def_width      )
        if (allocated(def_maxwidth )) deallocate(def_maxwidth   )
        if (allocated(def_tmsize   )) deallocate(def_tmsize     )
        if (allocated(def_bmsize   )) deallocate(def_bmsize     )
        if (allocated(def_unit     )) deallocate(def_unit       )
#if     setAborted_ENABLED
        if (allocated(def_help     )) deallocate(def_help       )
        if (allocated(def_stat     )) deallocate(def_stat       )
        if (allocated(def_renabled )) deallocate(def_renabled   )
#elif   !(setMarked_ENABLED || setNoted_ENABLED || setWarned_ENABLED)
#error  "Unrecognized interface."
#endif
#undef  SET_MARKED

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   mark_typer_ENABLED || note_typer_ENABLED || warn_typer_ENABLED || stop_typer_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(prefix      )) self%prefix      = prefix
        if (present(indent      )) self%indent      = indent
        if (present(break       )) self%break       = break
        if (present(newline     )) self%newline     = newline
        if (present(width       )) self%width       = width
        if (present(maxwidth    )) self%maxwidth    = maxwidth
        if (present(tmsize      )) self%tmsize      = tmsize
        if (present(bmsize      )) self%bmsize      = bmsize
        if (present(unit        )) self%unit        = unit
        if (present(sticky      )) self%sticky      = sticky
#if     stop_typer_ENABLED
        if (present(stat        )) self%stat        = stat
        if (present(help        )) self%help        = help
        if (present(renabled    )) self%renabled    = renabled
#elif   !(mark_typer_ENABLED || note_typer_ENABLED || warn_typer_ENABLED)
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif