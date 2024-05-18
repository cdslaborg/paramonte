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
!>  This module contains implementations of the procedures in [pm_sysShell](@ref pm_sysShell).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     isFailedGetShellShape_ENABLED || isFailedGetShellWidth_ENABLED || isFailedGetShellHeight_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_arraySplit, only: setSplit
        integer(IK)                     :: iostat, iline
        integer(IK)     , allocatable   :: sindex(:,:)
        character(:, SK), allocatable   :: string
        character(1, SK)                :: dumstr
        type(shell_type)                :: Shell
#if     isFailedGetShellHeight_ENABLED
        integer(IK) :: width
#elif   isFailedGetShellWidth_ENABLED
        integer(IK) :: height
#elif   !isFailedGetShellShape_ENABLED
#error  "Unrecognized interface."
#endif
        shell = shell_type(failed) ! Determine the shell type.
        if (failed) return ! LCOV_EXCL_LINE

        if (shell%is%cmd) then

            ! First try calling PowerShell.

            failed = isFailedGetOutput(SK_"powershell ^(Get-Host^).UI.RawUI.WindowSize.ToString()", string) ! `width`,`height`
            if (.not. failed) then
                read(string, *, iostat = iostat) width, height
                failed = logical(iostat /= 0_IK, LK)
            end if

            ! If failed, try a pure CMD approach.

            if (failed) then
                width = -1_IK
                height = -1_IK
                failed = isFailedGetOutput(SK_"mode con /status", string)
                ! string contains a table like the following:
                !   Lines:          30
                !   Columns:        120
                !   Keyboard rate:  31
                !   Keyboard delay: 0
                !   Code page:      437
                if (.not. failed) then
                    call setSplit(sindex, string, new_line(SK_"a"))
                    do iline = 1, size(sindex, 2, IK)
                        if (index(string(sindex(1,iline) : sindex(2,iline)), SK_"Lines", kind = IK) > 0_IK) then
                            read(string, *, iostat = iostat) dumstr, height
                            failed = logical(iostat /= 0_IK, LK)
                            if (failed) return
                        end if
                        if (index(string(sindex(1,iline) : sindex(2,iline)), SK_"Columns", kind = IK) > 0_IK) then
                            read(string, *, iostat = iostat) dumstr, width
                            failed = logical(iostat /= 0_IK, LK)
                            if (failed) return
                        end if
                        if (height /= -1_IK .and. width /= -1_IK) exit
                    end do
                    failed = logical(height == -1_IK .or. width == -1_IK, LK)
                end if
            end if

        elseif (shell%is%powershell) then

            failed = isFailedGetOutput("(Get-Host).UI.RawUI.WindowSize.ToString()", string) ! `width`,`height`
            if (.not. failed) then
                read(string, *, iostat = iostat) width, height
                failed = logical(iostat /= 0_IK, LK)
            end if

        else

#if         isFailedGetShellShape_ENABLED || isFailedGetShellWidth_ENABLED
            failed = isFailedGetEnvVar(SK_"COLUMNS", string, length = 3_IK)
            if (.not. failed) then
                read(string, *, iostat = iostat) width
                failed = logical(iostat /= 0_IK, LK)
            end if
#elif       isFailedGetShellShape_ENABLED || isFailedGetShellHeight_ENABLED
            failed = isFailedGetEnvVar(SK_"LINES", string, length = 3_IK)
            if (.not. failed) then
                read(string, *, iostat = iostat) height
                failed = logical(iostat /= 0_IK, LK)
            end if
#endif
            ! One last attempt through the `tput` command.

            if (failed) then
                failed = isFailedGetOutput(SK_"echo $(tput cols),$(tput lines)", string) ! `width`,`height`
                if (.not. failed) then
                    read(string, *, iostat = iostat) width, height
                    failed = logical(iostat /= 0_IK, LK)
                end if
            end if

        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif