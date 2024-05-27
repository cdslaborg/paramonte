character(*), parameter :: LICENSE(*) = &
[ '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' &
, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' &
, '!!!!                                                                                                                            !!!!' &
, '!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!' &
, '!!!!                                                                                                                            !!!!' &
, '!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!' &
, '!!!!                                                                                                                            !!!!' &
, '!!!!    This file is part of the ParaMonte library.                                                                             !!!!' &
, '!!!!                                                                                                                            !!!!' &
, '!!!!    LICENSE                                                                                                                 !!!!' &
, '!!!!                                                                                                                            !!!!' &
, '!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!' &
, '!!!!                                                                                                                            !!!!' &
, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' &
, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' &
]
integer :: narg, lenarg, inn_unit, out_unit, i, lenrec, iostat
character(:), allocatable :: inn_file, out_file, record
character(255) :: iomsg
logical :: valuable

narg = command_argument_count()
if (narg /= 2) error stop "The inn command arguments must be only the path to the FPP file followed by the path to the out cleaned file."

allocate(character(4096) :: inn_file, out_file)
call get_command_argument(number = 1, value = inn_file, length = lenarg); inn_file = trim(adjustl(inn_file(1 : lenarg)))
call get_command_argument(number = 2, value = out_file, length = lenarg); out_file = trim(adjustl(out_file(1 : lenarg)))
open(newunit = out_unit, file = out_file, status = "replace", action = "write", iostat = iostat, iomsg = iomsg)
if (isFailed(iostat, "Failed to open the output destination file ("//trim(out_file)//") with `write` access: "//trim(iomsg))) return
open(newunit = inn_unit, file = inn_file, status = "old", action = "read", iostat = iostat, iomsg = iomsg)
if (isFailed(iostat, "Failed to open the output destination file ("//trim(inn_file)//") with `write` access: "//trim(iomsg))) return
allocate(character(511) :: record)

! add license.
do i = 1, size(LICENSE)
    write(out_unit, "(A)") LICENSE(i)
end do
write(out_unit, "(A)")

block
    use iso_fortran_env, only: iostat_end
    do
        call setRecordFrom(inn_unit, record, lenrec, iostat, iomsg)
        if (iostat == iostat_end) exit
        if (isFailed(iostat, "Failed to read the source file contents: "//trim(iomsg))) return
        i = getLocNB(record(1 : lenrec))
        if (0 < i) then
            valuable = record(i : i) /= "#"
            if (valuable) then
                valuable = record(i : i) /= "!" &
                .or. record(i : min(lenrec, i + 3)) == "!DEC$" .or. record(i : min(lenrec, i + 4)) == "!$omp"
                !.or. (record(i : min(lenrec, i + 1)) == "!>" .and. index(inn_file, "@routines") == 0)
                if (valuable) then
                    write(out_unit, "(A)") record(1 : lenrec)
                    !print *, i, lenrec, """"//record(1 : lenrec)//""""
                    !if (record(i : min(lenrec, i + 2)) == "end") write(out_unit, "(A)")
                end if
            end if
        end if
    end do
end block

close(inn_unit)
close(out_unit)

contains

    function isFailed(iostat, iomsg) result(failed)
        character(*), intent(in) :: iomsg
        integer, intent(in) :: iostat
        logical :: failed
        failed = iostat /= 0
        if (failed) write(*, "(A)") trim(iomsg)
    end function

    ! Return the location of the first non-blank character, or 0 if str is empty or all blanks.
    function getLocNB(str) result(loc)
        character(*), intent(in) :: str
        integer :: loc
        do loc = 1, len(str)
            if (str(loc : loc) /= " ") return
        end do
        loc = 0
    end function

    subroutine setRecordFrom(unit, record, ub, iostat, iomsg)
        use iso_fortran_env, only: iostat_eor
        !integer, parameter :: LENLF = len(LF)
        !character(*), parameter :: LF = new_line("a")
        character(:), intent(inout), allocatable :: record
        character(*), intent(out) :: iomsg
        integer, intent(out) :: iostat
        integer, intent(in) :: unit
        integer, intent(out) :: ub
        integer :: lb, size
        integer :: lenrec
        if (.not. allocated(record)) error stop
        lenrec = len(record)
        lb = 0
        do
            read(unit, "(a)", advance = "no", size = size, iostat = iostat, iomsg = iomsg) record(lb + 1 : lenrec)
            if (iostat == iostat_eor) then ! Record reading is complete.
                ub = lb + size
                iostat = 0
                return
            elseif (iostat == 0) then ! There is still record to read.
                lb = lb + size
                record = record//repeat(" ", lenrec)
                lenrec = lenrec + lenrec
                cycle
            else
                return
            end if
        end do
    end subroutine

end