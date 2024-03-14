program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: isFailedList, isFailedMakeDir

    implicit none

    character(:, SK), allocatable   :: list
    integer(IK)     , allocatable   :: index(:,:)
    integer                         :: i

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isFailedMakeDir('./temp')")
    call disp%show( isFailedMakeDir('./temp') )
    call disp%show("if (isFailedList('.', list, index)) then; error stop; else; do i = 1, size(index,2); call disp%show(list(index(1,i):index(2,i)), deliml = SK_'""'); end do; end if")
                    if (isFailedList('.', list, index)) then; error stop; else; do i = 1, size(index,2); call disp%show(list(index(1,i):index(2,i)), deliml = SK_""""); end do; end if
    call disp%skip()

    call disp%skip()
    call disp%show("if (isFailedList('..', list, index)) then; error stop; else; do i = 1, size(index,2); call disp%show(list(index(1,i):index(2,i)), deliml = SK_'""'); end do; end if")
                    if (isFailedList('..', list, index)) then; error stop; else; do i = 1, size(index,2); call disp%show(list(index(1,i):index(2,i)), deliml = SK_""""); end do; end if
    call disp%skip()

end program example