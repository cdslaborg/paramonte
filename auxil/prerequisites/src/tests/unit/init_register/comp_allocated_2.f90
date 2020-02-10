program comp_allocated_2

    implicit none
    integer, parameter :: success = 0
    type :: subType
        real, allocatable :: r_comp
    end type

    type :: T
        type(subType), dimension(:), allocatable :: arr
    end type

    type(T), codimension[:], allocatable :: obj

    call assert(num_images() .GE. 2, 'Need at least two images.')

    associate(me => this_image())
        allocate(obj[*])

        if (me == 1) then
            call assert(.NOT. allocated(obj[2]%arr), 'obj%arr on image 2 allocated.')
        end if

        sync all

        if (me == 2) then
            allocate(obj%arr(3))
            allocate(obj%arr(2)%r_comp, source=13.7)
            print *, 'Image 2: memory allocated.'
        end if

        sync all

        if (me == 1) then
            call assert(allocated(obj[2]%arr), 'obj%arr on image 2 not allocated.')
            call assert(.NOT. allocated(obj[2]%arr(1)%r_comp), 'obj%arr(1)%r_comp should not be allocated')
            call assert(allocated(obj[2]%arr(2)%r_comp), 'obj%arr(2)%r_comp should be allocated')
            call assert(.NOT. allocated(obj[2]%arr(3)%r_comp), 'obj%arr(3)%r_comp should not be allocated')
            print *,'Test passed.'
        end if
        sync all
    end associate
contains
  subroutine assert(assertion,description,stat)
    logical, intent(in) :: assertion
    character(len=*), intent(in) :: description
    integer, intent(out), optional:: stat
    integer, parameter :: failure=1
    if (assertion) then
       if (present(stat)) stat=success
    else
      if (present(stat)) then
        stat=failure
      else
        error stop "Assertion "// description //" failed."
      end if
    end if
  end subroutine
end program

! vim:sw=4:ts=4:sts=4:
