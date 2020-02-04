program async_comp_alloc_2
  !! Author: Izaak Beekman, Andre Vehreschild
  !! Category: Regression
  !! Regression test for issue [#399](https://github.com/sourceryinstitute/OpenCoarrays/issues/399)
  implicit none
  type :: nonsymmetric
     real, dimension(:), allocatable :: arr
  end type

  type(nonsymmetric), codimension[*] :: parent_obj

  sync all

  if ( num_images() < 2 ) error stop "Test async_comp_alloc_2 requires at least 2 images to run"

  associate(me => this_image())
    if (me == 2) then
      allocate(parent_obj%arr(3))
      print *, 'Image 2: memory allocated.'
    end if

    sync all

    if (me == 1) then
      print *, 'Alloc status on [2]: ', allocated(parent_obj[2]%arr)
      print *, 'Image 2 has size ', size(parent_obj[2]%arr), ' asymmetric allocation'
      if (size(parent_obj[2]%arr) /= 3) error stop 'Test failed.'
      sync all
      print *, 'Test passed.'
    else
      sync all
    end if
  end associate
end program
