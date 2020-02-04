program main
  !! author: Damian Rouson
  !!
  !! Test co_broadcast with derived-type actual arguments
  implicit none

  integer, parameter :: sender=1 !! co_broadcast source_image
  character(len=*), parameter :: text="text" !! character message data

  type dynamic
    character(len=:), allocatable :: string
    character(len=len(text)), allocatable :: string_array(:)
    complex, allocatable :: scalar
    integer, allocatable :: vector(:)
    logical, allocatable :: matrix(:,:)
    real, allocatable ::  superstring(:,:,:, :,:,:, :,:,:, :,:,:, :,:,: )
  end type

  type(dynamic) alloc_message, alloc_content

  associate(me=>this_image())

    alloc_content = dynamic(                                               &
      string=text,                                                         &
      string_array=[text],                                                 &
      scalar=(0.,1.),                                                      &
      vector=reshape( [integer::], [0]),                                   &
      matrix=reshape( [.true.], [1,1]),                                    &
      superstring=reshape([1,2,3,4], [2,1,2, 1,1,1, 1,1,1, 1,1,1, 1,1,1 ]) &
      )

    alloc_message = alloc_content

    if(me /= sender) then
       alloc_message%vector = 0
       alloc_message%matrix = .false.
       alloc_message%superstring = 0
    endif

    sync all

   call co_broadcast(alloc_message,source_image=sender)

   associate( failures => [                                    &
     alloc_message%string /= alloc_content%string,             &
     alloc_message%string_array /= alloc_content%string_array, &
     alloc_message%scalar /= alloc_content%scalar,             &
     alloc_message%vector /= alloc_content%vector,             &
     alloc_message%matrix .neqv. alloc_content%matrix,         &
     alloc_message%superstring /= alloc_content%superstring    &
   ] )

      if ( any(failures) ) error stop "Test failed."

    end associate

    sync all  ! Wait for each image to pass the test
    if (me==sender) print *,"Test passed."

  end associate

end program main
