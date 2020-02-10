program main
  !! author: Damian Rouson
  !!
  !! Test co_broadcast with derived-type actual arguments
  implicit none

  integer, parameter :: sender=1 !! co_broadcast source_image
  character(len=*), parameter :: text="text" !! character message data

  interface
     function f(x) result(y)
       real x, y
     end function
  end interface

  type parent
    integer :: heritable=0
  end type

  type component
    integer :: subcomponent=0
  end type

  type, extends(parent) :: child

    ! Scalar and array derived-type components
    type(component) a, b(1,2,1, 1,1,1, 1)

    ! Scalar and array intrinsic-type components
    character(len=len(text)) :: c="", z(0)
    complex :: i=(0.,0.), j(1)=(0.,0.)
    integer :: k=0,       l(2,3)=0
    logical :: r=.false., s(1,2,3, 1,2,3, 1)=.false.
    real    :: t=0.,      u(3,2,1)=0.

    ! Scalar and array pointer components
    character(len=len(text)), pointer :: &
                        char_ptr=>null(), char_ptr_maxdim(:,:,:, :,:,:, :)=>null()
    complex, pointer :: cplx_ptr=>null(), cplx_ptr_maxdim(:,:,:, :,:,:, :)=>null()
    integer, pointer :: int_ptr =>null(), int_ptr_maxdim (:,:,:, :,:,:, :)=>null()
    logical, pointer :: bool_ptr=>null(), bool_ptr_maxdim(:,:,:, :,:,:, :)=>null()
    real, pointer    :: real_ptr=>null(), real_ptr_maxdim(:,:,:, :,:,:, :)=>null()
    procedure(f), pointer :: procedure_pointer=>null()
  end type

  type(child) message
  type(child) :: content = child( & ! define content using the insrinsic structure constructor
    parent=parent(heritable=-4),                                                                           & ! parent
    a=component(-3), b=reshape([component(-2),component(-1)], [1,2,1, 1,1,1, 1]),        & ! derived types
    c=text, z=[character(len=len(text))::], i=(0.,1.), j=(2.,3.), k=4, l=5, r=.true., s=.true., t=7., u=8. & ! intrinsic types
  )

  associate(me=>this_image())

    if (me==sender) then
      message = content
      allocate(message%char_ptr, message%char_ptr_maxdim(1,1,2, 1,1,1, 1), source=text   )
      allocate(message%cplx_ptr, message%cplx_ptr_maxdim(1,1,1, 1,1,2, 1), source=(0.,1.))
      allocate(message%int_ptr , message%int_ptr_maxdim (1,1,1, 1,1,1, 1), source=2      )
      allocate(message%bool_ptr, message%bool_ptr_maxdim(1,1,1, 1,2,1, 1), source=.true. )
      allocate(message%real_ptr, message%real_ptr_maxdim(1,1,1, 1,1,1, 1), source=3.     )
    end if

    call co_broadcast(message,source_image=sender)

    if (me==sender) then
      deallocate(message%char_ptr, message%char_ptr_maxdim)
      deallocate(message%cplx_ptr, message%cplx_ptr_maxdim)
      deallocate(message%int_ptr , message%int_ptr_maxdim )
      deallocate(message%bool_ptr, message%bool_ptr_maxdim)
      deallocate(message%real_ptr, message%real_ptr_maxdim)
    end if

    !! Verify correct broadcast of all non-pointer components (pointers become undefined on the receiving image).
    associate( failures => [                                &
      message%parent%heritable /= content%parent%heritable, &
      message%a%subcomponent /= content%a%subcomponent,     &
      message%c /= content%c,                               &
      message%z /= content%z,                               &
      message%i /= content%i,                               &
      message%j /= content%j,                               &
      message%k /= content%k,                               &
      message%l /= content%l,                               &
      message%r .neqv. content%r,                           &
      message%s .neqv. content%s,                           &
      message%t /= content%t,                               &
      any( message%u /= content%u )                         &
    ] )

      if ( any(failures) ) error stop "Test failed. "

    end associate

    sync all  ! Wait for each image to pass the test
    if (me==sender) print *,"Test passed."

  end associate

end program main
