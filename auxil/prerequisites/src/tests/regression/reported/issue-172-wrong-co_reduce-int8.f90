program co_reduce_factorial_int8
  !! author: Daniel Topa & Izaak Beekman
  !! category: regression
  !!
  !! [issue #172](https://github.com/sourceryinstitute/opencoarrays/issues/172)
  !! wherein co-reduce gets junk in the first image when binary
  !! operator's (pure function) arguments have `value` attribute
  !! instead of `intent(in)`

  implicit none
  integer(kind=1) :: value[ * ] !! Each image stores their image number here
  integer :: k
  integer(kind=1) :: np
  np = num_images ( )
  value = this_image ( )
  call co_reduce ( value, result_image = 1, operator = myProd )
  !! value[k /= 1] undefined, value[ k == 1 ] should equal $n!$ where $n$ is `num_images()`
  if ( this_image ( ) == 1 ) then
     write ( * , '( "Number of images = ", g0 )' ) num_images ( )
     do k = 1, num_images ( )
        write ( * , '( 2( a, i0 ) )' ) 'value [ ', k, ' ] is ', value [ k ]
        write ( * , '(a)' ) 'since RESULT_IMAGE is present, value on other images is undefined by the standard'
     end do
     write ( * , '( "Product  value = ", g0 )' ) value  !! should print num_images() factorial
     write ( * , 100 )
     if ( value == factorial( np ) ) then
        write ( * , '(a)' ) 'Test passed.'
     else
        write ( * , '(a, I0)') 'Answer should have been num_images()! = ', factorial( np )
        error stop 'Wrong answer for n! using co_reduce'
     end if
  end if
100 format ( "Expected value = num_images()!", /, " 2! = 2, 3! = 6, 4! = 24, ..." )

contains

  pure function myProd ( a, b ) result ( rslt )
    !! Product function to be used in `co_reduce` reduction for
    !! computing factorials. When `value` attribute is changed to
    !! `intent(in)` tests pass, and expected behavior is observed.
    integer(kind=1), value :: a, b
    !! multiply two inputs together.  If we change `value` to
    !! `intent(in)` the test passes and the issue goes away and
    !! according to C1276 of F2008:
    !!
    !! > C1276 The specification-part of a pure function subprogram
    !! > shall specify that all its nonpointer dummy data objects have
    !! > the INTENT (IN) or the VALUE attribute.
    !!
    !! Thanks to @LadaF for pointing this out.
    integer(kind=1)        :: rslt !! product of a*b
    rslt = a * b
  end function

  pure function factorial ( n ) result ( rslt )
    !! Compute $n!$
    integer(kind=1), intent(in) :: n
    integer(kind=1) :: rslt
    integer :: i
    rslt = 1
    do i = 1, n
      rslt = rslt*i
   end do
 end function
end program
