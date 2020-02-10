program co_reduce_res_im
  !! author: Daniel Topa & Izaak Beekman
  !! category: unit test
  !!
  !! This test is derived from
  !! [issue #172](https://github.com/sourceryinstitute/opencoarrays/issues/172)
  !! but tweaks the binary operator's (pure function) arguments have
  !! `intent(in)` which results in a working/passing test

  implicit none
  integer :: value[ * ] !! Each image stores their image number here
  integer :: k
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
     if ( value == factorial( num_images() ) ) then
        write ( * , '(a)' ) 'Test passed.'
     else
        write ( * , '(a, I0)') 'Answer should have been num_images()! = ', factorial( num_images() )
        error stop 'Wrong answer for n! using co_reduce'
     end if
  end if


contains

  pure function myProd ( a, b ) result ( rslt )
    !! Product function to be used in `co_reduce` reduction for
    !! computing factorials. When `intent(in)` attribute is changed
    !! to `value` tests fail
    integer, intent(in) :: a, b
    !! multiply two inputs together.  If we change `intent(in)` to
    !! `value` the test fails despite being correct according to C1276
    !! of F2008:
    !!
    !! > C1276 The specification-part of a pure function subprogram
    !! > shall specify that all its nonpointer dummy data objects have
    !! > the INTENT (IN) or the VALUE attribute.
    !!
    !! Thanks to @LadaF for pointing this out.
    integer        :: rslt !! product of a*b
    rslt = a * b
  end function

  pure function factorial ( n ) result ( rslt )
    !! Compute $n!$
    integer, intent(in) :: n
    integer :: rslt
    integer :: i
    rslt = 1
    do i = 1, n
      rslt = rslt*i
   end do
 end function
end program
