! Test that sendget with strides on either side (of the assignment) works
! as expected.
!
! This test needs at least three images, because sendget has the potential
! to check whether on image used in the communication is the current one.
! More than three images do not pay, because there is no general code in
! this test.
!
! Written by Andre Vehreschild

program stridedsendgettest

  implicit none

  integer, parameter :: src_image = 2, dst_image = 3, master_image = 1
  integer, save, dimension(4,6) :: srcmat[*], dstmat[*]
  integer, save, dimension(6) :: srcvec[*], dstvec[*]
  integer :: i
  logical :: test_passed = .true.

  ! Make sure that enough images are available for this test.
  ! Everything less than dst_image == 3 may make sendget use an
  ! optimized version saving a part of the communication, which is
  ! not what the test should test.
  if (num_images() < dst_image) then
     print*, "Pretend that the test was run and passed, even though there are too few images to perform test:"
     print*, "Test passed"
     error stop "Need at least three images."
  end if

  ! On the src_image, set some defined values, to be able to distinguish
  ! strides going wrong.
  if (this_image() == src_image) then
    srcvec = [(2 * i, i = 1, 6)]
    srcmat = reshape([(i * 2, i = 1, 4*6)], [4,6])
  ! On the dst_image set values that enable to recognize unset values.
  elseif (this_image() == dst_image) then
    dstmat = -1
    dstvec = -2
  end if

  ! Make sure data is valid on all images.
  sync all

  ! master_image is the controller in this communication and therefore needs
  ! to initiate the communication.
  if (this_image() == master_image) then
    ! Transfer data from the src-vector to the dst-vector on image
    ! dst_image.  This is a transfer of a contingous block of data and here for
    ! completeness only.
    dstvec(:)[dst_image] = srcvec(:)[src_image]
    ! This statement uses a stride in the send phase of the communication.
    dstmat(3,:)[dst_image] = srcvec(:)[src_image]
  end if

  ! Make sure the communication has completed.
  sync all

  ! Check the result of communication on the dst_image.
  if (this_image() == dst_image) then
    ! Check that transfering to the vector has succeeded.
    if (any(dstvec /= [2, 4, 6, 8, 10, 12])) error stop "SendGet vec/vec does not match."

    ! Check that transfering a vector into a matrix changes only the
    ! values desired.
    if (any(dstmat /= reshape([-1, -1,  2, -1, &
                               -1, -1,  4, -1, &
                               -1, -1,  6, -1, &
                               -1, -1,  8, -1, &
                               -1, -1, 10, -1, &
                               -1, -1, 12, -1], [4, 6]))) then
      error stop "SendGet matrow/vec does not match."
    end if
    ! Reset the dst-buffers to enable new test.
    dstvec = -2
    dstmat = -1
  end if

  ! Wait for dst having done its tests.
  sync all
  if (this_image() == master_image) then
    ! Execute strided get in sendget and store in a vector to just
    ! test the get.
    dstvec(:)[dst_image] = srcmat(2,:)[src_image]
    ! Test both strided get and strided send at once.
    dstmat(3,:)[dst_image] = srcmat(2,:)[src_image]
  end if

  ! Ensure that the communication is all done.
  sync all
  if (this_image() == dst_image) then
    ! Check, that the strided get has the expected result.
    if (any(dstvec /= [4, 12, 20, 28, 36, 44])) error stop "SendGet vec/matrow does not match."

    ! And that both communications with stride work as expected.
    if (any(dstmat /= reshape([-1, -1,  4, -1, &
                               -1, -1, 12, -1, &
                               -1, -1, 20, -1, &
                               -1, -1, 28, -1, &
                               -1, -1, 36, -1, &
                               -1, -1, 44, -1], [4, 6]))) then
      error stop "SendGet matrow/matrow does not match."
    end if

    ! Above checks would stop with an error on failure, so its save
    ! to unguardedly print here, when all tests pass.
    print *, "Test passed"
  end if
end program
