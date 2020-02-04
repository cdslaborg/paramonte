submodule(oc_assertions_interface) oc_assertions_implementation
  !! Define the assert procedure

  implicit none

contains

  module procedure assert

    character(len=:), allocatable :: message

    if (present(success)) success=assertion

    if (.not.assertion) then

      call set(message)
      if (present(error_message)) error_message = message

      if (.not. present(success)) then
#ifdef HAVE_VARIABLE_STOP_CODE
        error stop message
#else
        error stop "Assertion failed."
#endif
      end if
    end if

  contains

    pure subroutine set(msg)

      use opencoarrays_object_interface, only : object
        !! import abstract type with generic binding supporting user-defined derived type output

      character(len=:), intent(out),  allocatable :: msg
      integer, parameter :: max_image_num_digits=12, max_diagnostic_length=1024
      integer, parameter :: msg_len = max_diagnostic_length + max_diagnostic_length

      associate( base_msg_len=> len("Assertion '") +len_trim(adjustl(description)) +len("' failed on image ") +max_image_num_digits)
        allocate( character(len=msg_len) :: msg )
        write(msg,*) "Assertion '",trim(adjustl(description)),"' failed on image ",this_image()
      end associate

      if (present(diagnostic_data)) then
        block
          !! append diagnostic string
          character(len=max_diagnostic_length) diagnostic_string

          select type(diagnostic_data)
            type is(character(len=*))
              diagnostic_string = diagnostic_data
            type is(real)
              write(diagnostic_string,*) diagnostic_data
            type is(integer)
              write(diagnostic_string,*) diagnostic_data
            class is(object)
              write(diagnostic_string,*) diagnostic_data
            class default
              diagnostic_string = "of unrecognized type"
          end select

          msg = trim(adjustl(msg)) // " with diagnostic data " // trim(adjustl(diagnostic_string))
        end block
      end if

    end subroutine

  end procedure

end submodule oc_assertions_implementation
