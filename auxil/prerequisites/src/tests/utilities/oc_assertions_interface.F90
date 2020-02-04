#ifndef USE_ASSERTIONS
# define USE_ASSERTIONS .false.
#endif
module oc_assertions_interface
  !! author: Damian Rouson
  !!
  !! Utility for runtime checking of logical assertions.
  !!
  !! Instructions
  !! ------------
  !! Compile with -DUSE_ASSERTIONS=.false. to define the logical parameter named "assertions" and to thereby
  !! facilitate the elimination of assertions during the dead-code removal phase of optimizing compilers:
  !!
  !!    gfortran -cpp -DUSE_ASSERTIONS=.false. -c assertions_interface.f90
  !!
  !! or set the corresponding NO_ASSERTIONS variable defined in this directory's CMakeLists.txt:
  !!
  !!    FC=caf cmake <opencoarrays-source-path> -DNO_ASSERTIONS=ON
  !!
  !! Conditioning assertion calls on the "assertions" compile-time constant enables optimizing compilers
  !! to eliminate assertion calls via dead-code-removal optimiztion.
  !!
  !! Use case 1
  !! ----------
  !!    Pass the optional success argument & check for false return value as an indication of assertion failure:
  !!
  !!    use opencoarrays_assertions_interface, only : assert, assertions
  !!    if (assertions) call assert( 2 > 1, "always true inequality", success)
  !!    if (error_code/=0) call my_error_handler()
  !!
  !! Use case 2
  !! ----------
  !!    Error-terminate if the assertion fails:
  !!
  !!    use opencoarrays_assertions_interface, only : assert,assertions
  !!    if (assertions) call assert( 2 > 1, "always true inequality")
  !!
  implicit none
  private
  public :: assert
  public :: assertions

  logical, parameter :: assertions=USE_ASSERTIONS

  interface
#ifdef HAVE_ERROR_STOP_IN_PURE
    pure &
#endif
    module subroutine assert(assertion,description,diagnostic_data,success,error_message)
      !! On false assertion, error-terminate or, if present(success), copy assertion into success
      implicit none
      logical, intent(in) :: assertion
        !! Most assertions will be expressions, e.g., call assert( i>0, "positive i")
      character(len=*), intent(in) :: description
        !! Brief statement of what is being asserted
      class(*), intent(in), optional :: diagnostic_data
        !! Optional data to printed or added to error_message assertion evaluates to .false.
      logical, intent(out), optional :: success
        !! Optional copy of the assertion dummy argument
      character(len=:), intent(out), optional, allocatable :: error_message
        !! Optional informational message allocated only if assertion==.false. .and. present(success)
    end subroutine
  end interface
end module oc_assertions_interface
