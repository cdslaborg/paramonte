program main
  use iso_fortran_env, only : real64,int64,compiler_version,compiler_options
  use ieee_arithmetic, only : ieee_is_nan
  use global_field_module, only : global_field,initial_condition
  use ForTrilinos_assertion_utility, only : assert,error_message
  implicit none
  type(global_field) :: u,u_half,half_uu
  real(real64), parameter :: nu=1.,final_time=0.6_real64,tolerance=1.E-3_real64,safety_factor=0.1_real64
  real(real64) :: time=0.,dt,dx
  integer, parameter :: nodes=16
  procedure(initial_condition), pointer :: initial_u=>null()

  initial_u => ten_sin

#ifdef TAU
  call TAU_PROFILE_SET_NODE(this_image()-1) ! Start TAU (Cray or GNU compiler)
#else
#ifdef TAU_INTEL
  call TAU_PROFILE_SET_NODE(this_image())   ! Start TAU (Intel compiler)
#endif
#endif

  call u%set(initial_u,num_points=nodes)
  dx = u%grid_spacing()
  dt = safety_factor*diffusion_stability_limit(nu,dx,order_of_accuracy=2)
  do while(time<final_time)
    half_uu = 0.5_real64*(u*u)
    u_half = u + (dt/2._real64)*(nu*u%xx() - half_uu%x())
    half_uu = 0.5_real64*(u_half*u_half)
    u      = u + dt*(nu*u_half%xx() - half_uu%x())
    time = time + dt
  end do
  if (this_image()==1) print *,"Time =",time
  print *,"On image ",this_image(),"u =",u%state()
  call test(u)
  sync all
  if (this_image()==1) print *,"Test passed."

contains
  subroutine test(burgers_solution)
    type(global_field), intent(in) :: burgers_solution
    call assert(.not.any(ieee_is_nan(u%state())),error_message("Test failed: u is not a number."))
    call assert(sinusoid(u),error_message("Test failed: improper shape."))
  end subroutine

  function sinusoid(u_solution) result(is_sinusoid)
    type(global_field), intent(in) :: u_solution
    type(global_field) :: u_xx
    logical :: is_sinusoid
    real(real64), parameter :: threshold=-0.001,cap=0.001
    real(real64), allocatable :: u_xx_state(:)
    u_xx = u_solution%xx()
    u_xx_state = u_xx%state()
    if (num_images()/=1) then
      ! Ensure that the global midpoint is a local endpoint for whatever image contains the midpoint:
      call assert(mod(num_images(),2)==0,error_message("Test failed: uneven number of images."))
      ! Ensure that the left and right halves of the solution are concave down and up, respectively:
      if (this_image()<=num_images()/2) then
        call assert(all(u_xx_state<cap),error_message("Test failed: right half not concave up."))
      else
        call assert(all(u_xx_state>threshold),error_message("Test failed: left half not concave down."))
      end if
    else
      block
        integer :: size_u_xx
        size_u_xx = size(u_xx_state)
        call assert(all(u_xx_state(1:size_u_xx/2)<cap),error_message("Test failed: left half not concave down."))
        call assert(all(u_xx_state(size_u_xx/2+1:size_u_xx)>threshold),error_message("Test failed: right half not concave up."))
      end block
    end if
    is_sinusoid=.true.
  end function

  pure function diffusion_stability_limit(diffusivity,delta_x,order_of_accuracy)  result(stable_time_step)
    real(real64), intent(in) :: diffusivity,delta_x
    integer, intent(in) :: order_of_accuracy
    real(real64) :: stable_time_step
    real(real64), parameter, dimension(*) :: stability_limit=[2.,2.,2.5,2.79] ! third value needs to be checked
    ! See Moin, P. (2010) Fundamentals of Engineering Numerical Analysis, 2nd ed., pp. 111-116.
    stable_time_step = safety_factor*stability_limit(order_of_accuracy)*(delta_x**2)/(4._real64*diffusivity)
  end function

  pure function ten_sin(x) result(ten_sin_x)
    real(real64), intent(in) :: x
    real(real64) :: ten_sin_x
    ten_sin_x = 10._real64*sin(x)
  end function
end program
