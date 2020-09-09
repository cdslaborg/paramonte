!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a 
!!!!   copy of this software and associated documentation files (the "Software"), 
!!!!   to deal in the Software without restriction, including without limitation 
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!!!!   and/or sell copies of the Software, and to permit persons to whom the 
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be 
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of 
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
!!!!   work (education/research/industry/development/...) by citing the ParaMonte 
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   NOTE
!
!   This code has been modified and adapted from Introduction to Computational Economics Using Fortran by Hans Fehr and Fabian Kindermann
!
!   Currently this code is NOT USED IN ANY WAY IN ANY PARTS of the ParaMonte library and therefore the licensing of this file
!   does not apply to any parts of the ParaMonte library or its users.
!
!###################################################################################################################################
!
! This code is published under the GNU General Public License v3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
! Our program codes are distributed under the GNU General Public License v3. 
! This essentially means that you are allowed to copy and adapt the source code 
! of this toolbox and use it for your own purposes.
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
! It accompanies the textbook:
!     Fehr, H. & Kindermann, F. (2018). Introduction to Computational
!         Economics using Fortran. Oxford: Oxford University Press.
!
! #VC# VERSION: 1.3  (21 April 2020)
!
!##############################################################################
!##############################################################################

module EconomicsToolbox_mod

!##############################################################################
!##############################################################################
! Declaration of Variables
!##############################################################################
!##############################################################################

use Constants_mod, only: IK, RK; implicit none

private ! declare everything as private by default

character(*), parameter, public :: MODULE_NAME = "@EconomicsToolbox_mod"

public :: RK, IK ! needed for testing

! starting time for cpu timer
real(RK), private :: starttime_cpu

! should the random tbox_seed be set
logical, private :: tbox_seed = .true.

! Level of tolerance for all routines
real(RK),  private  :: tbox_gftol = 1d-8

! Maximum number of iterations
integer(IK), private  :: tbox_itermax_min = 200

! Maximum number of iterations for brent_pow
integer(IK), parameter, private  :: tbox_tbox_itermax_pow_b = 150

! Level of tolerance for all routines
real(RK),  private  :: tbox_gftol_root = 1d-8

! Maximum number of iterations for broydn
integer(IK), private  :: itermax_root = 200

! control variables for gnuplot
logical, private :: gnu_addtoplot = .false.
logical, private :: gnu_dolegend = .false.
logical, private :: gnu_histogram = .false.
integer(IK), private :: gnu_nmax, gnu_mmax

! plot data for gnuplot
real(RK), allocatable, private :: gnu_x(:, :), gnu_y(:, :)
real(RK), allocatable, private :: gnu_x_temp(:, :), gnu_y_temp(:, :)

! style data for gnuplot
character(LEN = 2000), dimension(1000), private :: gnu_definitions


!##############################################################################
!##############################################################################
! Define public access points
!##############################################################################
!##############################################################################

! matrix operations
public :: lu_solve, lu_invert, lu_dec
public :: cholesky

! rootfinding
public :: settol_root
public :: setiter_root
public :: fzero

! optimization
public :: settol_min
public :: setiter_min
public :: fminsearch, powell, brent

! linear programming with constraints
public :: solve_lin

! integration methods
public :: legendre

! discretization of normal distributions
public :: normal_discrete
public :: log_normal_discrete

! discretization of AR(1) process
public :: discretize_AR, discretize_log_AR

! probabilities and distributions
public :: uniformPDF, uniformCDF, uniformCDF_Inv
public :: normalPDF, normalCDF, normalCDF_Inv
public :: log_normalPDF, log_normalCDF, log_normalCDF_Inv
public :: GammaPDF, GammaCDF
public :: betaPDF, betaCDF
public :: bernoulliPDF, bernoulliCDF
public :: binomialPDF, binomialCDF
public :: binomial_coefficient

! simulation
public :: init_random_seed
public :: simulate_uniform
public :: simulate_normal
public :: simulate_log_normal
public :: simulate_Gamma
public :: simulate_beta
public :: simulate_bernoulli
public :: simulate_binomial
public :: simulate_AR

! discretization of continuous intervals
public :: grid_Cons_Equi, grid_Val_Equi, grid_Inv_Equi
public :: grid_Cons_Cheb, grid_Val_Cheb, grid_Inv_Cheb
public :: grid_Cons_Grow, grid_Val_Grow, grid_Inv_Grow

! interpolation
public :: poly_interpol
public :: linint_Equi, linint_Cheb, linint_Grow, linint_Gen
public :: spline_interp, spline_eval, spline

! plotting
public :: plot
public :: plot3d
public :: plot_hist
public :: execplot

! sorting
public :: sort

! the clock
public :: tic, toc

!##############################################################################
!##############################################################################
! Interface declarations
!##############################################################################
!##############################################################################


!##############################################################################
! INTERFACE assert_eq
!
! Interface for equality assertions by assert_eqx functions.
!##############################################################################
interface assert_eq

    module procedure assert_eq2, assert_eq3, assert_eq4, assert_eq5, assert_eqn

end interface


!##############################################################################
! INTERFACE normal_discrete
!
! Discretizes normal distribution.
!##############################################################################
interface normal_discrete

    module procedure normal_discrete_1, normal_discrete_2

end interface


!##############################################################################
! INTERFACE log_normal_discrete
!
! Discretizes log-normal distribution.
!##############################################################################
interface log_normal_discrete

    module procedure log_normal_discrete_1, log_normal_discrete_2

end interface


!##############################################################################
! INTERFACE simulate_uniform
!
! Simulates uniformly distributed random variables.
!##############################################################################
interface simulate_uniform

    module procedure simulate_uniform_1, simulate_uniform_n

end interface


!##############################################################################
! INTERFACE simulate_normal
!
! Simulates normallly distributed random variables.
!##############################################################################
interface simulate_normal

    module procedure simulate_normal_1, simulate_normal_n

end interface


!##############################################################################
! INTERFACE simulate_log-normal
!
! Simulates log-normallly distributed random variables.
!##############################################################################
interface simulate_log_normal

    module procedure simulate_log_normal_1, simulate_log_normal_n

end interface


!##############################################################################
! INTERFACE simulate_Gamma
!
! Simulates Gamma distributed random variables.
!##############################################################################
interface simulate_Gamma

    module procedure simulate_Gamma_1, simulate_Gamma_n

end interface


!##############################################################################
! INTERFACE simulate_beta
!
! Simulates beta distributed random variables.
!##############################################################################
interface simulate_beta

    module procedure simulate_beta_1, simulate_beta_n

end interface


!##############################################################################
! INTERFACE simulate_bernoulli
!
! Simulates bernoulli distributed random variables.
!##############################################################################
interface simulate_bernoulli

    module procedure simulate_bernoulli_1, simulate_bernoulli_n

end interface


!##############################################################################
! INTERFACE simulate_binomial
!
! Simulates binomial distributed random variables.
!##############################################################################
interface simulate_binomial

    module procedure simulate_binomial_1, simulate_binomial_n

end interface


!##############################################################################
! INTERFACE poly_interpol
!
! Interface for onedimensional polynomial interpolation.
!##############################################################################
interface poly_interpol

    module procedure poly_interpol_1, poly_interpol_m

end interface


!##############################################################################
! INTERFACE fminsearch
!
! Finds minimum of a function.
!##############################################################################
interface fminsearch

    module procedure brent, powell

end interface


!##############################################################################
! INTERFACE fzero
!
! Finds root of a function.
!##############################################################################
interface fzero

    module procedure newton_interpol, broydn

end interface


!##############################################################################
! INTERFACE grid_Inv_Equi
!
! Interface for inverting gridpoints.
!##############################################################################
interface grid_Inv_Equi

    module procedure grid_Inv_Equi_1, grid_Inv_Equi_m

end interface


!##############################################################################
! INTERFACE grid_Inv_Cheb
!
! Interface for inverting gridpoints.
!##############################################################################
interface grid_Inv_Cheb

    module procedure grid_Inv_Cheb_1, grid_Inv_Cheb_m

end interface


!##############################################################################
! INTERFACE grid_Inv_Grow
!
! Interface for inverting gridpoints.
!##############################################################################
interface grid_Inv_Grow

    module procedure grid_Inv_Grow_1, grid_Inv_Grow_m

end interface


!##############################################################################
! INTERFACE linint_Equi
!
! Interface for inverting gridpoints.
!##############################################################################
interface linint_Equi

    module procedure linint_Equi_1, linint_Equi_m

end interface


!##############################################################################
! INTERFACE linint_Cheb
!
! Interface for inverting gridpoints.
!##############################################################################
interface linint_Cheb

    module procedure linint_Cheb_1, linint_Cheb_m

end interface


!##############################################################################
! INTERFACE linint_Grow
!
! Interface for inverting gridpoints.
!##############################################################################
interface linint_Grow

    module procedure linint_Grow_1, linint_Grow_m

end interface


!##############################################################################
! INTERFACE spline_interp
!
! Interface for one- or multidimensional spline interpolation.
!##############################################################################
interface spline_interp

    module procedure spline_interp1, spline_interp2, spline_interp3, &
        spline_interp4, spline_interp5, spline_interp6, spline_interp7

end interface


!##############################################################################
! INTERFACE spline_eval
!
! Interface for evaluation of one- or multidimensional spline.
!##############################################################################
interface spline_eval

    module procedure &
        spline1, spline1_grid, &
        spline2, spline2_grid, &
        spline3, spline3_grid, &
        spline4, spline4_grid, &
        spline5, spline5_grid, &
        spline6, spline6_grid, &
        spline7, spline7_grid

end interface


!##############################################################################
! INTERFACE spline
!
! Interface for complete interpolation and evaluation of one- or
!     multidimensional spline.
!##############################################################################
interface spline

    module procedure &
        spline1_complete, spline1_complete_m , &
        spline2_complete, spline2_complete_m, &
        spline3_complete, spline3_complete_m , &
        spline4_complete, spline4_complete_m, &
        spline5_complete, spline5_complete_m, &
        spline6_complete, spline6_complete_m, &
        spline7_complete, spline7_complete_m

end interface


!##############################################################################
! INTERFACE plot_hist
!
! Creates histogram plot
!##############################################################################
interface plot_hist

    module procedure plot_hist_old, plot_hist_new

end interface



!##############################################################################
! INTERFACE plot3d
!
! Creates a two-dimensional plot.
!##############################################################################
interface plot3d

    module procedure plot3d_grid, plot3d_line

end interface


!##############################################################################
! INTERFACE sort
!
! Interface for sorting arrays.
!##############################################################################
interface sort

    module procedure sort_r, sort_r2, sort_i, sort_i2

end interface


contains















!##############################################################################
!##############################################################################
! MODULE errwarn
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE error
    !
    ! Throws error message and stops program.
    !##############################################################################
    subroutine error(routine, message)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! routine in which error occured
        character(len=*), intent(in) :: routine

        ! error message
        character(len=*), intent(in) :: message


        !##### ROUTINE CODE #######################################################

        ! write error message
        write(*,'(/a,a,a,a/)')'ERROR ',routine,': ',message

        ! stop program
        stop

    end subroutine error


    !##############################################################################
    ! SUBROUTINE warning
    !
    ! Throws warning message
    !##############################################################################
    subroutine warning(routine, message)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! routine in which warning occured
        character(len=*), intent(in) :: routine

        ! warning message
        character(len=*), intent(in) :: message


        !##### ROUTINE CODE #######################################################

        ! write warning message
        write(*,'(/a,a,a,a/)')'WARNING ',routine,': ',message

    end subroutine warning















!##############################################################################
!##############################################################################
! MODULE assertions
!##############################################################################
!##############################################################################


    !##############################################################################
    ! FUNCTION assert_eq2
    !
    ! Checks equality for two integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eq2(n1, n2, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer(IK), intent(in) :: n1, n2

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer(IK) :: assert_eq2


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (n1 == n2)then
            assert_eq2 = n1

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq2')
        end if

    end function assert_eq2


    !##############################################################################
    ! FUNCTION assert_eq3
    !
    ! Checks equality for three integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eq3(n1, n2, n3, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer(IK), intent(in) :: n1, n2, n3

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer(IK) :: assert_eq3


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (n1 == n2 .and. n2 == n3)then
            assert_eq3 = n1

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq3')
        end if

    end function assert_eq3


    !##############################################################################
    ! FUNCTION assert_eq4
    !
    ! Checks equality for four integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eq4(n1, n2, n3, n4, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer(IK), intent(in) :: n1, n2, n3, n4

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer(IK) :: assert_eq4


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4)then
            assert_eq4 = n1

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq4')
        end if

    end function assert_eq4


    !##############################################################################
    ! FUNCTION assert_eq5
    !
    ! Checks equality for five integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eq5(n1, n2, n3, n4, n5, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer(IK), intent(in) :: n1, n2, n3, n4, n5

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer(IK) :: assert_eq5


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4 .and. n4 == n5)then
            assert_eq5 = n1

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq5')
        end if

    end function assert_eq5


    !##############################################################################
    ! FUNCTION assert_eqn
    !
    ! Checks equality for n integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function assert_eqn(nn, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer(IK), intent(in) :: nn(:)

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer(IK) :: assert_eqn


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (all(nn(2:) == nn(1)))then
            assert_eqn = nn(1)

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eqn')
        end if

    end function assert_eqn















!##############################################################################
!##############################################################################
! MODULE clock
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE tic
    !
    ! Starts cpu timer.
    !##############################################################################
    subroutine tic()


        !##### ROUTINE CODE #######################################################

        ! get cpu time
        call cpu_time(starttime_cpu)

    end subroutine tic


    !##############################################################################
    ! SUBROUTINE toc
    !
    ! Stops cpu timer.
    !##############################################################################
    subroutine toc(file)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! optional file identifier
        integer(IK), intent(in), optional :: file


        !##### OTHER VARIABLES ####################################################

        real(RK) :: time
        integer(IK) :: outfile
        real(RK) :: times(4)


        !##### ROUTINE CODE #######################################################

        ! get output file identifier
        if(present(file))then
            outfile = file
        else
            outfile = 0
        endif

        ! get cpu time
        call cpu_time(time)

        ! calculate time difference
        time = time - starttime_cpu

        ! get number of days
        times(1) = floor(time/(24d0*60d0*60d0))
        time = time - times(1)*24d0*60d0*60d0

        ! get number of hours
        times(2) = floor(time/(60d0*60d0))
        time = time - times(2)*60d0*60d0

        ! get number of minutes
        times(3) = floor(time/60d0)
        time = time - times(3)*60d0

        ! get number of seconds
        times(4) = time

        call outTime(times, outfile)

    end subroutine toc


    !##############################################################################
    ! SUBROUTINE outTime
    !
    ! Writes time to file.
    !##############################################################################
    subroutine outTime(times, file)

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! time as integer(IK) array
        real(RK), intent(in) :: times(4)

        ! the output file identifier
        integer(IK), intent(in) :: file

        !##### OTHER VARIABLES ####################################################

        character(len=200) :: output1, output2


        !##### ROUTINE CODE #######################################################

        ! set up output
        write(output1, '(a)')'Time elapsed: '
        output2 = output1

        ! write time values
        if(times(1) > 0d0)then
            write(output1, '(a,1x,i3,a)') trim(output2), int(times(1)), ' d  '
            output2 = output1
        endif
        if(times(2) > 0d0)then
            write(output1, '(a,1x,i3,a)') trim(output2), int(times(2)), ' h  '
            output2 = output1
        endif
        if(times(3) > 0d0) then
            write(output1, '(a,1x,i3,a)')trim(output2), int(times(3)), ' min  '
            output2 = output1
        endif

        write(output1, '(a,1x,f7.3,a)')trim(output2), times(4), ' s  '

        if(file > 0) then
            write(file, '(/a/)')trim(output1)
        else
            write(*, '(/a/)')trim(output1)
        endif

    end subroutine outTime














!##############################################################################
!##############################################################################
! MODULE linint
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE grid_Cons_Equi
    !
    ! Constructs a whole equidistant grid on [left,right].
    !##############################################################################
    subroutine grid_Cons_Equi(a, left, right)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the array to fill
        real(RK), intent(out) :: a(1:)

        ! left and right interval point
        real(RK), intent(in) :: left, right


        !##### OTHER VARIABLES ####################################################

        real(RK) :: h
        integer(IK) :: j, n


        !##### ROUTINE CODE #######################################################

        ! get number of grid points
        n = size(a, 1)-1

        ! check for left <= right
        if(left >= right)call error('grid_Cons_Equi', &
            'left interval point greater than right point')

        ! calculate distance between grid points
        h = (right-left)/dble(n)

        ! calculate grid value
        a = h*(/(dble(j), j=0,n)/)+left

    end subroutine grid_Cons_Equi


    !##############################################################################
    ! SUBROUTINE grid_Cons_Cheb
    !
    ! Constructs a whole grid on [left,right] using chebychev nodes.
    !##############################################################################
    subroutine grid_Cons_Cheb(a, left, right)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the array to fill
        real(RK), intent(out) :: a(1:)

        ! left and right interval point
        real(RK), intent(in) :: left, right


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: j, n
        real(RK), parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! get number of grid points
        n = size(a, 1)-1

        ! check for left <= right
        if(left >= right)call error('grid_Cons_Cheb', &
            'left interval point greater than right point')

        ! calculate grid value
        a = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-(/(dble(j), j=0,n)/)+0.5d0)/dble(n+1)*pi)

    end subroutine grid_Cons_Cheb


    !##############################################################################
    ! SUBROUTINE grid_Cons_Grow
    !
    ! Constructs a growing grid on [left, right].
    !##############################################################################
    subroutine grid_Cons_Grow(a, left, right, growth)

        use Constants_mod, only: IK, RK; implicit none

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the array to fill
        real(RK), intent(out) :: a(1:)

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! the growth rate of the grid
        real(RK) :: growth


        !##### OTHER VARIABLES ####################################################

        real(RK) :: h
        integer(IK) :: j, n


        !##### ROUTINE CODE #######################################################

        ! get number of grid points
        n = size(a, 1)-1

        ! check for left <= right
        if(left >= right)call error('grid_Cons_Grow', &
            'left interval point greater than right point')

        ! check for growth
        if(growth <= 0d0)call error('grid_Cons_Grow', &
            'growth rate must be greater than zero')

        ! calculate factor
        h = (right-left)/((1d0+growth)**n-1d0)

        ! calculate grid value
        a = h*((1d0+growth)**(/(dble(j), j=0,n)/)-1d0)+left

    end subroutine grid_Cons_Grow


    !##############################################################################
    ! FUNCTION grid_Val_Equi
    !
    ! Calculates single gridpoint of an equidistant grid.
    !##############################################################################
    function grid_Val_Equi(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! value at the grid point x \in [0,n]
        real(RK) :: grid_Val_Equi


        !##### OTHER VARIABLES ####################################################

        real(RK) :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Val_Equi', &
            'left interval point greater than right point')

        ! calculate distance between grid points
        h = (right-left)/n

        ! calculate grid value
        grid_Val_Equi = h*x+left

    end function grid_Val_Equi


    !##############################################################################
    ! FUNCTION grid_Val_Cheb
    !
    ! Calculates single gridpoint of a Chebychev grid.
    !##############################################################################
    function grid_Val_Cheb(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! value at the grid point x \in [0,n]
        real(RK) :: grid_Val_Cheb


        !##### OTHER VARIABLES ####################################################

        real(RK), parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Val_Cheb', &
            'left interval point greater than right point')

        ! calculate grid value
        grid_Val_Cheb = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-x+0.5d0)/dble(n+1)*pi)

    end function grid_Val_Cheb


    !##############################################################################
    ! FUNCTION grid_Val_Grow
    !
    ! Calculates single gridpoint of a growing grid.
    !##############################################################################
    function grid_Val_Grow(x, left, right, growth, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! growth rate
        real(RK), intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! value at the grid point x \in [0,n]
        real(RK) :: grid_Val_Grow


        !##### OTHER VARIABLES ####################################################

        real(RK) :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Val', &
            'left interval point greater than right point')

        ! check for growth
        if(growth <= 0d0)call error('grid_Val_Grow', &
            'growth rate must be greater than zero')

        ! calculate factor
        h = (right-left)/((1+growth)**n-1)

        ! calculate grid value
        grid_Val_Grow = h*((1+growth)**x-1)+left

    end function grid_Val_Grow


    !##############################################################################
    ! FUNCTION grid_Inv_Equi_1
    !
    ! Calculates inverse of gridpoints of an equidistant grid.
    !##############################################################################
    function grid_Inv_Equi_1(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real(RK) :: grid_Inv_Equi_1


        !##### OTHER VARIABLES ####################################################

        real(RK) :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Equi', &
            'left interval point greater than right point')

        ! calculate distance between grid points
        h = (right-left)/n

        ! calculate grid value
        grid_Inv_Equi_1 = (x-left)/h

    end function grid_Inv_Equi_1


    !##############################################################################
    ! FUNCTION grid_Inv_Cheb_1
    !
    ! Calculates inverse of gridpoints of a Chebychev grid.
    !##############################################################################
    function grid_Inv_Cheb_1(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real(RK) :: grid_Inv_Cheb_1


        !##### OTHER VARIABLES ####################################################

        real(RK), parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Cheb', &
            'left interval point greater than right point')

        ! calculate grid value
        grid_Inv_Cheb_1 = dble(n) + 0.5d0 - acos((2d0*x-(left+right)) &
            /(right-left))*dble(n+1)/pi

    end function grid_Inv_Cheb_1


    !##############################################################################
    ! FUNCTION grid_Inv_Grow_1
    !
    ! Calculates inverse of gridpoints of a growing grid.
    !##############################################################################
    function grid_Inv_Grow_1(x, left, right, growth, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! growth rate
        real(RK), intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real(RK) :: grid_Inv_Grow_1


        !##### OTHER VARIABLES ####################################################

        real(RK) :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Grow', &
            'left interval point greater than right point')

        ! check for growth
        if(growth <= 0d0)call error('grid_Inv_Grow', &
            'growth rate must be greater than zero')

        ! calculate factor
        h = (right-left)/((1+growth)**n-1d0)

        ! calculate grid value
        grid_Inv_Grow_1 = log((x-left)/h+1d0)/log(1d0+growth)

    end function grid_Inv_Grow_1


    !##############################################################################
    ! FUNCTION grid_Inv_Equi_m
    !
    ! Calculates inverse of gridpoints of an equidistant grid.
    !##############################################################################
    function grid_Inv_Equi_m(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x(:)

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real(RK) :: grid_Inv_Equi_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK) :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Equi', &
            'left interval point greater than right point')

        ! calculate distance between grid points
        h = (right-left)/n

        ! calculate grid value
        grid_Inv_Equi_m = (x-left)/h

    end function grid_Inv_Equi_m


    !##############################################################################
    ! FUNCTION grid_Inv_Cheb_m
    !
    ! Calculates inverse of gridpoints of a Chebychev grid.
    !##############################################################################
    function grid_Inv_Cheb_m(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x(:)

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real(RK) :: grid_Inv_Cheb_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK), parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Cheb', &
            'left interval point greater than right point')

        ! calculate grid value
        grid_Inv_Cheb_m = dble(n) + 0.5d0 - acos((2d0*x-(left+right)) &
            /(right-left))*dble(n+1)/pi

    end function grid_Inv_Cheb_m


    !##############################################################################
    ! FUNCTION grid_Inv_Grow_m
    !
    ! Calculates inverse of gridpoints of a growing grid.
    !##############################################################################
    function grid_Inv_Grow_m(x, left, right, growth, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x(:)

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! growth rate
        real(RK), intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real(RK) :: grid_Inv_Grow_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK) :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Grow', &
            'left interval point greater than right point')

        ! check for growth
        if(growth <= 0d0)call error('grid_Inv_Grow', &
            'growth rate must be greater than zero')

        ! calculate factor
        h = (right-left)/((1d0+growth)**n-1d0)

        ! calculate grid value
        grid_Inv_Grow_m = log((x-left)/h+1d0)/log(1d0+growth)

    end function grid_Inv_Grow_m


    !##############################################################################
    ! subroutine linint_Equi_1
    !
    ! Calculates linear interpolant on an equidistant grid.
    !##############################################################################
    subroutine linint_Equi_1(x, left, right, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! left interpolation point
        integer(IK), intent(out) :: il

        ! right interpolation point
        integer(IK), intent(out) :: ir

        ! interpolation fraction
        real(RK), intent(out) :: phi

        !##### OTHER VARIABLES ####################################################

        real(RK) :: h, xinv, xl, xr


        !##### ROUTINE CODE #######################################################

        ! invert the grid to get point
        xinv = grid_Inv_Equi_1(min(max(x, left), right), left, right, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        h  = (right-left)/n
        xl = h*dble(il)+left
        xr = h*dble(ir)+left

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Equi_1


    !##############################################################################
    ! subroutine linint_Equi_m
    !
    ! Calculates linear interpolant on an equidistant grid.
    !##############################################################################
    subroutine linint_Equi_m(x, left, right, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x(:)

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! left interpolation point
        integer(IK), intent(out) :: il(:)

        ! right interpolation point
        integer(IK), intent(out) :: ir(:)

        ! interpolation fraction
        real(RK), intent(out) :: phi(:)

        !##### OTHER VARIABLES ####################################################

        integer(IK) :: m
        real(RK) :: h, xinv(size(x, 1)), xl(size(x, 1)), xr(size(x, 1))


        !##### ROUTINE CODE #######################################################

        ! check for sizes
        m = assert_eq(size(x, 1), size(il, 1), size(ir, 1), size(phi, 1), 'linint_Equi')
        m = m

        ! invert the grid to get point
        xinv = grid_Inv_Equi_m(min(max(x, left), right), left, right, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        h  = (right-left)/n
        xl = h*dble(il)+left
        xr = h*dble(ir)+left

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Equi_m


    !##############################################################################
    ! subroutine linint_Cheb_1
    !
    ! Calculates linear interpolant on an equidistant grid.
    !##############################################################################
    subroutine linint_Cheb_1(x, left, right, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! left interpolation point
        integer(IK), intent(out) :: il

        ! right interpolation point
        integer(IK), intent(out) :: ir

        ! interpolation fraction
        real(RK), intent(out) :: phi

        !##### OTHER VARIABLES ####################################################

        real(RK) :: xinv, xl, xr
        real(RK), parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! invert the grid to get point
        xinv = grid_Inv_Cheb_1(min(max(x, left), right), left, right, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        xl = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-dble(il)+0.5d0)/dble(n+1)*pi)
        xr = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-dble(ir)+0.5d0)/dble(n+1)*pi)

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Cheb_1


    !##############################################################################
    ! subroutine linint_Cheb_m
    !
    ! Calculates linear interpolant on a Chebychev grid.
    !##############################################################################
    subroutine linint_Cheb_m(x, left, right, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x(:)

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! left interpolation point
        integer(IK), intent(out) :: il(:)

        ! right interpolation point
        integer(IK), intent(out) :: ir(:)

        ! interpolation fraction
        real(RK), intent(out) :: phi(:)

        !##### OTHER VARIABLES ####################################################

        integer(IK) :: m
        real(RK) :: xinv(size(x, 1)), xl(size(x, 1)), xr(size(x, 1))
        real(RK), parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! check for sizes
        m = assert_eq(size(x, 1), size(il, 1), size(ir, 1), size(phi, 1), 'linint_Equi')
        m = m

        ! invert the grid to get point
        xinv = grid_Inv_Cheb_m(min(max(x, left), right), -1d0, 1d0, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        xl = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-dble(il)+0.5d0)/dble(n+1)*pi)
        xr = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-dble(ir)+0.5d0)/dble(n+1)*pi)

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Cheb_m


    !##############################################################################
    ! subroutine linint_Grow_1
    !
    ! Calculates linear interpolant on a growing grid.
    !##############################################################################
    subroutine linint_Grow_1(x, left, right, growth, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! growth rate
        real(RK), intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! left interpolation point
        integer(IK), intent(out) :: il

        ! right interpolation point
        integer(IK), intent(out) :: ir

        ! interpolation fraction
        real(RK), intent(out) :: phi

        !##### OTHER VARIABLES ####################################################

        real(RK) :: h, xinv, xl, xr


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('linint_Grow', &
            'left interval point greater than right point')

        ! invert the grid to get point
        xinv = grid_Inv_Grow_1(min(max(x, left), right), left, right, growth, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        h = (right-left)/((1+growth)**n-1)
        xl = h*((1+growth)**dble(il)-1d0)+left
        xr = h*((1+growth)**dble(ir)-1d0)+left

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Grow_1


    !##############################################################################
    ! subroutine linint_Grow_m
    !
    ! Calculates linear interpolant on an equidistant grid.
    !##############################################################################
    subroutine linint_Grow_m(x, left, right, growth, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real(RK), intent(in) :: x(:)

        ! left and right interval point
        real(RK), intent(in) :: left, right

        ! growth rate
        real(RK), intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer(IK), intent(in) :: n

        ! left interpolation point
        integer(IK), intent(out) :: il(:)

        ! right interpolation point
        integer(IK), intent(out) :: ir(:)

        ! interpolation fraction
        real(RK), intent(out) :: phi(:)

        !##### OTHER VARIABLES ####################################################

        integer(IK) :: m
        real(RK) :: h, xinv(size(x, 1)), xl(size(x, 1)), xr(size(x, 1))


        !##### ROUTINE CODE #######################################################

        ! check for sizes
        m = assert_eq(size(x, 1), size(il, 1), size(ir, 1), size(phi, 1), 'linint_Equi')
        m = m

        ! check for left <= right
        if(left >= right)call error('linint_Grow', &
            'left interval point greater than right point')

        ! invert the grid to get point
        xinv = grid_Inv_Grow_m(min(max(x, left), right), left, right, growth, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        h = (right-left)/((1+growth)**n-1)
        xl = h*((1+growth)**dble(il)-1d0)+left
        xr = h*((1+growth)**dble(ir)-1d0)+left

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Grow_m


    !##############################################################################
    ! function linint_Gen
    !
    ! For linear interpolation on irregular grids.
    !##############################################################################
    function linint_Gen(x, xi, yi, istart_in)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to evaluate
        real(RK), intent(in) :: x

        ! grid on which to evaluate
        real(RK), intent(in) :: xi(0:)

        ! data which to evaluate
        real(RK), intent(in) :: yi(0:)

        ! (optional) point where to start
        integer(IK), intent(in), optional :: istart_in

        ! output value
        real(RK) :: linint_Gen


        !##### OTHER VARIABLES ####################################################

        ! other variables
        integer(IK) :: ial, iar, n, istart
        real(RK) :: phi


        !##### ROUTINE CODE #######################################################

        n = assert_eq(size(xi,1), size(yi,1), 'linint_Gen')-1

        if(present(istart_in))then
            istart = min(max(istart_in, 0), n)
        else
            istart = n/2
        endif

        ! if grid value too large, search for first smaller point
        if(xi(istart) > x)then
            ial = istart
            do
                ial = ial - 1
                if(ial <= 0)exit
                if(xi(ial) <= x)exit
            enddo
            ial = max(ial, 0)
            ial = min(ial, n-1)
            iar = ial+1

        ! if grid value too small, search for first larger point
        else
            iar = istart
            do
                iar = iar + 1
                if(iar >= n)exit
                if(xi(iar) >= x)exit
            enddo
            iar = max(iar, 1)
            iar = min(iar, n)
            ial = iar-1
        endif

        ! linearly interpolate between the two points
        phi = 1d0 - (x-xi(ial))/(xi(iar)-xi(ial))
        linint_Gen = phi*yi(ial) + (1d0-phi)*yi(iar)

    end function





!##############################################################################
!##############################################################################
! MODULE matrixtools
!
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE lu_solve
    !
    ! Solves a linear equation system by lu-decomposition.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    subroutine lu_solve(a, b)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix of the system
        real(RK), intent(in) :: a(:, :)

        ! right side of equation and solution of the system
        real(RK), intent(inout) :: b(:)


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: indx(size(b))
        real(RK) :: worka(size(a, 1), size(a, 2))
        real(RK) :: d
        integer(IK) :: n


        !##### ROUTINE CODE #######################################################

        ! assert size equality
        n = assert_eq(size(a,1), size(a,2), size(b), 'lu_solve')
        n = n

        ! copy matrix to working matrix
        worka = a

        ! decompose matrix
        call lu_decomp(worka, indx, d)

        ! solve system
        call lu_back(worka, indx, b)

    end subroutine lu_solve


    !##############################################################################
    ! FUNCTION lu_invert
    !
    ! Inverts a matrix by lu-decomposition.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    function lu_invert(a)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix of the system
        real(RK), intent(in) :: a(:, :)

        ! the inverse of the matrix
        real(RK) :: lu_invert(size(a, 1), size(a, 2))


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: j, n


        !##### ROUTINE CODE #######################################################

        ! assert size equality
        n = assert_eq(size(a,1), size(a,2), 'lu_invert')

        ! set up unity matrix
        lu_invert = 0d0
        do j = 1, n
            lu_invert(j, j) = 1d0
        enddo

        ! succesively solve the system with unity matrix
        do j = 1, n
            call lu_solve(a, lu_invert(:, j))
        enddo

    end function lu_invert


    !##############################################################################
    ! SUBROUTINE lu_decomp
    !
    ! Calculates lu-decomposition of matrices.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    subroutine lu_decomp(a, indx, d)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix that shall be decomposed
        real(RK), intent(inout) :: a(:, :)

        ! row permutation indicator due to pivoting
        integer(IK), intent(out) :: indx(:)

        ! indicates whether number of row permutations was even or odd
        real(RK), intent(out) :: d


        !##### OTHER VARIABLES ####################################################

        real(RK) :: vv(size(a,1))
        real(RK), parameter :: tiny = 1.0e-20
        integer(IK) :: j, n, imax


        !##### ROUTINE CODE #######################################################

        ! check array sizes
        n = assert_eq(size(a,1), size(a,2), size(indx),'lu_decomp')

        ! initialize permutation indicator
        d = 1d0

        ! get maximum value in every row
        vv = maxval(abs(a), dim=2)

        ! if there is a zero row then matrix is singular
        if (any(abs(vv) <= 1d-100)) call error('lu_decomp', 'matrix is singular')

        ! invert v
        vv = 1d0/vv

        ! start lu-decomposition process
        do j = 1, n

            ! get index of pivot element
            imax = (j-1) + imaxloc(vv(j:n)*abs(a(j:n,j)))

            ! do pivoting if pivot element is not the first element
            if(j /= imax)then
                call swap(a(imax,:), a(j,:))
                d = -d
                vv(imax) = vv(j)
            endif

            ! indicate pivot element
            indx(j) = imax

            ! prevent division by 0
            if(abs(a(j,j))  <= 1d-100)call error('lu_decomp', 'matrix is singular')

            ! calculate new elements
            a(j+1:n, j) = a(j+1:n,j)/a(j,j)
            a(j+1:n, j+1:n) = a(j+1:n,j+1:n)-outerprod(a(j+1:n,j), a(j,j+1:n))
        enddo


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        function outerprod(a, b)

            real(RK), intent(in) :: a(:), b(:)
            real(RK) :: outerprod(size(a, 1),size(b, 1))

            outerprod = spread(a, dim=2, ncopies=size(b, 1)) * &
            spread(b, dim=1, ncopies=size(a, 1))

        end function outerprod


        subroutine swap(a, b)

            real(RK), intent(inout) :: a(:), b(:)
            real(RK) :: dum(size(a))

            dum = a
            a = b
            b = dum

        end subroutine swap


        function imaxloc(arr)

            real(RK), intent(in) :: arr(:)
            integer(IK) :: imaxloc
            integer(IK) :: imax(1)

            imax = maxloc(arr(:))
            imaxloc = imax(1)

        end function imaxloc

    end subroutine lu_decomp



    !##############################################################################
    ! SUBROUTINE lu_dec
    !
    ! Calculates lu-decomposition of matrices and returns L and U matrix.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    subroutine lu_dec(a, l, u)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix that shall be decomposed
        real(RK), intent(in) :: a(:, :)

        ! lower triangular matrix
        real(RK), intent(out) :: l(:, :)

        ! upper triangular martix
        real(RK), intent(out) :: u(:, :)


        !##### OTHER VARIABLES ####################################################

        real(RK) :: Awork(size(a,1), size(a,2))
        integer(IK) :: indx(size(a,1))
        real(RK) :: d
        integer(IK) :: n, j


        !##### ROUTINE CODE #######################################################

        ! check array sizes
        n = assert_eq((/size(a,1), size(a,2), size(l, 1), size(l, 2), &
        size(u, 1), size(u, 2)/), 'lu_dec')

        ! copy matrix
        Awork(:, :) = A(:, :)

        ! calculate decomposition
        call lu_decomp(Awork, indx, d)

        ! initialize matrices
        L(:, :) = 0d0
        U(:, :) = 0d0

        ! set up new matrices
        do j = 1, n

            ! diagonal element of L
            L(j, j) = 1d0

            ! other elements of L
            L(j, 1:j-1) = Awork(j, 1:j-1)

            ! elements of U
            U(j, j:n) = AWork(j, j:n)
        enddo

    end subroutine lu_dec


    !##############################################################################
    ! SUBROUTINE lu_back
    !
    ! Solves a lu decomposed linear equation system by backsubstitution.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    subroutine lu_back(a, indx, b)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! lu-decomposed matrix that defines system
        real(RK), intent(in) :: a(:, :)

        ! row permutation indicator due to pivoting
        integer(IK), intent(in) :: indx(:)

        ! right side of equation and solution of the system
        real(RK), intent(inout) :: b(:)


        !##### OTHER VARIABLES ####################################################

        real(RK) :: summ
        integer(IK) :: i, n, ii, ll


        !##### ROUTINE CODE #######################################################

        ! assert size equality
        n = assert_eq(size(a,1), size(a,2), size(indx), size(b), 'lu_back')

        ! start backward solving provess
        ii = 0
        do i = 1, n

            ll = indx(i)
            summ = b(ll)
            b(ll) = b(i)
            if(ii /= 0)then
                summ = summ-dot_product(a(i,ii:i-1), b(ii:i-1))
            elseif (abs(summ) >= 1d-100) then
                ii = i
            endif
            b(i)=summ
        enddo

        do i=n, 1, -1
            b(i) = (b(i)-dot_product(a(i,i+1:n), b(i+1:n)))/a(i,i)
        enddo

    end subroutine lu_back



    !##############################################################################
    ! SUBROUTINE cholesky
    !
    ! Calculates cholesky factorization of a symmetric matrix.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    subroutine cholesky(a, l)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the matrix that should be decomposed
        real(RK), intent(in) :: a(:, :)

        ! the cholesky factor
        real(RK), intent(out) :: l(:, :)


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: i, n
        real(RK) :: summ, p(size(a,1))


        !##### ROUTINE CODE #######################################################

        ! assert equalities
        n = assert_eq(size(a,1), size(a,2), size(l, 1), size(l, 2), &
            size(p), 'cholesky')

        ! copy matrix
        l = a

        ! decompose matrix
        do i = 1, n
            summ = l(i,i)-dot_product(l(i,1:i-1), l(i,1:i-1))
            if(summ <= 0d0)call error('cholesky', &
                'Cholesky decomposition failed')
            p(i) = sqrt(summ)
            l(i+1:n,i) = (l(i,i+1:n)-matmul(l(i+1:n,1:i-1),l(i,1:i-1)))/p(i)
        enddo

        ! copy matrix
        do i = 1, n
            l(i, i) = p(i)
            l(i, i+1:n) = 0d0
        enddo

    end subroutine cholesky















!##############################################################################
!##############################################################################
! MODULE probabilities
!##############################################################################
!##############################################################################


    !##############################################################################
    ! FUNCTION uniformPDF
    !
    ! Calculates uniform density at point x.
    !##############################################################################
    function uniformPDF(x, a, b)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! left end of the distribution
        real(RK), optional :: a

        ! right end of the distribution
        real(RK), optional :: b

        ! value of the uniform density
        real(RK) :: uniformPDF


        !##### OTHER VARIABLES ####################################################

        real(RK) :: a_c, b_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        a_c = 0d0
        if(present(a))a_c = a
        b_c = 1d0
        if(present(b))b_c = b

        if(b_c < a_c)then
            call error('uniformPDF','b is smaller than a')
        endif

        if(x < a_c .or. x > b_c)then
            uniformPDF = 0d0
        else
            uniformPDF = 1d0/(b_c-a_c)
        endif

    end function uniformPDF


    !##############################################################################
    ! FUNCTION uniformCDF
    !
    ! Calculates cumulated uniform distribution at point x.
    !##############################################################################
    function uniformCDF(x, a, b)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! left end of the distribution
        real(RK), optional :: a

        ! right end of the distribution
        real(RK), optional :: b

        ! value of the uniform distribution function
        real(RK) :: uniformCDF


        !##### OTHER VARIABLES ####################################################

        real(RK) :: a_c, b_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        a_c = 0d0
        if(present(a))a_c = a
        b_c = 1d0
        if(present(b))b_c = b

        if(b_c < a_c)then
            call error('uniformCDF','b is smaller than a')
        endif

        if(x <= a_c)then
            uniformCDF = 0d0
        elseif(x >= b_c)then
            uniformCDF = 1d0
        else
            uniformCDF = (x-a_c)/(b_c-a_c)
        endif

    end function uniformCDF


    !##############################################################################
    ! FUNCTION uniformCDF_Inv
    !
    ! Calculates inverse cumulated uniform distribution at point x.
    !##############################################################################
    function uniformCDF_Inv(x, a, b)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! left end of the distribution
        real(RK), optional :: a

        ! right end of the distribution
        real(RK), optional :: b

        ! value of the uniform distribution function
        real(RK) :: uniformCDF_Inv


        !##### OTHER VARIABLES ####################################################

        real(RK) :: a_c, b_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        a_c = 0d0
        if(present(a))a_c = a
        b_c = 1d0
        if(present(b))b_c = b

        if(b_c < a_c)then
            call error('uniformCDF_Inv','b is smaller than a')
        endif

        if(x <= 0d0)then
            uniformCDF_Inv = a_c
        elseif(x >= 1d0)then
            uniformCDF_Inv = b_c
        else
            uniformCDF_Inv = a_c + x*(b_c-a_c)
        endif

    end function uniformCDF_Inv


    !##############################################################################
    ! FUNCTION normalPDF
    !
    ! Calculates normal density functions at point x.
    !##############################################################################
    function normalPDF(x, mu, sigma)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! value of normal density at p
        real(RK) :: normalPDF


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c
        real(RK), parameter :: pi = 3.1415926535897d0


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sqrt(sigma)

        if(sigma_c <= 0d0)then
            call error('normalCDF','sigma has zero or negative value')
        endif

        normalPDF = 1d0/(sigma_c*sqrt(2d0*pi))*exp(-((x-mu_c)/sigma_c)**2/2d0)

    end function normalPDF


    !##############################################################################
    ! FUNCTION normalCDF
    !
    ! Calculates cumulated normal distribution at point x.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Fortran Code by John Burkardt available as Algorithm ASA066 from
    !     https://people.sc.fsu.edu/~jburkardt/f_src/asa066/asa066.html
    !
    !     REFERENCE: Hill, D. (1973). Algorithm AS 66: The Normal Integral.
    !                Applied Statistics, 22(3), 424-427.
    !##############################################################################
    function normalCDF(x, mu, sigma)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! value of the normal distribution at x
        real(RK) :: normalCDF


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c, xtrans, xabs, y
        real(RK), parameter :: a0  = 0.5d0
        real(RK), parameter :: a1  = 0.398942280444d0
        real(RK), parameter :: a2  = 0.399903438504d0
        real(RK), parameter :: a3  = 5.75885480458d0
        real(RK), parameter :: a4  = 29.8213557808d0
        real(RK), parameter :: a5  = 2.62433121679d0
        real(RK), parameter :: a6  = 48.6959930692d0
        real(RK), parameter :: a7  = 5.92885724438d0
        real(RK), parameter :: b0  = 0.398942280385d0
        real(RK), parameter :: b1  = 3.8052d-8
        real(RK), parameter :: b2  = 1.00000615302d0
        real(RK), parameter :: b3  = 3.98064794d-4
        real(RK), parameter :: b4  = 1.98615381364d0
        real(RK), parameter :: b5  = 0.151679116635d0
        real(RK), parameter :: b6  = 5.29330324926d0
        real(RK), parameter :: b7  = 4.8385912808d0
        real(RK), parameter :: b8  = 15.1508972451d0
        real(RK), parameter :: b9  = 0.742380924027d0
        real(RK), parameter :: b10 = 30.789933034d0
        real(RK), parameter :: b11 = 3.99019417011d0


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))then
            if(sigma_c > 0d0)sigma_c = sqrt(sigma)
        endif

        if(sigma_c <= 0d0)then
            call error('normalCDF','sigma has zero or negative value')
        endif

        ! standardize evaluation point
        xtrans = (x - mu_c)/sigma_c

        ! calculate absolute value and quadratic
        xabs = abs(xtrans)
        y = a0*xtrans**2

        ! choose the right interval for calculation
        if(xabs <= 1.28d0)then
            normalCDF = a0-xabs*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y+a7))))
        elseif(xabs <= 12.7d0)then
            normalCDF = b0*exp(-y)/(xabs-b1+b2/(xabs+b3+b4/(xabs-b5+b6/(xabs+b7-b8/ &
                    (xabs+b9+b10/(xabs+b11))))))
        else
            normalCDF = 0d0
        endif

        ! transform if other side of the bell
        if(xtrans > 0d0)normalCDF = 1d0-normalCDF

    end function normalCDF


    !##############################################################################
    ! FUNCTION normalCDF_Inv
    !
    ! Calculates inverse cumulated normal distribution at point x.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Fortran Code by John Burkardt available as Algorithm ASA241 from
    !     https://people.sc.fsu.edu/~jburkardt/f_src/asa241/asa241.html
    !
    !     REFERENCE: Wichura, M. (1988). Algorithm AS 241: The Percentage Points
    !                of the Normal Distribution. Applied Statistics, 37(3),
    !                477-484.
    !##############################################################################
    function normalCDF_Inv(x, mu, sigma)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK) :: x

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! value of the inverse of the normal distribution at x
        real(RK) :: normalCDF_Inv


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c, q, r
        real(RK), parameter :: a(8) = (/  3.3871328727963666080d0, &
                                        1.3314166789178437745d2, &
                                        1.9715909503065514427d3, &
                                        1.3731693765509461125d4, &
                                        4.5921953931549871457d4, &
                                        6.7265770927008700853d4, &
                                        3.3430575583588128105d4, &
                                        2.5090809287301226727d3 /)
        real(RK), parameter :: b(8) = (/  1.0d0, &
                                        4.2313330701600911252d1, &
                                        6.8718700749205790830d2, &
                                        5.3941960214247511077d3, &
                                        2.1213794301586595867d4, &
                                        3.9307895800092710610d4, &
                                        2.8729085735721942674d4, &
                                        5.2264952788528545610d3 /)
        real(RK), parameter :: c(8) = (/  1.42343711074968357734d0, &
                                        4.63033784615654529590d0, &
                                        5.76949722146069140550d0, &
                                        3.64784832476320460504d0, &
                                        1.27045825245236838258d0, &
                                        2.41780725177450611770d-1, &
                                        2.27238449892691845833d-2, &
                                        7.74545014278341407640d-4 /)
        real(RK), parameter :: const1 = 0.180625d0
        real(RK), parameter :: const2 = 1.6d0
        real(RK), parameter :: d(8) = (/  1.0d0, &
                                        2.05319162663775882187d0, &
                                        1.67638483018380384940d0, &
                                        6.89767334985100004550d-1, &
                                        1.48103976427480074590d-1, &
                                        1.51986665636164571966d-2, &
                                        5.47593808499534494600d-4, &
                                        1.05075007164441684324d-9 /)
        real(RK), parameter :: e(8) = (/  6.65790464350110377720d0, &
                                        5.46378491116411436990d0, &
                                        1.78482653991729133580d0, &
                                        2.96560571828504891230d-1, &
                                        2.65321895265761230930d-2, &
                                        1.24266094738807843860d-3, &
                                        2.71155556874348757815d-5, &
                                        2.01033439929228813265d-7 /)
        real(RK), parameter :: f(8) = (/   1.0d0, &
                                        5.99832206555887937690d-1, &
                                        1.36929880922735805310d-1, &
                                        1.48753612908506148525d-2, &
                                        7.86869131145613259100d-4, &
                                        1.84631831751005468180d-5, &
                                        1.42151175831644588870d-7, &
                                        2.04426310338993978564d-15 /)

        real(RK), parameter :: split1 = 0.425D+00
        real(RK), parameter :: split2 = 5.0D+00


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))then
            if(sigma_c > 0d0)sigma_c = sqrt(sigma)
        endif

        if(sigma_c <= 0d0)then
            call error('normalCDF_Inv','sigma has zero or negative value')
        endif

        ! check for very small or large values
        if(x <= 0d0)then
            normalcdf_Inv = - huge(x)
            return
        endif

        if (x >= 1d0)then
            normalcdf_Inv = huge(x)
            return
        endif

        ! calculate for reasonable values
        q = x - 0.5d0

        if(abs(q) <= split1)then
            r = const1 - q*q
            normalcdf_inv = q*r8poly_value(8, a, r)/r8poly_value(8, b, r)
        else

            if (q < 0d0) then
                r = x
            else
                r = 1d0 - x
            endif

            r = sqrt(-log(r))

            if(r <= split2)then
                r = r - const2
                normalcdf_Inv = r8poly_value(8, c, r)/r8poly_value(8, d, r)
            else
                r = r - split2
                normalcdf_Inv = r8poly_value(8, e, r)/r8poly_value (8, f, r)
            endif

            if(q < 0d0)then
                normalcdf_Inv = -normalcdf_inv
            endif
        endif

        ! transfer to mu and sigma
        normalCDF_Inv = mu_c + sigma_c*normalCDF_Inv


    contains

        ! stable evaluation of a polynomial
        function r8poly_value(n, a, x)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! degree of polynomial
            integer(IK), intent(in) :: n

            ! polynomial coefficients
            real(RK), intent(in) :: a(n)

            ! point where to evaluate the polynomial
            real(RK), intent(in) :: x

            ! the value of the polynomial
            real(RK) :: r8poly_value


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: i


            !##### ROUTINE CODE #######################################################

            r8poly_value = 0d0
            do i = n, 1, -1
                r8poly_value = r8poly_value * x + a(i)
            enddo

        end function

    end function normalCDF_Inv


    !##############################################################################
    ! FUNCTION log_normalPDF
    !
    ! Calculates log-normal density functions at point x.
    !##############################################################################
    function log_normalPDF(x, mu, sigma)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! value of normal density at p
        real(RK) :: log_normalPDF


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = exp(0.5d0)
        if(present(mu))mu_c = mu
        sigma_c = exp(1d0)*(exp(1d0)-1d0)
        if(present(sigma))sigma_c = sigma

        if(sigma_c <= 0d0)then
            call error('log_normalPDF','sigma has zero or negative value')
        endif
        if(mu_c <= 0d0)then
            call error('log_normalPDF','mu has zero or negative value')
        endif

        if(x < 0d0)then
            call error('log_normalPDF','x has negative value')
        endif

        ! get expectation and variance
        sigma_c = log(1d0+sigma_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sigma_c

        ! simulate normal and convert to log_normal
        if(x > 0d0)then
            log_normalPDF = normalPDF(log(x), mu_c, sigma_c)
        else
            log_normalPDF = 0d0
        endif

    end function log_normalPDF


    !##############################################################################
    ! FUNCTION log_normalCDF
    !
    ! Calculates log-normal distribution at point x.
    !##############################################################################
    function log_normalCDF(x, mu, sigma)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! value of normal density at p
        real(RK) :: log_normalCDF


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = exp(0.5d0)
        if(present(mu))mu_c = mu
        sigma_c = exp(1d0)*(exp(1d0)-1d0)
        if(present(sigma))sigma_c = sigma

        if(sigma_c <= 0d0)then
            call error('log_normalCDF','sigma has zero or negative value')
        endif
        if(mu_c <= 0d0)then
            call error('log_normalCDF','mu has zero or negative value')
        endif

        if(x < 0d0)then
            call error('log_normalCDF','x has negative value')
        endif

        ! get expectation and variance
        sigma_c = log(1d0+sigma_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sigma_c

        ! simulate normal and convert to log_normal
        if(x > 0d0)then
            log_normalCDF = normalCDF(log(x), mu_c, sigma_c)
        else
            log_normalCDF = 0d0
        endif

    end function log_normalCDF


    !##############################################################################
    ! FUNCTION log_normalCDF_Inv
    !
    ! Calculates cumulated log-normal distribution at point x.
    !##############################################################################
    function log_normalCDF_Inv(x, mu, sigma)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! value of normal density at p
        real(RK) :: log_normalCDF_Inv


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = exp(0.5d0)
        if(present(mu))mu_c = mu
        sigma_c = exp(1d0)*(exp(1d0)-1d0)
        if(present(sigma))sigma_c = sigma

        if(sigma_c <= 0d0)then
            call error('log_normalCDF_Inv','sigma has zero or negative value')
        endif
        if(mu_c <= 0d0)then
            call error('log_normalCDF_Inv','mu has zero or negative value')
        endif

        ! get expectation and variance
        sigma_c = log(1d0+sigma_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sigma_c

        ! simulate normal and convert to log_normal
        log_normalCDF_Inv = normalCDF_Inv(x, mu_c, sigma_c)
        log_normalCDF_Inv = exp(log_normalCDF_Inv)

    end function log_normalCDF_Inv


    !##############################################################################
    ! FUNCTION GammaPDF
    !
    ! Calculates Gamma density functions at point x.
    !##############################################################################
    function GammaPDF(x, alpha, beta)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! shape parameter of the distribution
        real(RK), optional :: alpha

        ! scale parameter of the distribution
        real(RK), optional :: beta

        ! value of Gamma density at p
        real(RK) :: gammaPDF


        !##### OTHER VARIABLES ####################################################

        real(RK) :: alpha_c, beta_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        alpha_c = 1d0
        if(present(alpha))alpha_c = alpha
        beta_c = 1d0
        if(present(beta))beta_c = beta

        ! check for validity of parameters
        if(alpha_c <= 0d0)then
            call error('GammaPDF','alpha has a non-positive value')
        endif
        if(beta_c <= 0d0)then
            call error('GammaPDF','beta has a non-positive value')
        endif

        if(x < 0d0 .or. alpha_c < 1d0 .and. x <= 0d0)then
            GammaPDF = 0d0
        else
            GammaPDF = beta_c**alpha_c/gamma_function(alpha_c)* &
                        x**(alpha_c-1d0)*exp(-beta_c*x)
        endif

    end function GammaPDF


    !##############################################################################
    ! FUNCTION GammaCDF
    !
    ! Calculates cumulated Gamma distribution at point x.
    !##############################################################################
    function GammaCDF(x, alpha, beta)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! shape parameter of the distribution
        real(RK), optional :: alpha

        ! scale parameter of the distribution
        real(RK), optional :: beta

        ! value of cumulated Gamma distribution at p
        real(RK) :: gammaCDF


        !##### OTHER VARIABLES ####################################################

        real(RK) :: alpha_c, beta_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        alpha_c = 1d0
        if(present(alpha))alpha_c = alpha
        beta_c = 1d0
        if(present(beta))beta_c = beta

        ! check for validity of parameters
        if(alpha_c <= 0d0)then
            call error('GammaCDF','alpha has a non-positive value')
        endif
        if(beta_c <= 0d0)then
            call error('GammaCDF','beta has a non-positive value')
        endif

        if(x <= 0d0)then
            GammaCDF = 0d0
        else
            GammaCDF = incomplete_gamma(beta_c*x, alpha_c)
        endif

    end function gammaCDF



    !##############################################################################
    ! FUNCTION betaPDF
    !
    ! Calculates beta density functions at point x.
    !##############################################################################
    function betaPDF(x, p, q)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! parameter of the distribution
        real(RK), optional :: p

        ! parameter of the distribution
        real(RK), optional :: q

        ! value of Gamma density at p
        real(RK) :: betaPDF


        !##### OTHER VARIABLES ####################################################

        real(RK) :: p_c, q_c, betanorm


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        p_c = 1d0
        if(present(p))p_c = p
        q_c = 1d0
        if(present(q))q_c = q

        ! check for validity of parameters
        if(p_c <= 0d0)then
            call error('betaPDF','p has a non-positive value')
        endif
        if(q_c <= 0d0)then
            call error('betaPDF','q has a non-positive value')
        endif

        if(x < 0d0 .or. x > 1d0)then
            betaPDF = 0d0
        else
            betanorm = exp(my_log_gamma(p_c) + my_log_gamma(q_c) - my_log_gamma(p_c+q_c))
            betaPDF = x**(p_c-1d0)*(1d0-x)**(q_c-1d0)/betanorm
        endif

    end function betaPDF


    !##############################################################################
    ! FUNCTION betaCDF
    !
    ! Calculates cumulated beta distribution at point x.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Fortran Code by John Burkardt available as Algorithm ASA063 from
    !     https://people.sc.fsu.edu/~jburkardt/f_src/asa063/asa063.html
    !
    !     REFERENCE: Majumder, K.L. & Bhattacharjee, G.P. (1973). Algorithm AS 63:
    !                The incomplete Beta Integral, Applied Statistics, 22(3),
    !                409-411.
    !##############################################################################
    function betaCDF(x, p_in, q_in)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! parameter of the distribution
        real(RK), optional :: p_in

        ! parameter of the distribution
        real(RK), optional :: q_in

        ! value of Gamma density at p
        real(RK) :: betaCDF


        !##### OTHER VARIABLES ####################################################

        real(RK), parameter :: acu = 0.1d-14
        real(RK) :: p, q, ai, beta_log, cx, pp, psq, qq, rx, temp, term, xx
        integer(IK) :: ns
        logical :: indx


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        p = 1d0
        if(present(p_in))p = p_in
        q = 1d0
        if(present(q_in))q = q_in

        ! check for validity of parameters
        if(p <= 0d0)then
            call error('betaPDF','p has a non-positive value')
        endif
        if(q <= 0d0)then
            call error('betaPDF','q has a non-positive value')
        endif

        ! check outside of range
        if(x <= 0d0)then
            betaCDF = 0d0
            return
        endif
        if(x >= 1d0)then
            betaCDF = 1d0
            return
        endif

        ! calculate logarithm of the complete beta function
        beta_log = my_log_gamma(p) + my_log_gamma(q) - my_log_gamma(p + q)

        psq = p + q
        cx = 1d0 - x

        if(p < psq*x)then
            xx = cx
            cx = x
            pp = q
            qq = p
            indx = .true.
        else
            xx = x
            pp = p
            qq = q
            indx = .false.
        endif

        term = 1d0
        ai = 1d0
        betaCDF = 1d0
        ns = int(qq + cx*psq)

        ! use Soper's reduction formula.
        rx = xx / cx
        temp = qq - ai
        if (ns == 0) then
            rx = xx
        endif

        do

            term = term*temp*rx / (pp + ai)
            betaCDF = betaCDF + term
            temp = abs (term)

            if(temp <= acu .and. temp <= acu*betaCDF) then

                betaCDF = betaCDF*exp (pp*log(xx) &
                    + (qq-1d0)*log(cx)-beta_log)/pp

                ! change if other tail
                if(indx)then
                    betaCDF = 1d0 - betaCDF
                endif
                exit
            endif

            ai = ai + 1d0
            ns = ns - 1

            if(0 <= ns)then
                temp = qq - ai
                if(ns == 0)then
                    rx = xx
                endif
            else
                temp = psq
                psq = psq + 1d0
            endif
        enddo

    end function betaCDF


    !##############################################################################
    ! FUNCTION incomplete_gamma
    !
    ! Calculates the incomplete gamma integral.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Fortran Code by John Burkardt available as Algorithm ASA239 from
    !     https://people.sc.fsu.edu/~jburkardt/f_src/asa239/asa239.html
    !
    !     REFERENCE: Shea, B. (1988). Algorithm AS 239: Chi-squared and Incomplete
    !                Gamma Integral, Applied Statistics, 37(3), 466-473.
    !##############################################################################
    function incomplete_gamma(x, p)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! parameter of the gamma function
        real(RK), intent(in) :: p

        ! value of normal density at p
        real(RK) :: incomplete_gamma


        !##### OTHER VARIABLES ####################################################

        real(RK) :: a, an, arg, b, c, pn1, pn2, pn3, pn4, pn5, pn6, rn
        real(RK), parameter :: elimit = - 88d0
        real(RK), parameter :: oflo = 1d37
        real(RK), parameter :: plimit = 1000d0
        real(RK), parameter :: tol = 1d-14
        real(RK), parameter :: xbig = 1d8


        !##### ROUTINE CODE #######################################################

        incomplete_gamma = 0d0

        ! check validity of inputs solution
        if(x < 0d0)then
            call error('incomeplete_gamma','x is smaller than zero')
        endif

        if(p <= 0d0)then
            call error('incomeplete_gamma','p is not positive')
        endif

        ! set incomeplete_gamma = 0 if x is zero
        if(x <= 0d0)then
            incomplete_gamma = 0d0
            return
        endif

        ! normal approximation for large values of p
        if(p > plimit)then
            pn1 = 3d0*sqrt(p)*((x/p)**(1d0/3d0) + 1d0/(9d0*p) - 1d0)
            incomplete_gamma = normalCDF(pn1)
            return
        endif

        ! set incomeplete_gamma = 1 if x is large
        if(x > xbig)then
            incomplete_gamma = 1d0
        endif

        ! Pearson series expansion
        if(x <= 1d0 .or. x < p)then
            arg = p*log(x) - x - my_log_gamma(p+1d0)
            c = 1d0
            incomplete_gamma = 1.0D+00
            a = p

            do
                a = a + 1d0
                c = c*x/a
                incomplete_gamma = incomplete_gamma + c

                if(c <= tol) then
                    exit
                endif

            enddo

            arg = arg + log (incomplete_gamma)

            if ( elimit <= arg ) then
              incomplete_gamma = exp ( arg )
            else
              incomplete_gamma = 0d0
            endif

        ! Use a continued fraction expansion.
        else

            arg = p*log (x) - x - my_log_gamma(p)
            a = 1d0 - p
            b = a + x + 1d0
            c = 0d0
            pn1 = 1d0
            pn2 = x
            pn3 = x + 1d0
            pn4 = x*b
            incomplete_gamma = pn3/pn4

            do

                a = a + 1d0
                b = b + 2d0
                c = c + 1d0
                an = a*c
                pn5 = b*pn3 - an*pn1
                pn6 = b*pn4 - an*pn2

                if ( abs(pn6) > 0d0 ) then

                    rn = pn5 / pn6

                    if ( abs(incomplete_gamma - rn) <= min(tol, tol*rn)) then
                        exit
                    endif
                    incomplete_gamma = rn
                endif

                pn1 = pn3
                pn2 = pn4
                pn3 = pn5
                pn4 = pn6

                !  Re-scale terms in continued fraction if terms are large.
                if (oflo <= abs(pn5) ) then
                    pn1 = pn1/oflo
                    pn2 = pn2/oflo
                    pn3 = pn3/oflo
                    pn4 = pn4/oflo
                endif

            enddo

            arg = arg + log (incomplete_gamma)

            if (arg >= elimit) then
              incomplete_gamma = 1d0 - exp (arg)
            else
              incomplete_gamma = 1d0
            endif
        endif

    end function incomplete_gamma


    !##############################################################################
    ! FUNCTION gamma_function
    !
    ! Calculates the gamma fuction.
    !##############################################################################
    function gamma_function(x)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x

        ! value of log of the gamma function
        real(RK) :: gamma_function


        !##### ROUTINE CODE #######################################################


        ! check validity of inputs solution
        if(x < 0d0)then
            call error('gamma_function','x is smaller than zero')
        endif

        gamma_function = my_log_gamma(x)
        gamma_function = exp(gamma_function)

    end function gamma_function


    !##############################################################################
    ! FUNCTION my_log_gamma
    !
    ! Calculates log of the gamma fuction.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Fortran Code by John Burkardt available as Algorithm AS245 from
    !     https://people.sc.fsu.edu/~jburkardt/f_src/asa245/asa245.html
    !
    !     REFERENCE: Macleod, A.J. (1989). Algorithm AS 245: A Robust and Reliable
    !                Algorithm for the  Logarithm of the Gamma Function. Applied
    !                Statistics, 38(2), 397-402.
    !##############################################################################
    function my_log_gamma(x_in)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real(RK), intent(in) :: x_in

        ! value of log of the gamma function
        real(RK) :: my_log_gamma


        !##### OTHER VARIABLES ####################################################

        real(RK), parameter :: alr2pi = 9.18938533204673d-1
        real(RK), parameter :: xlge = 5.10d6
        real(RK), parameter :: r1(9) = (/ -2.66685511495d0, &
                                        -2.44387534237d1, &
                                        -2.19698958928d1, &
                                         1.11667541262d1, &
                                         3.13060547623d0, &
                                         6.07771387771d-1, &
                                         1.19400905721d1, &
                                         3.14690115749d1, &
                                         1.52346874070d1 /)
        real(RK), parameter :: r2(9) = (/ -7.83359299449d1, &
                                        -1.42046296688d2, &
                                         1.37519416416d2, &
                                         7.86994924154d1, &
                                         4.16438922228d0, &
                                         4.70668766060d1, &
                                         3.13399215894d2, &
                                         2.63505074721d2, &
                                         4.33400022514d1 /)
        real(RK), parameter :: r3(9) = (/ -2.12159572323d5, &
                                         2.30661510616d5, &
                                         2.74647644705d4, &
                                        -4.02621119975d4, &
                                        -2.29660729780d3, &
                                        -1.16328495004d5, &
                                        -1.46025937511d5, &
                                        -2.42357409629d4, &
                                        -5.70691009324d2 /)
        real(RK), parameter :: r4(5) = (/  2.79195317918525d-1, &
                                         4.917317610505968d-1, &
                                         6.92910599291889d-2, &
                                         3.350343815022304d0, &
                                         6.012459259764103d0 /)
        real(RK) :: x, x1, x2, y


        !##### ROUTINE CODE #######################################################


        my_log_gamma = 0d0
        x = x_in

        ! check validity of inputs solution
        if(x < 0d0)then
            call error('my_log_gamma','x is smaller than zero')
        endif

        ! get solution for 0 < X < 0.5 and 0.5 <= x < 1.5
        if(x < 1.5d0)then
            if(x < 0.5d0)then
                my_log_gamma = -log(x)
                y = x + 1d0

                ! return if x is smaller than machine epsilon
                if(y <= 1d0)return
            else
                my_log_gamma = 0d0
                y = x
                x = (x-0.5d0) - 0.5d0
            endif
            my_log_gamma = my_log_gamma + x * ((((r1(5)*y + r1(4))*y + r1(3))*y &
                + r1(2))*y + r1(1)) / ((((y + r1(9))*y + r1(8))*y + r1(7))*y + r1(6))


        ! get solution for 1.5 <= x < 4.0
        elseif(x < 4.0d0)then
            y = (x - 1d0) - 1d0
            my_log_gamma = y * ((((r2(5)*x + r2(4))*x + r2(3))*x + r2(2))*x &
                + r2(1)) / ((((x + r2(9))*x + r2(8))*x + r2(7))*x + R2(6))

        ! get solution for 4.0 <= x < 12
        elseif(x < 12d0)then
            my_log_gamma = ((((r3(5)*x + r3(4))*x + r3(3))*x + r3(2))*x + r3(1)) / &
                ((((x + r3(9))*x + r3(8))*x + r3(7))*x + r3(6))
        else
            y = log(x)
            my_log_gamma = x * (y -1d0) - 0.5d0*y + alr2pi
            if(x <= XLGE)then
                x1 = 1d0/x
                x2 = x1 * x1
                my_log_gamma = my_log_gamma + x1 * ((r4(3)*x2 + r4(2))*x2 + r4(1)) / &
                    ((x2 + r4(5))*x2 + r4(4))
            endif
        endif

    end function my_log_gamma


    !##############################################################################
    ! FUNCTION bernoulliPDF
    !
    ! Calculates probabilities of the bernoulli distribution.
    !##############################################################################
    function bernoulliPDF(k, p)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! zero or 1 draw
        integer(IK), intent(in) :: k

        ! the probability of the positive draw
        real(RK), intent(in) :: p

        ! the probability of this to happen
        real(RK) :: bernoulliPDF


        !##### ROUTINE CODE #######################################################

        ! check for validity
        if(p < 0d0 .or. p > 1d0)then
            call error('bernoulliPDF', 'Probability p is out of [0,1]')
        endif

        ! outside of range
        if(k < 0)then
            bernoulliPDF = 0d0
            return
        endif

        if(k > 1)then
            bernoulliPDF = 1d0
            return
        endif

        if(k == 0)then
            bernoulliPDF = 1d0 - p
        else
            bernoulliPDF = p
        endif

    end function


    !##############################################################################
    ! FUNCTION bernoulliCDF
    !
    ! Calculates cumulated probabilities of the bernoulli distribution.
    !##############################################################################
    function bernoulliCDF(k, p)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! zero or 1 draw
        integer(IK), intent(in) :: k

        ! the probability of the positive draw
        real(RK), intent(in) :: p

        ! the cumulated probability of this to happen
        real(RK) :: bernoulliCDF


        !##### ROUTINE CODE #######################################################

        ! check for validity
        if(p < 0d0 .or. p > 1d0)then
            call error('bernoulliCDF', 'Probability p is out of [0,1]')
        endif

        ! outside of range
        if(k < 0)then
            bernoulliCDF = 0d0
            return
        endif

        if(k >= 1)then
            bernoulliCDF = 1d0
            return
        endif

        bernoulliCDF = 1d0 - p

    end function


    !##############################################################################
    ! FUNCTION binomialPDF
    !
    ! Calculates probabilities of the binomial distribution.
    !##############################################################################
    function binomialPDF(k, n, p)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the number of positive draws
        integer(IK), intent(in) :: k

        ! the total number of draws
        integer(IK), intent(in) :: n

        ! the probability of a positive draw
        real(RK), intent(in) :: p

        ! the probability of this to happen
        real(RK) :: binomialPDF


        !##### ROUTINE CODE #######################################################


        ! check for validity
        if(p < 0d0 .or. p > 1d0)then
            call error('binomialPDF', 'Probability p is out of [0,1]')
        endif

        ! outside of range
        if(k < 0 .or. k > n)then
            binomialPDF = 0d0
            return
        endif

        binomialPDF = binomial_coefficient(n, k)*p**k*(1d0-p)**(n-k)

    end function


    !##############################################################################
    ! FUNCTION binomialCDF
    !
    ! Calculates cumulated probabilities of the binomial distribution.
    !##############################################################################
    function binomialCDF(k, n, p)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the number of positive draws
        integer(IK), intent(in) :: k

        ! the total number of draws
        integer(IK), intent(in) :: n

        ! the probability of a positive draw
        real(RK), intent(in) :: p

        ! the cumulated probability of this to happen
        real(RK) :: binomialCDF


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: ii


        !##### ROUTINE CODE #######################################################

        ! check for validity
        if(p < 0d0 .or. p > 1d0)then
            call error('binomialCDF', 'Probability p is out of [0,1]')
        endif

        ! outside of range
        if(k < 0)then
            binomialCDF = 0d0
            return
        endif

        if(k > n)then
            binomialCDF = 1d0
            return
        endif

        binomialCDF = 0d0
        do ii = 0, k
            binomialCDF = binomialCDF +  binomial_coefficient(n, ii)*p**ii*(1d0-p)**(n-ii)
        enddo

    end function


    !##############################################################################
    ! FUNCTION binomial_coefficient
    !
    ! Calculate binomial coefficients efficiently
    !##############################################################################
    function binomial_coefficient(n, k)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! out of number
        integer(IK), intent(in) :: n

        ! number
        integer(IK), intent(in) :: k

        ! the probability of this to happen
        real(RK) :: binomial_coefficient


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n0, k0, ii


        !##### ROUTINE CODE #######################################################


        if(k == 0 .or. k == n)then
            binomial_coefficient = 1d0
        endif

        ! copy values
        n0 = n
        k0 = k

        ! use symmetry of binomial coefficients
        if(k0 > n0 - k0)then
            k0 = n0 - k0
        endif

        ! use factorial formula for stability
        binomial_coefficient = 1d0
        do ii = 1, k
            binomial_coefficient = binomial_coefficient*dble(n+1-ii)
            binomial_coefficient = binomial_coefficient/dble(ii)
        enddo

    end function



    !##############################################################################
    ! SUBROUTINE init_random_seed
    !
    ! To ensure that each random draw is a new sequence
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !    The GNU GCC Fortran compiler library documentation
    !    https://gcc.gnu.org/onlinedocs/gcc-4.8.5/gfortran/RANDOM_005fSEED.html
    !##############################################################################
    subroutine init_random_seed(fixed)

        use Constants_mod, only: IK, RK; implicit none
        integer(IK), allocatable :: seed(:)
        integer(IK) :: i, n, dt(8)
        integer(IK), parameter :: int64 = selected_int_kind(16)
        integer(kind=int64) :: t
        logical, optional :: fixed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(t)
        if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24_int64 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
        endif

        if(present(fixed))then
            if(fixed)t = 0
        endif
        do i = 1, n
            seed(i) = lcg(t)
        enddo
        call random_seed(put=seed)

    contains

        ! This simple PRNG might not be good enough for real work, but is
        ! sufficient for seeding a better PRNG.
        function lcg(s)

            use Constants_mod, only: IK, RK; implicit none
            integer(IK) :: lcg
            integer(int64) :: s

            if (s == 0) then
                s = 104729
            else
                s = mod(s, 4294967296_int64)
            endif
            s = mod(s * 279470273_int64, 4294967291_int64)
            lcg = int(mod(s, int(huge(0), int64)), kind(0))

        end function lcg

    end subroutine init_random_seed


    !##############################################################################
    ! SUBROUTINE simulate_uniform_1
    !
    ! Simulates one draw from a uniform distribution.
    !##############################################################################
    subroutine simulate_uniform_1(x, a, b, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x

        ! left end of the distribution
        real(RK), optional :: a

        ! right end of the distribution
        real(RK), optional :: b

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: a_c, b_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        a_c = 0d0
        if(present(a))a_c = a
        b_c = 1d0
        if(present(b))b_c = b

        ! initialize the random seed
        if(tbox_seed)then
            if(present(fixed))then
                call init_random_seed(fixed)
            else
                call init_random_seed()
            endif
            tbox_seed = .false.
        endif

        ! draw the random number
        call random_number(x)

        x = a_c + (b_c-a_c)*x

    end subroutine simulate_uniform_1


    !##############################################################################
    ! SUBROUTINE simulate_uniform_n
    !
    ! Simulates a series draw from a uniform distribution.
    !##############################################################################
    subroutine simulate_uniform_n(x, a, b, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x(:)

        ! left end of the distribution
        real(RK), optional :: a

        ! right end of the distribution
        real(RK), optional :: b

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: a_c, b_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        a_c = 0d0
        if(present(a))a_c = a
        b_c = 1d0
        if(present(b))b_c = b

        ! initialize the random seed
        if(tbox_seed)then
            if(present(fixed))then
                call init_random_seed(fixed)
            else
                call init_random_seed()
            endif
            tbox_seed = .false.
        endif

        call random_number(x)

        x = a_c + (b_c-a_c)*x

    end subroutine simulate_uniform_n


    !##############################################################################
    ! Subroutine simulate_normal_1
    !
    ! Simulates one draw from a normal distribution using
    !     Box-Muller tranformation.
    !
    ! REFERENCE: Box, G.E.P., Muller, M.E. (1958). A Note on the Generation of
    !            Random Normal Deviates, The Annals of Mathematical Statistics,
    !            29(2), 610–611.
    !##############################################################################
    subroutine simulate_normal_1(x, mu, sigma, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: uni1, uni2, mu_c, sigma_c
        real(RK) :: pi = 3.141592653589793d0


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sigma

        if(sigma_c <= 0d0)then
            call error('simulate_normal','sigma has zero or negative value')
        endif

        ! simulate a uniform variable draw
        if(present(fixed))then
            call simulate_uniform_1(uni1, fixed=fixed)
            call simulate_uniform_1(uni2, fixed=fixed)
        else
            call simulate_uniform_1(uni1)
            call simulate_uniform_1(uni2)
        endif

        ! transform by Box-Muller transformation
        x = sqrt(-2d0*log(uni1))*cos(2d0*pi*uni2)

        ! transform to mean and variance
        x = mu_c + sqrt(sigma_c)*x

    end subroutine simulate_normal_1


    !##############################################################################
    ! SUBROUTINE simulate_normal_n
    !
    ! Simulates draws from a normal distribution using
    !     Box-Muller tranformation.
    !
    ! REFERENCE: Box, G.E.P., Muller, M.E. (1958). A Note on the Generation of
    !            Random Normal Deviates, The Annals of Mathematical Statistics,
    !            29(2), 610–611.
    !##############################################################################
    subroutine simulate_normal_n(x, mu, sigma, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x(:)

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: uni1(size(x, 1)/2), uni2(size(x, 1)/2), mu_c, sigma_c
        integer(IK) :: n, in, n2
        real(RK) :: pi = 3.141592653589793d0


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sigma

        if(sigma_c <= 0d0)then
            call error('simulate_normal','sigma has zero or negative value')
        endif

        ! get size of x
        n = size(x, 1)
        n2 = n/2

        ! simulate a uniform variable draw
        if(present(fixed))then
            call simulate_uniform_n(uni1, fixed=fixed)
            call simulate_uniform_n(uni2, fixed=fixed)
        else
            call simulate_uniform_n(uni1)
            call simulate_uniform_n(uni2)
        endif

        ! transform by Box-Muller transformation
        do in = 1, n2
            x(2*(in-1)+1) = sqrt(-2d0*log(uni1(in)))*cos(2d0*pi*uni2(in))
            x(2*(in-1)+2) = sqrt(-2d0*log(uni1(in)))*sin(2d0*pi*uni2(in))
        enddo

        if(2*n2 /= n)then
            call simulate_normal_1(x(n), 0d0, 1d0)
        endif

        ! transform to mean and variance
        x = mu_c + sqrt(sigma_c)*x

    end subroutine simulate_normal_n


    !##############################################################################
    ! SUBROUTINE simulate_log_normal_1
    !
    ! Simulates one draw from a log normal distribution.
    !##############################################################################
    subroutine simulate_log_normal_1(x, mu, sigma, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = exp(0.5d0)
        if(present(mu))mu_c = mu
        sigma_c = exp(1d0)*(exp(1d0)-1d0)
        if(present(sigma))sigma_c = sigma

        if(sigma_c <= 0d0)then
            call error('simulate_log_normal','sigma has zero or negative value')
        endif
        if(mu_c <= 0d0)then
            call error('simulate_log_normal','mu has zero or negative value')
        endif

        ! get expectation and variance
        sigma_c = log(1d0+sigma_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sigma_c

        ! simulate normal and convert to log_normal
        if(present(fixed))then
            call simulate_normal_1(x, mu_c, sigma_c, fixed)
        else
            call simulate_normal_1(x, mu_c, sigma_c)
        endif

        ! transform to log normal
        x = exp(x)

    end subroutine simulate_log_normal_1


    !##############################################################################
    ! SUBROUTINE simulate_log_normal_n
    !
    ! Simulates draws from a log normal distribution.
    !##############################################################################
    subroutine simulate_log_normal_n(x, mu, sigma, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x(:)

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = exp(0.5d0)
        if(present(mu))mu_c = mu
        sigma_c = exp(1d0)*(exp(1d0)-1d0)
        if(present(sigma))sigma_c = sigma

        if(sigma_c <= 0d0)then
            call error('simulate_log_normal','sigma has zero or negative value')
        endif
        if(mu_c <= 0d0)then
            call error('simulate_log_normal','mu has zero or negative value')
        endif

        ! get expectation and variance
        sigma_c = log(1d0+sigma_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sigma_c

        ! simulate normal and convert to log_normal
        if(present(fixed))then
            call simulate_normal_n(x, mu_c, sigma_c, fixed)
        else
            call simulate_normal_n(x, mu_c, sigma_c)
        endif

        ! transform to log normal
        x = exp(x)

    end subroutine simulate_log_normal_n


    !##############################################################################
    ! SUBROUTINE simulate_Gamma_1
    !
    ! Simulates draws from a gamma distribution using the
    !     Marsaglia/Tsang method.
    !
    ! REFERENCE: Marsaglia, G., Tsang, W.W. (2000). A simple method for generating
    !            gamma variables, ACM Transactions on Mathematical Software, 26(2),
    !            363-372.
    !##############################################################################
    subroutine simulate_Gamma_1(x, alpha, beta, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x

        ! shape parameter
        real(RK), optional :: alpha

        ! rate parameter
        real(RK), optional :: beta

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: alpha_c, beta_c, d, c, uni


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        alpha_c = 1d0
        if(present(alpha))alpha_c = alpha
        beta_c = 1d0
        if(present(beta))beta_c = beta

        ! check for validity of parameters
        if(alpha_c <= 0d0)then
            call error('simulate_Gamma','alpha has a non-positive value')
        endif
        if(beta_c <= 0d0)then
            call error('simulate_Gamma','beta has a non-positive value')
        endif

        ! check whether alpha >= 1
        if(alpha_c >= 1d0)then

            ! calculate c and d
            d = alpha_c - 1d0/3d0
            c = 1d0/sqrt(9d0*d)

            ! simulate the gammas step by step
            if(present(fixed))then
                x = simulate_Gamma_plain(c, d, fixed)
            else
                x = simulate_Gamma_plain(c, d)
            endif
        else

            ! generate alpha+1 variables
            d = alpha_c + 2d0/3d0
            c = 1d0/sqrt(9d0*d)

            ! simulate the gammas step by step
            if(present(fixed))then
                x = simulate_Gamma_plain(c, d, fixed)
            else
                x = simulate_Gamma_plain(c, d)
            endif

            ! add uniform transformation
            if(present(fixed))then
                call simulate_uniform_1(uni, fixed=fixed)
            else
                call simulate_uniform_1(uni)
            endif
            x = x*uni**(1d0/alpha_c)

        endif

        ! apply scaling parameter
        x = x/beta_c


    end subroutine simulate_Gamma_1


    !##############################################################################
    ! SUBROUTINE simulate_Gamma_n
    !
    ! Simulates draws from a gamma distribution using the
    !     Marsaglia/Tsang method.
    !
    ! REFERENCE: Marsaglia, G., Tsang, W.W. (2000). A simple method for generating
    !            gamma variables, ACM Transactions on Mathematical Software, 26(2),
    !            363-372.
    !##############################################################################
    subroutine simulate_Gamma_n(x, alpha, beta, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x(:)

        ! shape parameter
        real(RK), optional :: alpha

        ! rate parameter
        real(RK), optional :: beta

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: alpha_c, beta_c, d, c, uni(size(x, 1))
        integer(IK) :: n, in


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        alpha_c = 1d0
        if(present(alpha))alpha_c = alpha
        beta_c = 1d0
        if(present(beta))beta_c = beta

        ! check for validity of beta
        if(beta_c <= 0d0)then
            call error('simulate_Gamma','beta has a non-positive value')
        endif

        ! get size of x
        n = size(x, 1)

        ! check whether alpha >= 1
        if(alpha_c >= 1d0)then

            ! calculate c and d
            d = alpha_c - 1d0/3d0
            c = 1d0/sqrt(9d0*d)

            ! simulate the gammas step by step
            do in = 1, n
                if(present(fixed))then
                    x(in) = simulate_Gamma_plain(c, d, fixed)
                else
                    x(in) = simulate_Gamma_plain(c, d)
                endif
            enddo
        else

            ! generate alpha+1 variables
            d = alpha_c + 2d0/3d0
            c = 1d0/sqrt(9d0*d)

            ! simulate the gammas step by step
            do in = 1, n
                if(present(fixed))then
                    x(in) = simulate_Gamma_plain(c, d, fixed)
                else
                    x(in) = simulate_Gamma_plain(c, d)
                endif
            enddo

            ! add uniform transformation
            if(present(fixed))then
                call simulate_uniform_n(uni, fixed=fixed)
            else
                call simulate_uniform_n(uni)
            endif
            x = x*uni**(1d0/alpha_c)

        endif

        ! apply scaling parameter
        x = x/beta_c


    end subroutine simulate_Gamma_n


    !##############################################################################
    ! FUNCTION simulate_Gamma_plain
    !
    ! Simulates one draw from gamma distribution without any testing.
    !
    ! REFERENCE: Marsaglia, G., Tsang, W.W., A simple method for generating gamma
    !            variables, ACM Transactions on Mathematical Software, Vol. 26,
    !            No. 2, 2000, 363-372.
    !##############################################################################
    function simulate_Gamma_plain(c, d, fixed) result(res)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! c parameter
        real(RK) :: c

        ! d parameter
        real(RK) :: d

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed

        ! a gamma realization
        real(RK) :: res


        !##### OTHER VARIABLES ####################################################

        real(RK) :: v, x, u


        !##### ROUTINE CODE #######################################################

        do
            ! generate v
            if(present(fixed))then
                call simulate_normal_1(x, fixed=fixed)
            else
                call simulate_normal_1(x)
            endif
            v = (1d0+c*x)
            v = v*v*v

            ! test v
            if(v <= 0d0)cycle

            ! generate u
            if(present(fixed))then
                call simulate_uniform_1(u, fixed=fixed)
            else
                call simulate_uniform_1(u)
            endif

            ! check squeeze
            if(u < 1d0 - 0.0331d0*(x*x)*(x*x))then
                res = d*v
                return
            endif

            ! check real condition
            if(log(u) < 0.5d0*x*x + d*(1d0 - v + log(v)))then
                res = d*v
                return
            endif
        enddo

    end function simulate_Gamma_plain


    !##############################################################################
    ! SUBROUTINE simulate_beta_1
    !
    ! Simulates one draw from a beta distribution.
    !##############################################################################
    subroutine simulate_beta_1(x, p, q, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x

        ! shape parameter
        real(RK), optional :: p

        ! rate parameter
        real(RK), optional :: q

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: p_c, q_c, gamma1, gamma2, bern


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        p_c = 1d0
        if(present(p))p_c = p
        q_c = 1d0
        if(present(q))q_c = q

        ! check for validity of parameters
        if(p_c <= 0d0)then
            call error('simulate_beta','p has a non-positive value')
        endif
        if(q_c <= 0d0)then
            call error('simulate_beta','q has a non-positive value')
        endif

        ! get bernoulli parametzer for small p and q
        bern = p_c/(p_c+q_c)

        if(present(fixed))then
            call simulate_Gamma_1(gamma1, p_c, fixed=fixed)
            call simulate_Gamma_1(gamma2, q_c, fixed=fixed)
        else
            call simulate_Gamma_1(gamma1, p_c)
            call simulate_Gamma_1(gamma2, q_c)
        endif

        ! check whether gammas both greaters 0
        if(gamma1 > 0d0 .or. gamma2 > 0d0)then
            x = gamma1/(gamma1+gamma2)
        else
            if(present(fixed))then
                call simulate_bernoulli_1(x, bern, fixed=fixed)
            else
                call simulate_bernoulli_1(x, bern)
            endif
        endif

    end subroutine simulate_beta_1


    !##############################################################################
    ! SUBROUTINE simulate_beta_n
    !
    ! Simulates draws from a beta distribution.
    !##############################################################################
    subroutine simulate_beta_n(x, p, q, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x(:)

        ! shape parameter
        real(RK), optional :: p

        ! rate parameter
        real(RK), optional :: q

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        real(RK) :: p_c, q_c, gamma1(size(x,1)), gamma2(size(x,1)), bern
        integer(IK) :: n, in


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        p_c = 1d0
        if(present(p))p_c = p
        q_c = 1d0
        if(present(q))q_c = q

        ! check for validity of parameters
        if(p_c <= 0d0)then
            call error('simulate_beta','p has a non-positive value')
        endif
        if(q_c <= 0d0)then
            call error('simulate_beta','q has a non-positive value')
        endif

        ! get size of x
        n = size(x, 1)

        ! get bernoulli parametzer for small p and q
        bern = p_c/(p_c+q_c)

        if(present(fixed))then
            call simulate_Gamma_n(gamma1, p_c, fixed=fixed)
            call simulate_Gamma_n(gamma2, q_c, fixed=fixed)
        else
            call simulate_Gamma_n(gamma1, p_c)
            call simulate_Gamma_n(gamma2, q_c)
        endif

        do in = 1, n

            ! check whether gammas both greaters 0
            if(gamma1(in) > 0d0 .or. gamma2(in) > 0d0)then
                x(in) = gamma1(in)/(gamma1(in)+gamma2(in))
            else
                if(present(fixed))then
                    call simulate_bernoulli_1(x(in), bern, fixed=fixed)
                else
                    call simulate_bernoulli_1(x(in), bern)
                endif
            endif
        enddo

    end subroutine simulate_beta_n


    !##############################################################################
    ! SUBROUTINE simulate_bernoulli_1
    !
    ! Simulates one draw from bernoulli distribution
    !##############################################################################
    subroutine simulate_bernoulli_1(x, p, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x

        ! probability to draw 1
        real(RK), intent(in) :: p

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### ROUTINE CODE #######################################################

        if(p < 0d0)then
            call error('simulate_bernoulli', 'p has a negative value')
        endif

        if(present(fixed))then
            call simulate_uniform_1(x, fixed=fixed)
        else
            call simulate_uniform_1(x)
        endif

        if(x <= p)then
            x = 1d0
        else
            x = 0d0
        endif

    end subroutine simulate_bernoulli_1


    !##############################################################################
    ! SUBROUTINE simulate_bernoulli_n
    !
    ! Simulates draws from bernoulli distribution
    !##############################################################################
    subroutine simulate_bernoulli_n(x, p, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x(:)

        ! probability to draw 1
        real(RK), intent(in) :: p

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: ii


        !##### ROUTINE CODE #######################################################

        if(p < 0d0)then
            call error('simulate_bernoulli', 'p has a negative value')
        endif

        if(present(fixed))then
            call simulate_uniform_n(x, fixed=fixed)
        else
            call simulate_uniform_n(x)
        endif

        do ii = 1, size(x, 1)
            if(x(ii) <= p)then
                x(ii) = 1d0
            else
                x(ii) = 0d0
            endif
        enddo

    end subroutine simulate_bernoulli_n


    !##############################################################################
    ! SUBROUTINE simulate_binomial_1
    !
    ! Simulates one draw from binomial distribution
    !##############################################################################
    subroutine simulate_binomial_1(x, n, p, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x

        ! total number of draws
        integer(IK), intent(in) :: n

        ! probability to draw 1
        real(RK), intent(in) :: p

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: ii
        real(RK) :: cumprob


        !##### ROUTINE CODE #######################################################

        if(p < 0d0)then
            call error('simulate_binomial', 'p has a negative value')
        endif

        if(n < 0)then
            call error('simulate_binomial', 'n must be positive')
        endif

        if(present(fixed))then
            call simulate_uniform_1(x, fixed=fixed)
        else
            call simulate_uniform_1(x)
        endif

        ! derive bernoulli probability schedule
        cumprob = binomialPDF(0, n, p)
        if(x <= cumprob)then
            x = 0d0
        else
            do ii = 1, n
                cumprob = cumprob + binomialPDF(ii, n, p)
                if(x <= cumprob)then
                    x = dble(ii)
                    exit
                endif
            enddo
        endif

    end subroutine simulate_binomial_1


    !##############################################################################
    ! SUBROUTINE simulate_binomial_n
    !
    ! Simulates draws from binomial distribution
    !##############################################################################
    subroutine simulate_binomial_n(x, n, p, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real(RK), intent(out) :: x(:)

        ! total number of draws
        integer(IK), intent(in) :: n

        ! probability to draw 1
        real(RK), intent(in) :: p

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: ii, it
        real(RK) :: cumprob(0:n)


        !##### ROUTINE CODE #######################################################

        if(p < 0d0)then
            call error('simulate_binomial', 'p has a negative value')
        endif

        if(n < 0)then
            call error('simulate_binomial', 'n must be positive')
        endif

        ! simulate
        if(present(fixed))then
            call simulate_uniform_n(x, fixed=fixed)
        else
            call simulate_uniform_n(x)
        endif

        ! derive bernoulli probability schedule
        cumprob(0) = binomialPDF(0, n, p)
        do ii = 1, n
            cumprob(ii) = cumprob(ii-1) + binomialPDF(ii, n, p)
        enddo

        ! get draws
        do it = 1, size(x, 1)
            do ii = 0, n
                if(x(it) <= cumprob(ii))then
                    x(it) = dble(ii)
                    exit
                endif
            enddo
        enddo

    end subroutine simulate_binomial_n


    !##############################################################################
    ! SUBROUTINE normal_discrete_1
    !
    ! Creates n points and probabilities for a normal distribution.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     CompEcon toolbox of Mario Miranda and Paul Fackler available under
    !     http://www4.ncsu.edu/~pfackler/compecon/toolbox.html
    !
    ! REFERENCE: Miranda, M. & Fackler, P. (2002). Applied Computational Economics
    !            and Finance. Cambridge: MIT Press.
    !##############################################################################
    subroutine normal_discrete_1(x, prob, mu, sigma)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! discrete points of normal distribution
        real(RK), intent(out) :: x(:)

        ! probability weights
        real(RK), intent(out) :: prob(:)

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c, pim4, z=0d0, z1, p1, p2, p3, pp
        integer(IK) :: n, m, i, j, its
        integer(IK), parameter :: maxit = 200
        real(RK), parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sqrt(sigma)

        if(sigma_c < 0d0)then
            call error('normal_discrete','sigma has negative value')
        endif

        ! check for right array sizes
        n = assert_eq(size(x,1), size(prob,1), 'normal_discrete')

        ! calculate 1/pi^0.25
        pim4 = 1d0/pi**0.25d0

        ! get number of points
        m = (n+1)/2

        ! initialize x and prob
        x = 0d0
        prob = 0d0

        ! start iteration
        do i = 1, m

            ! set reasonable starting values
            if(i == 1)then
                z = sqrt(dble(2*n+1))-1.85575d0*(dble(2*n+1)**(-1d0/6d0))
            elseif(i == 2)then
                z = z - 1.14d0*(dble(n)**0.426d0)/z
            elseif(i == 3)then
                z = 1.86d0*z+0.86d0*x(1)
            elseif(i == 4)then
                z = 1.91d0*z+0.91d0*x(2);
            else
                z = 2d0*z+x(i-2);
            endif

            ! root finding iterations
            its = 0
            do while(its < maxit)
                its = its+1
                p1 = pim4
                p2 = 0d0
                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = z*sqrt(2d0/dble(j))*p2-sqrt(dble(j-1)/dble(j))*p3
                enddo
                pp = sqrt(2d0*dble(n))*p2
                z1 = z
                z  = z1-p1/pp
                if(abs(z-z1) < 1e-14)exit
            enddo
            if(its >= maxit)then
                call error('normal_discrete', &
                    'Could not discretize normal distribution')
            endif
            x(n+1-i) = z
            x(i) = -z
            prob(i) = 2d0/pp**2
            prob(n+1-i) = prob(i)
        enddo

        ! set output data
        prob = prob/sqrt(pi)
        x = x*sqrt(2d0)*sigma_c + mu_c

    end subroutine normal_discrete_1


    !##############################################################################
    ! SUBROUTINE normal_discrete_2
    !
    ! Creates n1*n2 points and probabilities for a two-dimensional normal
    !     distribution.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     CompEcon toolbox of Mario Miranda and Paul Fackler available under
    !     http://www4.ncsu.edu/~pfackler/compecon/toolbox.html
    !
    ! REFERENCE: Miranda, M. & Fackler, P. (2002). Applied Computational Economics
    !            and Finance. Cambridge: MIT Press.
    !##############################################################################
    subroutine normal_discrete_2(n, x, prob, mu, sigma, rho)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! number of points in every direction
        integer(IK), intent(in) :: n(2)

        ! discrete points of normal distribution
        real(RK), intent(out) :: x(:, :)

        ! probability weights
        real(RK), intent(out) :: prob(:)

        ! expectation of distribution
        real(RK), optional :: mu(2)

        ! variance of distribution
        real(RK), optional :: sigma(2)

        ! correlation of distribution
        real(RK), optional :: rho


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c(2), sig_c(2), rho_c, sigma_c(2,2), l(2,2)
        real(RK) :: x1(n(1)), x2(n(2)), p1(n(1)), p2(n(2))
        integer(IK) :: m, j, k


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sig_c(1) = 1d0
        if(present(sigma))sig_c = sigma
        rho_c = 0d0
        if(present(rho))rho_c = rho

        ! set up variance covariance matrix
        sigma_c(1, 1) = sig_c(1)
        sigma_c(2, 2) = sig_c(2)
        sigma_c(1, 2) = rho_c*sqrt(sig_c(1)*sig_c(2))
        sigma_c(2, 1) = sigma_c(1, 2)

        ! check for right array sizes
        m = assert_eq(size(x,1), size(prob,1), n(1)*n(2), 'normal_discrete')
        m = assert_eq(size(x,2), 2, 'normal_discrete')

        ! check whether sigma is symmetric
        if(any(abs(transpose(sigma_c) - sigma_c) > 1d-20)) &
            call error('normal_discrete', &
            'Variance-Covariance matrix is not symmetric')

        ! get standard normal distributed random variables
        call normal_discrete(x1, p1, 0d0, 1d0)
        call normal_discrete(x2, p2, 0d0, 1d0)

        ! get joint distribution
        m = 1
        do k = 1, n(2)
            do j = 1, n(1)
                prob(m) = p1(j)*p2(k)
                x(m, :) = (/x1(j), x2(k)/)
                m = m+1
            enddo
        enddo

        ! decompose var-cov matrix
        if(.not.any(abs(sig_c)  <= 1d-100))then
            call cholesky(sigma_c, l)
        else
            l = 0d0
            l(1,1) = sqrt(sig_c(1))
            l(2,2) = sqrt(sig_c(2))
        endif

        ! calculate distribution
        x = matmul(x, transpose(l))
        x(:, 1) = x(:, 1) + mu_c(1)
        x(:, 2) = x(:, 2) + mu_c(2)

    end subroutine normal_discrete_2


    !##############################################################################
    ! SUBROUTINE log_normal_discrete_1
    !
    ! Creates n points and probabilities for a log-normal distribution.
    !
    ! REFERENCE: Miranda, M. & Fackler, P. (2002). Applied Computational Economics
    !            and Finance. Cambridge: MIT Press.
    !##############################################################################
    subroutine log_normal_discrete_1(x, prob, mu, sigma)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! discrete points of normal distribution
        real(RK), intent(out) :: x(:)

        ! probability weights
        real(RK), intent(out) :: prob(:)

        ! expectation of distribution
        real(RK), optional :: mu

        ! variance of distribution
        real(RK), optional :: sigma


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c, sigma_c
        integer(IK) :: n

        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = exp(0.5d0)
        if(present(mu))mu_c = mu
        sigma_c = exp(1d0)*(exp(1d0)-1d0)
        if(present(sigma))sigma_c = sigma

        if(sigma_c < 0d0)then
            call error('log_normal_discrete','sigma has negative value')
        endif
        if(mu_c <= 0d0)then
            call error('log_normal_discrete','mu has zero or negative value')
        endif

        ! get expectation and variance
        sigma_c = log(1d0+sigma_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sigma_c

        ! check for right array sizes
        n = assert_eq(size(x,1), size(prob,1), 'normal_discrete')
        n = n

        call normal_discrete(x, prob, mu_c, sigma_c)

        x = exp(x)

    end subroutine log_normal_discrete_1


    !##############################################################################
    ! SUBROUTINE log_normal_discrete_2
    !
    ! Creates n1*n2 points and probabilities for a two-dimensional log-normal
    !     distribution.
    !
    ! REFERENCE: Miranda, M. & Fackler, P. (2002). Applied Computational Economics
    !            and Finance. Cambridge: MIT Press.
    !##############################################################################
    subroutine log_normal_discrete_2(n, x, prob, mu, sigma, rho)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! number of points in every direction
        integer(IK), intent(in) :: n(2)

        ! discrete points of normal distribution
        real(RK), intent(out) :: x(:, :)

        ! probability weights
        real(RK), intent(out) :: prob(:)

        ! expectation of distribution
        real(RK), optional :: mu(2)

        ! variance of distribution
        real(RK), optional :: sigma(2)

        ! correlation of distribution
        real(RK), optional :: rho


        !##### OTHER VARIABLES ####################################################

        real(RK) :: mu_c(2), sig_c(2), rho_c, sigma_c(2,2), l(2,2)
        real(RK) :: x1(n(1)), x2(n(2)), p1(n(1)), p2(n(2))
        integer(IK) :: m, j, k


        !##### ROUTINE CODE #######################################################

        ! initialize parameters
        mu_c = exp(0.5d0)
        if(present(mu))mu_c = mu
        sig_c = exp(1d0)*(exp(1d0)-1d0)
        if(present(sigma))sig_c = sigma
        rho_c = 0d0
        if(present(rho))rho_c = rho

        if(any(sig_c < 0d0))then
            call error('log_normal_discrete','sigma has negative value')
        endif
        if(any(mu_c <= 0d0))then
            call error('log_normal_discrete','mu has zero or negative value')
        endif
        if(rho_c < -1d0 .or. rho_c > 1d0)then
            call error('log_normal_discrete','rho is outside -1 to 1')
        endif

        ! get expectation and variance
        sig_c = log(1d0+sig_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sig_c

        ! set up covariance matrix
        sigma_c(1, 1) = sig_c(1)
        sigma_c(2, 2) = sig_c(2)
        sigma_c(1, 2) = log(rho_c*sqrt(exp(sig_c(1))-1d0)* &
            sqrt(exp(sig_c(2))-1d0)+1d0)
        sigma_c(2, 1) = sigma_c(1, 2)

        ! check for right array sizes
        m = assert_eq(size(x,1), size(prob,1), n(1)*n(2), 'normal_discrete')
        m = assert_eq(size(x,2), 2, 'normal_discrete')

        ! check whether sigma is symmetric
        if(any(abs(transpose(sigma_c) - sigma_c) > 1d-20)) &
            call error('normal_discrete', &
            'Variance-Covariance matrix is not symmetric')

        ! get standard normal distributed random variables
        call normal_discrete(x1, p1, 0d0, 1d0)
        call normal_discrete(x2, p2, 0d0, 1d0)

        ! get joint distribution
        m = 1
        do k = 1, n(2)
            do j = 1, n(1)
                prob(m) = p1(j)*p2(k)
                x(m, :) = (/x1(j), x2(k)/)
                m = m+1
            enddo
        enddo

        ! decompose var-cov matrix
        if(.not.any(abs(sig_c) <= 1d-100))then
            call cholesky(sigma_c, l)
        else
            l = 0d0
            l(1,1) = sqrt(sig_c(1))
            l(2,2) = sqrt(sig_c(2))
        endif

        ! calculate distribution
        x = matmul(x, transpose(l))
        x(:, 1) = x(:, 1) + mu_c(1)
        x(:, 2) = x(:, 2) + mu_c(2)
        x = exp(x)

    end subroutine log_normal_discrete_2














!##############################################################################
!##############################################################################
! MODULE polynomial
!##############################################################################
!##############################################################################


    !##############################################################################
    ! FUNCTION poly_interpol_1
    !
    ! Constructs interpolating polynomial given nodes xi and data yi.
    !##############################################################################
    function poly_interpol_1(x, xi, yi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate polynomial
        real(RK), intent(in) :: x

        ! nodes of interpolation
        real(RK), intent(in) :: xi(0:)

        ! data of interpolation
        real(RK), intent(in) :: yi(0:)

        ! return value
        real(RK) :: poly_interpol_1


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: j, n, i
        real(RK) :: lagrange_poly(0:size(xi, 1)-1)

        !##### ROUTINE CODE #######################################################

        ! get number of interpolation nodes
        n = assert_eq(size(xi, 1), size(yi, 1), 'poly_interpol') - 1

        ! initialize lagrange basis polynomials
        lagrange_poly(:) = 1d0

        ! span polynomials
        do j = 0, n
            do i = 0, n
                if(j /= i)then
                    lagrange_poly(j) = lagrange_poly(j) * (x-xi(i))/(xi(j)-xi(i))
                endif
            enddo
        enddo

        poly_interpol_1 = sum(yi*lagrange_poly, 1)

    end function poly_interpol_1


    !##############################################################################
    ! FUNCTION poly_interpol_m
    !
    ! Constructs interpolating polynomial given nodes xi and data yi and several
    !     points x.
    !##############################################################################
    function poly_interpol_m(x, xi, yi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate polynomial
        real(RK), intent(in) :: x(1:)

        ! nodes of interpolation
        real(RK), intent(in) :: xi(0:)

        ! data of interpolation
        real(RK), intent(in) :: yi(0:)

        ! return value
        real(RK) :: poly_interpol_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: j, n, i
        real(RK) :: lagrange_poly(size(x, 1), 0:size(xi, 1))


        !##### ROUTINE CODE #######################################################

        ! get number of interpolation nodes
        n = assert_eq(size(xi, 1), size(yi, 1), 'poly_interpol') - 1

        ! initialize lagrange basis polynomials
        lagrange_poly(:, :) = 1d0

        ! span polynomials
        do j = 0, n
            do i = 0, n
                if(j /= i)then
                    lagrange_poly(:, j) = lagrange_poly(:, j) * &
                        (x-xi(i))/(xi(j)-xi(i))
                endif
            enddo
        enddo

        do j = 1, size(x, 1)
            poly_interpol_m(j) = sum(yi*lagrange_poly(j, :), 1)
        enddo

    end function poly_interpol_m















!##############################################################################
!##############################################################################
! MODULE minimization
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE settol_min
    !
    ! For setting global tolerance level.
    !##############################################################################
    subroutine settol_min(tol)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! tolerance level
        real(RK), intent(in) :: tol


        !##### ROUTINE CODE #######################################################

        ! check whether tolerance level is valid
        if(tol > 0d0)then
            tbox_gftol = tol
        else
            call warning('settol_min', 'tolerance level is not valid')
        endif

    end subroutine settol_min


    !##############################################################################
    ! SUBROUTINE setiter_min
    !
    ! For setting maximum number of iterations.
    !##############################################################################
    subroutine setiter_min(iter)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! tolerance level
        integer(IK), intent(in) :: iter


        !##### ROUTINE CODE #######################################################

        ! check whether tolerance level is valid
        if(iter > 0)then
            tbox_itermax_min = iter
        else
            call warning('setiter_min', 'number of iterations is not valid')
        endif

    end subroutine setiter_min


    !##############################################################################
    ! SUBROUTINE brent
    !
    ! Minimizes a one dimensional function.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    subroutine brent(xmin, fret, minimum, maximum, func)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! minimum value found
        real(RK), intent(inout) :: xmin

        ! function value at minimum
        real(RK), intent(out) :: fret

        ! left, middle and right interval points
        real(RK), intent(in) :: minimum, maximum


        !##### OTHER VARIABLES ####################################################

        real(RK) :: tol
        real(RK), parameter :: cgold = 0.3819660d0
        real(RK), parameter :: zeps = 1.0e-3*epsilon(xmin)
        real(RK) :: a=0d0, b=0d0, d=0d0, e=0d0, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm, ax, bx, cx
        integer(IK) :: iter


        !##### INTERFACES #########################################################

        ! interface for the function
        interface
            function func(p)
                use Constants_mod, only: IK, RK; implicit none
                real(RK), intent(in) :: p
                real(RK) :: func
            end function func
        end interface


        !##### ROUTINE CODE #######################################################

        ! set tolerance level
        tol =  tbox_gftol

        ! set ax, bx and cx
        ax = minimum
        cx = maximum

        a = min(ax, cx)
        b = max(ax, cx)

        if(abs(xmin-a) <= 1d-6)then
            bx = a + 1d-6
        elseif(abs(xmin-b) <= 1d-6)then
            bx = b - 1d-6
        elseif(xmin > a .and. xmin < b)then
            bx = xmin
        else
            bx = (ax+cx)/2d0
        endif

        v = bx
        w = v
        x = v
        e = 0d0
        fx = func(x)
        fv = fx
        fw = fx

        do iter = 1,tbox_itermax_min
            xm = 0.5d0*(a+b)
            tol1 = tol*abs(x)+zeps
            tol2 = 2.0d0*tol1

            if(abs(x-xm) <= (tol2-0.5d0*(b-a)))then
                xmin = x
                fret = fx
                return
            endif

            if(abs(e) > tol1)then
                r = (x-w)*(fx-fv)
                q = (x-v)*(fx-fw)
                p = (x-v)*q-(x-w)*r
                q = 2.0d0*(q-r)
                if (q > 0.0d0) p = -p
                q = abs(q)
                etemp = e
                e = d
                if(abs(p) >= abs(0.5d0*q*etemp) .or. &
                        p <= q*(a-x) .or. p >= q*(b-x))then
                    e = merge(a-x, b-x, x >= xm )
                    d = CGOLD*e
                else
                    d = p/q
                    u = x+d
                    if(u-a < tol2 .or. b-u < tol2) d = sign(tol1, xm-x)
                endif

            else
                e = merge(a-x, b-x, x >= xm )
                d = CGOLD*e
            endif

            u = merge(x+d, x+sign(tol1, d), abs(d) >= tol1 )
            fu = func(u)
            if(fu <= fx)then
                if(u >= x)then
                    a = x
                else
                b = x
                endif
                call shft(v, w, x, u)
                call shft(fv, fw, fx, fu)
            else
                if(u < x)then
                    a = u
                else
                    b = u
                endif
                if(fu <= fw .or. abs(w-x)  <= 1d-100)then
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                elseif(fu <= fv .or. abs(v-x) <= 1d-100 .or. abs(v-w) <= 1d-100)then
                    v = u
                    fv = fu
                endif
            endif
        enddo

        call warning('fminsearch', 'maximum iterations exceeded')


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! SUBROUTINE shft
        !
        ! Shifts b to a, c to b and d to c.
        !##########################################################################
        subroutine shft(a, b, c, d)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            real(RK), intent(out)   :: a
            real(RK), intent(inout) :: b, c
            real(RK), intent(in   ) :: d


            !##### ROUTINE CODE ###################################################
            a = b
            b = c
            c = d
        end subroutine shft

    end subroutine brent



    !##############################################################################
    ! SUBROUTINE powell
    !
    ! Powell is a multidimensional function minimizer.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    subroutine powell(p, fret, minimum, maximum, func)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! starting and ending point
        real(RK), intent(inout) :: p(:)

        ! value of function in minimum
        real(RK), intent(out) :: fret

        ! minimum optimization interval point
        real(RK), intent(in) :: minimum(:)

        ! maximum optimization interval point
        real(RK), intent(in) :: maximum(:)


        !##### OTHER VARIABLES ####################################################

        real(RK) :: xi(size(p, 1), size(p, 1))
        real(RK) :: ftol
        real(RK), parameter :: tiny = 1.0e-25
        integer(IK) :: i, ibig, n, iter
        real(RK) :: del, fp, fptt, t
        real(RK), dimension(size(p)) :: pt, ptt, xit
        real(RK) :: xicom(size(p,1))


        !##### INTERFACES #########################################################

        ! interface for the function
        interface
            function func(p)
                use Constants_mod, only: IK, RK; implicit none
                real(RK), intent(in) :: p(:)
                real(RK) :: func
            end function func
        end interface


        !##### ROUTINE CODE #######################################################

        ! set tolerance level
        ftol = tbox_gftol

        ! Set number of points
        n = assert_eq(size(p), size(minimum, 1), size(maximum,1), 'fminsearch')

        ! initialize direction set
        xi = 0d0
        do i = 1, n
            xi(i, i) = 1d0
        enddo

        ! calculate function value
        fret = func(p)

        ! store old p
        pt(:) = p(:)

        ! start iteration
        iter = 0
        do
            ! step counter
            iter = iter+1

            ! save old function value
            fp = fret

            ! ibig will be direction of steepest decline
            ibig = 0
            del = 0.0d0

            ! iterate over all dimensions
            do i = 1, n

                ! copy direction i and store old function value
                xit(:) = xi(:,i)
                fptt = fret

                ! minimize along this direction
                call linmin(p, xit, n, fret, func)

                ! store i into i big if i is the direction of steepest decline
                if (fptt-fret > del) then
                    del=fptt-fret
                    ibig=i
                endif
            enddo

            ! termination criterion
            if (2d0*(fp - fret) <= ftol*(abs(fp) + abs(fret)) + tiny) return

            ! quit if maximum iterations reached
            if(iter == tbox_itermax_min)then
                call warning('fminsearch', 'maximum iterations exceeded')
            endif

            ! construct extrapolated point
            ptt(:) = 2d0*p(:) - pt(:)
            xit(:) = p(:) - pt(:)
            pt(:) = p(:)

            ! calculate function value at extrapolated point
            fptt = func(ptt)

            ! if function value greater than actual value -> no change of directions
            if (fptt >= fp) cycle
            t = 2d0*(fp - 2d0*fret + fptt) * (fp - fret - del)**2 - del*(fp - fptt)**2

            ! if t > 0 -> no change of directions
            if (t >= 0d0) cycle

            ! else minimize along new direction and start new iteration
            call linmin(p, xit, n, fret, func)
            xi(:, ibig) = xi(:, n)
            xi(:,n) = xit(:)
        enddo


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! SUBROUTINE linmin
        !
        ! Minimizes multidimensional function along a given direction.
        !##########################################################################
        subroutine linmin(p, xi, n, fret, func)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! number of dimensions
            integer(IK), intent(in) :: n

            ! point where to start and minimum if minimization is done
            real(RK), intent(inout) :: p(n)

            ! direction in which to optimize
            real(RK), intent(inout) :: xi(n)

            ! value of function at minimum
            real(RK), intent(out)  :: fret


            !##### OTHER VARIABLES ################################################

            real(RK) :: tol
            real(RK)  :: ax, bx, cx, xmin


            !##### INTERFACES #####################################################

                ! interface for the function
            interface
                function func(p)
                    use Constants_mod, only: IK, RK; implicit none
                    real(RK), intent(in) :: p(:)
                    real(RK) :: func
                end function func
            end interface


            !##### ROUTINE CODE ###################################################

            ! set tolerance level
            tol = tbox_gftol

            xicom(:) = xi(:)

            ! get optimization interval
            call getintervall(ax, bx, cx)

            ! minimize function using one dimensional optimizer
            fret = brent_pow(ax, bx, cx, tol, xmin, func)

            ! calculate new direction and endpoint
            xi(:) = xmin*xi(:)
            p(:) = p(:) + xi(:)

        end subroutine linmin


        !##########################################################################
        ! SUBROUTINE getinterval
        !
        ! Calculates optimization interval along a given direction.
        !##########################################################################
        subroutine getintervall(ax, bx, cx)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! left, middle and right interval point
            real(RK), intent(out) :: ax, bx, cx


            !##### OTHER VARIABLES ################################################

            integer(IK) :: i
            real(RK)  :: w(0:1,n)


            !##### ROUTINE CODE ###################################################

            ! calculate right interval point
            cx = -1.e20
            do i = 1, n
                if(abs(xicom(i)) >= 1d-100)then
                    w(0,i) = (extr(i, 0, xicom(i))-p(i))/xicom(i)
                else
                    w(0,i) = -1.e20
                endif
                if(w(0,i) > cx)cx = w(0, i)
            enddo

            ! calculate left interval point
            ax = 1.e20
            do i=1, n
                if(abs(xicom(i)) >= 1d-100)then
                    w(1,i) = (extr(i, 1, xicom(i))-p(i))/xicom(i)
                else
                    w(1,i) = 1.e20
                endif
                if(w(1,i) < ax)ax = w(1,i)
            enddo

            ! calculate point in between [ax, cx]
            bx = 0d0

        end subroutine getintervall


        !##########################################################################
        ! FUNCTION extr
        !
        ! Calculates interval endpoint in a given dimension.
        !##########################################################################
        real(RK) function extr(dim, maxxx, richt)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! dimension in which to search
            integer(IK), intent(in) :: dim

            ! do you want the maximum or minimum endpoint
            integer(IK), intent(in) :: maxxx

            ! what is the optimization direction
            real(RK) :: richt


            !##### ROUTINE CODE ###################################################

            if(richt > 0d0 .and. maxxx == 1 .or. richt < 0d0 .and. maxxx == 0)then
                extr = maximum(dim)
            else
                extr = minimum(dim)
            endif

        end function extr


        !##########################################################################
        ! FUNCTION f1dim
        !
        ! Maps multidimensional function and point/direction combo into
        !     one-dimensional function.
        !##########################################################################
        function f1dim(x, func)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! point where to evaluate the multidimensional function
            real(RK), intent(in)  :: x

            ! function value
            real(RK) :: f1dim


            !##### OTHER VARIABLES ################################################

            real(RK) :: xt(n)


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function func(p)
                    use Constants_mod, only: IK, RK; implicit none
                    real(RK), intent(in) :: p(:)
                    real(RK) :: func
                end function func
            end interface


            !##### ROUTINE CODE ###################################################

            ! create point where to evaluate func
            xt(:) = p(:)+x*xicom(:)

            ! evaluate func at this point
            f1dim = func(xt)

        end function f1dim


        !##########################################################################
        ! FUNCTION brent_pow
        !
        ! Minimizes a one dimensional function.
        !##########################################################################
        function brent_pow(ax, bx, cx, tol, xmin, func)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! left, middle and right interval points
            real(RK), intent(in) :: ax, bx, cx

            ! level of tolerance
            real(RK), intent(in) :: tol

            ! minimum value found
            real(RK), intent(out) :: xmin

            ! function value at minimum
            real(RK) :: brent_pow


            !##### OTHER VARIABLES ################################################

            real(RK), parameter :: cgold = 0.3819660d0
            real(RK), parameter :: zeps=1.0e-3*epsilon(ax)
            integer(IK) :: iter
            real(RK) :: a=0d0, b=0d0, d=0d0, e=0d0, etemp=0d0
            real(RK) :: fu, fv, fw, fx, p, q, r, tol1, tol2, &
                u, v, w, x, xm


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function func(p)
                    use Constants_mod, only: IK, RK; implicit none
                    real(RK), intent(in) :: p(:)
                    real(RK) :: func
                end function func
            end interface


            !##### ROUTINE CODE ###################################################

            a = min(ax, cx)
            b = max(ax, cx)
            v = bx
            w = v
            x = v
            e = 0d0
            fx = f1dim(x, func)
            fv = fx
            fw = fx

            do iter = 1,tbox_tbox_itermax_pow_b

                xm = 0.5d0*(a+b)
                tol1 = tol*abs(x)+zeps
                tol2 = 2.0d0*tol1

                if(abs(x-xm) <= (tol2-0.5d0*(b-a)))then
                    xmin = x
                    brent_pow = fx
                    return
                endif

                if(abs(e) > tol1)then
                    r = (x-w)*(fx-fv)
                    q = (x-v)*(fx-fw)
                    p = (x-v)*q-(x-w)*r
                    q = 2.0d0*(q-r)
                    if (q > 0.0d0) p = -p
                    q = abs(q)
                    etemp = e
                    e = d
                    if(abs(p) >= abs(0.5d0*q*etemp ) .or. &
                            p <= q*(a-x) .or. p >= q*(b-x))then
                        e = merge(a-x, b-x, x >= xm )
                        d = CGOLD*e
                    else
                        d = p/q
                        u = x+d
                        if(u-a < tol2 .or. b-u < tol2) d = sign(tol1, xm-x)
                    endif

                else
                    e = merge(a-x, b-x, x >= xm )
                    d = CGOLD*e
                endif
                u = merge(x+d, x+sign(tol1, d), abs(d) >= tol1 )
                fu = f1dim(u, func)
                if(fu <= fx)then
                    if(u >= x)then
                        a = x
                    else
                        b = x
                    endif
                    call shft(v, w, x, u)
                    call shft(fv, fw, fx, fu)
                else
                    if(u < x)then
                        a = u
                    else
                        b = u
                    endif
                    if(fu <= fw .or. abs(w-x) <= 1d-100)then
                        v = w
                        fv = fw
                        w = u
                        fw = fu
                    elseif(fu <= fv .or. abs(v-x) <= 1d-100 .or. abs(v-w) <= 1d-100)then
                        v = u
                        fv = fu
                    endif
               endif
            enddo

            xmin = x
            brent_pow = fx
            call warning('fminsearch', 'maximum iterations exceeded')

        end function brent_pow


        !##########################################################################
        ! SUBROUTINE shft
        !
        ! Shifts b to a, c to b and d to c.
        !##########################################################################
        subroutine shft(a, b, c, d)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            real(RK), intent(out)   :: a
            real(RK), intent(inout) :: b, c
            real(RK), intent(in   ) :: d


            !##### ROUTINE CODE ###################################################
            a = b
            b = c
            c = d
        end subroutine shft

    end subroutine powell















!##############################################################################
!##############################################################################
! MODULE simplex
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE solve_lin
    !
    ! For solving a linear program in normal form by means of the simplex
    !   algorithm.
    !##############################################################################
    subroutine solve_lin(x, c, A, b, numle, numge, numeq)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! solution of the linear program
        real(RK), intent(inout) :: x(:)

        ! coefficients in the function to minimize
        real(RK), intent(in) :: c(:)

        ! constraint matrix
        real(RK), intent(in) :: A(:, :)

        ! target vectors of constraint
        real(RK), intent(in) :: b(:)

        ! number of lower equal constraints
        integer(IK), intent(in) :: numle

        ! number of greater equal constraints
        integer(IK), intent(in) :: numge

        ! number of equality constraints
        integer(IK), intent(in) :: numeq


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: m, n, i1
        real(RK) :: A_h(size(A, 1), size(A, 2)+numle+numge)
        real(RK) :: c_h(size(c, 1)+numle+numge), b_h(size(b, 1))
        real(RK) :: x_h(size(x, 1)+numle+numge)
        real(RK) :: newA(size(A, 1), size(A, 2)+numle+numge)


        !##### ROUTINE CODE #######################################################

        ! check for sizes
        n = assert_eq(size(x, 1), size(c, 1), size(A, 2), 'solve_lin')
        m = assert_eq(size(A, 1), size(b, 1), 'solve_lin')

        ! check for correct inputs
        if(numle < 0)then
           call error('solve_lin', 'Number of lower equal constraints must '// &
               'not be negative')
        elseif(numge < 0)then
           call error('solve_lin', 'Number of greater equal constraints must '// &
               'not be negative')
        elseif(numeq < 0)then
           call error('solve_lin', 'Number of equality constraints must '// &
               'not be negative')
        elseif(numle+numge+numeq /= size(b,1))then
           call error('solve_lin', 'Number of equations does not match size of b')
        endif

        ! set up optimization problem
        A_h = 0d0
        A_h(1:m, 1:n) = A(:, :)
        do i1 = 1, numle
            A_h(i1, n+i1) = 1d0
        enddo
        do i1 = 1, numge
            A_h(numle+i1, n+numle+i1) = -1d0
        enddo

        ! check for negative bs
        b_h = b
        do i1 = 1, m
            if(b(i1) < 0d0)then
                A_h(i1, :) = -A_h(i1, :)
                b_h(i1) = -b_h(i1)
            endif
        enddo

        ! initialize c
        c_h = 0d0
        c_h(1:n) = c(:)

        call get_starting_value(x_h, A_h, b_h, newA)

        call solve_simplex(x_h, c_h, newA)

        x = x_h(1:n)

    contains


        !##############################################################################
        ! SUBROUTINE get_starting_value
        !
        ! Calculates a starting value for the linear program.
        !##############################################################################
        subroutine get_starting_value(x0, A, b, newA)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! starting value of the linear program
            real(RK), intent(out) :: x0(:)

            ! constraint matrix
            real(RK), intent(in) :: A(:, :)

            ! target vectors of constraint
            real(RK), intent(in) :: b(:)

            ! new matrix for simplex
            real(RK), intent(out) :: newA(:, :)


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: m, n, j
            real(RK) :: Astart(size(A,1), size(A,1)+size(A,2))
            real(RK) :: cstart(size(A,1)+size(A,2)), xstart(size(A,1)+size(A,2))

            !##### ROUTINE CODE #######################################################

            ! get sizes
            n = size(A, 2)
            m = size(A, 1)

            ! set up help problem
            cstart(1:n) = 0d0
            cstart(n+1:n+m) = 1d0

            ! set up help matrix
            Astart(1:m, 1:n) = A
            Astart(1:m, n+1:n+m) = 0d0
            do j = 1, m
                Astart(j,n+j) = 1d0
            enddo

            ! get initial guess
            xstart(1:n) = 0d0
            xstart(n+1:n+m) = b

            ! solve linear program
            call solve_simplex(xstart, cstart, Astart, newA)

            ! set starting value
            x0 = xstart(1:n)

        end subroutine get_starting_value


        !##############################################################################
        ! SUBROUTINE solve_simplex
        !
        ! Solves a linear program in canonic form.
        !##############################################################################
        subroutine solve_simplex(x, c, A, newA)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! starting value of the linear program and result
            real(RK), intent(inout) :: x(:)

            ! coefficients of the program
            real(RK), intent(in) :: c(:)

            ! constraint matrix
            real(RK), intent(in) :: A(:, :)

            ! new tableau if starting value calculation
            real(RK), intent(out), optional :: newA(:, :)


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: k, m, n, j, ibas, inot, piv(2), ihelp, i1, i2, n1, bhelp
            integer(IK) :: bas(size(A, 1)), nbas(size(A, 2)-size(A, 1))
            real(RK) :: alpha(size(A, 1), size(A,2)-size(A, 1))
            real(RK) :: gamma(size(A, 2)-size(A,1)), x0(size(A, 1))
            real(RK) :: pivcheck(size(A,1)), phelp, alpha_h(size(alpha, 2))

            !##### ROUTINE CODE #######################################################

            ! get sizes
            n = size(A, 2)
            m = size(A, 1)
            k = size(A, 2)-size(A,1)

            ! set up basis and non basis elements
            ibas = 1
            inot = 1
            do j = 1, n
                if(abs(x(j)) >= 1d-100)then
                    bas(ibas) = j
                    ibas = ibas + 1
                else
                    nbas(inot) = j
                    inot = inot + 1
                endif
            enddo

            ! set up x0 and c'x
            do ibas = 1, m
                x0(ibas) = x(bas(ibas))
            enddo

            ! set up alphas
            do inot = 1, k
                alpha(:, inot) = -A(:, nbas(inot))
            enddo

            ! set up gammas
            do inot = 1, k
                gamma(inot) = 0d0
                do ibas = 1, m
                    gamma(inot) = gamma(inot) + alpha(ibas, inot)*c(bas(ibas))
                enddo
                gamma(inot) = gamma(inot) + c(nbas(inot))
            enddo

            ! start algorithm
            do

                ! choose pivot column
                piv = 0
                do inot = 1, k
                    if(gamma(inot) < 0d0)then
                        piv(2) = inot
                        exit
                    endif
                enddo

                ! algorithm ends of no gamma < 0
                if(piv(2) == 0)exit

                ! else choose pivot row
                do ibas = 1, m
                    if(abs(alpha(ibas, piv(2))) >= 1d-100)then
                        pivcheck(ibas) = x0(ibas)/alpha(ibas, piv(2))
                    else
                        pivcheck(ibas) = -1d300
                    endif
                enddo
                phelp = -1d300
                do ibas = 1, m
                    if(alpha(ibas, piv(2)) < 0d0 .and. pivcheck(ibas) > phelp)then
                        phelp = pivcheck(ibas)
                        piv(1) = ibas
                    endif
                enddo

                ! no solution in piv(1) == 0
                if(piv(1) == 0)then
                    call error('solve_lin','Problem has no solution')
                endif

                ! Apply basis change
                Ihelp = nbas(piv(2))
                nbas(piv(2)) = bas(piv(1))
                bas(piv(1)) = Ihelp

                ! change pivot element
                alpha(piv(1), piv(2)) = 1d0/alpha(piv(1), piv(2))

                ! change pivot column
                do ibas = 1, m
                    if(ibas /= piv(1))then
                        alpha(ibas, piv(2)) = alpha(ibas, piv(2))* &
                            alpha(piv(1), piv(2))
                    endif
                enddo

                ! change pivot row
                do inot = 1, k
                    if(inot /= piv(2))then
                        alpha(piv(1), inot) = -alpha(piv(1), inot)* &
                            alpha(piv(1), piv(2))
                    endif
                enddo

                ! change other elements of alpha
                do ibas = 1, m
                    do inot = 1, k
                        if(ibas /= piv(1) .and. inot /= piv(2))then
                            alpha(ibas, inot) = alpha(ibas, inot) + &
                                alpha(ibas, piv(2))*alpha(piv(1), inot)/ &
                                alpha(piv(1), piv(2))
                        endif
                    enddo
                enddo

                ! change x0
                x0(piv(1)) = -x0(piv(1))*alpha(piv(1), piv(2))
                do ibas = 1, m
                    if(ibas /= piv(1))then
                        x0(ibas) = x0(ibas) + alpha(ibas, piv(2))*x0(piv(1))/ &
                            alpha(piv(1), piv(2))
                    endif
                enddo

                ! change gammas
                gamma(piv(2)) = gamma(piv(2))*alpha(piv(1), piv(2))
                do inot = 1, k
                    if(inot /= piv(2))then
                        gamma(inot) = gamma(inot) + alpha(piv(1), inot)* &
                            gamma(piv(2))/alpha(piv(1), piv(2))
                    endif
                enddo

            enddo

            ! get solution
            x = 0d0
            do ibas = 1, m
                x(bas(ibas)) = x0(ibas)
            enddo

            ! set up new tableau if needed
            if(present(newA))then

                ! check for existence of a solution
                n1 = 1
                do i1 = 1, n
                    if(c(i1) > 0d0)then
                        n1 = i1
                        exit
                    endif
                enddo
                if(any(x(n1:n) > 0d0))then
                    call error('solve_lin', &
                        'Linear Program does not have a solution')
                endif

                ! sort tableau in ascending order
                do i1 = m-1,1,-1
                    do i2 = 1, i1
                        if(bas(i2) > bas(i2+1))then
                            bhelp = bas(i2)
                            bas(i2) = bas(i2+1)
                            bas(i2+1) = bhelp

                            alpha_h = alpha(i2, :)
                            alpha(i2, :) = alpha(i2+1, :)
                            alpha(i2+1, :) = alpha_h
                        endif
                    enddo
                enddo

                ! get new matrix
                newA = 0d0
                do i1 = 1, k
                    if(nbas(i1) <= n1-1)then
                        newA(:, nbas(i1)) = -alpha(:, i1)
                    endif
                enddo
            endif

        end subroutine solve_simplex

    end subroutine solve_lin















!##############################################################################
!##############################################################################
! MODULE rootfinding
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE settol_root
    !
    ! For setting global tolerance level.
    !##############################################################################
    subroutine settol_root(tol)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! tolerance level
        real(RK), intent(in) :: tol


        !##### ROUTINE CODE #######################################################

        ! check whether tolerance level is valid
        if(tol > 0d0)then
            tbox_gftol_root = tol
        else
            call warning('settol_root', 'tolerance level is not valid')
        endif

    end subroutine settol_root


    !##############################################################################
    ! SUBROUTINE setiter_root
    !
    ! For setting maximum number of iterations.
    !##############################################################################
    subroutine setiter_root(iter)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! number of iterations
        integer(IK), intent(in) :: iter


        !##### ROUTINE CODE #######################################################

        ! check whether number of iterations is valid
        if(iter > 0)then
            itermax_root = iter
        else
            call warning('setiter_root', 'number of iterations is not valid')
        endif

    end subroutine setiter_root


    !##############################################################################
    ! SUBROUTINE newton_interpol
    !
    ! Find root of one-dimensional function by interpolatory newton method.
    !##############################################################################
    subroutine newton_interpol(x, funcv, check_return)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! initial guess and root of the function
        real(RK), intent(inout) :: x

        ! check is true if newton_interpol converged to local minimum or can make no
        !     further progress
        logical, intent(out), optional :: check_return

        !##### OTHER VARIABLES ####################################################

        real(RK), parameter :: eps = epsilon(x)
        real(RK) :: tolf, tolmin
        real(RK), parameter :: tolx = eps
        real(RK), parameter :: stpmx = 100d0
        real(RK) :: x1, x2, f1, f2, xnew, fnew, h
        integer(IK) :: its


        !##### INTERFACES #########################################################

        ! interface for the function
        interface
            function funcv(p)
                use Constants_mod, only: IK, RK; implicit none
                real(RK), intent(in) :: p
                real(RK) :: funcv
            end function funcv
        end interface


        !##### ROUTINE CODE #######################################################

        ! set tolerance levels
        tolf = tbox_gftol_root
        tolmin = tbox_gftol_root
        if(present(check_return))check_return = .false.

        ! initialize values
        x1 = x
        h = 1d-6*max(abs(x), 0.01d0)
        x2 = x + h

        ! calculate function values at x1, x2
        f1 = funcv(x1)
        f2 = funcv(x2)

        ! check if already in zero
        if(abs(f1) < 0.01d0*tolf)then
            x = x1
            if(present(check_return))check_return = .false.
            return
        endif

        if(abs(f2) < 0.01d0*tolf)then
            x = x2
            if(present(check_return))check_return = .false.
            return
        endif

        if(abs((1d0-f1/f2)/(x2-x1)) < tolmin .or. &
            abs((1d0-f1/f2)) < epsilon(1d0))then
            x = x1
            if(present(check_return))check_return = .true.
            return
        endif

        ! start iteration
        do its = 1, itermax_root

            ! calculate new point xnew
            xnew = x2 - (x2-x1)/(1d0-f1/f2)

            ! calculate new function value
            fnew = funcv(xnew)

            ! check wether function is small enough
            if(abs(fnew) < tolf)then
                x = xnew
                if(present(check_return))check_return = .false.
                return
            endif

            ! check whether you are in a minimum or cannot proceed further
            if(2d0*abs(xnew-x2) < tolx*abs(xnew+x2))then
                x = x2
                if(present(check_return))check_return = .true.
                return
            endif

            ! else set new data and repeat step
            x1 = x2
            f1 = f2
            x2 = xnew
            f2 = fnew

            if(abs((1d0-f1/f2)/(x2-x1)) < tolmin .or. &
               abs((1d0-f1/f2)) < epsilon(1d0))then
                x = x1
                if(present(check_return))check_return = .true.
                return
            endif
        enddo

        ! throw warning if newton didn't converge
        if(present(check_return))check_return = .true.

        x = xnew

    end subroutine newton_interpol


    !##############################################################################
    ! SUBROUTINE broydn
    !
    ! Find root of multidimensional function of.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    subroutine broydn(x, funcv, check_return)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! initial guess and root of the function
        real(RK), intent(inout) :: x(:)

        ! check is true if broydn converged to local minimum or can make no
        !     further progress
        logical, intent(out), optional :: check_return


        !##### OTHER VARIABLES ####################################################

        real(RK), parameter :: eps = epsilon(x)
        real(RK) :: tolf, tolmin
        real(RK), parameter :: tolx = eps
        real(RK), parameter :: stpmx = 100d0
        integer(IK) :: i, j, k, its, n
        real(RK) :: f, fold, stpmax
        real(RK), dimension(size(x)) :: c, d, fvcold, g, p, s, t, w, xold
        real(RK), dimension(size(x),size(x)) :: qt, r
        logical :: restrt, sing, check
        real(RK) :: fvec(size(x,1))


        !##### INTERFACES #########################################################

        ! interface for the function
        interface
            function funcv(p)
                use Constants_mod, only: IK, RK; implicit none
                real(RK), intent(in) :: p(:)
                real(RK) :: funcv(size(p, 1))
            end function funcv
        end interface


        !##### ROUTINE CODE #######################################################

        ! set tolerance levels
        tolf = tbox_gftol_root
        tolmin = tbox_gftol_root
        if(present(check_return))check_return = .false.

        ! get size of x
        n = size(x)

        ! calculate function euklidean norm at starting point
        f = fmin(x, funcv)

        ! check if root has been found
        if (maxval(abs(fvec(:))) < 0.01d0*tolf) then
            if(present(check_return))check_return = .false.
            return
        endif

        stpmax = stpmx*max(sqrt(dot_product(x(:),x(:))),dble(n))
        restrt = .true.

        ! iterate broydn steps
        do its=1,itermax_root

            ! If restart then calculate jacobian of function
            if (restrt) then

                ! calculate jacobian of func at x
                call fdjac(x, fvec, r, funcv)

                ! make q-r-decomposition of jacobian
                call qrdcmp(r, c, d, sing)

                ! throw error if jacobian is singular
                if(sing)then
                    call warning('fzero', 'singular jacobian')
                    if(present(check_return))check_return = .true.
                    return
                endif

                ! create unity matrix
                qt(:,:) = 0d0
                    do j = 1, n
                            qt(j, j) = 1d0
                    enddo

                ! for Q^T explicitly
                do k = 1, n-1
                    if (abs(c(k)) >= 1d-100) then
                        qt(k:n,:) = qt(k:n, :)-outerprod(r(k:n, k), &
                            matmul(r(k:n, k), qt(k:n, :)))/c(k)
                    endif
                enddo
                where(lower_triangle(n,n))r(:, :) = 0d0

                ! puts diagonal elements of R matrix to r
                do j = 1, n
                    r(j, j) = d(j)
                enddo

            ! else do Broydn update step
            else

                ! set up s as delta x
                s(:) = x(:)-xold(:)

                ! t = R*delta x
                do i = 1, n
                    t(i) = dot_product(r(i,i:n), s(i:n))
                enddo

                ! w = delta f - B*s = delta f - R*s*Q^T
                w(:) = fvec(:)-fvcold(:)-matmul(t(:), qt(:,:))

                ! if w entries are small enough, set them to zero
                where(abs(w(:)) < EPS*(abs(fvec(:))+abs(fvcold(:)))) w(:) = 0d0

                ! update for non-noisy components of w
                if(any(abs(w(:)) >= 1d-100))then

                    ! update t and s
                    t(:) = matmul(qt(:,:),w(:))
                    s(:)=s(:)/dot_product(s,s)

                    ! update R and Q^T
                    call qrupdt(r,qt,t,s)

                    ! get diagonal of matrix r
                    do j = 1, size(r,1)
                        d(j) = r(j,j)
                    enddo

                    ! if any diagonal value of r is 0, then jacobian is singular
                    if(any(abs(d(:)) <= 1d-100))then
                        call warning('fzero', 'singular jacobian')
                        if(present(check_return))check_return = .true.
                        return
                    endif
                endif
            endif

            ! perform the newton step by inverting jacobian
            p(:) = -matmul(qt(:,:), fvec(:))
            do i = 1, n
                g(i) = -dot_product(r(1:i,i), p(1:i))
            enddo

            ! store old x, function value and function norm
            xold(:) = x(:)
            fvcold(:) = fvec(:)
            fold = f

            ! solve linear equation with upper triangular matrix r
            call rsolv(r, d, p)

            ! searches along the new gradient direction for new x and f
            call lnsrch(xold, fold, g, p, x, f, stpmax, check, funcv)

            ! check whether root was found
            if(maxval(abs(fvec(:))) < tolf)then
                if(present(check_return))check_return = .false.
                return
            endif

            ! if check is true
            if(check)then

                ! check if improvement can be made, if not, return
                if(restrt .or. maxval(abs(g(:))*max(abs(x(:)), &
                        1d0)/max(f, 0.5d0*n)) < tolmin)then
                    if(present(check_return))check_return = check
                    return
                endif

                ! else calculate new jacobian
                restrt=.true.

            ! if check is false
            else

                ! do broydn step
                restrt=.false.

                ! check for convergence
                if (maxval((abs(x(:)-xold(:)))/max(abs(x(:)), &
                    1.0d0)) < tolx)then
                    if(present(check_return))check_return = check
                    return
                endif
            endif
        enddo

        ! throw warning if broydn didn't converge
        if(present(check_return))check_return = .true.


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! FUNCTION fdjac
        !
        ! Calculates finite difference jacobian.
        !##########################################################################
        subroutine fdjac(x, fvec, df, funcv)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! value where to calculate finite difference jacobian
            real(RK), intent(inout) :: x(:)

            ! function value at x
            real(RK), intent(in) :: fvec(:)

            ! resulting finite difference jacobian
            real(RK), intent(out) :: df(:, :)


            !##### OTHER VARIABLES ################################################

            real(RK), parameter :: eps = 1.0e-6
            integer(IK) :: j, n
            real(RK), dimension(size(x)) :: xsav, xph, h


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function funcv(p)
                    use Constants_mod, only: IK, RK; implicit none
                    real(RK), intent(in) :: p(:)
                    real(RK) :: funcv(size(p, 1))
                end function funcv
            end interface


            !##### ROUTINE CODE ###################################################

            ! check equality of sizes
            n = assert_eq(size(x), size(fvec), size(df,1), size(df,2), 'fdjac')

            ! store old x
            xsav = x

            ! calculate difference
            h = eps*abs(xsav)
            where(abs(h) <= 1d-100)h = EPS

            ! calculate x + h
            xph = xsav + h
            h = xph - xsav

            ! itertate over dimensions and calculate difference
            do j = 1, n
                x(j) = xph(j)
                df(:,j) = (funcv(x)-fvec(:))/h(j)
                x(j) = xsav(j)
            enddo

        end subroutine fdjac


        !##########################################################################
        ! FUNCTION lnsrch
        !
        ! Finds point along a line, given function value and gradient, where
        !     function has decreased sufficiently (for one dimensional function).
        !##########################################################################
        subroutine lnsrch(xold, fold, g, p, x, f, stpmax, check, funcv)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! point where to start line search
            real(RK), intent(in) :: xold(:)

            ! the old function value
            real(RK), intent(in) :: fold

            ! gradient at this point
            real(RK), intent(in) :: g(:)

            ! a line search direction
            real(RK), intent(inout) :: p(:)

            ! new value along the search line
            real(RK), intent(out) :: x(:)

            ! function value at new x
            real(RK), intent(out) :: f

            ! maximum size of steps such that lnsrch does not search un undefined
            !     areas
            real(RK), intent(in) :: stpmax

            ! is true if x is too close at xold
            logical, intent(out) :: check


            !##### OTHER VARIABLES ################################################

            real(RK), parameter :: alf = 1.0e-4
            real(RK), parameter :: tolx = epsilon(x)
            integer(IK) :: ndum
            real(RK) :: a, alam, alam2=0d0, alamin, b, disc, f2=0d0, pabs, rhs1, rhs2, &
                slope, tmplam


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function funcv(p)
                    use Constants_mod, only: IK, RK; implicit none
                    real(RK), intent(in) :: p(:)
                    real(RK) :: funcv(size(p, 1))
                end function funcv
            end interface


            !##### ROUTINE CODE ###################################################

            ! assert sizes or arrays
            ndum = assert_eq(size(g), size(p), size(x), size(xold), 'lnsrch')
            ndum = ndum

            ! set check's default value
            check=.false.

            ! calculate norm of p
            pabs = sqrt(dot_product(p, p))

            ! restrict p to maximum stepsize
            if(pabs > stpmax)p(:) = p(:)*stpmax/pabs

            ! calculate slope
            slope = dot_product(g, p)

            ! throw error if you would go uphill
            if(slope >= 0d0)then
                call warning('lnsrch', 'roundoff problem, I cannot go uphill')
                return
            endif

            ! calculate newton stepsize
            alamin = tolx/maxval(abs(p(:))/max(abs(xold(:)),1d0))
            alam = 1d0

            ! start iteration
            do
                ! calculate calculate new x
                x(:) = xold(:)+alam*p(:)

                ! calculate new function value at x
                f = fmin(x, funcv)

                ! if new x is not away enough return with check=true
                if(alam < alamin)then
                    x(:) = xold(:)
                    check = .true.
                    return

                ! if optimal value found return with false
                elseif(f <= fold+alf*alam*slope)then
                    return

                ! else do backtracking
                else
                    if(abs(alam -1d0) <= 1d-100)then
                        tmplam = -slope/(2d0*(f-fold-slope))
                    else
                        rhs1 = f-fold-alam*slope
                        rhs2 = f2-fold-alam2*slope
                        a = (rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
                        b = (-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
                        if(abs(a) <= 1d-100)then
                            tmplam = -slope/(2d0*b)
                        else
                            disc = b*b-3d0*a*slope
                            if(disc < 0d0)then
                                tmplam = 0.5d0*alam
                            elseif(b <= 0d0)then
                                tmplam = (-b+sqrt(disc))/(3d0*a)
                            else
                                tmplam = -slope/(b+sqrt(disc))
                            endif
                        endif
                        if(tmplam > 0.5d0*alam)tmplam = 0.5d0*alam
                    endif
                endif
                alam2 = alam
                f2 = f
                alam = max(tmplam,0.1d0*alam)
            enddo

        end subroutine lnsrch


        !##########################################################################
        ! FUNCTION fmin
        !
        ! Calculates vector norm of multidimensional function.
        !##########################################################################
        function fmin(x, funcv)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! value where to evaluate function
            real(RK), intent(in) :: x(:)

            ! euklidean square norm of function at x
            real(RK) :: fmin


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function funcv(p)
                    use Constants_mod, only: IK, RK; implicit none
                    real(RK), intent(in) :: p(:)
                    real(RK) :: funcv(size(p, 1))
                end function funcv
            end interface

            ! calculate function value
            fvec = funcv(x)

            ! calculate squared norm
            fmin = 0.5d0*dot_product(fvec, fvec)

        end function fmin


        !##########################################################################
        ! SUBROUTINE qrdcmp
        !
        ! Calculates QR decomposition of a matrix.
        !##########################################################################
        subroutine qrdcmp(a, c, d, sing)

            use Constants_mod, only: IK, RK; implicit none

            real(RK), intent(inout) :: a(:, :)
            real(RK), intent(out) :: c(:), d(:)
            logical, intent(out) :: sing
            integer(IK) :: k, n
            real(RK) :: scale, sigma

            n = assert_eq(size(a,1), size(a,2), size(c), size(d), 'qrdcmp')
            sing = .false.
            do k = 1, n-1
                scale = maxval(abs(a(k:n, k)))
                if(abs(scale) <= 1d-100)then
                        sing = .true.
                        c(k) = 0d0
                        d(k) = 0d0
                else
                        a(k:n, k) = a(k:n, k)/scale
                        sigma = sign(sqrt(dot_product(a(k:n, k),a(k:n, k))),a(k, k))
                        a(k,k) = a(k, k)+sigma
                        c(k) = sigma*a(k, k)
                        d(k) = -scale*sigma
                        a(k:n, k+1:n) = a(k:n, k+1:n)-outerprod(a(k:n, k),&
                                matmul(a(k:n, k),a(k:n, k+1:n)))/c(k)
                endif
            enddo
            d(n) = a(n, n)
            if (abs(d(n)) <= 1d-100) sing = .true.

        end subroutine qrdcmp


        !##########################################################################
        ! SUBROUTINE qrupdt
        !
        ! Updates qr-matrices.
        !##########################################################################
        subroutine qrupdt(r,qt,u,v)

            use Constants_mod, only: IK, RK; implicit none

            real(RK), intent(inout) :: r(:, :), qt(:, :)
            real(RK), intent(inout) :: u(:)
            real(RK), intent(in) :: v(:)
            integer(IK) :: i, k, n

            n = assert_eq((/ size(r,1), size(r,2), size(qt,1), size(qt,2), &
                size(u), size(v)/), 'qrupdt')
            k = n+1-ifirstloc(abs(u(n:1:-1)) >= 1d-100)
            if(k < 1)k=1
            do i = k-1, 1, -1
                call rotate(r,qt,i,u(i),-u(i+1))
                u(i) = pythag(u(i),u(i+1))
            enddo
            r(1,:) = r(1,:)+u(1)*v
            do i = 1,k-1
                call rotate(r,qt,i,r(i,i),-r(i+1,i))
            enddo
        end subroutine qrupdt

        !##########################################################################
        ! SUBROUTINE rsolv
        !
        ! Solves upper diagonal system.
        !##########################################################################
        subroutine rsolv(a, d, b)

            use Constants_mod, only: IK, RK; implicit none

            real(RK), intent(in) :: a(:, :), d(:)
            real(RK), intent(inout) :: b(:)
            integer(IK) :: i, n

            n = assert_eq(size(a,1), size(a,2), size(b), size(d), 'rsolv')
            b(n) = b(n)/d(n)
            do i = n-1, 1, -1
                    b(i) =( b(i)-dot_product(a(i, i+1:n),b(i+1:n)))/d(i)
            enddo

        end subroutine rsolv


        subroutine rotate(r, qt, i, a, b)

            use Constants_mod, only: IK, RK; implicit none

            real(RK), intent(inout) :: r(:, :), qt(:, :)
            integer(IK), intent(in) :: i
            real(RK), intent(in) :: a, b
            integer(IK) :: n
            real(RK) :: c, fact, s, temp(size(r,1))

            n = assert_eq(size(r,1), size(r,2), size(qt,1), size(qt,2), 'rotate')
            if(abs(a) <= 1d-100)then
                c = 0d0
                s = sign(1d0, b)
            elseif(abs(a) > abs(b))then
                fact = b/a
                c = sign(1d0/sqrt(1d0+fact**2), a)
                s = fact*c
            else
                fact = a/b
                s = sign(1d0/sqrt(1d0+fact**2), b)
                c=fact*s
            endif
            temp(i:n) = r(i, i:n)
            r(i, i:n) = c*temp(i:n)-s*r(i+1, i:n)
            r(i+1, i:n) = s*temp(i:n)+c*r(i+1, i:n)
            temp = qt(i, :)
            qt(i, :) = c*temp-s*qt(i+1, :)
            qt(i+1, :) = s*temp+c*qt(i+1, :)

        end subroutine rotate


        function pythag(a, b)

            use Constants_mod, only: IK, RK; implicit none

            real(RK), intent(in) :: a, b
            real(RK) :: pythag
            real(RK) :: absa, absb

            absa = abs(a)
            absb = abs(b)
            if(absa > absb)then
                pythag = absa*sqrt(1d0+(absb/absa)**2)
            else
                if(abs(absb) <= 1d-100)then
                    pythag = 0d0
                else
                    pythag = absb*sqrt(1d0+(absa/absb)**2)
                endif
            endif

        end function pythag


        function ifirstloc(mask)

            logical, intent(in) :: mask(:)
            integer(IK) :: ifirstloc, loca(1)

            loca = maxloc(merge(1, 0, mask))
            ifirstloc = loca(1)
            if(.not. mask(ifirstloc))ifirstloc = size(mask)+1

        end function ifirstloc


        function lower_triangle(j, k, extra)

            integer(IK), intent(in) :: j, k
            integer(IK), intent(in), optional :: extra
            logical :: lower_triangle(j, k)
            integer(IK) :: n

            n = 0
            if(present(extra))n = extra

            lower_triangle = (outerdiff(arth_i(1, 1, j), arth_i(1, 1, k)) > -n)

        end function lower_triangle


        function outerdiff(a, b)

            integer(IK), intent(in) :: a(:), b(:)
            integer(IK) :: outerdiff(size(a, 1),size(b, 1))

            outerdiff = spread(a, dim=2, ncopies=size(b, 1)) - &
                spread(b, dim=1, ncopies=size(a, 1))

        end function outerdiff


        function outerprod(a, b)

            real(RK), intent(in) :: a(:), b(:)
            real(RK) :: outerprod(size(a, 1),size(b, 1))

            outerprod = spread(a, dim=2, ncopies=size(b, 1)) * &
                spread(b, dim=1, ncopies=size(a, 1))

        end function outerprod


        function arth_i(first, increment, n)

            integer(IK), intent(in) :: first, increment, n
            integer(IK), parameter :: npar_arth = 16
            integer(IK), parameter :: npar2_arth = 8
            integer(IK) :: arth_i(n)
            integer(IK) :: k, k2, temp

            if(n > 0)arth_i(1) = first
            if(n <= npar_arth) then
                do k = 2, n
                    arth_i(k) = arth_i(k-1) + increment
                enddo
            else
                do k = 2, npar2_arth
                    arth_i(k) = arth_i(k-1) + increment
                enddo
                temp = increment*npar2_arth
                k = npar2_arth
                do
                    if(k >= n)exit
                    k2 = k+k
                    arth_i(k+1:min(k2,n)) = temp+arth_i(1:min(k,n-k))
                    temp = temp + temp
                    k = k2
                enddo
            endif
        end function arth_i

    end subroutine broydn















!##############################################################################
!##############################################################################
! MODULE splines
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE spline_interp1
    !
    ! Subroutine for one-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp1(yi, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real(RK), intent(in) :: yi(0:)

        ! coefficients for spline interpolation
        real(RK), intent(out) :: c(1:)


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j
        real(RK), allocatable :: r(:), d(:)


        !##### ROUTINE CODE #######################################################

        ! assert sizes for the two arrays do fit
        n = assert_eq(size(yi,1)+2, size(c,1), 'spline_interp')

        ! deallocate help arrays
        if(allocated(r))deallocate(r)
        if(allocated(d))deallocate(d)

        ! allocate help arrays
        allocate(r(n))
        allocate(d(n))

        ! calculate real n
        n = n-3

        ! calculate numerical derivatives at end points
        r(1) = (2d0*yi(0)-5d0*yi(1)+4d0*yi(2)-yi(3))/6d0
        r(n+3) = (2d0*yi(n)-5d0*yi(n-1)+4d0*yi(n-2)-yi(n-3))/6d0

        ! set rest of right side of equation system
        r(2:n+2) = yi(0:n)

        ! solve the spline interpolation equation system
        c(2) = (yi(0)-r(1))/6d0
        c(n+2) = (yi(n)-r(n+3))/6d0

        d(3) = 4d0
        r(3) = yi(1)-c(2)
        r(n+1) = yi(n-1)-c(n+2)

        do j = 4, n+1
            d(j) = 4d0-1d0/d(j-1)
            r(j) = r(j)-r(j-1)/d(j-1)
        enddo

        c(n+1) = r(n+1)/d(n+1)

        do j = n, 3, -1
            c(j) = (r(j)-c(j+1))/d(j)
        enddo

        c(1) = r(1)+2d0*c(2)-c(3)
        c(n+3) = r(n+3)+2d0*c(n+2)-c(n+1)

    end subroutine spline_interp1


    !##############################################################################
    ! SUBROUTINE spline_interp2
    !
    ! Subroutine for two-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp2(yi, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real(RK), intent(in) :: yi(0:, 0:)

        ! coefficients for spline interpolation
        real(RK), intent(out) :: c(1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n(2), j
        real(RK), allocatable :: tempc(:, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 2
            n(j) = assert_eq(size(yi,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, 0:n(2)))

        ! calculate temporary coefficients
        do j = 0, n(2)
            call spline_interp1(yi(:, j), tempc(:, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            call spline_interp1(tempc(j, :), c(j, :))
        enddo

    end subroutine spline_interp2


    !##############################################################################
    ! SUBROUTINE spline_interp3
    !
    ! Subroutine for three-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp3(yi, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real(RK), intent(in) :: yi(0:, 0:, 0:)

        ! coefficients for spline interpolation
        real(RK), intent(out) :: c(1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n(3), j, j2
        real(RK), allocatable :: tempc(:, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 3
            n(j) = assert_eq(size(yi,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, 0:n(3)))

        ! calculate temporary coefficients
        do j = 0, n(3)
            call spline_interp2(yi(:, :, j), tempc(:, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                call spline_interp1(tempc(j, j2, :), c(j, j2, :))
            enddo
        enddo

    end subroutine spline_interp3


    !##############################################################################
    ! SUBROUTINE spline_interp4
    !
    ! Subroutine for four-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp4(yi, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:)

        ! coefficients for spline interpolation
        real(RK), intent(out) :: c(1:, 1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n(4), j, j2, j3
        real(RK), allocatable :: tempc(:, :, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 4
            n(j) = assert_eq(size(yi,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, n(3)+3, 0:n(4)))

        ! calculate temporary coefficients
        do j = 0, n(4)
            call spline_interp3(yi(:, :, :, j), tempc(:, :, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                do j3 = 1, n(3)+3
                    call spline_interp1(tempc(j, j2, j3, :), c(j, j2, j3, :))
                enddo
            enddo
        enddo

    end subroutine spline_interp4


    !##############################################################################
    ! SUBROUTINE spline_interp5
    !
    ! Subroutine for five-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp5(yi, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:, 0:)

        ! coefficients for spline interpolation
        real(RK), intent(out) :: c(1:, 1:, 1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n(5), j, j2, j3, j4
        real(RK), allocatable :: tempc(:, :, :, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 5
            n(j) = assert_eq(size(yi,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, n(3)+3, n(4)+3, 0:n(5)))

        ! calculate temporary coefficients
        do j = 0, n(5)
            call spline_interp4(yi(:, :, :, :, j), tempc(:, :, :, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                do j3 = 1, n(3)+3
                    do j4 = 1, n(4)+3
                        call spline_interp1(tempc(j, j2, j3, j4, :), &
                            c(j, j2, j3, j4, :))
                    enddo
                enddo
            enddo
        enddo

    end subroutine spline_interp5


    !##############################################################################
    ! SUBROUTINE spline_interp6
    !
    ! Subroutine for six-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp6(yi, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:, 0:, 0:)

        ! coefficients for spline interpolation
        real(RK), intent(out) :: c(1:, 1:, 1:, 1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n(6), j, j2, j3, j4, j5
        real(RK), allocatable :: tempc(:, :, :, :, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 6
            n(j) = assert_eq(size(yi,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, n(3)+3, n(4)+3, n(5)+3, 0:n(6)))

        ! calculate temporary coefficients
        do j = 0, n(6)
            call spline_interp5(yi(:, :, :, :, :, j), tempc(:, :, :, :, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                do j3 = 1, n(3)+3
                    do j4 = 1, n(4)+3
                        do j5 = 1, n(5)+3
                            call spline_interp1(tempc(j, j2, j3, j4, j5, :), &
                                c(j, j2, j3, j4, j5, :))
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine spline_interp6


    !##############################################################################
    ! SUBROUTINE spline_interp7
    !
    ! Subroutine for seven-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp7(yi, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:, 0:, 0:, 0:)

        ! coefficients for spline interpolation
        real(RK), intent(out) :: c(1:, 1:, 1:, 1:, 1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n(7), j, j2, j3, j4, j5, j6
        real(RK), allocatable :: tempc(:, :, :, :, :, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 7
            n(j) = assert_eq(size(yi,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, n(3)+3, n(4)+3, n(5)+3, n(6)+3, 0:n(7)))

        ! calculate temporary coefficients
        do j = 0, n(7)
            call spline_interp6(yi(:, :, :, :, :, :, j), tempc(:, :, :, :, :, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                do j3 = 1, n(3)+3
                    do j4 = 1, n(4)+3
                        do j5 = 1, n(5)+3
                            do j6 = 1, n(6)+3
                                call spline_interp1(tempc(j, j2, j3, &
                                    j4, j5, j6, :),c(j, j2, j3, j4, j5, j6, :))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine spline_interp7


    !##############################################################################
    ! FUNCTION spline1
    !
    ! Function for evaluation of one-dimensional spline.
    !##############################################################################
    function spline1(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:)

        ! value of spline function
        real(RK) :: spline1


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n1, j1, p1, q1
        real(RK) :: phi1, xtemp1


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n1 = size(c, 1)

        ! calculate left and right summation end point
        p1 = max(floor(x)+1, 1)
        q1 = min(p1+3, n1)

        spline1 = 0d0

        do j1 = p1, q1

            ! calculate value where to evaluate basis function
            xtemp1 = abs(x-j1+2)

            ! calculate basis function
            if(xtemp1 <= 1d0)then
                phi1 = 4d0+xtemp1**2*(3d0*xtemp1-6d0)
            elseif(xtemp1 <= 2d0)then
                phi1 = (2d0-xtemp1)**3
            else
                phi1 = 0d0
            endif

            ! calculate spline value
            spline1 = spline1+c(j1)*phi1
        enddo

    end function spline1


    !##############################################################################
    ! FUNCTION spline1_grid
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline1_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:)

        ! left interval endpoint
        real(RK), intent(in) :: left

        ! right interval endpoint
        real(RK), intent(in) :: right

        ! growth rate of grid
        real(RK), intent(in), optional :: growth

        ! value of spline function
        real(RK) :: spline1_grid


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n
        real(RK) :: xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        n = size(c, 1)-3

        ! invert grid
        if(present(growth))then
            xtemp = grid_Inv_Grow(x, left, right, growth, n)
        else
            xtemp = grid_Inv_Equi(x, left, right, n)
        endif

        ! calculate spline value
        spline1_grid = spline1(xtemp, c)

    end function spline1_grid


    !##############################################################################
    ! FUNCTION spline1_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline1_complete(x, yi, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:)

        ! left interval endpoint
        real(RK), intent(in) :: left

        ! right interval endpoint
        real(RK), intent(in) :: right

        ! growth rate of grid
        real(RK), intent(in), optional :: growth

        ! value of spline function
        real(RK) :: spline1_complete


        !##### OTHER VARIABLES ####################################################

        real(RK) :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline1_complete_m((/x/), yi, left, right, growth)
        else
            spline_temp = spline1_complete_m((/x/), yi, left, right)
        endif

        ! paste data
        spline1_complete = spline_temp(1)

    end function spline1_complete


    !##############################################################################
    ! FUNCTION spline1_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline1_complete_m(x, yi, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(1:)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:)

        ! left interval endpoint
        real(RK), intent(in) :: left

        ! right interval endpoint
        real(RK), intent(in) :: right

        ! growth rate of grid
        real(RK), intent(in), optional :: growth

        ! value of spline function
        real(RK) :: spline1_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK) :: c(1:size(yi, 1)+2)
        real(RK) :: xtemp(1:size(x, 1))
        integer(IK) :: n, m, j


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        n = size(yi, 1)-1

        ! calculate number of evaluation points
        m = size(x, 1)

        ! invert grid for every evaluation point
        if(present(growth))then
            xtemp(:) = grid_Inv_Grow(x(:), left, right, growth, n)
        else
            xtemp(:) = grid_Inv_Equi(x(:), left, right, n)
        endif

        ! interpolate data
        call spline_interp1(yi, c)

        ! calculate spline values at point
        do j = 1, m
            spline1_complete_m(j) = spline1(xtemp(j), c)
        enddo

    end function spline1_complete_m


    !##############################################################################
    ! FUNCTION spline2
    !
    ! Function for evaluation of two-dimensional spline.
    !##############################################################################
    function spline2(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(2)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:)

        ! value of spline function
        real(RK) :: spline2


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n(2), p(2), q(2)
        integer(IK) :: j1, j2
        real(RK) :: phi1, xtemp1, phi2, xtemp2
        real(RK) :: s2


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n(1) = size(c, 1)
        n(2) = size(c, 2)

        ! calculate left and right summation end point
        p = max(floor(x)+1, 1)
        q = min(p+3, n)

        spline2 = 0d0

        do j1 = p(1), q(1)

            ! calculate value where to evaluate basis function
            xtemp1 = abs(x(1)-j1+2)

            ! calculate basis function
            if(xtemp1 <= 1d0)then
                phi1 = 4d0+xtemp1**2*(3d0*xtemp1-6d0)
            elseif(xtemp1 <= 2d0)then
                phi1 = (2d0-xtemp1)**3
            else
                phi1 = 0d0
            endif


            !#### calculate spline for second dimension ###########################

            s2 = 0d0

            do j2 = p(2), q(2)

                ! calculate value where to evaluate basis function
                xtemp2 = abs(x(2)-j2+2)

                ! calculate basis function
                if(xtemp2 <= 1d0)then
                    phi2 = 4d0+xtemp2**2*(3d0*xtemp2-6d0)
                elseif(xtemp2 <= 2d0)then
                    phi2 = (2d0-xtemp2)**3
                else
                    phi2 = 0d0
                endif

                ! calculate spline value
                s2 = s2+c(j1, j2)*phi2
            enddo

            ! calculate spline value
            spline2 = spline2+s2*phi1
        enddo

    end function spline2


    !##############################################################################
    ! FUNCTION spline2_grid
    !
    ! Function for evaluation of two-dimensional spline, includ inverting grid.
    !##############################################################################
    function spline2_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(2)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:)

        ! left interval endpoint
        real(RK), intent(in) :: left(2)

        ! right interval endpoint
        real(RK), intent(in) :: right(2)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(2)

        ! value of spline function
        real(RK) :: spline2_grid


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j
        real(RK) :: xtemp(2)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline2_grid = spline2(xtemp, c)

    end function spline2_grid


    !##############################################################################
    ! FUNCTION spline2_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline2_complete(x, yi, left, right, growth)


        integer(IK), parameter :: dim = 2

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(dim)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline2_complete


        !##### OTHER VARIABLES ####################################################

        real(RK) :: xtemp(1, dim)
        real(RK) :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline2_complete_m(xtemp, yi, left, right, growth)
        else
            spline_temp = spline2_complete_m(xtemp, yi, left, right)
        endif

        ! paste data
        spline2_complete = spline_temp(1)

    end function spline2_complete


    !##############################################################################
    ! FUNCTION spline2_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline2_complete_m(x, yi, left, right, growth)

        integer(IK), parameter :: dim = 2

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline2_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK) :: c(1:size(yi, 1)+2, 1:size(yi, 2)+2)
        real(RK) :: xtemp(1:size(x, 1), dim)
        integer(IK) :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp2(yi, c)

        ! calculate spline values at point
        do j = 1, m
            spline2_complete_m(j) = spline2(xtemp(j, :), c)
        enddo

    end function spline2_complete_m


    !##############################################################################
    ! FUNCTION spline3
    !
    ! Function for evaluation of three-dimensional spline.
    !##############################################################################
    function spline3(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(3)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:)

        ! value of spline function
        real(RK) :: spline3


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n(3), p(3), q(3)
        integer(IK) :: j1, j2, j3
        real(RK) :: phi1, xtemp1, phi2, xtemp2, phi3, xtemp3
        real(RK) :: s2, s3


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n(1) = size(c, 1)
        n(2) = size(c, 2)
        n(3) = size(c, 3)

        ! calculate left and right summation end point
        p = max(floor(x)+1, 1)
        q = min(p+3, n)

        spline3 = 0d0

        do j1 = p(1), q(1)

            ! calculate value where to evaluate basis function
            xtemp1 = abs(x(1)-j1+2)

            ! calculate basis function
            if(xtemp1 <= 1d0)then
                phi1 = 4d0+xtemp1**2*(3d0*xtemp1-6d0)
            elseif(xtemp1 <= 2d0)then
                phi1 = (2d0-xtemp1)**3
            else
                phi1 = 0d0
            endif


            !#### calculate spline for second dimension ###########################

            s2 = 0d0

            do j2 = p(2), q(2)

                ! calculate value where to evaluate basis function
                xtemp2 = abs(x(2)-j2+2)

                ! calculate basis function
                if(xtemp2 <= 1d0)then
                    phi2 = 4d0+xtemp2**2*(3d0*xtemp2-6d0)
                elseif(xtemp2 <= 2d0)then
                    phi2 = (2d0-xtemp2)**3
                else
                    phi2 = 0d0
                endif


                !#### calculate spline for second dimension #######################

                s3 = 0d0

                do j3 = p(3), q(3)

                    ! calculate value where to evaluate basis function
                    xtemp3 = abs(x(3)-j3+2)

                    ! calculate basis function
                    if(xtemp3 <= 1d0)then
                        phi3 = 4d0+xtemp3**2*(3d0*xtemp3-6d0)
                    elseif(xtemp3 <= 2d0)then
                        phi3 = (2d0-xtemp3)**3
                    else
                        phi3 = 0d0
                    endif

                    ! calculate spline value
                    s3 = s3+c(j1, j2, j3)*phi3
                enddo

                ! calculate spline value
                s2 = s2+s3*phi2
            enddo

            ! calculate spline value
            spline3 = spline3+s2*phi1
        enddo

    end function spline3


    !##############################################################################
    ! FUNCTION spline3_grid
    !
    ! Function for evaluation of three-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline3_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(3)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:)

        ! left interval endpoint
        real(RK), intent(in) :: left(3)

        ! right interval endpoint
        real(RK), intent(in) :: right(3)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(3)

        ! value of spline function
        real(RK) :: spline3_grid


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j
        real(RK) :: xtemp(3)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline3_grid = spline3(xtemp, c)

    end function spline3_grid


    !##############################################################################
    ! FUNCTION spline3_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline3_complete(x, yi, left, right, growth)


        integer(IK), parameter :: dim = 3

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(dim)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline3_complete


        !##### OTHER VARIABLES ####################################################

        real(RK) :: xtemp(1, dim)
        real(RK) :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline3_complete_m(xtemp, yi, left, right, growth)
        else
            spline_temp = spline3_complete_m(xtemp, yi, left, right)
        endif

        ! paste data
        spline3_complete = spline_temp(1)

    end function spline3_complete


    !##############################################################################
    ! FUNCTION spline3_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline3_complete_m(x, yi, left, right, growth)

        integer(IK), parameter :: dim = 3

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline3_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK) :: c(1:size(yi, 1)+2, 1:size(yi, 2)+2, 1:size(yi, 3)+2)
        real(RK) :: xtemp(1:size(x, 1), dim)
        integer(IK) :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp3(yi, c)

        ! calculate spline values at point
        do j = 1, m
            spline3_complete_m(j) = spline3(xtemp(j, :), c)
        enddo

    end function spline3_complete_m


    !##############################################################################
    ! FUNCTION spline4
    !
    ! Function for evaluation of four-dimensional spline.
    !##############################################################################
    function spline4(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(4)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:, 1:)

        ! value of spline function
        real(RK) :: spline4


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j, p, q
        real(RK) :: phi, xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n = size(c, 1)

        ! calculate left and right summation end point
        p = max(floor(x(1))+1, 1)
        q = min(p+3, n)

        spline4 = 0d0

        do j = p, q

            ! calculate value where to evaluate basis function
            xtemp = abs(x(1)-j+2)

            ! calculate basis function
            if(xtemp <= 1d0)then
                phi = 4d0+xtemp**2*(3d0*xtemp-6d0)
            elseif(xtemp <= 2d0)then
                phi = (2d0-xtemp)**3
            else
                phi = 0d0
            endif

            ! calculate spline value
            spline4 = spline4+spline3(x(2:4), c(j, :, :, :))*phi
        enddo

    end function spline4


    !##############################################################################
    ! FUNCTION spline4_grid
    !
    ! Function for evaluation of four-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline4_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(4)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:, 1:)

        ! left interval endpoint
        real(RK), intent(in) :: left(4)

        ! right interval endpoint
        real(RK), intent(in) :: right(4)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(4)

        ! value of spline function
        real(RK) :: spline4_grid


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j
        real(RK) :: xtemp(4)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline4_grid = spline4(xtemp, c)

    end function spline4_grid


    !##############################################################################
    ! FUNCTION spline4_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline4_complete(x, yi, left, right, growth)


        integer(IK), parameter :: dim = 4

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(dim)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline4_complete


        !##### OTHER VARIABLES ####################################################

        real(RK) :: xtemp(1, dim)
        real(RK) :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline4_complete_m(xtemp, yi, left, right, growth)
        else
            spline_temp = spline4_complete_m(xtemp, yi, left, right)
        endif

        ! paste data
        spline4_complete = spline_temp(1)

    end function spline4_complete


    !##############################################################################
    ! FUNCTION spline4_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline4_complete_m(x, yi, left, right, growth)

        integer(IK), parameter :: dim = 4

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline4_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK) :: c(1:size(yi, 1)+2, 1:size(yi, 2)+2, 1:size(yi, 3)+2, &
            1:size(yi, 4)+2)
        real(RK) :: xtemp(1:size(x, 1), dim)
        integer(IK) :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp4(yi, c)

        ! calculate spline values at point
        do j = 1, m
            spline4_complete_m(j) = spline4(xtemp(j, :), c)
        enddo

    end function spline4_complete_m


    !##############################################################################
    ! FUNCTION spline5
    !
    ! Function for evaluation of five-dimensional spline.
    !##############################################################################
    function spline5(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(5)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:, 1:, 1:)

        ! value of spline function
        real(RK) :: spline5


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j, p, q
        real(RK) :: phi, xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n = size(c, 1)

        ! calculate left and right summation end point
        p = max(floor(x(1))+1, 1)
        q = min(p+3, n)

        spline5 = 0d0

        do j = p, q

            ! calculate value where to evaluate basis function
            xtemp = abs(x(1)-j+2)

            ! calculate basis function
            if(xtemp <= 1d0)then
                phi = 4d0+xtemp**2*(3d0*xtemp-6d0)
            elseif(xtemp <= 2d0)then
                phi = (2d0-xtemp)**3
            else
                phi = 0d0
            endif

            ! calculate spline value
            spline5 = spline5+spline4(x(2:5), c(j, :, :, :, :))*phi
        enddo

    end function spline5


    !##############################################################################
    ! FUNCTION spline5_grid
    !
    ! Function for evaluation of five-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline5_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(5)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:, 1:, 1:)

        ! left interval endpoint
        real(RK), intent(in) :: left(5)

        ! right interval endpoint
        real(RK), intent(in) :: right(5)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(5)

        ! value of spline function
        real(RK) :: spline5_grid


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j
        real(RK) :: xtemp(5)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline5_grid = spline5(xtemp, c)

    end function spline5_grid


    !##############################################################################
    ! FUNCTION spline5_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline5_complete(x, yi, left, right, growth)


        integer(IK), parameter :: dim = 5

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(dim)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline5_complete


        !##### OTHER VARIABLES ####################################################

        real(RK) :: xtemp(1, dim)
        real(RK) :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline5_complete_m(xtemp, yi, left, right, growth)
        else
            spline_temp = spline5_complete_m(xtemp, yi, left, right)
        endif

        ! paste data
        spline5_complete = spline_temp(1)

    end function spline5_complete


    !##############################################################################
    ! FUNCTION spline5_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline5_complete_m(x, yi, left, right, growth)

        integer(IK), parameter :: dim = 5

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline5_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK) :: c(1:size(yi, 1)+2, 1:size(yi, 2)+2, 1:size(yi, 3)+2, &
            1:size(yi, 4)+2, 1:size(yi, 5)+2)
        real(RK) :: xtemp(1:size(x, 1), dim)
        integer(IK) :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp5(yi, c)

        ! calculate spline values at point
        do j = 1, m
            spline5_complete_m(j) = spline5(xtemp(j, :), c)
        enddo

    end function spline5_complete_m


    !##############################################################################
    ! FUNCTION spline6
    !
    ! Function for evaluation of six-dimensional spline.
    !##############################################################################
    function spline6(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(6)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:, 1:, 1:, 1:)

        ! value of spline function
        real(RK) :: spline6


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j, p, q
        real(RK) :: phi, xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n = size(c, 1)

        ! calculate left and right summation end point
        p = max(floor(x(1))+1, 1)
        q = min(p+3, n)

        spline6 = 0d0

        do j = p, q

            ! calculate value where to evaluate basis function
            xtemp = abs(x(1)-j+2)

            ! calculate basis function
            if(xtemp <= 1d0)then
                phi = 4d0+xtemp**2*(3d0*xtemp-6d0)
            elseif(xtemp <= 2d0)then
                phi = (2d0-xtemp)**3
            else
                phi = 0d0
            endif

            ! calculate spline value
            spline6 = spline6+spline5(x(2:6), c(j, :, :, :, :, :))*phi
        enddo

    end function spline6


    !##############################################################################
    ! FUNCTION spline6_grid
    !
    ! Function for evaluation of six-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline6_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(6)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:, 1:, 1:, 1:)

        ! left interval endpoint
        real(RK), intent(in) :: left(6)

        ! right interval endpoint
        real(RK), intent(in) :: right(6)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(6)

        ! value of spline function
        real(RK) :: spline6_grid


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j
        real(RK) :: xtemp(6)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline6_grid = spline6(xtemp, c)

    end function spline6_grid


    !##############################################################################
    ! FUNCTION spline6_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline6_complete(x, yi, left, right, growth)


        integer(IK), parameter :: dim = 6

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(dim)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline6_complete


        !##### OTHER VARIABLES ####################################################

        real(RK) :: xtemp(1, dim)
        real(RK) :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline6_complete_m(xtemp, yi, left, right, growth)
        else
            spline_temp = spline6_complete_m(xtemp, yi, left, right)
        endif

        ! paste data
        spline6_complete = spline_temp(1)

    end function spline6_complete


    !##############################################################################
    ! FUNCTION spline6_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline6_complete_m(x, yi, left, right, growth)

        integer(IK), parameter :: dim = 6

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline6_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK) :: c(1:size(yi, 1)+2, 1:size(yi, 2)+2, 1:size(yi, 3)+2, &
            1:size(yi, 4)+2, 1:size(yi, 5)+2, 1:size(yi, 6)+2)
        real(RK) :: xtemp(1:size(x, 1), dim)
        integer(IK) :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp6(yi, c)

        ! calculate spline values at point
        do j = 1, m
            spline6_complete_m(j) = spline6(xtemp(j, :), c)
        enddo

    end function spline6_complete_m


    !##############################################################################
    ! FUNCTION spline7
    !
    ! Function for evaluation of seven-dimensional spline.
    !##############################################################################
    function spline7(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(7)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:, 1:, 1:, 1:, 1:)

        ! value of spline function
        real(RK) :: spline7


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j, p, q
        real(RK) :: phi, xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n = size(c, 1)

        ! calculate left and right summation end point
        p = max(floor(x(1))+1, 1)
        q = min(p+3, n)

        spline7 = 0d0

        do j = p, q

            ! calculate value where to evaluate basis function
            xtemp = abs(x(1)-j+2)

            ! calculate basis function
            if(xtemp <= 1d0)then
                phi = 4d0+xtemp**2*(3d0*xtemp-6d0)
            elseif(xtemp <= 2d0)then
                phi = (2d0-xtemp)**3
            else
                phi = 0d0
            endif

            ! calculate spline value
            spline7 = spline7+spline6(x(2:7), c(j, :, :, :, :, :, :))*phi
        enddo

    end function spline7


    !##############################################################################
    ! FUNCTION spline7_grid
    !
    ! Function for evaluation of seven-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline7_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(7)

        ! coefficients for spline interpolation
        real(RK), intent(in) :: c(1:, 1:, 1:, 1:, 1:, 1:, 1:)

        ! left interval endpoint
        real(RK), intent(in) :: left(7)

        ! right interval endpoint
        real(RK), intent(in) :: right(7)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(7)

        ! value of spline function
        real(RK) :: spline7_grid


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, j
        real(RK) :: xtemp(7)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline7_grid = spline7(xtemp, c)

    end function spline7_grid


    !##############################################################################
    ! FUNCTION spline7_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline7_complete(x, yi, left, right, growth)


        integer(IK), parameter :: dim = 7

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(dim)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline7_complete


        !##### OTHER VARIABLES ####################################################

        real(RK) :: xtemp(1, dim)
        real(RK) :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline7_complete_m(xtemp, yi, left, right, growth)
        else
            spline_temp = spline7_complete_m(xtemp, yi, left, right)
        endif

        ! paste data
        spline7_complete = spline_temp(1)

    end function spline7_complete


    !##############################################################################
    ! FUNCTION spline7_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline7_complete_m(x, yi, left, right, growth)

        integer(IK), parameter :: dim = 7

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real(RK), intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real(RK), intent(in) :: yi(0:, 0:, 0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real(RK), intent(in) :: left(dim)

        ! right interval endpoint
        real(RK), intent(in) :: right(dim)

        ! growth rate of grid
        real(RK), intent(in), optional :: growth(dim)

        ! value of spline function
        real(RK) :: spline7_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real(RK) :: c(1:size(yi, 1)+2, 1:size(yi, 2)+2, 1:size(yi, 3)+2, &
            1:size(yi, 4)+2, 1:size(yi, 5)+2, 1:size(yi, 6)+2, 1:size(yi, 7)+2)
        real(RK) :: xtemp(1:size(x, 1), dim)
        integer(IK) :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(yi, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp7(yi, c)

        ! calculate spline values at point
        do j = 1, m
            spline7_complete_m(j) = spline7(xtemp(j, :), c)
        enddo

    end function spline7_complete_m















!##############################################################################
!##############################################################################
! MODULE gaussian_int
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE legendre
    !
    ! Calculates Gauss-Legendre abscissas and weights on [x1, x2].
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992).
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
    subroutine legendre(x1, x2, x, w)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! left interval point on which integral should be calculated
        real(RK), intent(in) :: x1

        ! left interval point on which integral should be calculated
        real(RK), intent(in) :: x2

        ! abscissas of gaussian integration formula
        real(RK), intent(out) :: x(:)

        ! weights of gaussian integration formula
        real(RK), intent(out) :: w(:)


        !##### OTHER VARIABLES ####################################################

        real(RK), parameter :: eps = 3.0e-14
        integer(IK), parameter :: maxits = 10
        real(RK), parameter :: pi = 3.14159265358979d0
        integer(IK) :: its, j, m, n
        real(RK) :: xl, xm
        real(RK), dimension((size(x)+1)/2) :: p1, p2, p3, pp, z, z1
        logical, dimension((size(x)+1)/2) :: unfinished


        !##### ROUTINE CODE #######################################################

        ! assert size equality
        n = assert_eq(size(x), size(w), 'legendre')

        ! calculate only up to (n+1)/2 due to symmetry
        m = (n+1)/2

        ! calculate interval midpoint
        xm = 0.5d0*(x2+x1)

        ! calculate half of interval length
        xl = 0.5d0*(x2-x1)

        ! set initial guess for the roots
        z = cos(pi*(arth(1,1,m)-0.25d0)/(n+0.5d0))

        ! initialized unfinished
        unfinished = .true.

        ! iterate Newton steps up to maximum iterations
        do its = 1, maxits

            ! calculate Legendre polynomial at z where root has not yet been found

            ! initialize p1 and p2
            where (unfinished)
                p1 = 1d0
                p2 = 0d0
            endwhere

            ! calculate polynomial value at z by recursive formula
            do j = 1, n

                ! only where root has not yet been found
                where (unfinished)

                    ! the polynomial of order n - 2
                    p3 = p2

                    ! the polynomial of order n - 1
                    p2 = p1

                    ! the legendre polynomial
                    p1 = ((2d0*j-1d0)*z*p2-(j-1d0)*p3)/j
                endwhere
            enddo

            ! calculate derivative of polynomial p1 at z
            where (unfinished)

                ! derivative
                pp = n*(z*p1-p2)/(z*z-1d0)

                ! store old z
                z1 = z

                ! perform the newton step
                z = z1-p1/pp

                ! check for difference between old and new guess being small enough
                unfinished=(abs(z-z1) > EPS)
            endwhere

            ! if all values have sufficiently converged, stop iteration
            if (.not. any(unfinished)) exit
        end do

        ! throw error message if not sufficiently convergerd
        if(its == maxits+1)call error('legendre', 'too many iterations')

        ! else calculate abscissas
        x(1:m) = xm-xl*z

        ! symmetry for abscissas
        x(n:n-m+1:-1) = xm+xl*z

        ! calculate weights
        w(1:m) = 2d0*xl/((1d0-z**2)*pp**2)

        ! symmetry for weights
        w(n:n-m+1:-1) = w(1:m)

    end subroutine legendre


    !##############################################################################
    ! FUNCTION arth
    !
    ! Calculates incremented array from first with n entries.
    !##############################################################################
    function arth(first, increment, n)

        integer(IK), intent(in) :: first, increment, n
        integer(IK), parameter :: npar_arth = 16
        integer(IK), parameter :: npar2_arth = 8
        integer(IK) :: arth(n)
        integer(IK) :: k, k2, temp

        ! initialize first element
        if(n > 0)arth(1) = first

        ! calculate by hand if n <= 16
        if(n <= npar_arth) then
            do k = 2, n
                arth(k) = arth(k-1) + increment
            enddo

        ! else set entries stepwise by 8 steps
        else
            do k = 2, npar2_arth
                arth(k) = arth(k-1) + increment
            enddo
            temp = increment*npar2_arth
            k = npar2_arth
            do
                if(k >= n)exit
                k2 = k+k
                arth(k+1:min(k2,n)) = temp+arth(1:min(k,n-k))
                temp = temp + temp
                k = k2
            enddo
        endif

    end function arth















!##############################################################################
!##############################################################################
! MODULE AR_discrete
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE discretize_AR
    !
    ! Discretizes an AR(1) process of the form z_j = \rho*z_{j-1} + eps using
    !     the Rouwenhorst method.
    !
    ! REFERENCE: Kopecky, K.A., Suen, R.M.H., Finite state Markov-chain
    !            approximations to highly persistent processes, Review of Economic
    !            Dynamics, Vol. 13, No. 3, 2010, 701-714.
    !##############################################################################
    subroutine discretize_AR(rho, mu, sigma_eps, z, pi, w)

       use Constants_mod, only: IK, RK; implicit none


       !##### INPUT/OUTPUT VARIABLES #############################################

       ! autoregression parameter
       real(RK), intent(in) :: rho

       ! unconditional mean of the process
       real(RK), intent(in) :: mu

       ! variance of the shock
       real(RK), intent(in) :: sigma_eps

       ! discrete shock values
       real(RK), intent(out) :: z(:)

       ! transition matrix
       real(RK), intent(out) :: pi(:, :)

       ! the stationary distribution
       real(RK), intent(out), optional :: w(:)


       !##### OTHER VARIABLES ####################################################

       integer(IK) :: n, in
       real(RK) :: psi, sigma_eta


       !##### ROUTINE CODE #######################################################

       ! assert size equality and get approximation points
       n = assert_eq(size(z), size(pi,1), size(pi,2), 'discretize_AR')

       ! calculate variance of the overall process
       sigma_eta = sigma_eps/(1d0-rho**2)

       ! determine the transition matrix
       call rouwenhorst_matrix(rho, pi)

       ! determine the nodes
       psi = sqrt(dble(n-1))*sqrt(sigma_eta)
       do in = 1, n
           z(in) = -psi + 2d0*psi*dble(in-1)/dble(n-1)
       enddo
       z = z + mu

       if(present(w))then
           w = 1d0/dble(n)
           do in = 1, 10000
               w = matmul(transpose(pi), w)
           enddo
       endif

       !##########################################################################
       ! Subroutines and functions
       !##########################################################################

       contains


       !##########################################################################
       ! subroutine rouwenhorst_matrix
       !
       ! Calculates value of function that should be integrated for pis.
       !##########################################################################
       recursive subroutine rouwenhorst_matrix(rho, pi_new)

           use Constants_mod, only: IK, RK; implicit none
           real(RK), intent(in) :: rho
           real(RK), intent(out) :: pi_new(:, :)
           integer(IK) :: n
           real(RK) :: p, pi_old(size(pi_new,1)-1, size(pi_new,1)-1)

           n = size(pi_new, 1)
           p = (1d0 + rho)/2d0

           if(n == 2)then
               pi_new(1, :) = (/p, 1d0-p/)
               pi_new(2, :) = (/1d0-p, p/)
           else
               call rouwenhorst_matrix(rho, pi_old)
               pi_new = 0d0

               pi_new(1:n-1, 1:n-1) = pi_new(1:n-1, 1:n-1) + p*pi_old
               pi_new(1:n-1, 2:n  ) = pi_new(1:n-1, 2:n  ) + (1d0-p)*pi_old
               pi_new(2:n  , 1:n-1) = pi_new(2:n  , 1:n-1) + (1d0-p)*pi_old
               pi_new(2:n  , 2:n  ) = pi_new(2:n  , 2:n  ) + p*pi_old

               pi_new(2:n-1, :) = pi_new(2:n-1, :)/2d0
           endif
       end subroutine

    end subroutine discretize_AR


    !##############################################################################
    ! SUBROUTINE discretize_log_AR
    !
    ! Discretizes a log-AR(1) process using the Rouwenhorst method.
    !
    ! REFERENCE: Kopecky, K.A., Suen, R.M.H., Finite state Markov-chain
    !            approximations to highly persistent processes, Review of Economic
    !            Dynamics, Vol. 13, No. 3, 2010, 701-714.
    !##############################################################################
    subroutine discretize_log_AR(rho, mu, sigma_eps, z, pi, w)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! autoregression parameter
        real(RK), intent(in) :: rho

        ! unconditional mean of the process
        real(RK), intent(in) :: mu

        ! variance of the shock
        real(RK), intent(in) :: sigma_eps

        ! discrete shock values
        real(RK), intent(out) :: z(:)

        ! transition matrix
        real(RK), intent(out) :: pi(:, :)

        ! the stationary distribution
        real(RK), intent(out), optional :: w(:)


        !##### OTHER VARIABLES ####################################################

        real(RK) :: sigma_eta, mu_c, sigma_c


        !##### ROUTINE CODE #######################################################

        ! calculate variance of the overall process
        sigma_eta = sigma_eps/(1d0-rho**2)

        ! get the transformed variance and expectation
        sigma_c = log(1d0+sigma_eta/mu**2)
        mu_c  = log(mu)-0.5d0*sigma_c
        sigma_c = sigma_c*(1d0-rho**2)

        ! discretize the log distribution
        if(present(w))then
            call discretize_AR(rho, mu_c, sigma_c, z, pi, w)
        else
            call discretize_AR(rho, mu_c, sigma_c, z, pi)
        endif

        ! take exponentials
        z = exp(z)

    end subroutine discretize_log_AR


    !##############################################################################
    ! SUBROUTINE simulate_AR
    !
    ! Simulates a discrete AR(1) process.
    !##############################################################################
    subroutine simulate_AR(pi, shocks, fixed)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! transition matrix
        real(RK), intent(in) :: pi(:, :)

        ! simulated schocks
        integer(IK), intent(out) :: shocks(:)

        ! should the random seed be initialized at a fixed values
        logical, optional :: fixed


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: T, n, j


        !##### ROUTINE CODE #######################################################

        ! assert size equality and get number of simulated schocks
        n = assert_eq(size(pi,1), size(pi,2), 'tauchen')
        T = size(shocks)

        ! initialize the random seed
        if(tbox_seed)then
            if(present(fixed))then
                call init_random_seed(fixed)
            else
                call init_random_seed()
            endif
            tbox_seed = .false.
        endif

        ! get first entry
        shocks(1) = n/2+1

        ! now calculate other shocks
        do j = 2, T
            shocks(j) = get_tomorrow(pi(shocks(j-1), :))
        enddo


    !##########################################################################
    ! Subroutines and functions
    !##########################################################################

    contains


        !##########################################################################
        ! FUNCTION get_tomorrow
        !
        ! Calculates value of function that should be integrated for pis.
        !##########################################################################
        function get_tomorrow(pi)

            use Constants_mod, only: IK, RK; implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! transition probabilities
            real(RK), intent(in) :: pi(:)

            ! tomorrows shock
            integer(IK) :: get_tomorrow


            !##### OTHER VARIABLES ################################################

            real(RK) :: rand
            integer(IK) :: i1


            !##### ROUTINE CODE ###################################################

            ! get random number
            call random_number(rand)

            ! get tomorrows value
            do i1 = 1, size(pi, 1)-1

                if(rand <= sum(pi(1:i1), 1))then
                    get_tomorrow = i1
                    return
                endif
            enddo

            ! else choose last value
            get_tomorrow = i1
            return

        end function

    end subroutine simulate_AR












!##############################################################################
!##############################################################################
! MODULE gnuplot
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE execplot
    !
    ! Actually creates the plot files.
    !##############################################################################
    subroutine execplot(xlim, xticks, xlabel, ylim, yticks, ylabel, title, &
            legend, filename, filetype, output)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! x axis definitial
        real(RK), optional :: xlim(2)

        ! x axis tick definitions
        real(RK), optional :: xticks

        ! y axis definitial
        real(RK), optional :: ylim(2)

        ! y axis tick definitions
        real(RK), optional :: yticks

        ! output file name
        character(LEN=*), optional :: title

        ! output file name
        character(LEN=*), optional :: xlabel

        ! output file name
        character(LEN=*), optional :: ylabel

        ! legend position
        character(LEN=2), optional :: legend

        ! output file name
        character(LEN=*), optional :: filename

        ! file type
        character(LEN=*), optional :: filetype

        ! output file name
        character(LEN=*), optional :: output

        !##### OTHER VARIABLES ####################################################

        integer(IK) :: i1, i2
        character(LEN=3) :: ft
        character(LEN=150) :: cfile, dfile


        !##### ROUTINE CODE #######################################################

        if(present(output))then
            cfile = output//'_c.dat'
            dfile = output//'_d.dat'
        else
            cfile = 'command13545431.dat'
            dfile = 'plotdata13545431.dat'
        endif

        ! write the output data file
        open(213659,file=trim(dfile))
        do i1 = 1, gnu_nmax
            write(213659, '(2000e21.10e4)')(gnu_x(i1, i2), gnu_y(i1, i2), i2=1, gnu_mmax)
        enddo
        close(213659)

        ! write the command file
        open(213659,file=trim(cfile))

        ! set terminal and grid
        write(213659, '(a)')'if (strstrt(GPVAL_TERMINALS, "wxt") > 0) {'
        write(213659, '(a)')'    set terminal wxt title "Gnuplot"'
        write(213659, '(a)')'} else {'
        write(213659, '(a)')'    if (strstrt(GPVAL_TERMINALS, "x11") > 0) {'
        write(213659, '(a)')'        set terminal x11'
        write(213659, '(a)')'    } else {'
        write(213659, '(a)')'        if (strstrt(GPVAL_TERMINALS, "qt") > 0) {'
        write(213659, '(a)')'            set terminal qt'
        write(213659, '(a)')'        } else {'
        write(213659, '(a)')'            print ""'
        write(213659, '(a)')'            print "ATTENTION: There seems to be NO VALID TERMINAL installed in "'
        write(213659, '(a)')'            print "           your version of gnuplot. It might be a good idea to "'
        write(213659, '(a)')'            print "           completely uninstall your Fortran/gnuplot/geany installation"'
        write(213659, '(a)')'            print "           using the uninstallation files on www.ce-fortran.com. "'
        write(213659, '(a)')'            print "           Afterwards you should reinstall the Fortran system again."'
        write(213659, '(a)')'            print "           If this does not help, please refer to the forum."'
        write(213659, '(a)')'            print ""'
        write(213659, '(a)')'            q'
        write(213659, '(a)')'        }'
        write(213659, '(a)')'    }'
        write(213659, '(a)')'}'
        write(213659, '(a)')'set grid'
        if(gnu_histogram)then
            write(213659, '(a)')'set style data histograms'
            write(213659, '(a)')'set style fill solid border -1'
        endif

        ! set x axis
        if(present(xlim))write(213659, '(a,e13.5,a,e13.5,a)')'set xrange [',minval(xlim),':',maxval(xlim),']'
        if(present(xticks))write(213659, '(a,e13.5)')'set xtics ',xticks
        if(present(xlabel))write(213659, '(a)')'set xlabel "'//xlabel//'"font ",12"'

        ! set y axis
        if(present(ylim) .and. .not. gnu_histogram) &
            write(213659, '(a,e13.5,a,e13.5,a)')'set yrange [',minval(ylim),':',maxval(ylim),']'
        if(gnu_histogram)then
            if(present(ylim))then
                write(213659, '(a,e13.5,a,e13.5,a)') &
                    'set yrange [',min(minval(ylim),0d0),':',maxval(ylim),']'
            else
                write(213659, '(a)') 'set yrange [0:]'
            endif
        endif

        if(present(yticks))write(213659, '(a,e13.5)')'set ytics ',yticks
        if(present(ylabel))write(213659, '(a)')'set ylabel "'//ylabel//'"font ",12"'

        ! set title
        if(present(title))write(213659, '(a)')'set title "'//title//'" font ",16"'

        ! legend statement
        if(gnu_dolegend)then
            if(present(legend))then
                select case (legend(1:2))
                    case("ln")
                        write(213659, '(a)')'set key inside left top'
                    case("ls")
                        write(213659, '(a)')'set key inside left bottom'
                    case("lo")
                        write(213659, '(a)')'set key outside vert center left'
                    case("lb")
                        write(213659, '(a)')'set key outside left below'
                    case("rn")
                        write(213659, '(a)')'set key inside right top'
                    case("rs")
                        write(213659, '(a)')'set key inside right bottom'
                    case("ro")
                        write(213659, '(a)')'set key outside vert center right'
                    case("rb")
                        write(213659, '(a)')'set key outside right below'
                    case("cn")
                        write(213659, '(a)')'set key inside center top'
                    case("cs")
                        write(213659, '(a)')'set key inside center bottom'
                    case("co")
                        write(213659, '(a)')'set key outside vert top center'
                    case("cb")
                        write(213659, '(a)')'set key outside center below'
                    case default
                        write(213659, '(a)')'set key outside vert bottom center'
                end select
            else
                write(213659, '(a)')'set key center below'
            endif
        else
            write(213659, '(a)')'unset key'
        endif

        ! write plot lines
        if(gnu_mmax == 1)then
            write(213659, '(a)')'plot "'//trim(dfile)//'" using 1:2 '//trim(gnu_definitions(1))
        else
            write(213659, '(a)')'plot "'//trim(dfile)//'" using 1:2 '//trim(gnu_definitions(1))//',\'
            do i1 = 2, gnu_mmax-1
                write(213659, '(a,i4,a,i4,a)')'     "'//trim(dfile)//'" using '&
                    ,2*(i1-1)+1,':',2*(i1-1)+2,' '//trim(gnu_definitions(i1))//',\'
            enddo
            write(213659, '(a,i4,a,i4,a)')'     "'//trim(dfile)//'" using ',  &
                2*(gnu_mmax-1)+1,':',2*(gnu_mmax-1)+2,' '//trim(gnu_definitions(i1))
        endif

        write(213659, '(a)')'pause -1 "Press RETURN to continue..."'

        ! write graph to file
        if(present(filename))then
            ft = "eps"
            if(present(filetype))then
                if(filetype(1:3) == "png")ft= "png"
            endif
            write(213659, *)
            write(213659, '(a)')'set terminal '//ft
            write(213659, '(a)')'set output "'//filename//'.'//ft//'"'
            write(213659, '(a)')'replot'
        endif

        write(213659, '(a)')'q'
        close(213659)

        call system('gnuplot "'//trim(cfile)//'"')
        if(.not.present(output))then
            open(213659,file=trim(dfile))
            close(213659, status='delete')
            open(213659,file=trim(cfile))
            close(213659, status='delete')
        endif

        gnu_addtoplot = .false.
        gnu_dolegend = .false.
        gnu_histogram = .false.

    end subroutine


    !##############################################################################
    ! SUBROUTINE plot
    !
    ! Plots a x-y-data column to the output file.
    !##############################################################################
    subroutine plot(xin, yin, color, linewidth, marker, markersize, noline, legend)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        real(RK), intent(in) :: xin(:), yin(:)
        character(LEN=*), optional :: color
        real(RK), optional :: linewidth
        integer(IK), optional :: marker
        real(RK), optional :: markersize
        logical, optional :: noline
        character(LEN=*), optional :: legend


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, i1
        logical :: lines, points


        !##### ROUTINE CODE #######################################################

        n = assert_eq(size(xin, 1), size(yin, 1), 'plot')

        ! generate the respective plot data
        if(.not.gnu_addtoplot)then

            ! allocate new x array
            if(allocated(gnu_x))deallocate(gnu_x)
            allocate(gnu_x(n, 1))
            gnu_x(:, 1) = xin

            ! allocate new y array
            if(allocated(gnu_y))deallocate(gnu_y)
            allocate(gnu_y(n, 1))
            gnu_y(:, 1) = yin

            ! allocate new x_temp array
            if(allocated(gnu_x_temp))deallocate(gnu_x_temp)
            allocate(gnu_x_temp(n, 1))
            gnu_x_temp(:, 1) = xin

            ! allocate new y_temp array
            if(allocated(gnu_y_temp))deallocate(gnu_y_temp)
            allocate(gnu_y_temp(n, 1))
            gnu_y_temp(:, 1) = yin

            gnu_addtoplot = .true.
            gnu_nmax = n
            gnu_mmax = 1

        else

            ! get new number of lines
            gnu_mmax = gnu_mmax+1
            if(gnu_mmax > 1000)then
                write(*,'(/a/)')'SORRY: I CANNOT PLOT MORE THAN 1000 LINES'
                return
            endif

            ! deallocate arrays
            deallocate(gnu_x)
            deallocate(gnu_y)

            ! if the new array is the longer one
            if(n > gnu_nmax)then
                allocate(gnu_x(n, gnu_mmax))
                allocate(gnu_y(n, gnu_mmax))

                gnu_x(1:gnu_nmax, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(1:gnu_nmax, 1:gnu_mmax-1) = gnu_y_temp

                ! fill up with the same values
                do i1 = gnu_nmax+1, n
                    gnu_x(i1, 1:gnu_mmax-1) = gnu_x(gnu_nmax, 1:gnu_mmax-1)
                    gnu_y(i1, 1:gnu_mmax-1) = gnu_y(gnu_nmax, 1:gnu_mmax-1)
                enddo

                gnu_x(:, gnu_mmax) = xin
                gnu_y(:, gnu_mmax) = yin

                ! set new nmax
                gnu_nmax = n

            ! if the old array is the longer one
            else

                allocate(gnu_x(gnu_nmax, gnu_mmax))
                allocate(gnu_y(gnu_nmax, gnu_mmax))

                gnu_x(:, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(:, 1:gnu_mmax-1) = gnu_y_temp

                gnu_x(1:n, gnu_mmax) = xin
                gnu_y(1:n, gnu_mmax) = yin

                ! fill up with same values
                do i1 = n+1, gnu_nmax
                    gnu_x(i1, gnu_mmax) = gnu_x(n, gnu_mmax)
                    gnu_y(i1, gnu_mmax) = gnu_y(n, gnu_mmax)
                enddo
            endif

            deallocate(gnu_x_temp)
            allocate(gnu_x_temp(size(gnu_x,1), size(gnu_x,2)))
            gnu_x_temp = gnu_x

            deallocate(gnu_y_temp)
            allocate(gnu_y_temp(size(gnu_y,1), size(gnu_y,2)))
            gnu_y_temp = gnu_y

        endif

        ! check for lines and points
        lines = .true.
        if(present(noline))then
            if(noline)lines = .false.
        endif

        if(.not.lines)then
            points = .true.
        else
            points = .false.
        endif
        if(present(marker))points = .true.

        ! set up definitions
        gnu_definitions(gnu_mmax) = 'with'

        ! get lines and points
        if(lines .and. points)then
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' linespoints'
        elseif(lines)then
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' lines'
        else
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' points'
        endif

        ! get the line color
        if(present(color))then
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' lc rgb "'//adjustl(trim(color))//'"'
        else
            if(gnu_mmax == 1)then
                gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' lc rgb "blue"'
            else
                write(gnu_definitions(gnu_mmax), '(a,i4)')trim(gnu_definitions(gnu_mmax))//' lc ', gnu_mmax-1
            endif
        endif

        ! get line width
        if(lines)then
            if(present(linewidth)) then
                write(gnu_definitions(gnu_mmax), '(a,f8.2)')trim(gnu_definitions(gnu_mmax))//' lw ', max(linewidth, 0d0)
            else
                write(gnu_definitions(gnu_mmax), '(a,f8.2)')trim(gnu_definitions(gnu_mmax))//' lw ', 2d0
            endif
        endif

        ! get marker definition
        if(points)then

            if(present(marker)) then
                write(gnu_definitions(gnu_mmax), '(a,i2)')trim(gnu_definitions(gnu_mmax))//' pt ', min(marker, 13)
            else
                write(gnu_definitions(gnu_mmax), '(a,i2)')trim(gnu_definitions(gnu_mmax))//' pt ', 1
            endif

            if(present(markersize))then
                write(gnu_definitions(gnu_mmax), '(a,f8.2)')trim(gnu_definitions(gnu_mmax))//' ps ', max(markersize, 0d0)
            else
                write(gnu_definitions(gnu_mmax), '(a,f8.2)')trim(gnu_definitions(gnu_mmax))//' ps ', 1d0
            endif
        endif

        ! set the legend
        if(present(legend))then
            gnu_dolegend = .true.
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' title "'//adjustl(trim(legend))//'"'
        else
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' notitle '
        endif

    end subroutine


    !##############################################################################
    ! SUBROUTINE plot_hist_old
    !
    ! Plots a histogram.
    !##############################################################################
    subroutine plot_hist_old(xin, yin, color, legend)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        real(RK), intent(in) :: xin(:), yin(:)
        character(LEN=*), optional :: color
        character(LEN=*), optional :: legend


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, i1


        !##### ROUTINE CODE #######################################################

        n = assert_eq(size(xin, 1), size(yin, 1), 'plot')

        ! generate the respective plot data
        if(.not.gnu_addtoplot)then

            ! allocate new x array
            if(allocated(gnu_x))deallocate(gnu_x)
            allocate(gnu_x(n, 1))
            gnu_x(:, 1) = xin

            ! allocate new y array
            if(allocated(gnu_y))deallocate(gnu_y)
            allocate(gnu_y(n, 1))
            gnu_y(:, 1) = yin

            ! allocate new x_temp array
            if(allocated(gnu_x_temp))deallocate(gnu_x_temp)
            allocate(gnu_x_temp(n, 1))
            gnu_x_temp(:, 1) = xin

            ! allocate new y_temp array
            if(allocated(gnu_y_temp))deallocate(gnu_y_temp)
            allocate(gnu_y_temp(n, 1))
            gnu_y_temp(:, 1) = yin

            gnu_addtoplot = .true.
            gnu_nmax = n
            gnu_mmax = 1

        else

            ! get new number of lines
            gnu_mmax = gnu_mmax+1
            if(gnu_mmax > 1000)then
                write(*,'(/a/)')'SORRY: I CANNOT PLOT MORE THAN 1000 LINES'
                return
            endif

            ! deallocate arrays
            deallocate(gnu_x)
            deallocate(gnu_y)

            ! if the new array is the longer one
            if(n > gnu_nmax)then
                allocate(gnu_x(n, gnu_mmax))
                allocate(gnu_y(n, gnu_mmax))

                gnu_x(1:gnu_nmax, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(1:gnu_nmax, 1:gnu_mmax-1) = gnu_y_temp

                ! fill up with the same values
                do i1 = gnu_nmax+1, n
                    gnu_x(i1, 1:gnu_mmax-1) = gnu_x(gnu_nmax, 1:gnu_mmax-1)
                    gnu_y(i1, 1:gnu_mmax-1) = gnu_y(gnu_nmax, 1:gnu_mmax-1)
                enddo

                gnu_x(:, gnu_mmax) = xin
                gnu_y(:, gnu_mmax) = yin

                ! set new nmax
                gnu_nmax = n

            ! if the old array is the longer one
            else

                allocate(gnu_x(gnu_nmax, gnu_mmax))
                allocate(gnu_y(gnu_nmax, gnu_mmax))

                gnu_x(:, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(:, 1:gnu_mmax-1) = gnu_y_temp

                gnu_x(1:n, gnu_mmax) = xin
                gnu_y(1:n, gnu_mmax) = yin

                ! fill up with same values
                do i1 = n+1, gnu_nmax
                    gnu_x(i1, gnu_mmax) = gnu_x(n, gnu_mmax)
                    gnu_y(i1, gnu_mmax) = gnu_y(n, gnu_mmax)
                enddo
            endif

            deallocate(gnu_x_temp)
            allocate(gnu_x_temp(size(gnu_x,1), size(gnu_x,2)))
            gnu_x_temp = gnu_x

            deallocate(gnu_y_temp)
            allocate(gnu_y_temp(size(gnu_y,1), size(gnu_y,2)))
            gnu_y_temp = gnu_y

        endif

        ! set up definitions
        gnu_definitions(gnu_mmax) = 'with boxes'

        ! get the line color
        if(present(color))then
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' lc rgb "'//adjustl(trim(color))//'"'
        else
            write(gnu_definitions(gnu_mmax), '(a,i4)')trim(gnu_definitions(gnu_mmax))//' lc ', gnu_mmax
        endif

        ! set the legend
        if(present(legend))then
            gnu_dolegend = .true.
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' title "'//adjustl(trim(legend))//'"'
        else
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' notitle '
        endif

        gnu_histogram = .true.

    end subroutine


    !##############################################################################
    ! SUBROUTINE plot_hist_new
    !
    ! Plots a histogram.
    !##############################################################################
    subroutine plot_hist_new(xvalues, nbins, left, right, incl_outside, absolute, color, legend)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        real(RK), intent(in) :: xvalues(1:)
        integer(IK), intent(in) :: nbins
        real(RK), optional :: left
        real(RK), optional :: right
        logical, optional :: incl_outside
        logical, optional :: absolute
        character(LEN=*), optional :: color
        character(LEN=*), optional :: legend


        !##### OTHER VARIABLES ####################################################

        integer(IK) :: n, i1
        real(RK) :: lower, upper, total
        real(RK) :: thresholds(0:nbins), freq(0:nbins)
        real(RK) :: xin(nbins), yin(nbins)



        !##### ROUTINE CODE #######################################################

        n = nbins

        ! get lower and upper endpoint as well as width
        lower = minval(xvalues)
        if(present(left))lower = left

        upper = maxval(xvalues)
        if(present(right))upper = right

        ! get thresholds
        call grid_Cons_Equi(thresholds, lower, upper)

        ! get frequency
        freq(0) = dble(count(xvalues < lower))
        do i1 = 1, nbins
            freq(i1) = dble(count(xvalues <= thresholds(i1)))
        enddo
        total = freq(nbins) - freq(0)

        if(present(incl_outside))then
            if(incl_outside)then
                freq(0) = 0d0
                freq(nbins) = freq(nbins) + dble(count(xvalues >= upper))
                total = freq(nbins)
            endif
        endif

        ! get histogram data
        do i1 = 1, nbins
            xin(i1) = (thresholds(i1-1) + thresholds(i1))/2d0
            yin(i1) = (freq(i1) - freq(i1-1))
        enddo

        if(present(absolute))then
            if(.not. absolute)then
                yin = yin/total
            endif
        else
            yin = yin/total
        endif

        ! generate the respective plot data
        if(.not.gnu_addtoplot)then

            ! allocate new x array
            if(allocated(gnu_x))deallocate(gnu_x)
            allocate(gnu_x(n, 1))
            gnu_x(:, 1) = xin

            ! allocate new y array
            if(allocated(gnu_y))deallocate(gnu_y)
            allocate(gnu_y(n, 1))
            gnu_y(:, 1) = yin

            ! allocate new x_temp array
            if(allocated(gnu_x_temp))deallocate(gnu_x_temp)
            allocate(gnu_x_temp(n, 1))
            gnu_x_temp(:, 1) = xin

            ! allocate new y_temp array
            if(allocated(gnu_y_temp))deallocate(gnu_y_temp)
            allocate(gnu_y_temp(n, 1))
            gnu_y_temp(:, 1) = yin

            gnu_addtoplot = .true.
            gnu_nmax = n
            gnu_mmax = 1

        else

            ! get new number of lines
            gnu_mmax = gnu_mmax+1
            if(gnu_mmax > 1000)then
                write(*,'(/a/)')'SORRY: I CANNOT PLOT MORE THAN 1000 LINES'
                return
            endif

            ! deallocate arrays
            deallocate(gnu_x)
            deallocate(gnu_y)

            ! if the new array is the longer one
            if(n > gnu_nmax)then
                allocate(gnu_x(n, gnu_mmax))
                allocate(gnu_y(n, gnu_mmax))

                gnu_x(1:gnu_nmax, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(1:gnu_nmax, 1:gnu_mmax-1) = gnu_y_temp

                ! fill up with the same values
                do i1 = gnu_nmax+1, n
                    gnu_x(i1, 1:gnu_mmax-1) = gnu_x(gnu_nmax, 1:gnu_mmax-1)
                    gnu_y(i1, 1:gnu_mmax-1) = gnu_y(gnu_nmax, 1:gnu_mmax-1)
                enddo

                gnu_x(:, gnu_mmax) = xin
                gnu_y(:, gnu_mmax) = yin

                ! set new nmax
                gnu_nmax = n

            ! if the old array is the longer one
            else

                allocate(gnu_x(gnu_nmax, gnu_mmax))
                allocate(gnu_y(gnu_nmax, gnu_mmax))

                gnu_x(:, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(:, 1:gnu_mmax-1) = gnu_y_temp

                gnu_x(1:n, gnu_mmax) = xin
                gnu_y(1:n, gnu_mmax) = yin

                ! fill up with same values
                do i1 = n+1, gnu_nmax
                    gnu_x(i1, gnu_mmax) = gnu_x(n, gnu_mmax)
                    gnu_y(i1, gnu_mmax) = gnu_y(n, gnu_mmax)
                enddo
            endif

            deallocate(gnu_x_temp)
            allocate(gnu_x_temp(size(gnu_x,1), size(gnu_x,2)))
            gnu_x_temp = gnu_x

            deallocate(gnu_y_temp)
            allocate(gnu_y_temp(size(gnu_y,1), size(gnu_y,2)))
            gnu_y_temp = gnu_y

        endif

        ! set up definitions
        gnu_definitions(gnu_mmax) = 'with boxes'

        ! get the line color
        if(present(color))then
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' lc rgb "'//adjustl(trim(color))//'"'
        else
            write(gnu_definitions(gnu_mmax), '(a,i4)')trim(gnu_definitions(gnu_mmax))//' lc ', gnu_mmax
        endif

        ! set the legend
        if(present(legend))then
            gnu_dolegend = .true.
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' title "'//adjustl(trim(legend))//'"'
        else
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' notitle '
        endif

        gnu_histogram = .true.

    end subroutine


    !##############################################################################
    ! SUBROUTINE plot3d_grid
    !
    ! Plots a x-y-z-data column on a grid. Plot will be executed immediately.
    !
    ! Credits to Patrick Wiesmann for providing a first version of a 3d plotting
    !     subroutine during his research assistantship.
    !##############################################################################
    subroutine plot3d_grid(xin, yin, zin, color, linewidth, marker, markersize, noline, &
                xlim, xticks, xlabel, ylim, yticks, ylabel, zlim, zticks, zlevel, zlabel, &
                surf, surf_color, transparent, view, title, filename, filetype, output)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! input variables
        real(RK), intent(in) :: xin(:), yin(:), zin(:,:)

        ! line and marker color
        character(LEN=*), optional :: color

        ! width of the plotting line
        real(RK), optional :: linewidth

        ! marker definition
        integer(IK), optional :: marker

        ! size of the marker
        real(RK), optional :: markersize

        ! if no line should be plotted
        logical, optional :: noline

        ! x axis definitial
        real(RK), optional :: xlim(2)

        ! x axis tick definitions
        real(RK), optional :: xticks

        ! y axis definitial
        real(RK), optional :: ylim(2)

        ! y axis tick definitions
        real(RK), optional :: yticks

           ! z axis definitial
        real(RK), optional :: zlim(2)

        ! z axis tick definitions
        real(RK), optional :: zticks

        ! point where the z axis is placed
        real(RK), optional :: zlevel

        ! rotates the graph
        real(RK), optional :: view(2)

        ! colored surface or not
        logical, optional :: surf

        ! color of the surface
        integer(IK), optional :: surf_color

        ! see through surface or not
        logical, optional :: transparent

        ! output file name
        character(LEN=*), optional :: xlabel

        ! output file name
        character(LEN=*), optional :: ylabel

        ! output file name
        character(LEN=*), optional :: zlabel

        ! output file name
        character(LEN=*), optional :: title

        ! output file name
        character(LEN=*), optional :: filename

        ! file type
        character(LEN=*), optional :: filetype

        ! output file name
        character(LEN=*), optional :: output


        !##### OTHER VARIABLES ####################################################

        logical :: lines, points
        character(LEN = 2000) :: definition
        integer(IK) :: i1, i2, n1, n2
        character(LEN=3) :: ft
        character(LEN=150) :: cfile, dfile


        !##### ROUTINE CODE #######################################################

        ! assert that the sizes are coorect
        n1 = assert_eq(size(xin, 1), size(zin, 1), 'plot3d')
        n2 = assert_eq(size(yin, 1), size(zin, 2), 'plot3d')

        ! check for lines and points
        lines = .true.
        if(present(noline))then
            if(noline)lines = .false.
        endif

        if(.not.lines)then
            points = .true.
        else
            points = .false.
        endif
        if(present(marker))points = .true.

        ! set up definitions
        definition = 'with'

        ! get lines and points
        if(lines .and. points)then
            definition = trim(definition)//' linespoints'
        elseif(lines)then
            definition = trim(definition)//' lines'
        else
            definition = trim(definition)//' points'
        endif

        ! get the line color
        if(present(color))then
            definition = trim(definition)//' lc rgb "'//adjustl(trim(color))//'"'
        else
            definition = trim(definition)//' lc rgb "blue"'
        endif

        ! get line width
        if(lines)then
            if(present(linewidth)) then
                write(definition, '(a,f8.2)')trim(definition)//' lw ', max(linewidth, 0d0)
            else
                if(present(surf))then
                    write(definition, '(a,f8.2)')trim(definition)//' lw ', 0.3d0
                else
                    write(definition, '(a,f8.2)')trim(definition)//' lw ', 1d0
                endif
            endif
        endif

        ! get marker definition
        if(points)then

            if(present(marker)) then
                write(definition, '(a,i2)')trim(definition)//' pt ', min(marker, 13)
            else
                write(definition, '(a,i2)')trim(definition)//' pt ', 1
            endif

            if(present(markersize))then
                write(definition, '(a,f8.2)')trim(definition)//' ps ', max(markersize, 0d0)
            else
                write(definition, '(a,f8.2)')trim(definition)//' ps ', 1d0
            endif
        endif

        ! set the legend
        definition = trim(definition)//' notitle '

        ! now directly create output
        if(present(output))then
            cfile = output//'_c.dat'
            dfile = output//'_d.dat'
        else
            cfile = 'command135454313d.dat'
            dfile = 'plotdata135454313d.dat'
        endif

        ! write the output data file
        open(213659,file=trim(dfile))
        do i1 = 1, n1
            do i2 = 1, n2
                write(213659, '(3e21.10e4)')xin(i1), yin(i2), zin(i1, i2)
             enddo
            write(213659,*)
        enddo
        close(213659)

        ! write the command file
        open(213659,file=trim(cfile))

        ! set terminal and grid
        write(213659, '(a)')'if (strstrt(GPVAL_TERMINALS, "wxt") > 0) {'
        write(213659, '(a)')'    set terminal wxt title "Gnuplot"'
        write(213659, '(a)')'} else {'
        write(213659, '(a)')'    if (strstrt(GPVAL_TERMINALS, "x11") > 0) {'
        write(213659, '(a)')'        set terminal x11'
        write(213659, '(a)')'    } else {'
        write(213659, '(a)')'        if (strstrt(GPVAL_TERMINALS, "qt") > 0) {'
        write(213659, '(a)')'            set terminal qt'
        write(213659, '(a)')'        } else {'
        write(213659, '(a)')'            print ""'
        write(213659, '(a)')'            print "ATTENTION: There seems to be NO VALID TERMINAL installed in "'
        write(213659, '(a)')'            print "           your version of gnuplot. It might be a good idea to "'
        write(213659, '(a)')'            print "           completely uninstall your Fortran/gnuplot/geany installation"'
        write(213659, '(a)')'            print "           using the uninstallation files on www.ce-fortran.com. "'
        write(213659, '(a)')'            print "           Afterwards you should reinstall the Fortran system again."'
        write(213659, '(a)')'            print "           If this does not help, please refer to the forum."'
        write(213659, '(a)')'            print ""'
        write(213659, '(a)')'            q'
        write(213659, '(a)')'        }'
        write(213659, '(a)')'    }'
        write(213659, '(a)')'}'
        write(213659, '(a)')'set grid'

        if(present(zlevel)) then
            write(213659, '(a,e13.5)')'set ticslevel ', -min(max(zlevel, 0d0), 1d0)
        else
            write(213659, '(a,e13.5)')'set ticslevel 0'
        endif

        if(present(surf) .and. surf) then
            write(213659, '(a)')'set pm3d'
            if(present(surf_color)) then
                select case (surf_color)
                    case(1)
                         write(213659,'(a)')'set palette gray'
                    case(2)
                        write(213659,'(a)')'set palette rgbformulae 33,13,10'
                    case(3)
                        write(213659,'(a)')'set palette rgbformulae 3,11,6'
                    case(4)
                        write(213659,'(a)')'set palette rgbformulae 23,28,3'
                    case(5)
                        write(213659, '(a)')'set palette rgbformulae 21,22,23'
                    case default
                        write(213659, '(a)')'set palette rgbformulae 7,5,15'
                    end select
            endif
        endif

        if(present(transparent) .and. .not. transparent) then
            write(213659, '(a)')'set hidden3d'
        endif

        if(present(view)) then
            write(213659, '(a,e13.5,a,e13.5)')'set view ', min(max(view(1), 0d0), 360d0), &
                ',',min(max(view(2), 0d0), 360d0)
        endif

        ! set x axis
        if(present(xlim))write(213659, '(a,e13.5,a,e13.5,a)')'set xrange [',minval(xlim),':',maxval(xlim),']'
        if(present(xticks))write(213659, '(a,e13.5)')'set xtics ',xticks
        if(present(xlabel))write(213659, '(a)')'set xlabel "'//xlabel//'"font ",12"'

        ! set y axis
        if(present(ylim))write(213659, '(a,e13.5,a,e13.5,a)')'set yrange [',minval(ylim),':',maxval(ylim),']'
        if(present(yticks))write(213659, '(a,e13.5)')'set ytics ',yticks
        if(present(ylabel))write(213659, '(a)')'set ylabel "'//ylabel//'"font ",12"'

        ! set z axis
        if(present(zlim))write(213659, '(a,e13.5,a,e13.5,a)')'set zrange [',minval(zlim),':',maxval(zlim),']'
        if(present(zticks))write(213659, '(a,e13.5)')'set ztics ',zticks
        if(present(zlabel))write(213659, '(a)')'set zlabel "'//zlabel//'"font ",12"'

        ! set title
        if(present(title))write(213659, '(a)')'set title "'//title//'" font ",16"'

        ! legend statement
        write(213659, '(a)')'unset key'

        ! write plot lines for 3D-data
        write(213659, '(a)')'splot "'//trim(dfile)//'" using 1:2:3 '//trim(definition)

        write(213659, '(a)')'pause -1 "Press RETURN to continue..."'

        ! write graph to file
        if(present(filename))then
            ft = "eps"
            if(present(filetype))then
                if(filetype(1:3) == "png")ft= "png"
            endif
            write(213659, *)
            write(213659, '(a)')'set terminal '//ft
            write(213659, '(a)')'set output "'//filename//'.'//ft//'"'
            write(213659, '(a)')'replot'
        endif

        write(213659, '(a)')'q'
        close(213659)

        call system('gnuplot "'//trim(cfile)//'"')
        if(.not.present(output))then
            open(213659,file=trim(dfile))
            close(213659, status='delete')
            open(213659,file=trim(cfile))
            close(213659, status='delete')
        endif

    end subroutine


    !##############################################################################
    ! SUBROUTINE plot3d_line
    !
    ! Plots a x-y-z-data column with random points. Plot will be executed immediately.
    !
    ! Credits to Patrick Wiesmann for providing a first version of a 3d plotting
    !     subroutine during his research assistantship.
    !##############################################################################
    subroutine plot3d_line(xin, yin, zin, color, linewidth, marker, markersize, noline, &
                xlim, xticks, xlabel, ylim, yticks, ylabel, zlim, zticks, zlevel, zlabel, &
                view, title, filename, filetype, output)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! input variables
        real(RK), intent(in) :: xin(:), yin(:), zin(:)

        ! line and marker color
        character(LEN=*), optional :: color

        ! width of the plotting line
        real(RK), optional :: linewidth

        ! marker definition
        integer(IK), optional :: marker

        ! size of the marker
        real(RK), optional :: markersize

        ! if no line should be plotted
        logical, optional :: noline

        ! x axis definitial
        real(RK), optional :: xlim(2)

        ! x axis tick definitions
        real(RK), optional :: xticks

        ! y axis definitial
        real(RK), optional :: ylim(2)

        ! y axis tick definitions
        real(RK), optional :: yticks

           ! z axis definitial
        real(RK), optional :: zlim(2)

        ! z axis tick definitions
        real(RK), optional :: zticks

        ! point where the z axis is placed
        real(RK), optional :: zlevel

        ! rotates the graph
        real(RK), optional :: view(2)

        ! output file name
        character(LEN=*), optional :: xlabel

        ! output file name
        character(LEN=*), optional :: ylabel

        ! output file name
        character(LEN=*), optional :: zlabel

        ! output file name
        character(LEN=*), optional :: title

        ! output file name
        character(LEN=*), optional :: filename

        ! file type
        character(LEN=*), optional :: filetype

        ! output file name
        character(LEN=*), optional :: output


        !##### OTHER VARIABLES ####################################################

        logical :: lines, points
        character(LEN = 2000) :: definition
        integer(IK) :: i1, n
        character(LEN=3) :: ft
        character(LEN=150) :: cfile, dfile


        !##### ROUTINE CODE #######################################################

        ! assert that the sizes are coorect
        n = assert_eq(size(xin, 1), size(yin, 1), size(zin, 1), 'plot3d')

        ! check for lines and points
        lines = .true.
        if(present(noline))then
            if(noline)lines = .false.
        endif

        if(.not.lines)then
            points = .true.
        else
            points = .false.
        endif
        if(present(marker))points = .true.

        ! set up definitions
        definition = 'with'

        ! get lines and points
        if(lines .and. points)then
            definition = trim(definition)//' linespoints'
        elseif(lines)then
            definition = trim(definition)//' lines'
        else
            definition = trim(definition)//' points'
        endif

        ! get the line color
        if(present(color))then
            definition = trim(definition)//' lc rgb "'//adjustl(trim(color))//'"'
        else
            definition = trim(definition)//' lc rgb "blue"'
        endif

        ! get line width
        if(lines)then
            if(present(linewidth)) then
                write(definition, '(a,f8.2)')trim(definition)//' lw ', max(linewidth, 0d0)
            else
                write(definition, '(a,f8.2)')trim(definition)//' lw ', 1d0
            endif
        endif

        ! get marker definition
        if(points)then

            if(present(marker)) then
                write(definition, '(a,i2)')trim(definition)//' pt ', min(marker, 13)
            else
                write(definition, '(a,i2)')trim(definition)//' pt ', 1
            endif

            if(present(markersize))then
                write(definition, '(a,f8.2)')trim(definition)//' ps ', max(markersize, 0d0)
            else
                write(definition, '(a,f8.2)')trim(definition)//' ps ', 1d0
            endif
        endif

        ! set the legend
        definition = trim(definition)//' notitle '

        ! now directly create output
        if(present(output))then
            cfile = output//'_c.dat'
            dfile = output//'_d.dat'
        else
            cfile = 'command135454313d.dat'
            dfile = 'plotdata135454313d.dat'
        endif

        ! write the output data file
        open(213659,file=trim(dfile))
        do i1 = 1, n
            write(213659, '(3e21.10e4)')xin(i1), yin(i1), zin(i1)
        enddo
        close(213659)

        ! write the command file
        open(213659,file=trim(cfile))

        ! set terminal and grid
        write(213659, '(a)')'if (strstrt(GPVAL_TERMINALS, "wxt") > 0) {'
        write(213659, '(a)')'    set terminal wxt title "Gnuplot"'
        write(213659, '(a)')'} else {'
        write(213659, '(a)')'    if (strstrt(GPVAL_TERMINALS, "x11") > 0) {'
        write(213659, '(a)')'        set terminal x11'
        write(213659, '(a)')'    } else {'
        write(213659, '(a)')'        if (strstrt(GPVAL_TERMINALS, "qt") > 0) {'
        write(213659, '(a)')'            set terminal qt'
        write(213659, '(a)')'        } else {'
        write(213659, '(a)')'            print ""'
        write(213659, '(a)')'            print "ATTENTION: There seems to be NO VALID TERMINAL installed in "'
        write(213659, '(a)')'            print "           your version of gnuplot. It might be a good idea to "'
        write(213659, '(a)')'            print "           completely uninstall your Fortran/gnuplot/geany installation"'
        write(213659, '(a)')'            print "           using the uninstallation files on www.ce-fortran.com. "'
        write(213659, '(a)')'            print "           Afterwards you should reinstall the Fortran system again."'
        write(213659, '(a)')'            print "           If this does not help, please refer to the forum."'
        write(213659, '(a)')'            print ""'
        write(213659, '(a)')'            q'
        write(213659, '(a)')'        }'
        write(213659, '(a)')'    }'
        write(213659, '(a)')'}'
        write(213659, '(a)')'set grid'

        if(present(zlevel)) then
            write(213659, '(a,e13.5)')'set ticslevel ', -min(max(zlevel, 0d0), 1d0)
        else
            write(213659, '(a,e13.5)')'set ticslevel 0'
        endif

        if(present(view)) then
            write(213659, '(a,e13.5,a,e13.5)')'set view ', min(max(view(1), 0d0), 360d0), &
                ',',min(max(view(2), 0d0), 360d0)
        endif

        ! set x axis
        if(present(xlim))write(213659, '(a,e13.5,a,e13.5,a)')'set xrange [',minval(xlim),':',maxval(xlim),']'
        if(present(xticks))write(213659, '(a,e13.5)')'set xtics ',xticks
        if(present(xlabel))write(213659, '(a)')'set xlabel "'//xlabel//'"font ",12"'

        ! set y axis
        if(present(ylim))write(213659, '(a,e13.5,a,e13.5,a)')'set yrange [',minval(ylim),':',maxval(ylim),']'
        if(present(yticks))write(213659, '(a,e13.5)')'set ytics ',yticks
        if(present(ylabel))write(213659, '(a)')'set ylabel "'//ylabel//'"font ",12"'

        ! set z axis
        if(present(zlim))write(213659, '(a,e13.5,a,e13.5,a)')'set zrange [',minval(zlim),':',maxval(zlim),']'
        if(present(zticks))write(213659, '(a,e13.5)')'set ztics ',zticks
        if(present(zlabel))write(213659, '(a)')'set zlabel "'//zlabel//'"font ",12"'

        ! set title
        if(present(title))write(213659, '(a)')'set title "'//title//'" font ",16"'

        ! legend statement
        write(213659, '(a)')'unset key'

        ! write plot lines for 3D-data
        write(213659, '(a)')'splot "'//trim(dfile)//'" using 1:2:3 '//trim(definition)

        write(213659, '(a)')'pause -1 "Press RETURN to continue..."'

        ! write graph to file
        if(present(filename))then
            ft = "eps"
            if(present(filetype))then
                if(filetype(1:3) == "png")ft= "png"
            endif
            write(213659, *)
            write(213659, '(a)')'set terminal '//ft
            write(213659, '(a)')'set output "'//filename//'.'//ft//'"'
            write(213659, '(a)')'replot'
        endif

        write(213659, '(a)')'q'
        close(213659)

        call system('gnuplot "'//trim(cfile)//'"')
        if(.not.present(output))then
            open(213659,file=trim(dfile))
            close(213659, status='delete')
            open(213659,file=trim(cfile))
            close(213659, status='delete')
        endif

    end subroutine











!##############################################################################
!##############################################################################
! MODULE sorting
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE sort_r
    !
    ! Sorts an array of type real(RK) in ascending order.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     The Wikibook "Algorithm Implementation" available at
    !     https://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Quicksort#FORTRAN_90.2F95
    !
    !     and follows closely the Qsort implementation found in
    !     "A FORTRAN 90 Numerical Library" (AFNL), which is available at
    !     https://sourceforge.net/projects/afnl/
    !
    !     REFERENCE: Ramos, A. (2006). A FORTRAN 90 numerical library.
    !##############################################################################
    subroutine sort_r(x)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the array that should be sorted
        real(RK), intent(inout) :: x(:)


        !##### OTHER VARIABLES ####################################################

        ! from which list size should insertion sort be used
        integer(IK), parameter :: Isw = 10

        ! type for the left and right bounds
        type Limits
           integer(IK) :: Ileft, Iright
        end type Limits

        ! other variables
        integer(IK) :: Ipvn, Ileft, Iright, ISpos, ISmax
        type(Limits), allocatable :: Stack(:)


        !##### ROUTINE CODE #######################################################

        allocate(Stack(Size(X)*2))

        Stack(:)%Ileft = 0

        ! Iniitialize the stack
        Ispos = 1
        Ismax = 1
        Stack(ISpos)%Ileft  = 1
        Stack(ISpos)%Iright = size(X)

        do while (Stack(ISpos)%Ileft /= 0)

           Ileft = Stack(ISPos)%Ileft
           Iright = Stack(ISPos)%Iright

           ! choose between inseration and quick sort
           if (Iright-Ileft <= Isw) then
              call InsrtLC(X, Ileft, Iright)
              ISpos = ISPos + 1
           else
              Ipvn = ChoosePiv(X, Ileft, Iright)
              Ipvn = Partition(X, Ileft, Iright, Ipvn)

              Stack(ISmax+1)%Ileft = Ileft
              Stack(ISmax+1) %Iright = Ipvn-1
              Stack(ISmax+2)%Ileft = Ipvn + 1
              Stack(ISmax+2)%Iright = Iright
              ISpos = ISpos + 1
              ISmax = ISmax + 2
           endif
        enddo

        ! deallocate all arrays
        deallocate(Stack)


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! FUNCTION ChoosePiv
        !
        ! Determines the pivotal element for the quicksort algorithm.
        !##########################################################################
        function ChoosePiv(XX, IIleft, IIright) result (IIpv)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            real(RK), intent(in) :: XX(:)

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIleft, IIright

            ! the pivotal element that is returned
            integer(IK) :: IIpv


            !##### OTHER VARIABLES ####################################################

            real(RK) :: XXcp(3)
            integer(IK) :: IIpt(3), IImd


            !##### ROUTINE CODE #######################################################

            IImd = Int((IIleft+IIright)/2)
            IIpv = IImd
            XXcp(1) = XX(IIleft)
            XXcp(2) = XX(IImd)
            XXcp(3) = XX(IIright)
            IIpt = (/1,2,3/)

            call InsrtLC(XXcp, 1, 3, IIpt)

            select case (IIpt(2))
            case (1)
                IIpv = IIleft
            case (2)
                IIpv = IImd
            case (3)
                IIpv = IIright
            End Select

        end function


        !##########################################################################
        ! SUBROUTINE InsrtLC
        !
        ! Perform an insertion sort of the list XX(:) between index
        !     values IIl and IIr.
        !##########################################################################
        subroutine InsrtLC(XX, IIl, IIr, IIpt)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            real(RK), intent(inout) :: XX(:)

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIl, IIr

            ! the permutations
            integer(IK), intent(inout), optional :: IIpt(:)


            !##### OTHER VARIABLES ####################################################

            real(RK) :: RRtmp
            integer(IK) :: II, JJ


            !##### ROUTINE CODE #######################################################

            do II = IIl+1, IIr
                RRtmp = XX(II)
                do JJ = II-1, 1, -1
                    if (RRtmp < XX(JJ)) then
                        XX(JJ+1) = XX(JJ)
                        if(present(IIpt))call Swap_IN(IIpt, JJ, JJ+1)
                    else
                        Exit
                    endif
                enddo
                XX(JJ+1) = RRtmp
            enddo

        end subroutine InsrtLC


        !##########################################################################
        ! FUNCTION Partition
        !
        ! Arranges the array X between the index values Ileft and Iright
        !     positioning elements smallers than X(Ipv) at the left and the others
        !     at the right.
        !##########################################################################
        function Partition(X, Ileft, Iright, Ipv) result(Ipvfn)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            real(RK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: Ileft, Iright, Ipv

            ! return value
            integer(IK) :: Ipvfn


            !##### OTHER VARIABLES ####################################################

            real(RK) :: Rpv
            integer(IK) :: I


            !##### ROUTINE CODE #######################################################

            Rpv = X(Ipv)
            call Swap(X, Ipv, Iright)
            Ipvfn = Ileft

            do I = Ileft, Iright-1
                if (X(I) <= Rpv) then
                    call Swap(X, I, Ipvfn)
                    Ipvfn = Ipvfn + 1
                endif
            enddo

            call Swap(X, Ipvfn, Iright)

        end function Partition


        !##########################################################################
        ! SUBROUTINE Swap
        !
        ! Swaps elements i and j of array x
        !##########################################################################
        subroutine Swap(X, I, J)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            real(RK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: I, J


            !##### OTHER VARIABLES ####################################################

            real(RK) :: Xtmp


            !##### ROUTINE CODE #######################################################

            Xtmp = X(I)
            X(I) = X(J)
            X(J) = Xtmp

        end subroutine Swap


        !##########################################################################
        ! SUBROUTINE Swap_IN
        !
        ! Swaps elements i and j of an integer(IK) array
        !##########################################################################
        subroutine Swap_IN(X, I, J)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: I, J


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: Xtmp


            !##### ROUTINE CODE #######################################################

            Xtmp = X(I)
            X(I) = X(J)
            X(J) = Xtmp

        end subroutine Swap_IN

     end subroutine sort_r



    !##############################################################################
    ! SUBROUTINE sort_r2
    !
    ! Sorts an array of type real(RK) in ascending order and returns new ordering.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     The Wikibook "Algorithm Implementation" available at
    !     https://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Quicksort#FORTRAN_90.2F95
    !
    !     and follows closely the Qsort implementation found in
    !     "A FORTRAN 90 Numerical Library" (AFNL), which is available at
    !     https://sourceforge.net/projects/afnl/
    !
    !     REFERENCE: Ramos, A. (2006). A FORTRAN 90 numerical library.
    !##############################################################################
    subroutine sort_r2(x, iorder)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the array that should be sorted
        real(RK), intent(inout) :: x(:)

        ! an array that will contain the new sorting order
        integer(IK), intent(out) :: iorder(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        ! from which list size should insertion sort be used
        integer(IK), parameter :: Isw = 10

        ! type for the left and right bounds
        type Limits
           integer(IK) :: Ileft, Iright
        end type Limits

        ! other variables
        integer(IK) :: Ipvn, Ileft, Iright, ISpos, ISmax, ii
        type(Limits), allocatable :: Stack(:)


        !##### ROUTINE CODE #######################################################

        ! initialize the sorting order array
        iorder = (/(ii, ii = 1, size(x, 1))/)

        allocate(Stack(Size(X)*2))

        Stack(:)%Ileft = 0

        ! Iniitialize the stack
        Ispos = 1
        Ismax = 1
        Stack(ISpos)%Ileft  = 1
        Stack(ISpos)%Iright = size(X)

        do while (Stack(ISpos)%Ileft /= 0)

           Ileft = Stack(ISPos)%Ileft
           Iright = Stack(ISPos)%Iright

           ! choose between inseration and quick sort
           if (Iright-Ileft <= Isw) then
              call InsrtLC(X, iorder, Ileft, Iright)
              ISpos = ISPos + 1
           else
              Ipvn = ChoosePiv(X, Ileft, Iright)
              Ipvn = Partition(X, iorder, Ileft, Iright, Ipvn)

              Stack(ISmax+1)%Ileft = Ileft
              Stack(ISmax+1) %Iright = Ipvn-1
              Stack(ISmax+2)%Ileft = Ipvn + 1
              Stack(ISmax+2)%Iright = Iright
              ISpos = ISpos + 1
              ISmax = ISmax + 2
           endif
        enddo

        ! deallocate all arrays
        deallocate(Stack)


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! FUNCTION ChoosePiv
        !
        ! Determines the pivotal element for the quicksort algorithm.
        !##########################################################################
        function ChoosePiv(XX, IIleft, IIright) result (IIpv)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            real(RK), intent(in) :: XX(:)

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIleft, IIright

            ! the pivotal element that is returned
            integer(IK) :: IIpv


            !##### OTHER VARIABLES ####################################################

            real(RK) :: XXcp(3)
            integer(IK) :: IIpt(3), IImd


            !##### ROUTINE CODE #######################################################

            IImd = Int((IIleft+IIright)/2)
            IIpv = IImd
            XXcp(1) = XX(IIleft)
            XXcp(2) = XX(IImd)
            XXcp(3) = XX(IIright)
            IIpt = (/1,2,3/)

            call InsrtLC_help(XXcp, 1, 3, IIpt)

            select case (IIpt(2))
            case (1)
                IIpv = IIleft
            case (2)
                IIpv = IImd
            case (3)
                IIpv = IIright
            End Select

        end function

        !##########################################################################
        ! SUBROUTINE InsrtLC_help
        !
        ! Just a helping routine. Same as below but without iorder.
        !##########################################################################
        subroutine InsrtLC_help(XX, IIl, IIr, IIpt)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            real(RK), intent(inout) :: XX(:)

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIl, IIr

            ! the permutations
            integer(IK), intent(inout) :: IIpt(:)


            !##### OTHER VARIABLES ####################################################

            real(RK) :: RRtmp
            integer(IK) :: II, JJ


            !##### ROUTINE CODE #######################################################

            do II = IIl+1, IIr
                RRtmp = XX(II)
                do JJ = II-1, 1, -1
                    if (RRtmp < XX(JJ)) then
                        XX(JJ+1) = XX(JJ)
                        call Swap_IN(IIpt, JJ, JJ+1)
                    else
                        Exit
                    endif
                enddo
                XX(JJ+1) = RRtmp
            enddo

        end subroutine InsrtLC_help


        !##########################################################################
        ! SUBROUTINE InsrtLC
        !
        ! Perform an insertion sort of the list XX(:) between index
        !     values IIl and IIr.
        !##########################################################################
        subroutine InsrtLC(XX, iorder, IIl, IIr)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            real(RK), intent(inout) :: XX(:)

            ! an array that will contain the new sorting order
            integer(IK), intent(inout) :: iorder(size(XX, 1))

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIl, IIr


            !##### OTHER VARIABLES ####################################################

            real(RK) :: RRtmp
            integer(IK) :: IItmp
            integer(IK) :: II, JJ


            !##### ROUTINE CODE #######################################################

            do II = IIl+1, IIr
                RRtmp = XX(II)
                IItmp = iorder(II)
                do JJ = II-1, 1, -1
                    if (RRtmp < XX(JJ)) then
                        XX(JJ+1) = XX(JJ)
                        iorder(JJ+1) = iorder(JJ)
                    else
                        Exit
                    endif
                enddo
                XX(JJ+1) = RRtmp
                iorder(JJ+1) = IItmp
            enddo

        end subroutine InsrtLC


        !##########################################################################
        ! FUNCTION Partition
        !
        ! Arranges the array X between the index values Ileft and Iright
        !     positioning elements smallers than X(Ipv) at the left and the others
        !     at the right.
        !##########################################################################
        function Partition(X, iorder, Ileft, Iright, Ipv) result(Ipvfn)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            real(RK), intent(inout) :: X(:)

            ! an array that will contain the new sorting order
            integer(IK), intent(inout) :: iorder(size(X, 1))

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: Ileft, Iright, Ipv

            ! return value
            integer(IK) :: Ipvfn


            !##### OTHER VARIABLES ####################################################

            real(RK) :: Rpv
            integer(IK) :: I


            !##### ROUTINE CODE #######################################################

            Rpv = X(Ipv)
            call Swap(X, Ipv, Iright)
            call Swap_IN(iorder, Ipv, Iright)
            Ipvfn = Ileft

            do I = Ileft, Iright-1
                if (X(I) <= Rpv) then
                    call Swap(X, I, Ipvfn)
                    call Swap_IN(iorder, I, Ipvfn)
                    Ipvfn = Ipvfn + 1
                endif
            enddo

            call Swap(X, Ipvfn, Iright)
            call Swap_IN(iorder, Ipvfn, Iright)

        end function Partition


        !##########################################################################
        ! SUBROUTINE Swap
        !
        ! Swaps elements i and j of array x
        !##########################################################################
        subroutine Swap(X, I, J)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            real(RK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: I, J


            !##### OTHER VARIABLES ####################################################

            real(RK) :: Xtmp


            !##### ROUTINE CODE #######################################################

            Xtmp = X(I)
            X(I) = X(J)
            X(J) = Xtmp

        end subroutine Swap


        !##########################################################################
        ! SUBROUTINE Swap_IN
        !
        ! Swaps elements i and j of an integer(IK) array
        !##########################################################################
        subroutine Swap_IN(X, I, J)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: I, J


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: Xtmp


            !##### ROUTINE CODE #######################################################

            Xtmp = X(I)
            X(I) = X(J)
            X(J) = Xtmp

        end subroutine Swap_IN

    end subroutine sort_r2




    !##############################################################################
    ! SUBROUTINE sort_i
    !
    ! Sorts an array of type integer(IK) in ascending order.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     The Wikibook "Algorithm Implementation" available at
    !     https://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Quicksort#FORTRAN_90.2F95
    !
    !     and follows closely the Qsort implementation found in
    !     "A FORTRAN 90 Numerical Library" (AFNL), which is available at
    !     https://sourceforge.net/projects/afnl/
    !
    !     REFERENCE: Ramos, A. (2006). A FORTRAN 90 numerical library.
    !##############################################################################
    subroutine sort_i(x)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the array that should be sorted
        integer(IK), intent(inout) :: x(:)


        !##### OTHER VARIABLES ####################################################

        ! from which list size should insertion sort be used
        integer(IK), parameter :: Isw = 10

        ! type for the left and right bounds
        type Limits
           integer(IK) :: Ileft, Iright
        end type Limits

        ! other variables
        integer(IK) :: Ipvn, Ileft, Iright, ISpos, ISmax
        type(Limits), allocatable :: Stack(:)


        !##### ROUTINE CODE #######################################################

        allocate(Stack(Size(X)*2))

        Stack(:)%Ileft = 0

        ! Iniitialize the stack
        Ispos = 1
        Ismax = 1
        Stack(ISpos)%Ileft  = 1
        Stack(ISpos)%Iright = size(X)

        do while (Stack(ISpos)%Ileft /= 0)

           Ileft = Stack(ISPos)%Ileft
           Iright = Stack(ISPos)%Iright

           ! choose between inseration and quick sort
           if (Iright-Ileft <= Isw) then
              call InsrtLC(X, Ileft, Iright)
              ISpos = ISPos + 1
           else
              Ipvn = ChoosePiv(X, Ileft, Iright)
              Ipvn = Partition(X, Ileft, Iright, Ipvn)

              Stack(ISmax+1)%Ileft = Ileft
              Stack(ISmax+1) %Iright = Ipvn-1
              Stack(ISmax+2)%Ileft = Ipvn + 1
              Stack(ISmax+2)%Iright = Iright
              ISpos = ISpos + 1
              ISmax = ISmax + 2
           endif
        enddo

        ! deallocate all arrays
        deallocate(Stack)


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! FUNCTION ChoosePiv
        !
        ! Determines the pivotal element for the quicksort algorithm.
        !##########################################################################
        function ChoosePiv(XX, IIleft, IIright) result (IIpv)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(in) :: XX(:)

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIleft, IIright

            ! the pivotal element that is returned
            integer(IK) :: IIpv


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: XXcp(3)
            integer(IK) :: IIpt(3), IImd


            !##### ROUTINE CODE #######################################################

            IImd = Int((IIleft+IIright)/2)
            IIpv = IImd
            XXcp(1) = XX(IIleft)
            XXcp(2) = XX(IImd)
            XXcp(3) = XX(IIright)
            IIpt = (/1,2,3/)

            call InsrtLC(XXcp, 1, 3, IIpt)

            select case (IIpt(2))
            case (1)
                IIpv = IIleft
            case (2)
                IIpv = IImd
            case (3)
                IIpv = IIright
            End Select

        end function


        !##########################################################################
        ! SUBROUTINE InsrtLC
        !
        ! Perform an insertion sort of the list XX(:) between index
        !     values IIl and IIr.
        !##########################################################################
        subroutine InsrtLC(XX, IIl, IIr, IIpt)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: XX(:)

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIl, IIr

            ! the permutations
            integer(IK), intent(inout), optional :: IIpt(:)


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: RRtmp
            integer(IK) :: II, JJ


            !##### ROUTINE CODE #######################################################

            do II = IIl+1, IIr
                RRtmp = XX(II)
                do JJ = II-1, 1, -1
                    if (RRtmp < XX(JJ)) then
                        XX(JJ+1) = XX(JJ)
                        if(present(IIpt))call Swap_IN(IIpt, JJ, JJ+1)
                    else
                        Exit
                    endif
                enddo
                XX(JJ+1) = RRtmp
            enddo

        end subroutine InsrtLC


        !##########################################################################
        ! FUNCTION Partition
        !
        ! Arranges the array X between the index values Ileft and Iright
        !     positioning elements smallers than X(Ipv) at the left and the others
        !     at the right.
        !##########################################################################
        function Partition(X, Ileft, Iright, Ipv) result(Ipvfn)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: Ileft, Iright, Ipv

            ! return value
            integer(IK) :: Ipvfn


            !##### OTHER VARIABLES ####################################################

            real(RK) :: Rpv
            integer(IK) :: I


            !##### ROUTINE CODE #######################################################

            Rpv = X(Ipv)
            call Swap(X, Ipv, Iright)
            Ipvfn = Ileft

            do I = Ileft, Iright-1
                if (X(I) <= Rpv) then
                    call Swap(X, I, Ipvfn)
                    Ipvfn = Ipvfn + 1
                endif
            enddo

            call Swap(X, Ipvfn, Iright)

        end function Partition


        !##########################################################################
        ! SUBROUTINE Swap
        !
        ! Swaps elements i and j of array x
        !##########################################################################
        subroutine Swap(X, I, J)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: I, J


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: Xtmp


            !##### ROUTINE CODE #######################################################

            Xtmp = X(I)
            X(I) = X(J)
            X(J) = Xtmp

        end subroutine Swap


        !##########################################################################
        ! SUBROUTINE Swap_IN
        !
        ! Swaps elements i and j of an integer(IK) array
        !##########################################################################
        subroutine Swap_IN(X, I, J)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: I, J


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: Xtmp


            !##### ROUTINE CODE #######################################################

            Xtmp = X(I)
            X(I) = X(J)
            X(J) = Xtmp

        end subroutine Swap_IN

     end subroutine sort_i



    !##############################################################################
    ! SUBROUTINE sort_i2
    !
    ! Sorts an array of type integer(IK) in ascending order and returns new ordering.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     The Wikibook "Algorithm Implementation" available at
    !     https://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Quicksort#FORTRAN_90.2F95
    !
    !     and follows closely the Qsort implementation found in
    !     "A FORTRAN 90 Numerical Library" (AFNL), which is available at
    !     https://sourceforge.net/projects/afnl/
    !
    !     REFERENCE: Ramos, A. (2006). A FORTRAN 90 numerical library.
    !##############################################################################
    subroutine sort_i2(x, iorder)

        use Constants_mod, only: IK, RK; implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the array that should be sorted
        integer(IK), intent(inout) :: x(:)

        ! an array that will contain the new sorting order
        integer(IK), intent(out) :: iorder(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        ! from which list size should insertion sort be used
        integer(IK), parameter :: Isw = 10

        ! type for the left and right bounds
        type Limits
           integer(IK) :: Ileft, Iright
        end type Limits

        ! other variables
        integer(IK) :: Ipvn, Ileft, Iright, ISpos, ISmax, ii
        type(Limits), allocatable :: Stack(:)


        !##### ROUTINE CODE #######################################################

        ! initialize the sorting order array
        iorder = (/(ii, ii = 1, size(x, 1))/)

        allocate(Stack(Size(X)*2))

        Stack(:)%Ileft = 0

        ! Iniitialize the stack
        Ispos = 1
        Ismax = 1
        Stack(ISpos)%Ileft  = 1
        Stack(ISpos)%Iright = size(X)

        do while (Stack(ISpos)%Ileft /= 0)

           Ileft = Stack(ISPos)%Ileft
           Iright = Stack(ISPos)%Iright

           ! choose between inseration and quick sort
           if (Iright-Ileft <= Isw) then
              call InsrtLC(X, iorder, Ileft, Iright)
              ISpos = ISPos + 1
           else
              Ipvn = ChoosePiv(X, Ileft, Iright)
              Ipvn = Partition(X, iorder, Ileft, Iright, Ipvn)

              Stack(ISmax+1)%Ileft = Ileft
              Stack(ISmax+1) %Iright = Ipvn-1
              Stack(ISmax+2)%Ileft = Ipvn + 1
              Stack(ISmax+2)%Iright = Iright
              ISpos = ISpos + 1
              ISmax = ISmax + 2
           endif
        enddo

        ! deallocate all arrays
        deallocate(Stack)


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! FUNCTION ChoosePiv
        !
        ! Determines the pivotal element for the quicksort algorithm.
        !##########################################################################
        function ChoosePiv(XX, IIleft, IIright) result (IIpv)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(in) :: XX(:)

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIleft, IIright

            ! the pivotal element that is returned
            integer(IK) :: IIpv


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: XXcp(3)
            integer(IK) :: IIpt(3), IImd


            !##### ROUTINE CODE #######################################################

            IImd = Int((IIleft+IIright)/2)
            IIpv = IImd
            XXcp(1) = XX(IIleft)
            XXcp(2) = XX(IImd)
            XXcp(3) = XX(IIright)
            IIpt = (/1,2,3/)

            call InsrtLC_help(XXcp, 1, 3, IIpt)

            select case (IIpt(2))
            case (1)
                IIpv = IIleft
            case (2)
                IIpv = IImd
            case (3)
                IIpv = IIright
            End Select

        end function

        !##########################################################################
        ! SUBROUTINE InsrtLC_help
        !
        ! Just a helping routine. Same as below but without iorder.
        !##########################################################################
        subroutine InsrtLC_help(XX, IIl, IIr, IIpt)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: XX(:)

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIl, IIr

            ! the permutations
            integer(IK), intent(inout) :: IIpt(:)


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: RRtmp
            integer(IK) :: II, JJ


            !##### ROUTINE CODE #######################################################

            do II = IIl+1, IIr
                RRtmp = XX(II)
                do JJ = II-1, 1, -1
                    if (RRtmp < XX(JJ)) then
                        XX(JJ+1) = XX(JJ)
                        call Swap_IN(IIpt, JJ, JJ+1)
                    else
                        Exit
                    endif
                enddo
                XX(JJ+1) = RRtmp
            enddo

        end subroutine InsrtLC_help


        !##########################################################################
        ! SUBROUTINE InsrtLC
        !
        ! Perform an insertion sort of the list XX(:) between index
        !     values IIl and IIr.
        !##########################################################################
        subroutine InsrtLC(XX, iorder, IIl, IIr)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: XX(:)

            ! an array that will contain the new sorting order
            integer(IK), intent(inout) :: iorder(size(XX, 1))

            ! the left and right definition of the sub-array
            integer(IK), intent(in) :: IIl, IIr


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: RRtmp
            integer(IK) :: IItmp
            integer(IK) :: II, JJ


            !##### ROUTINE CODE #######################################################

            do II = IIl+1, IIr
                RRtmp = XX(II)
                IItmp = iorder(II)
                do JJ = II-1, 1, -1
                    if (RRtmp < XX(JJ)) then
                        XX(JJ+1) = XX(JJ)
                        iorder(JJ+1) = iorder(JJ)
                    else
                        Exit
                    endif
                enddo
                XX(JJ+1) = RRtmp
                iorder(JJ+1) = IItmp
            enddo

        end subroutine InsrtLC


        !##########################################################################
        ! FUNCTION Partition
        !
        ! Arranges the array X between the index values Ileft and Iright
        !     positioning elements smallers than X(Ipv) at the left and the others
        !     at the right.
        !##########################################################################
        function Partition(X, iorder, Ileft, Iright, Ipv) result(Ipvfn)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: X(:)

            ! an array that will contain the new sorting order
            integer(IK), intent(inout) :: iorder(size(X, 1))

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: Ileft, Iright, Ipv

            ! return value
            integer(IK) :: Ipvfn


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: Rpv
            integer(IK) :: I


            !##### ROUTINE CODE #######################################################

            Rpv = X(Ipv)
            call Swap(X, Ipv, Iright)
            call Swap_IN(iorder, Ipv, Iright)
            Ipvfn = Ileft

            do I = Ileft, Iright-1
                if (X(I) <= Rpv) then
                    call Swap(X, I, Ipvfn)
                    call Swap_IN(iorder, I, Ipvfn)
                    Ipvfn = Ipvfn + 1
                endif
            enddo

            call Swap(X, Ipvfn, Iright)
            call Swap_IN(iorder, Ipvfn, Iright)

        end function Partition


        !##########################################################################
        ! SUBROUTINE Swap
        !
        ! Swaps elements i and j of array x
        !##########################################################################
        subroutine Swap(X, I, J)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: I, J


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: Xtmp


            !##### ROUTINE CODE #######################################################

            Xtmp = X(I)
            X(I) = X(J)
            X(J) = Xtmp

        end subroutine Swap


        !##########################################################################
        ! SUBROUTINE Swap_IN
        !
        ! Swaps elements i and j of an integer(IK) array
        !##########################################################################
        subroutine Swap_IN(X, I, J)


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! the array to work on
            integer(IK), intent(inout) :: X(:)

            ! the left and right definition of the sub-array and pivotal element
            integer(IK), intent(in) :: I, J


            !##### OTHER VARIABLES ####################################################

            integer(IK) :: Xtmp


            !##### ROUTINE CODE #######################################################

            Xtmp = X(I)
            X(I) = X(J)
            X(J) = Xtmp

        end subroutine Swap_IN

    end subroutine sort_i2


end module EconomicsToolbox_mod
