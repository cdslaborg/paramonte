! Coarray 3D Navier-Stokes Solver Test
!
! Copyright (c) 2012-2014, Sourcery, Inc.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!     * Neither the name of Sourcery, Inc., nor the
!       names of its contributors may be used to endorse or promote products
!       derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL SOURCERY, INC., BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!

!(*----------------------------------------------------------------------------------------------------------------------
!         basic in-core shear code ( 7 words/node, not threaded, in-core, no file read/write )
!------------------------------------------------------------------------------------------------------------------------*)

! Define universal constants:
! In the case of exactly representable numbers, the definitions are useful
! to ensure subprogram argument type/kind/rank matching without having to
! repeat kind specifiers everywhere.
module constants_module
  use iso_fortran_env, only : int64
  implicit none
  private
  public :: one,zero
  integer(int64), parameter :: one=1_int64,zero=0_int64
end module

! Initialize the random seed with a varying seed to ensure a different
! random number sequence for each invocation of subroutine, e.g. for
! invocations on different images of a coarray parallel program.
! Setting any seed values to zero is deprecated because it can result
! in low-quality random number sequences.
! (Source: https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html)
module random_module
  implicit none
  private
  public :: init_random_seed
contains
  subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t

            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               if (this_image()==1) print *,"OS provides random number generator"
               read(un) seed
               close(un)
            else
               if (this_image()==1) print *,"OS does not provide random number generator"
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
          end subroutine init_random_seed
end module random_module

module run_size
    use iso_fortran_env, only : int64,real64 ! 64-bit integer and real kind parameters
    use constants_module, only : one         ! 64-bit unit to ensure argument kind match
    implicit none
        real, codimension[*] ::  viscos, shear, b11, b22, b33, b12, velmax
        integer(int64), codimension[*] ::  nx, ny, nz, nsteps, output_step
        integer(int64), codimension[*] :: my, mx, first_y, last_y, first_x, last_x
        real(real64), codimension[*] :: cpu_time, tran_time, sync_time, total_time
        real(real64), codimension[*] :: max_cpu_time, max_tran_time, max_sync_time, max_total_time
        real(real64), codimension[*] :: min_cpu_time, min_tran_time, min_sync_time, min_total_time

        real ::  time, cfl, dt
        integer(int64) :: my_node, num_nodes
        real, parameter :: pi = 3.141592653589793

contains

    subroutine max_velmax()
        integer(int64) :: i

        sync all
        if( my_node == 1) then
            do i = 2, num_nodes;     velmax = max( velmax, velmax[i] );    end do
        end if
        sync all
        if (my_node>1) velmax = velmax[1]
        sync all
    end subroutine max_velmax

    subroutine global_times()
        integer(int64) :: i, stage

        max_cpu_time = cpu_time
        max_tran_time = tran_time
        max_total_time = sync_time
        max_total_time = total_time
        min_cpu_time = cpu_time
        min_tran_time = tran_time
        min_total_time = sync_time
        min_total_time = total_time

        do stage = 1, num_nodes-1
            i = 1 + mod( my_node-1+stage, num_nodes )
            max_cpu_time = max( max_cpu_time, cpu_time[i] )
            min_cpu_time = min( min_cpu_time, cpu_time[i] )
            max_tran_time = max( max_tran_time, tran_time[i] )
            min_tran_time = min( min_tran_time, tran_time[i] )
            max_sync_time = max( max_sync_time, sync_time[i] )
            min_sync_time = min( min_sync_time, sync_time[i] )
            max_total_time = max( max_total_time, total_time[i] )
            min_total_time = min( min_total_time, total_time[i] )
        end do
        sync all
    end subroutine global_times

subroutine copy3( A,B, n1, sA1, sB1, n2, sA2, sB2, n3, sA3, sB3 )
  implicit none
  complex, intent(in)  :: A(0:*)
  complex, intent(out) :: B(0:*)
  integer(int64), intent(in) :: n1, sA1, sB1
  integer(int64), intent(in) :: n2, sA2, sB2
  integer(int64), intent(in) :: n3, sA3, sB3
  integer(int64) i,j,k

  do k=0,n3-1
     do j=0,n2-1
        do i=0,n1-1
           B(i*sB1+j*sB2+k*sB3) = A(i*sA1+j*sA2+k*sA3)
        end do
     end do
  end do
end subroutine copy3

end module run_size

program cshear

  !(***********************************************************************************************************
  !                   m a i n   p r o g r a m
  !***********************************************************************************************************)
      use iso_fortran_env, only : int64,real64 ! 64-bit integer and real kind parameters
      use run_size
      implicit none

      interface
        subroutine solve_navier_stokes
        end subroutine solve_navier_stokes
      end interface

nx=128;ny=128;nz=128
viscos=0.; shear=0.
b11=1.; b22=1.; b33=1.; b12 =0.
nsteps=10;output_step=1

      num_nodes = num_images()
      my_node = this_image()

      if( my_node == 1 ) then

!       write(6,*) "nx,ny,nz : ";               read(5,*) nx, ny, nz
        if ( mod(nx/2,num_nodes) /= 0) then;    write(6,*) "nx/2 not multiple of num_nodes";    stop;   end if
        if ( mod(ny,num_nodes) /= 0) then;      write(6,*) " ny not multiple of num_nodes";     stop;   end if

!       write(6,*) "viscos, shear : ";          read(5,*) viscos, shear
!       write(6,*) "b11 b22 b33 b12 : ";        read(5,*) b11, b22, b33, b12
!       write(6,*) "nsteps, output_step : ";    read(5,*) nsteps, output_step

        write(6,fmt="(3(A,i4))") "nx =",nx,   "   ny =",ny,   "   nz =",nz
        write(6,fmt="(2(A,f7.3))") "viscos = ", viscos, "      shear = ", shear
        write(6,fmt="(A,4f7.3)") "b11 b22 b33 b12 = ", b11, b22, b33, b12
        write(6,fmt="(2(A,i6))") "nsteps = ", nsteps, "       output_step = ", output_step

        write(6,fmt="(A,i4,A)") "----------------- running on ", num_nodes, " images -------------------"

      end if

      sync all  !--- images > 1 wait on inputs from image = 1 !

      if( my_node > 1 ) then
        nx = nx[1];     ny = ny[1];     nz = nz[1]
        viscos = viscos[1];     shear = shear[1]
        b11 = b11[1];   b22 = b22[1];   b33 = b33[1];   b12 = b12[1]
        nsteps = nsteps[1];     output_step = output_step[1]
      end if

 	    mx = nx/2 / num_nodes;  first_x = (my_node-1)*mx + 1;   last_x  = (my_node-1)*mx + mx
        my = ny / num_nodes;    first_y = (my_node-1)*my + 1;   last_y  = (my_node-1)*my + my

      if(my_node == 1 ) write(6,fmt="(A, f6.2)") "message size (MB) = ", real(nz*4*mx*my*8)/real(1024*1024)

    call solve_navier_stokes

end program cshear

!  (***********************************************************************************************************
!             n a v i e r - s t o k e s   s o l v e r
!  ************************************************************************************************************)

  subroutine solve_navier_stokes
      use run_size
      implicit none

  !(*****************************   declarations     ****************************************)

       integer(int64) ::  stop, rflag, oflag, step, rkstep, nshells
       real ::  k1(nx/2), k2(ny), k3(nz), mk1(nx/2), mk2(ny), mk3(nz) &
              , kx(nx/2), ky_(nx/2,ny), ky(nx/2,ny), kz(nz)
       complex :: sx(nx/2,3), sy(ny,3), sz(nz,3)
       integer(int64) :: trigx, trigy, trigz, trigxy

       complex, allocatable ::  u(:,:,:,:)[:]    ! u(nz,4,first_x:last_x,ny)[*]    !(*-- x-y planes --*)
       complex, allocatable ::  ur(:,:,:,:)[:]   !ur(nz,4,first_y:last_y,nx/2)[*]  !(*-- x-z planes --*)
       complex, allocatable ::  un(:,:,:,:)      !un(nz,3,first_x:last_x,ny)[*]    !(*-- x-y planes --*)
       complex, allocatable :: bufr_X_Y(:,:,:,:)
       complex, allocatable :: bufr_Y_X(:,:,:,:)

interface
!--------    note: integer(int64)'s required for FFT's and other assembly-coded externals   ------

   function ctrig( len ) bind(C)               !(*-- define complex FFT trig table --*)
        import int64
        integer(int64), value, intent(in) :: len
        integer(int64) :: ctrig     !-- C pointer!
   end function ctrig

   function rtrig( len ) bind(C)              !(*-- define real FFT trig table --*)
        import int64
        integer(int64), value, intent(in):: len
        integer(int64) :: rtrig     !-- C pointer!
   end function rtrig

   subroutine cfft( len, lot, data, inc, jmp, ctrig, isign ) bind(C)    !(*-- complex FFT --*)
        import int64
        integer(int64), value, intent(in) :: len, lot, inc, jmp, ctrig, isign
        complex, dimension(0:0), intent(in) :: data
   end subroutine cfft

   subroutine rfft( len, lot, data, inc, jmp, rtrig, isign ) bind(C)   !(*-- real FFT --*)
        import int64
        integer(int64), value, intent(in) :: len, lot, inc, jmp, rtrig, isign
        complex, dimension(0:0), intent(in) :: data
   end subroutine rfft

   function WALLTIME() bind(C, name = "WALLTIME")
       import real64
       real(real64) :: WALLTIME
   end function WALLTIME


end interface

       trigx = rtrig( nx )
       trigy = ctrig( ny )
       trigz = ctrig( nz )
	   trigxy = ctrig( nx+ny )

       allocate (  u(nz , 4 , first_x:last_x , ny)[*] )     !(*-- y-z planes --*)
       allocate ( ur(nz , 4 , first_y:last_y , nx/2)[*] )   !(*-- x-z planes --*)
       allocate ( un(nz , 3 , first_x:last_x , ny) )        !(*-- y-z planes --*)
       allocate ( bufr_X_Y(nz,4,mx,my) )
       allocate ( bufr_Y_X(nz,4,my,mx) )


       stop = 0;   step = 0;  rkstep = 2;  rflag = 0;  cfl = 1;   dt = 0
       nshells = max( nx,ny,nz )

                        call define_kspace
                        call define_field
                        call enforce_conjugate_symmetry
						call copy_n_s
                        call define_shifts

       total_time =  -WALLTIME()      !-- start the clock

       tran_time = 0;       cpu_time = -WALLTIME()

  !(*********************************   begin execution loop   *****************************************)

       do while (stop == 0)

                        call phase1
            rkstep = 1
                        call transpose_X_Y
                        call phase2
                        call transpose_Y_X
                        call define_step
                        call define_shifts
                        call phase3
                        call pressure
		if (oflag /= 0) call spectra
                        call advance
                        call phase1
            rkstep = 2
                        call transpose_X_Y
                        call phase2
                        call transpose_Y_X
                        call phase3
                        call advance
                        call pressure
        if (rflag /= 0) call remesh
                        call copy_s_n

                		step = step + 1
                        time = time + dt
       end do

  !(*********************************   end execution loop   ***********************************************)

       deallocate ( u, ur, un )
       deallocate ( bufr_X_Y );    deallocate ( bufr_Y_X )
       sync all     !-- wait for all images to finish!

       total_time =  total_time + WALLTIME()    !-- stop the clock
       cpu_time =  cpu_time + WALLTIME()    !-- stop the clock
       call global_times

       if (my_node == 1 )   write(6,fmt="(3(10X,A,2f7.2))")  &
                            , "total_time ", min_total_time/step, max_total_time/step &
                            , "cpu_time ", min_cpu_time/step, max_cpu_time/step &
                            , "tran_time ", min_tran_time/step, max_tran_time/step


       write(6,fmt="(A,i4,3f7.2)")  "image ", my_node, total_time/step, cpu_time/step, tran_time/step


contains

 !(***********************************************************************************************************
 !                          transpose the Y and Z planes
 !***********************************************************************************************************)

!-----                   u(nz,4,mx,my*num_nodes) [num_nodes]
!-----                  ur(nz,4,my,mx*num_nodes) [num_nodes]
!-----                bufr(nz,4,my,mx) or bufr(nz,4,mx,my)

!-------------   out-of-place transpose data_s --> data_r  ----------------------------

 subroutine transpose_X_Y

    use constants_module, only : one
    use run_size
    implicit none

    integer(int64) :: i,stage

    cpu_time = cpu_time + WALLTIME()
    sync all   !--  wait for other nodes to finish compute
    tran_time = tran_time - WALLTIME()

    call copy3 (    u(1,1,first_x,1+(my_node-1)*my) &                   !-- intra-node transpose
                ,  ur(1,1,first_y,1+(my_node-1)*mx) &                   !-- no inter-node transpose needed
                ,   nz*3, one, one        &                                 !-- note: only 3 of 4 words needed
                ,   mx, nz*4, nz*4*my &
                ,   my, nz*4*mx, nz*4 )

    do stage = 1, num_nodes-1
        i = 1 + mod( my_node-1+stage, num_nodes )
        bufr_X_Y(:,:,:,:) = u(:,:,:,1+(my_node-1)*my:my_node*my)[i]         !-- inter-node transpose to buffer
        call copy3 ( bufr_X_Y, ur(1,1,first_y,1+(i-1)*mx)  &                !-- intra-node transpose from buffer
                        ,   nz*3, one, one        &                             !-- note: only 3 of 4 words needed
                        ,   mx, nz*4, nz*4*my &
                        ,   my, nz*4*mx, nz*4 )
    end do

    sync all     !--  wait for other nodes to finish transpose
    tran_time = tran_time + WALLTIME()
    cpu_time = cpu_time - WALLTIME()

 end  subroutine transpose_X_Y

!-------------   out-of-place transpose data_r --> data_s  ----------------------------

subroutine transpose_Y_X
    use run_size
    implicit none

    integer(int64) :: i, stage

    cpu_time = cpu_time + WALLTIME()
    sync all   !--  wait for other nodes to finish compute
    tran_time = tran_time - WALLTIME()

    call copy3 (   ur(1,1,first_y,1+(my_node-1)*mx) &                   !-- intra-node transpose
                ,   u(1,1,first_x,1+(my_node-1)*my) &                   !-- no inter-node transpose needed
                ,   nz*4, one, one        &                                 !-- note: all 4 words needed
                ,   my, nz*4, nz*4*mx &
                ,   mx, nz*4*my, nz*4 )

    do stage = 1, num_nodes-1
        i = 1 + mod( my_node-1+stage, num_nodes )
        bufr_Y_X(:,:,:,:) = ur(:,:,:,1+(my_node-1)*mx:my_node*mx)[i]        !-- inter-node transpose to buffer
        call copy3 ( bufr_Y_X, u(1,1,first_x,1+(i-1)*my)  &                 !-- intra-node transpose from buffer
                    ,   nz*4, one, one        &
                    ,   my, nz*4, nz*4*mx &
                    ,   mx, nz*4*my, nz*4 )
    end do

    sync all     !--  wait for other nodes to finish transpose
    tran_time = tran_time + WALLTIME()
    cpu_time = cpu_time - WALLTIME()

 end  subroutine transpose_Y_X


!(*************************************************************************************************************
!           enforce conjugate symmetry for plane kx=0 of wavespace  (half of this plane is redundant)
!***************************************************************************************************************)

  subroutine enforce_conjugate_symmetry

    integer(int64) ::  i, x, y, z

!(*------------------------- un( K ) = conjg( un( -K ) ) ---------------------*)

    if (my_node == 1 ) then     !-- x=1 is in node=1
        x = 1
        do i = 1, 3
            z = 1;              y = 1;              un(z,i,x,y) = 0
            z = 1;              do y = 2, ny/2;     un(z,i,x,y) = conjg( un(z,i,x,ny+2-y) );        end do
            do z = 2, nz/2;     y = 1;              un(z,i,x,y) = conjg( un(nz+2-z,i,x,y) );        end do
            do z = 2, nz/2;     do y = 2, ny;       un(z,i,x,y) = conjg( un(nz+2-z,i,x,ny+2-y) );   end do;	end do
        end do
    end if
end  subroutine enforce_conjugate_symmetry

 !(***********************************************************************************************************
 !                spectra :  accumulate spectra and other statistics over flow field
 !***********************************************************************************************************)

 subroutine  spectra

      use run_size
     implicit none

    integer(int64) :: k, x, y, z
    real :: kk, ww, uw, uu, uv, duu, factor   &
          , ek(nshells), dk(nshells), hk(nshells), tk(nshells), sample(nshells)
    real, save, codimension[*] ::  sum_ek, sum_dk, sum_hk, sum_tk

      total_time =  total_time + WALLTIME()     !-- stop the clock!  time/step does not include spectra time

	  oflag = 0
      ek = 0;       dk = 0;       hk = 0;     tk = 0;   sample = 0

 !(*---------------------   three dimensional spectra  -----------------------*)

      do x = first_x, last_x;   do y = 1, ny;      do  z = 1, nz

            if( mk1(x)+mk2(y)+mk3(z) > 2./9. ) &
                              then;   factor = 0
            else if (x == 1)  then;   factor = 1
                              else;   factor = 2
            end if

            kk = kx(x)**2 + ky(x,y)**2 + kz(z)**2
            k = 1 + int( sqrt( kk ) + 0.5 )

			  uu = factor * real( un(z,1,x,y) * conjg( un(z,1,x,y) ) &
                                + un(z,2,x,y) * conjg( un(z,2,x,y) ) &
                                + un(z,3,x,y) * conjg( un(z,3,x,y) ) )
              ww = kk * uu
              uv = factor * real( un(z,1,x,y) * conjg( un(z,2,x,y) ) )

              uw = factor * 2 * aimag( kx(x) *   un(z,2,x,y) * conjg( un(z,3,x,y) ) &
                                     + ky(x,y) * un(z,3,x,y) * conjg( un(z,1,x,y) ) &
                                     + kz(z) *   un(z,1,x,y) * conjg( un(z,2,x,y) ) )

			  duu = factor * real( un(z,1,x,y) * conjg( u(z,1,x,y) ) &
                                 + un(z,2,x,y) * conjg( u(z,2,x,y) ) &
								 + un(z,3,x,y) * conjg( u(z,3,x,y) ) ) / (dt/2) + shear * uv

              sample(k) = sample(k) + factor        !(*-- shell sample --*)
              ek(k) = ek(k) + uu                    !(*-- 2 * energy sum --*)
              dk(k) = dk(k) + ww                    !(*-- enstrophy sum --*)
              hk(k) = hk(k) + uw                    !(*-- helicity sum --*)
              tk(k) = tk(k) + duu                   !(*-- transfer sum --*)

      end do;   end do;  end do

 !(************************     finished accumulation :  compute final statistics     *************************)

    sum_ek = 0;     sum_dk = 0;     sum_hk = 0;     sum_tk = 0
    do k = nshells, 1, -1
        sum_ek  = sum_ek + ek(k)
        sum_dk  = sum_dk + dk(k)
        sum_hk  = sum_hk + hk(k)
        sum_tk  = sum_tk + tk(k)
    end do

    sync all
    if (my_node == 1)  then
        do k = 2, num_nodes
            sum_ek = sum_ek + sum_ek[k]
            sum_dk = sum_dk + sum_dk[k]
            sum_hk = sum_hk + sum_hk[k]
            sum_tk = sum_tk + sum_tk[k]
        end do

        if (step == 0)   write(6,*) "step   time     energy    enstrophy   helicity   transfer"
        write(6,fmt="(i3, 5e11.3)")  step,  time,    sum_ek,   sum_dk,     sum_hk,    sum_tk
    end if

     total_time =  total_time - WALLTIME()     !-- restart the clock!
 end  subroutine  spectra

  !(************************************************************************************************************
  !        define_field  :   define initial flow field from scratch
  !************************************************************************************************************)

 subroutine define_field

    use constants_module, only : zero
    use run_size
    use random_module
    implicit none

    real ::     k, k12, f, phi, theta1, theta2
    complex ::  alpha, beta
    integer(int64) ::  x, y, z
    real, parameter :: klo=8, khi=16

   call init_random_seed !(* seed a different pseudo-random number sequence for each image *)
   time = 0

   do x = first_x, last_x
       do  y = 1, ny
            do z = 1, nz
                 call random_number(theta1 )
                 call random_number(theta2 )
                 call random_number(phi    )
                 k   = sqrt( kx(x)**2 + ky(x,y)**2 + kz(z)**2 )
                 k12 = sqrt( kx(x)**2 + ky(x,y)**2 )

                 if ( k == 0  .or.  mk1(x)+mk2(y)+mk3(z)>2./9.  .or.  k < klo  .or.  k > khi ) &
                     then;   f = 0
                     else;   f = sqrt( 1./(2*pi) ) / k
                 end if

                 alpha = f * exp( (0,2) * pi * theta1 ) * cos( 2*pi * phi )
                 beta  = f * exp( (0,2) * pi * theta2 ) * sin( 2*pi * phi )

                 if (k12 == 0) &
                 then; un(z,1,x,y) = alpha
                       un(z,2,x,y) = beta
                       un(z,3,x,y) = 0

                 else; un(z,1,x,y) = ( beta * kz(z) * kx(x)   + alpha * k * ky(x,y) ) / ( k * k12 )
                       un(z,2,x,y) = ( beta * kz(z) * ky(x,y) - alpha * k * kx(x)   ) / ( k * k12 )
                       un(z,3,x,y) = - beta * k12 / k
                 end if

   end do;  end do; end do
 end  subroutine define_field

 !(***********************************************************************************************************
 !          define_shifts  :    define coordinate shifts for control of 1-d alias errors
 ! ***********************************************************************************************************)

    subroutine  define_shifts
    use constants_module, only : zero
    use run_size
    implicit none

           integer(int64) ::  x, y, z
           integer(int64), save ::  init = 0
           real :: delta_x, delta_y, delta_z
           integer :: i,seed_size

           if (init == 0) &     !-- Note: delta's not carried over from previous run
           then;
                init = 1
                call random_seed(size=seed_size)
                call random_seed(put=[(1234567,i=1,seed_size)])!(* same random numbers for each image! *)
                do  x = 1, nx/2;  sx(x,3) = exp (  (0,1) * ( pi / nx ) * k1(x) ); end do
                do  y = 1, ny  ;  sy(y,3) = exp (  (0,1) * ( pi / ny ) * k2(y) ); end do
                do  z = 1, nz  ;  sz(z,3) = exp (  (0,1) * ( pi / nz ) * k3(z) ); end do
            else;
                call random_number(delta_x); delta_x = 2*pi / nx * delta_x
                do  x = 1, nx/2;  sx(x,1) = sx(x,3)
                                  sx(x,2) = exp (  (0,1) * delta_x * k1(x) )
                                  sx(x,3) = exp (  (0,1) * ( delta_x + pi / nx ) * k1(x) ); end do

                call random_number(delta_y); delta_y = 2*pi / ny * delta_y
                do  y = 1, ny  ;  sy(y,1) = sy(y,3)
                                  sy(y,2) = exp (  (0,1) * delta_y * k2(y) )
                                  sy(y,3) = exp (  (0,1) * ( delta_y + pi / ny ) * k2(y) ); end do

                call random_number(delta_z); delta_z = 2*pi / nz * delta_z
                do  z = 1, nz  ;  sz(z,1) = sz(z,3)
                                  sz(z,2) = exp (  (0,1) * delta_z * k3(z) )
                                  sz(z,3) = exp (  (0,1) * ( delta_z + pi / nz ) * k3(z) ); end do
           end if

 end  subroutine  define_shifts

  !(***********************************************************************************************************
  !       define_step  :   update time, metric, shifts for the next step
  !**********************************************************************************************************)

  subroutine  define_step
      use run_size
      implicit none

    sync all

    if (cfl /= 0) then
cpu_time = cpu_time + WALLTIME()
        call max_velmax
cpu_time = cpu_time - WALLTIME()
        dt = cfl / velmax
    end if

    if (        shear > 0               &
        .and.  .01*b11*shear*dt < b12   &
        .and.  b12 <= b11*shear*dt )    then
            dt = b12 / ( b11 * shear )          !(* limit dt, hit the orthognal mesh *)
            oflag = 1
    else if ( mod (step,output_step) == 0 ) then
            oflag = 1
    end if

    b12 = b12 - b11 * shear * dt

    if ( b12 < -b22/2 )     rflag = 1                                   !(* remesh at the end of the step? *)
    if ( step == nsteps )   stop = 1                                    !(* last step? *)

 end   subroutine    define_step

   !(***********************************************************************************************************
   !      define_kspace  :   define physical wavespace from computational wavespace and metric
   !**********************************************************************************************************)

  subroutine    define_kspace
    use run_size
    implicit none

     integer(int64) ::  x, y, z

       do  x = 1, nx/2   ;   k1(x) = x - 1;     end do
       do  y = 1, ny/2+1 ;   k2(y) = y - 1;     end do
       do  z = 1, nz/2+1 ;   k3(z) = z - 1;     end do

       do  y = ny/2+2, ny;   k2(y) = y - 1 - ny;    end do
       do  z = nz/2+2, nz;   k3(z) = z - 1 - nz;    end do

       do  x = 1, nx/2 ;     mk1(x) = (k1(x)/nx)**2;    kx(x) =   b11 * k1(x);  end do
       do  z = 1, nz   ;     mk3(z) = (k3(z)/nz)**2;    kz(z) =   b33 * k3(z);  end do
       do  y = 1, ny   ;     mk2(y) = (k2(y)/ny)**2
                             do  x = 1, nx/2 ;    ky(x,y) = b22 * k2(y) + b12 * k1(x);  end do; end do

end   subroutine    define_kspace

  !(***********************************************************************************************************
  !   phase 1 :  on entry, data-plane contains velocity in wave space.  interpolate database, shifted mesh,
  !              and proceed to physical y space .
  !************************************************************************************************************)

  subroutine  phase1

      use run_size
      implicit none

    complex :: shift
    integer(int64) :: i, x, y, z

   	do  x = first_x, last_x

   		do  y = 1, ny;    do  z = 1, nz
             shift = sz(z,rkstep+1) * sy(y,rkstep+1) * sx(x,rkstep+1)
             u(z,1,x,y) = shift * u(z,1,x,y)
             u(z,2,x,y) = shift * u(z,2,x,y)
             u(z,3,x,y) = shift * u(z,3,x,y)
        end do; end do

!(*---------------------------   LEAVING FOURIER WAVE SPACE  --------------------------*)

        do i = 1, 3
            call cfft ( ny, nz, u(1,i,x,1), nz*4*mx, one, trigy, one );    end do
	end do

 end   subroutine  phase1

 !(**********************************************************************************************************
 !     phase 2 :  on entry, data-plane contains velocity in physical y space, and wave x,z space on shifted
 !                mesh.  Proceed to physical x,z space,  form nonlinear terms, and return to wave x,z space.
 !***********************************************************************************************************)

 subroutine  phase2

     use run_size
     implicit none

     complex :: s2(nz,nx/2), vs(nz,nx/2)
     integer(int64) ::  i, x, y, z
     real :: v2r, v2i, s2r, s2i, u1r, u1i, u2r, u2i, u3r, u3i, u4r, u4i

  velmax = 0

  do  y = first_y, last_y

      do  x = 1, nx/2 ;  do  z = 1, nz ;   vs(z,x) = ur(z,2,y,x);    end do; end do

      do  i = 1, 3
           call cfft ( nz, nx/2, ur(1,i,y,1), one, nz*4*my, trigz, one )
           call rfft ( nx, nz,   ur(1,i,y,1), nz*4*my, one, trigx, one )
      end do

!(*----------------------------  WELCOME TO PHYSICAL SPACE  --------------------------*)

      do  x = 1, nx/2;  do  z = 1, nz
           u1r = real(ur(z,1,y,x));  u1i = aimag(ur(z,1,y,x))
           u2r = real(ur(z,2,y,x));  u2i = aimag(ur(z,2,y,x))
           u3r = real(ur(z,3,y,x));  u3i = aimag(ur(z,3,y,x))

      if ( rkstep == 1 ) velmax = max( velmax &
                                     , b11*nx*abs(u1r) + b22*ny*abs(u2r) + b33*nz*abs(u3r) &
                                     , b11*nx*abs(u1i) + b22*ny*abs(u2i) + b33*nz*abs(u3i) )

           v2r = u2r * u2r;             v2i = u2i * u2i
           s2r = u1r * u3r;             s2i = u1i * u3i
           u4r = u2r * u3r;             u4i = u2i * u3i
           u3r = u3r * u3r  -  v2r;     u3i = u3i * u3i  -  v2i
           u2r = u1r * u2r;             u2i = u1i * u2i
           u1r = u1r * u1r  -  v2r;     u1i = u1i * u1i  -  v2i

           s2(z,x)     = cmplx(s2r, s2i)
           ur(z,1,y,x) = cmplx(u1r, u1i)
           ur(z,2,y,x) = cmplx(u2r, u2i)
           ur(z,3,y,x) = cmplx(u3r, u3i)
           ur(z,4,y,x) = cmplx(u4r, u4i)
      end do;   end do

!(*----------------------------  LEAVING PHYSICAL SPACE  --------------------------*)

      do  i = 1,  4
           call rfft ( nx, nz,   ur(1,i,y,1), nz*4*my, one, trigx, -one )
           do  z = 1, nz ;   ur(z,i,y,1) = cmplx(real(ur(z,i,y,1)),0);   end do
           call cfft ( nz, nx/2, ur(1,i,y,1), one, nz*4*my, trigz, -one )
      end do

      call rfft ( nx, nz, s2, nz, one, trigx, -one )
      do  z = 1, nz ;   s2(z,1) = cmplx(real(s2(z,1)),0);    end do
      call cfft ( nz, nx/2, s2, one, nz, trigz, -one )

      do  x = 1, nx/2;  do  z = 1, nz
          ur(z,1,y,x) = kx(x) * ur(z,1,y,x) + kz(z) * s2(z,x)   - (0,1) * 2*nx*nz*shear * vs(z,x)
          ur(z,3,y,x) = kx(x) * s2(z,x)     + kz(z) * ur(z,3,y,x)
      end do;   end do
  end do

 end  subroutine  phase2

  !(***********************************************************************************************************
  !     phase 3 :  on entry, the data-plane contains the four stresses on a shifted mesh in physical y space,
  !                wave x,z space.   Return to y  wave space on unshifted mesh and complete time derivative of
  !                velocity ( not divergence free yet )
  !***********************************************************************************************************)

  subroutine  phase3

      use run_size
      implicit none

    integer(int64) :: i, x, y, z
    complex :: shift

    do  x = first_x, last_x

       do i = 1, 4
           call cfft ( ny, nz, u(1,i,x,1), nz*4*mx, one, trigy, -one )
       end do

!(*---------------------------   WELCOME TO FOURIER WAVE SPACE  --------------------------*)

       do  y = 1, ny ;   do  z = 1, nz
               shift = -dt / (4*nx*ny*nz) * (0,1)*conjg( sy(y,rkstep) * sz(z,rkstep) * sx(x,rkstep) )
               u(z,1,x,y) = shift * (         u(z,1,x,y) + ky(x,y) * u(z,2,x,y) )
               u(z,2,x,y) = shift * ( kx(x) * u(z,2,x,y) + kz(z)   * u(z,4,x,y) )
               u(z,3,x,y) = shift * (         u(z,3,x,y) + ky(x,y) * u(z,4,x,y) )
       end do;  end do
    end do

 end   subroutine  phase3

  !(***********************************************************************************************************
  !   pressure :  add the gradient of a scalar, enforce continuity ( zero divergence )
  !***********************************************************************************************************)

  subroutine  pressure

      use run_size
      implicit none

     complex :: psi
     integer(int64) :: x, y, z

      do x = first_x, last_x ;     do  y = 1, ny

            if ( x /= 1 )  then
                  do  z = 1, nz
                        psi = ( kx(x) * u(z,1,x,y) + ky(x,y) * u(z,2,x,y) + kz(z) * u(z,3,x,y) ) &
                              / ( kx(x)**2 + ky(x,y)**2 + kz(z)**2 )
                        u(z,1,x,y) = u(z,1,x,y) - kx(x) * psi
                        u(z,2,x,y) = u(z,2,x,y) - ky(x,y) * psi
                        u(z,3,x,y) = u(z,3,x,y) - kz(z) * psi
                   end do
			else if ( y /= 1 )  then
                  do  z = 1, nz
                        psi = ( ky(1,y) * u(z,2,1,y) + kz(z) * u(z,3,1,y) ) &
                              / ( ky(1,y)**2 + kz(z)**2 )
                        u(z,2,1,y) = u(z,2,1,y) - ky(1,y) * psi
                        u(z,3,1,y) = u(z,3,1,y) - kz(z) * psi
                  end do
            else
                  do  z = 1, nz ;    u(z,3,1,1) = 0;     end do
            end if
	 end do;    end do

end   subroutine  pressure

!(*****************************************************************************************************************
!                                remesh  :  remesh the sheared coordinate system
!*****************************************************************************************************************)

subroutine   remesh

    use constants_module, only : one
    use run_size
    implicit none

    complex :: u2(nx+ny,nz), shift(nx+ny)
    integer(int64) :: i, x, y, z

    write(6,fmt="(A,i4)") "remesh image ", my_node

    total_time =  total_time + WALLTIME()     !-- stop the clock!

    do x = first_x, last_x

        do  y = 1, nx+ny ;   shift(y) =  exp( (0,-2) * pi / (nx+ny) * k1(x) * (y - 1) ) / (nx+ny);    end do

        do  i = 1,  3
            do  z = 1, nz
                do  y = 1, ny/2           ;   u2(y,z) = u(z,i,x,y);     end do
                do  y = ny/2+1, nx+ny/2+1 ;   u2(y,z) = 0;              end do
                do  y = nx+ny/2+2, nx+ny  ;   u2(y,z) = u(z,i,x,y-nx);  end do
            end do

            call cfft ( nx+ny, nz, u2, one , nx+ny, trigxy, one )

            do  z = 1, nz ;  do  y = 1, nx+ny ;     u2(y,z) = u2(y,z) * shift(y);  end do;  end do

            call cfft ( nx+ny, nz, u2, one, nx+ny, trigxy, -one )

            do  z = 1, nz
                do  y = 1, ny/2
                    if (mk1(x)+mk2(y)+mk3(z) > 2./9.) &
                        then;  u(z,i,x,y) = 0
                        else;  u(z,i,x,y) = u2(y,z)
                    end if
                end do
                do  y = ny/2+1, ny
                    if (mk1(x)+mk2(y)+mk3(z) > 2./9.) &
                        then;  u(z,i,x,y) = 0
                        else;  u(z,i,x,y) = u2(y+nx,z)
                    end if
                end do
            end do
        end do

        do  y = 1, ny ;    ky(x,y) = ky(x,y) + b22 * k1(x);     end do           !(* update ky for this x *)

    end do

    b12 = b12 + b22;     rflag = 0                                               !(*   update metric, account for remesh    *)

    total_time =  total_time - WALLTIME()     !-- restart the clock!
 end subroutine   remesh

 !(***********************************************************************************************************
 !         copy_n_s,   copy_s_n :  copy data between data_s and data_n
 !***********************************************************************************************************)

 subroutine  copy_n_s

     use run_size
     implicit none

	integer(int64) ::  x, y, z

	do y = 1, ny;    do x = first_x, last_x;   do z = 1, nz
		u(z,1,x,y) = un(z,1,x,y)
		u(z,2,x,y) = un(z,2,x,y)
    	u(z,3,x,y) = un(z,3,x,y)
	end do; end do; end do
 end  subroutine  copy_n_s

 subroutine  copy_s_n

     use run_size
     implicit none

	integer(int64) ::  x, y, z

	do y = 1, ny;	do x = first_x, last_x;     do z = 1, nz
		un(z,1,x,y) = u(z,1,x,y)
		un(z,2,x,y) = u(z,2,x,y)
		un(z,3,x,y) = u(z,3,x,y)
    end do; end do; end do

 end  subroutine  copy_s_n

 !(***********************************************************************************************************
 !                         advance :     second-order runge-kutta time step algorithm
 !***********************************************************************************************************)

 subroutine  advance

     use run_size
     implicit none

    integer(int64) ::  x, y, z
    real :: factor, xyfac, zfac(nz)	  !(* viscous integrating factors *)

    if (rkstep == 1) then
        do z = 1, nz;     zfac(z) = exp( - viscos * dt * kz(z)**2 ); end do

        do x = first_x, last_x
            do  y = 1, ny
                ky_(x,y) = ky(x,y)
                ky(x,y) = b22 * k2(y) + b12 * k1(x)

                do z = 1, nz
                    if (mk1(x)+mk2(y)+mk3(z) > 2./9.) then
                        u(z,1,x,y) = 0;     u(z,2,x,y) = 0;     u(z,3,x,y) = 0
                    else
                        factor = zfac(z) * exp( - viscos * dt * ( kx(x)**2 + ( ky_(x,y)**2 + ky_(x,y)*ky(x,y) + ky(x,y)**2 )/3 ) )

                        un(z,1,x,y) = factor * ( un(z,1,x,y) + u(z,1,x,y) )
                        u(z,1,x,y)  = un(z,1,x,y) + factor * u(z,1,x,y)

                        un(z,2,x,y) = factor * ( un(z,2,x,y) + u(z,2,x,y) )
                        u(z,2,x,y)  = un(z,2,x,y) + factor * u(z,2,x,y)

                        un(z,3,x,y) = factor * ( un(z,3,x,y) + u(z,3,x,y) )
                        u(z,3,x,y)  = un(z,3,x,y) + factor * u(z,3,x,y)
                    end if
        end do;  end do; end do

    else if (rkstep == 2) then

        do x = first_x, last_x
            do  y = 1, ny
                do z = 1, nz
                    if (mk1(x)+mk2(y)+mk3(z) > 2./9.) then
                        u(z,1,x,y) = 0;     u(z,2,x,y) = 0;     u(z,3,x,y) = 0
                    else
                        u(z,1,x,y) = un(z,1,x,y) + u(z,1,x,y)
                        u(z,2,x,y) = un(z,2,x,y) + u(z,2,x,y)
                        u(z,3,x,y) = un(z,3,x,y) + u(z,3,x,y)
                    end if
        end do; end do; end do

    end if

 end  subroutine  advance

 end   subroutine solve_navier_stokes
