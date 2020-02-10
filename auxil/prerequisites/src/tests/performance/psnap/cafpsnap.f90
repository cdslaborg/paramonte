! P-SNAP test
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
!     * Neither the name of the Sourcery, Inc., nor the
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

!  Fortran translation of

!/*
! * P-SNAP v1.2 -- PAL System Noise Activity Program -- LA-CC-06-025
! *        <http://www.c3.lanl.gov/pal/software/psnap/>
! *
! * Copyright (C) 2006, The Regents of the University of California
! *
! *                PAL -- Performance and Architecture Laboratory
! *                  <http://www.c3.lanl.gov/pal/>
! *                Los Alamos National Laboratory
! *                  <http://www.lanl.gov/>
! */

!  by Dan Nagle

! * This program is free software; you can redistribute it and/or modify
! * it under the terms of the GNU General Public License as published by
! * the Free Software Foundation; either version 2 of the License, or
! * (at your option) any later version.
! *
! * This program is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! * GNU General Public License for more details.
! *
! * You should have received a copy of the GNU General Public License
! * along with this program; if not, write to the Free Software
! * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
! * 02110-1301 USA.

program psnap

!  #define VERSION_STR "v1.2"

use, intrinsic :: iso_fortran_env, only: error_unit, output_unit, int64, real64, int32

use, intrinsic :: iso_c_binding

character( len= *), parameter :: psnap_rcs_id = &
   '$Id$'

character( len= *), parameter :: string_fmt = '( a)'

type, bind( c) :: counters

   integer( c_long) :: val
   integer( c_int) :: index

end type counters

!  int rank, np; // globals

integer :: rank
integer :: np
integer :: np_half

!integer( int64), save :: n = 100000
integer( int64), codimension[ *], save :: n = 1000
integer( int64), codimension[ *], save :: w = 1000
integer( int64), dimension( :), allocatable :: r
integer( int64) :: i
integer( int64), codimension[ *], save :: iteration_count = 0
integer( int64), codimension[ *] :: localmax, globalmax
integer( int64), codimension[ *] :: localsum
integer( int64), dimension( :), allocatable :: sum_all
integer( int64), dimension( :), allocatable, codimension[ :] :: localhist
integer( int64), codimension[ *], save :: granularity = 1000
integer( int64), codimension[ *], save :: barrier = 0

!character( kind= c_char, len= 1024), codimension[ *] :: hostname

interface
  subroutine start_timer() bind(C, name="start_timer")
    use iso_c_binding
  end subroutine
  subroutine stop_timer() bind(C, name="stop_timer")
    use iso_c_binding
  end subroutine
  function elapsed_time() bind(c,name="elapsed_time") result(res)
    use iso_c_binding
    !use, intrinsic :: iso_fortran_env, only: int64
    integer(c_int) :: res
  end function
end interface

type( counters), codimension[ *] :: count_loc, globalminloc, globalmaxloc

   integer( int64) :: j
   integer :: astat

   character( len= 32) :: cl_arg

! ----------------------------------------------------------------------

!  psnap text

continue

   rank = this_image()
   np = num_images()
   np_half = np / 2

   j = 0

   if( rank == 1 )then

      cl_args: do

         j = j + 1

         if( j >= command_argument_count() ) exit cl_args

         call get_command_argument( number= int( j, int32), value= cl_arg)

         which_arg: select case( cl_arg( 2: 2))

         case( 'b')

            j = j + 1
            call get_command_argument( number= int( j, int32), value= cl_arg)
            write( unit= barrier, fmt= *) cl_arg
            do i = 2, np
               barrier[ i] = barrier
            end do

         case( 'n')

            j = j + 1
            call get_command_argument( number= int( j, int32), value= cl_arg)
            write( unit= n, fmt= *) cl_arg
            w = n / 10
            do i = 2, np
               n[ i] = n
               w[ i] = w
            end do

         case( 'w')

            j = j + 1
            call get_command_argument( number= int( j, int32), value= cl_arg)
            write( unit= w, fmt= *) cl_arg
            do i = 2, np
               w[ i] = w
            end do

         case( 'c')

            j = j + 1
            call get_command_argument( number= int( j, int32), value= cl_arg)
            write( unit= iteration_count, fmt= *) cl_arg
            do i = 2, np
               iteration_count[ i] = iteration_count
            end do

         case( 'g')

            j = j + 1
            call get_command_argument( number= int( j, int32), value= cl_arg)
            write( unit= granularity, fmt= *) cl_arg
            do i = 2, np
               granularity[ i] = granularity
            end do

         case( 'h')

            call usage()
            stop 'normal exit in psnap'

         case default

            call usage()
            stop 'normal exit in psnap'

         end select which_arg

      end do cl_args

   end if

!  distribute sizes before allocating

   sync all

   !call get_environment_variable( name= 'HOSTNAME', value= hostname)

   allocate( r( 1: n + w), stat= astat)

   alloc_r_error: if( astat > 0 )then

      stop 'error allocating r'

   end if alloc_r_error

   allocate( sum_all( 1: np), stat= astat)

   alloc_sum_error: if( astat > 0 )then

      stop 'error allocating sum_all'

   end if alloc_sum_error

   only_rank_0: if( rank == 1 )then

      call print_banner()

   end if only_rank_0

! ----------------------------------------------------------------------

!  warmup loop here; calibration follows

   if( w > 0 ) call warmup_loop( w)

   if( iteration_count == 0 )then

      iteration_count = calibrate_loop( granularity)

	!write(*,*) 'Iteration after calibrate', iteration_count,'proc',rank

      count_loc% val = iteration_count
      count_loc% index = rank

      globalminloc = count_loc
      globalmaxloc = count_loc

!  compute global counts before communicating them

      sync all

      do i = 2, np
         if( globalminloc[ i]% val < globalminloc% val )then
            globalminloc = globalminloc[ i]
         end if
         if( globalmaxloc[ i]% val > globalmaxloc% val )then
            globalmaxloc = globalmaxloc[ i]
         end if
      end do

      if( rank == 1 )then
         write( unit= output_unit, fmt= '( a, i0/ a, i0, a, i0/ a, i0, a, i0)') "my_count= ", iteration_count, &
                                    "global_min= ", globalminloc% val, " min_loc= ", globalminloc% index, &
                                    "global_max= ", globalmaxloc% val, " max_loc= ", globalmaxloc% index
         write( unit= output_unit, fmt= string_fmt) "Using Global max for calibration"
      end if

      iteration_count = globalmaxloc% val

   end if

   r = 0

! ----------------------------------------------------------------------

!  measurement loop

   sync all

   do i = 1, n + w
      r( i) = loop( iteration_count)
      if( barrier /= 0 )then
         if(  mod( i, barrier) == 0 ) sync all
      end if
   end do

   sync all

! ----------------------------------------------------------------------

!   build histograms

   localsum = sum( r( w+1: ) )

   if( rank == 1 )then

      sum_all( 1) = localsum

      do i = 2, np
         sync images( i)
         sum_all( i) = localsum[ i]
      end do

   else
      sync images( 1)
   end if

   localmax = maxval( r( w+1: ))

   if( rank == 1 )then

      globalmax = localmax

      do i = 2, np
         sync images( i)
         globalmax = max( globalmax, localmax[ i])
      end do

   else
      sync images( 1)
   end if

   if( rank == 1 )then

      do i = 2, np
         globalmax[ i] = globalmax
      end do

   end if

   sync all

   allocate( localhist( 0: globalmax)[ *], stat= astat )

   alloc_localhist_error: if( astat > 0 )then

      stop 'error allocating localhist'

   end if alloc_localhist_error

   localhist = 0

   make_hist: do i = 1+w, n+w

      localhist( r( i)) = localhist( r( i)) + 1

   end do make_hist

   final_print: if( rank == 1 )then

!  print rank 0's histogram

      if( n > 0 )then
         write( unit= output_unit, fmt= '(a, i9, 9x, i9, 3x)') "#", 1, sum_all( 1)!, trim( hostname)
         do i = 0, globalmax
            if( localhist( i) > 0 )then
               write( unit= output_unit, fmt= '(1x, i9, i9, i9, 3x)') rank, i, localhist( i)!, trim( hostname)
            end if
         end do

      end if

!  print rank i's histogram

      do i = 2, np

         sync images( i)
         !localhist(:) = localhist(:)[ i]
         do j=lbound(localhist,dim=1),ubound(localhist,dim=1)
           localhist(j) = localhist(j)[ i]
         end do

         !hostname = hostname[ i]
         if( n > 0 )then
            write( unit= output_unit, fmt= '(a, i9, 9x, i9, 3x)') "#", i, sum_all( i)!, trim( hostname)
            do j = 0, globalmax
               if( localhist( j) > 0 )then
                  write( unit= output_unit, fmt= '(1x, i9, i9, i9, 3x)') i, j, localhist( j)!, trim( hostname)
               end if
            end do

         end if
      end do

   else final_print

      sync images( 1)

   end if final_print

stop 'normal exit in psnap'

contains

! ---------------------------------------------------------------------

function get_usecs() result( usecs)
!integer( int64) :: usecs

!  usec per sec

!integer( int64), parameter :: c = 1000000

!   integer( int64) :: t

!   integer( int64) :: r

integer(kind=8) :: t,r
real(real64) :: usecs
integer(kind=8),parameter :: c = 1000000

continue

   !call system_clock( count= t, count_rate= r)
   call cpu_time(usecs)
   usecs = usecs*1.d6
   !if( r /= c ) usecs = int( real( t, 8) / real( r, 8) * real( c, 8))
   !if( r /= c ) usecs = int( t/r * c)

return

end function get_usecs

! ---------------------------------------------------------------------

function loop( iterations) result( dt)
integer(int64) :: dt
integer( int64), intent( in) :: iterations

   integer( int64) :: i

   integer(int64) :: usecs_init, usecs_final

   integer :: next_rank, prev_rank

   integer, codimension[ *], save :: coarray

continue

!   usecs_init = get_usecs()

   call start_timer()

   next_rank = to_upper_half( rank)
   prev_rank = from_lower_half( rank)

!   write(*,*) 'Proc',rank,'Next rank',next_rank,'Prev rank',prev_rank

   counter: do i = 1, iterations

      even_odd: if( sending_half( rank) )then

!  send rank to next then fetch rank

         coarray[ next_rank] = coarray

         sync images( next_rank)

      else even_odd

!  stay calm

         coarray[ prev_rank] = coarray

         sync images( prev_rank)

      end if even_odd

   end do counter

!   usecs_final = get_usecs()
   call stop_timer()

!   write(*,*) 'usec_init',usec_init,'usec_final',usec_final

   dt = elapsed_time()

   !write(*,*) 'usec_init',usec_init,'usec_final',usec_final,'dt',dt

return

end function loop

! ---------------------------------------------------------------------

subroutine warmup_loop( wa)
integer( int64), intent( in) :: wa

!integer( int64), parameter :: counter = 1000000
integer( int64), parameter :: counter = 10000

   integer( int64) :: min_time_usecs
   integer( int64) :: loop_time

   integer :: i

continue

   min_time_usecs = huge( 0_int64)

   reloop: do i = 1, wa

      loop_time = loop( counter)

      min_time_usecs = min( loop_time, min_time_usecs)

   end do reloop

return

end subroutine warmup_loop

! ---------------------------------------------------------------------

function calibrate_loop( usecs) result( cl)
integer( int64) :: cl
integer( int64), intent( in) :: usecs

integer( int64), parameter :: calibrate_useconds = 100000000
!real( real64), parameter :: preset_tolerance = 0.001_real64
real( real64), parameter :: preset_tolerance = 1.0_real64

integer( int64), parameter :: initial_ntrial = 1000

!integer( int64), parameter :: initial_counter = 1000000
integer( int64), parameter :: initial_counter = 100000

   integer( int64) :: counter
   integer( int64) :: min_time_usecs
   integer( int64) :: tolerance
   integer( int64) :: difference
   integer( int64) :: total_time

   integer( int64) :: ntrial
   integer( int64) :: i

   integer( int64) :: loop_time

continue

   counter = initial_counter

   !write(*,*) 'Counter after initial counter',counter

!  if usecs / granularity is less than 1/preset_tolerance then use zero

   tolerance = int( real( usecs, real64) * preset_tolerance, int64)

   total_time = 0

   trials: do

      ntrial = initial_ntrial
      min_time_usecs = huge( 0_int64)

      get_min: do i = 1, ntrial

         loop_time = loop( counter)
	 !write(*,*) 'loop_time',loop_time
         min_time_usecs = min( min_time_usecs, loop_time)

      end do get_min

!  keep an estimate of total calibration time

      total_time = total_time + min_time_usecs * ntrial

      counter = int( real( counter, real64) * real( usecs, real64) / real( min_time_usecs, real64), int64)

      !write(*,*) 'Counter after assignment',counter

      difference = abs( min_time_usecs - usecs)

      if( difference <= tolerance .or. total_time >= calibrate_useconds ) exit trials

   end do trials

   cl = counter

   write( unit= output_unit, fmt= '( a, i2, a, i10, a, i10, a, i10, a, i10)' ) "#rank= ", rank, &
                                                                           " count= ", counter, " time= ", min_time_usecs, &
                                                                           " difference= ", difference, " tolerance= ", tolerance

   time_out: if( total_time > calibrate_useconds )then

      write( unit= output_unit, fmt= '( a, i2, a, f10.4, a, f10.4, a, f10.4, a, i0)' ) "PSNAP: WARNING rank ", rank, &
                                                                  " didn't converge in 10 seconds tolerance = ",  &
                                                                  real( difference) / real( usecs), &
                                                                  " should be ", preset_tolerance, " approx ", &
                                                                  preset_tolerance * 100.0, " percent, granularity= ", usecs

   end if time_out

return

end function calibrate_loop

! ---------------------------------------------------------------------

subroutine print_banner()

continue

   write( unit= output_unit, fmt= string_fmt) '########'
   write( unit= output_unit, fmt= string_fmt) '##P-SNAP: PAL System Noise Activity Program'
   write( unit= output_unit, fmt= string_fmt) '##' // psnap_rcs_id
   write( unit= output_unit, fmt= string_fmt) '##This is a Fortran translation of P-SNAP v 1.2 from'
   write( unit= output_unit, fmt= string_fmt) '##http://www.c3.lanl.gov/pal/software/psnap/'
   write( unit= output_unit, fmt= string_fmt) '##This program is the coarray ping-pong version'
   write( unit= output_unit, fmt= string_fmt) '########'

return

end subroutine print_banner

! ---------------------------------------------------------------------

subroutine usage()

character( len= *), dimension( 17), parameter :: msg = &
 [ "Usage: psnap [OPTIONS]                                             ", &
   "                                                                   ", &
   "  -n <reps>   number of repetitions                                ", &
   "                default: 100000                                    ", &
   "  -w <reps>   number of warm-up repetitions                        ", &
   "                default: 10%% of the number of reps                ", &
   "  -c <count>  calibration count                                    ", &
   "                default: perform a calibration to match granularity", &
   "  -g <usecs>  granularity of the test in microseconds              ", &
   "                default: 1000                                      ", &
   "  -b <N>      perform a barrier between every N loops              ", &
   "                default: no                                        ", &
   "  -h          this message                                         ", &
   "                                                                   ", &
   "  Example: psnap -n 1000000 -w 10 > psnap.out                      ", &
   "    runs a test with 1000000 repetitions and 10 warm-up reps.      ", &
   "                                                                   " ]

   integer :: i

continue

   write( unit= error_unit, fmt= string_fmt) ( trim( msg( i)), i = 1, size( msg, 1))

stop 'normal exit in usage'

end subroutine usage

! ---------------------------------------------------------------------

function sending_half( i) result( l)

integer, intent( in) :: i
logical :: l

continue

!  this must process np == even values only

   l = i <= np_half

return

end function sending_half

! ---------------------------------------------------------------------

function to_upper_half( i) result( l)

integer, intent( in) :: i
integer :: l

continue

!  this must process lower half ranks only

   l = i + np_half

return

end function to_upper_half

! ---------------------------------------------------------------------

function from_lower_half( i) result( l)

integer, intent( in) :: i
integer :: l

continue

!  this must process upper half ranks only

   l = i - np_half

return

end function from_lower_half

! ---------------------------------------------------------------------

end program psnap
