!*********************************************************************
! ForTrilinos: Object-Oriented Fortran 2003 interface to Trilinos
!                Copyright 2010 Sandia Corporation
!
! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
! the U.S. Government retains certain rights in this software.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the Corporation nor the names of the
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Questions? Contact Karla Morris  (knmorri@sandia.gov) or
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************

module ForTrilinos_error
#include "compiler_capabilities.txt"
  implicit none
  private
  public :: error

#ifdef ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS
  integer ,parameter :: max_string_length=256
#endif /* ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS */

  type :: error
    private
    integer code
    class(*), allocatable :: data_(:)
#ifdef ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS
    character(len=max_string_length) message ! gfortran 4.7.0 workaround
#else
    character(:) ,allocatable :: message
#endif /* ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS */
  contains
    procedure :: error_code
    procedure :: define_error
    generic :: error=>define_error
  end type

contains

  subroutine define_error(this,new_code,new_message,new_data)
    class(error), intent(out) :: this
    integer ,intent(in) :: new_code
    character(len=*) ,intent(in) :: new_message
    class(*) ,intent(in), optional :: new_data(:)
    this%code = new_code
    this%message = new_message
    if (present(new_data)) allocate(this%data_(lbound(new_data,1):ubound(new_data,1)),source=new_data)
  end subroutine

  integer function error_code(this)
    class(error) ,intent(in) :: this
    error_code = this%code
  end function

  function error_message(this)
    class(error) ,intent(in) :: this
    character(:), allocatable :: error_message
    error_message = this%message
  end function

end module ForTrilinos_error
