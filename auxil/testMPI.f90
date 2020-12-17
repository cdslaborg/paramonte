!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   MIT License
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   Permission is hereby granted, free of charge, to any person obtaining a 
!   copy of this software and associated documentation files (the "Software"), 
!   to deal in the Software without restriction, including without limitation 
!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!   and/or sell copies of the Software, and to permit persons to whom the 
!   Software is furnished to do so, subject to the following conditions:
!
!   The above copyright notice and this permission notice shall be 
!   included in all copies or substantial portions of the Software.
!
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!   ACKNOWLEDGMENT
!
!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!   As per the ParaMonte library license agreement terms, if you use any parts of 
!   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
!   work (education/research/industry/development/...) by citing the ParaMonte 
!   library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program testMPI

    use iso_fortran_env, only: compiler_version, IK => int32
    use mpi

    implicit none

    integer(IK) :: ierrMPI, imageID, imageCount
    logical     :: isInitialized
    logical     :: isFinalized

    call mpi_initialized( isInitialized, ierrMPI )
    if (.not. isInitialized) call mpi_init(ierrMPI)
    call mpi_comm_rank(mpi_comm_world, imageID, ierrMPI)
    call mpi_comm_size(mpi_comm_world, imageCount, ierrMPI)
    call mpi_finalized( isFinalized, ierrMPI )
    if (.not. isFinalized) then
        call mpi_barrier(mpi_comm_world,ierrMPI)
        call mpi_finalize(ierrMPI)
    end if

end program testMPI
