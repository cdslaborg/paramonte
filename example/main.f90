!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

! Keep in mind that Fortran is case-insensitive, except for character values and string values. 
! So, feel free to use upper-case or lower-case for the Fortran syntax and entity names.
! The ParaMonte library uses camelCase convention for variable naming.

#if defined IS_COMPATIBLE_COMPILER

    program main

    ! This is the Object-Oriented-Programming (OOP) style interface to the ParaMonte routines.
    ! This is more flexible but less portable, as its compilation requires the same compiler 
    ! brand and version with which the ParaMonte library has been built.

    use paramonte, only: ParaDRAM
    use LogFunc_mod, only: getLogFunc, NDIM

    implicit none
    type(ParaDRAM) :: pd

    call pd%runSampler  ( ndim = NDIM &
                        , getLogFunc = getLogFunc &
                        , inputFile = "./paramonte.in" &    ! this is optional argument
                        ! You can also specify simulation specifications as input arguments, like 
                        ! the following. This is possible only from the OOP interface to ParaDRAM.
                        , greedyAdaptationCount = 0 &       ! this is optional argument
                        , description = "an example run" &  ! this is optional argument
                        ! More optional arguments can appear here. 
                        ! See the ParaDRAM routine's list of input arguments.
                        )

    end program main

#else

    ! This is the default simple procedural interface to the ParaMonte routines.
    ! The first two arguments to the sampler routines are mandatory.
    ! The inputFile argument is optional.

    program main

    use paramonte, only: runParaDRAM
    use LogFunc_mod, only: getLogFunc, NDIM

    implicit none

    call runParaDRAM( NDIM &
                    , getLogFunc &
                    , "./paramonte.in" & ! this is optional argument
                    )

    end program main

#endif