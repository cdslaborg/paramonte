!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%  
!%  Description:
!%      +   Run the Monte Carlo sampler of the ParaMonte library given the input log-target density function `getLogFunc()`.
!%  Output:
!%      +   The simulation output files.
!%  Author:
!%      +   Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
!%  Visit:
!%      +   https://www.cdslab.org/paramonte
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!   Keep in mind that Fortran is case-insensitive, except for character values and string values. 
!   So, feel free to use upper-case or lower-case for the Fortran syntax and entity names.
!   The ParaMonte library uses camelCase convention for variable naming.

#if defined IS_COMPATIBLE_COMPILER

    program main

    ! This is the Object-Oriented-Programming (OOP) style interface to the ParaMonte routines.
    ! This is more flexible but less portable, as its compilation requires the same compiler 
    ! brand and version with which the ParaMonte library has been built.

    use paramonte, only: ParaDRAM, IK
    use LogFunc_mod, only: getLogFunc, NDIM

    implicit none
    type(ParaDRAM) :: pmpd

    call pmpd%runSampler( ndim = NDIM &
                        , getLogFunc = getLogFunc &
                        , inputFile = "./paramonte.in" &    ! this is optional argument
                        ! You can also specify simulation specifications as input arguments, like 
                        ! the following. This is possible only from the OOP interface to ParaDRAM.
                        , greedyAdaptationCount = 0_IK &    ! this is optional argument
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%