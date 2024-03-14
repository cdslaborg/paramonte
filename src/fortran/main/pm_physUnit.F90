!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This module contains relevant physical constants.
!>
!>  \note
!>  The constants of this module are saved with the highest available `real` precision kind.<br>
!>  To use the constants at expressions involving lower-precision `real` kinds, simply convert the numbers to the desired kind
!>  via the Fortran intrinsic `real(x, kind = RKC)` where `RKC` refers to the target kind parameter used in the expression.<br>
!>
!>  \devnote
!>  This is an experimental module still under development.<br>
!>  Other similar work includes but is not limited to *Kim et al. 2017, Fortran 90 Programming With Physical Unit Annotations.*<br>
!>
!>  \finmain
!>
!>  \todo
!>  \pmed
!>  Th contents of this module should be substantially expanded.<br>
!>  Conversion procedures for converting physical quantities in various units must be added.<br>
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_physUnit

    use pm_kind, only: SK, RKB

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_physUnit"

    !>  \brief
    !>  This the derived type for constructing physical constants or unit conversions that are uncertain.
    !>
    !>  \brief
    !>  This derived minimally contains two components for holding the physical value and its associated standard deviation.
    !>
    !>  \finmain{uncertain_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type uncertain_type
        real(RKB) :: val    !<  \public The scalar component of type `real` of kind \RKB, containing the value of the physical quantity.
        real(RKB) :: std    !<  \public The scalar component of type `real` of kind \RKB, containing the standard deviation of the associated uncertainty with the reported value.
    end type

    real(RKB)               , parameter :: PI = acos(-1._RKB)                                           !<  \public The scalar constant of type `real` of kind \RKB representing the irrational number \f$\pi\f$.<br>
    real(RKB)               , parameter :: NAPIER = exp(1._RKB)                                         !<  \public The scalar constant of type `real` of kind \RKB representing the Napier constant (a.k.a. Euler number) \f$e = \exp(1)\f$.
    real(RKB)               , parameter :: AU2METER = 149597870700._RKB                                 !<  \public The scalar constant of type `real` of kind \RKB representing one **astronomical unit** unit of length in meters.<br>
    real(RKB)               , parameter :: DEGREE2RADIAN = PI / 180._RKB                                !<  \public The scalar constant of type `real` of kind \RKB representing one degree in radians.<br>
    real(RKB)               , parameter :: ARCMIN2RADIAN = PI / 10800._RKB                              !<  \public The scalar constant of type `real` of kind \RKB representing one arcminute in radians.<br>
    real(RKB)               , parameter :: ARCSEC2RADIAN = PI / 648000._RKB                             !<  \public The scalar constant of type `real` of kind \RKB representing one arcsecond in radians.<br>
    type(uncertain_type)    , parameter :: DALTON2KG = uncertain_type(1.660539040e-27_RKB, 2.e-35_RKB)  !<  \public The scalar constant of type `real` of kind \RKB representing one Dalton in kilograms.<br>
    real(RKB)               , parameter :: EV2JOULES = 1.602176634e-19_RKB                              !<  \public The scalar constant of type `real` of kind \RKB representing one electronvolt in Joules.<br>
                                                                                                        !!  An electronvolt is the amount of kinetic energy gained or lost by a single electron accelerating from rest
                                                                                                        !!  through an electric potential difference of one volt in vacuum. Hence, it has a value of one volt, \f$1 J/C\f$,
                                                                                                        !!  multiplied by the elementary charge \f$e = 1.602176634\times10−19 C\f$.<br>
                                                                                                        !!  Therefore, one electronvolt is equal to \f$1.602176634×10−19 J\f$.<br>
                                                                                                        !!  The electronvolt (\f$\ms{eV}\f$) is a unit of energy, but is not an SI unit.<br>
    real(RKB)               , parameter :: JOULES2ERGS = 1.e+7_RKB                                      !<  \public The scalar constant of type `real` of kind \RKB representing one Joules of energy in ergs.<br>
    real(RKB)               , parameter :: ERGS2JOULES = 1.e-7_RKB                                      !<  \public The scalar constant of type `real` of kind \RKB representing one ergs of energy in Joules.<br>
    real(RKB)               , parameter :: KEV2ERGS = 1000 * EV2JOULES * JOULES2ERGS                    !<  \public The scalar constant of type `real` of kind \RKB representing one kilo-electronvolts of energy in ergs.<br>
    real(RKB)               , parameter :: ERGS2KEV = 1._RKB / KEV2ERGS                                 !<  \public The scalar constant of type `real` of kind \RKB representing one ergs of energy in kilo-electronvolts.<br>

end module pm_physUnit ! LCOV_EXCL_LINE