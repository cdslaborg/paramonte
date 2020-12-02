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

!> \brief This module contains the mathematical and programming constants.
!> @author Amir Shahmoradi

module Constants_mod

    use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
    use, intrinsic :: iso_c_binding, only: CIK => c_int32_t, CRK => c_double
#if defined CFI_ENABLED
    use, intrinsic :: iso_c_binding, only: IK => c_int32_t, RK => c_double
#else
    use, intrinsic :: iso_fortran_env, only: IK => int32, RK => real64
#endif

    implicit none

    character(*), parameter :: MODULE_NAME = "@Constants_mod"

    ! Constants for computational accuracy

    integer     , parameter :: SPR = real32                                             !< @public single precision real kind
    integer     , parameter :: DPR = real64                                             !< @public double precision real kind
    integer     , parameter :: SPI = int32                                              !< @public single precision integer kind
    integer     , parameter :: DPI = int64                                              !< @public double precision integer kind
    integer     , parameter :: SPC = kind((1._SPR,1._SPR))                              !< @public single-precision complex kind
    integer     , parameter :: DPC = kind((1._DPR,1._DPR))                              !< @public double-precision complex kind
    integer     , parameter :: CK = kind((1._RK,1._RK))                                 !< @public complex kind
    integer     , parameter :: RKP = precision(1._RK)                                   !< @public real kind precision
    integer(IK) , parameter :: MAX_REC_LEN = 9999                                       !< @public maximum string record length

    ! Mathematical constants

    real(RK)    , parameter :: PI = 3.141592653589793238462643383279502884197_DPR       !< @public = acos(-1._RK) : The irrational number Pi.
    real(RK)    , parameter :: INVPI = 1._DPR / PI                                      !< @public = inverse Pi.
    real(RK)    , PARAMETER :: TWOPI = 6.283185307179586476925286766559005768394_DPR    !< @public 2*PI
    real(RK)    , parameter :: LN2 = log(2._RK)                                         !< @public Natural Log of 2 (= 0.693147180559945_RK).
    real(RK)    , parameter :: INVLN2 = 1._RK / LN2                                     !< @public Inverse of the natural Log of 2 (= 0.693147180559945_RK).
    real(RK)    , parameter :: LN10 = log(1.e1_RK)                                      !< @public Natural Log of 10 (= 2.302585092994046_RK).
    real(RK)    , parameter :: LNPI = log(PI)                                           !< @public ln(2pi) (= 1.144729885849400_RK)
    real(RK)    , parameter :: LN2PI = log(2._RK*PI)                                    !< @public ln(2pi) (= 1.837877066409345_RK)
    real(RK)    , parameter :: SQRT2 = sqrt(2._RK)                                      !< @public Square root of 2.
    real(RK)    , parameter :: NAPIER = exp(1._RK)                                      !< @public Napier number e.
    real(RK)    , parameter :: SQRTPI = sqrt(PI)                                        !< @public Square root of Pi.
    real(RK)    , parameter :: SQRT2PI = sqrt(2._RK*acos(-1._RK))                       !< @public Square root of 2Pi.
    real(RK)    , parameter :: INVSQRT2 = 1._RK / SQRT2PI                               !< @public Square root of 2.
    real(RK)    , parameter :: HALFLN2PI = 0.5_RK*LN2PI                                 !< @public ln(sqrt(2pi))
    real(RK)    , parameter :: INVSQRTPI = sqrt(INVPI)                                  !< @public = inverse of the square root of Pi.
    real(RK)    , parameter :: INVSQRT2PI = 1._RK / SQRT2PI                             !< @public 1/sqrt(2*Pi) (= 0.398942280401432_RK)
    real(RK)    , parameter :: LOGINVSQRT2PI = log(INVSQRT2PI)                          !< @public Log(1/sqrt(2Pi)), used in Gaussian distribution.
    real(RK)    , parameter :: SQRT_HALF_PI = sqrt(0.5_RK*PI)                           !< @public Square root of PI/2 (= 1.2533141373155_RK)
    real(RK)    , parameter :: LOG10NAPIER = log10(NAPIER)                              !< @public Log10 of Napier constant (= 0.434294481903259_RK).
    real(RK)    , parameter :: EPS_RK = epsilon(1._RK)                                  !< @public the smallest representable real increment (highest precision) by the machine
    real(RK)    , parameter :: HUGE_IK = huge(1_IK)                                     !< @public largest number of kind RK
    real(RK)    , parameter :: HUGE_RK = huge(1._RK)                                    !< @public largest number of kind RK
    real(RK)    , parameter :: TINY_RK = tiny(1._RK)                                    !< @public tiniest number of kind RK
    real(RK)    , parameter :: LOGHUGE_RK = log(HUGE_RK)                                !< @public log of the largest number of kind RK
    real(RK)    , parameter :: LOGTINY_RK = log(TINY_RK)                                !< @public log of the smallest number of kind RK
    real(RK)    , parameter :: POSINF_RK =  HUGE_RK / 1.e1_RK                           !< @public the division is done to avoid overflow in output
    real(RK)    , parameter :: POSINF_IK =  HUGE_IK / 2_IK                              !< @public the division is done to avoid overflow in output
    real(RK)    , parameter :: LOGINF_RK =  log(POSINF_RK)
    real(RK)    , parameter :: NEGLOGINF_RK = -LOGINF_RK
    real(RK)    , parameter :: LOGINF_IK =  log(POSINF_IK)
    real(RK)    , parameter :: NEGINF_RK = -POSINF_RK
    real(RK)    , parameter :: NEGINF_IK = -POSINF_IK
    real(RK)    , parameter :: NULL_RK = -HUGE_RK
    integer(IK) , parameter :: NULL_IK = -HUGE_IK
    character(1), parameter :: NULL_SK = achar(30)                                      ! This must remain a single character as it is assumed in multiple routines: Record separator
    character(1), parameter :: NLC = achar(10)                                          ! the New Line Character
    character(1), parameter :: TAB = achar(9)                                           ! the TAB Character
    character(*), parameter :: UNDEFINED = "UNDEFINED"

    ! null values

    type, private  :: NullType
        real(RK)     :: RK = NULL_RK
        integer(IK)  :: IK = NULL_IK
        character(1) :: SK = NULL_SK
    end type NullType
  
    type(NullType), protected :: NullVal

    ! Physical constants

    real(RK), parameter :: ERG2KEV = 6.241509125883258e8_RK                      !< @public 1 (erg) = ERG2KEV (keV)
    real(RK), parameter :: KEV2ERG = 1.60217662080000e-9_RK                      !< @public 1 (keV) = KEV2ERG (erg)
    real(RK), parameter :: LOG_ERG2KEV = log(ERG2KEV)
    real(RK), parameter :: LOG_KEV2ERG = log(KEV2ERG)

    ! Cosmological constants

    !real(RK), parameter :: LIGHT_SPEED = 3.e5_RK                                ! LIGHT_SPEED is the speed of light (Km/s).
    !real(RK), parameter :: HUBBLE_TIME_GYRS = 13.8_RK		                     ! hubble time (liddle 2003, page 57) in units of gyrs
    !real(RK), parameter :: HUBBLE_CONST = 7.1e1_RK                              ! HUBBLE_CONST is the Hubble constant in units of km/s/MPc.
    !real(RK), parameter :: LS2HC = LIGHT_SPEED / HUBBLE_CONST                   ! the speed of light in units of km/s divided by the Hubble constant.
    !real(RK), parameter :: MPC2CM = 3.09e24_RK                                  ! 1 Mega Parsec = MPC2CM centimeters.
    !real(RK), parameter :: LOG10MPC2CMSQ4PI = log10(4._RK*PI) + 2*log10(MPC2CM) ! Log10(MPC2CM centimeters.
    !real(RK), parameter :: OMEGA_DE = 0.7_RK                                    ! Dark Energy density.
    !real(RK), parameter :: OMEGA_DM = 0.3_RK                                    ! Dark Matter density.

    character(len=1), parameter :: CARRIAGE_RETURN = achar(13)
    character(len=1), parameter :: ESCAPE = achar(27)
    character(len=1), parameter :: ESC = achar(27)
    character(len=1), parameter :: CLOCK_TICK(4) = [ "|" , "/" , "-" , "\" ]

    interface getPosInf
        module procedure :: getPosInf_RK
    end interface getPosInf

    interface getNegInf
        module procedure :: getNegInf_RK
    end interface getNegInf

    ! file extentions

    type, private       :: FileType_type
        character(6)    :: binary  = "binary"
        character(6)    :: matlab  = "MATLAB"
        character(6)    :: python  = "Python"
        character(5)    :: julia   = "Julia"
        character(5)    :: ascii   = "ASCII"
        character(1)    :: rlang   = "R"
    end type FileType_type

    type, private       :: FileExt_type
        character(4)    :: binary  = ".bin"
        character(2)    :: matlab  = ".m"
        character(3)    :: python  = ".py"
        character(3)    :: julia   = ".jl"
        character(4)    :: ascii   = ".txt"
        character(2)    :: r       = ".R"
    end type FileExt_type   

    type(FileExt_type)  , parameter :: FILE_EXT = FileExt_type()
    type(FileType_type) , parameter :: FILE_TYPE = FileType_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ParaMonte methods
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, private       :: ParaMonteSamplingMethod_type
        integer(IK)     :: count = 4_IK
        character(8)    :: ParaDRAM = "ParaDRAM"
        character(8)    :: ParaDISE = "ParaDISE"
        character(8)    :: ParaHDMC = "ParaHDMC"
        character(8)    :: ParaTemp = "ParaTemp"
        character(8)    :: ParaNest = "ParaNest"
    end type ParaMonteSamplingMethod_type

    type(ParaMonteSamplingMethod_type), parameter :: PMSM = ParaMonteSamplingMethod_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getPosInf_RK() result(posInf)
        !> \brief Return IEEE positive infinity.
        !> @param[out] posInf The positive infinity.
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
        implicit none
        real(RK) :: posInf
        posInf = ieee_value(0._RK, ieee_positive_inf)
    end function getPosInf_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getNegInf_RK() result(negInf)
        !> \brief Return IEEE negative infinity.
        !> @param[out] posInf The negative infinity.
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_negative_inf
        implicit none
        real(RK) :: negInf
        negInf = ieee_value(0._RK, ieee_negative_inf)
    end function getNegInf_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Constants_mod ! LCOV_EXCL_LINE