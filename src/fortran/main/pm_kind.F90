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
!>  This module defines the relevant Fortran kind type-parameters frequently used in the ParaMonte library
!>  for the two standard supported Fortran and C-Fortran Interoperation (**CFI**) modes.<br>
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_kind

    use, intrinsic :: iso_fortran_env   , only: character_kinds
    use, intrinsic :: iso_fortran_env   , only: integer_kinds
    use, intrinsic :: iso_fortran_env   , only: logical_kinds
    use, intrinsic :: iso_fortran_env   , only: real_kinds

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Specific (interoperable) kind type parameters of the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CFI_ENABLED

    !>  \cond excluded
    use, intrinsic :: iso_c_binding     , only:  IK64_DEF => c_int64_t
    use, intrinsic :: iso_c_binding     , only:  IK32_DEF => c_int32_t
    use, intrinsic :: iso_c_binding     , only:  IK16_DEF => c_int16_t
    use, intrinsic :: iso_c_binding     , only:   IK8_DEF => c_int8_t
    use, intrinsic :: iso_c_binding     , only: RK128_DEF => c_long_double
    use, intrinsic :: iso_c_binding     , only:  RK64_DEF => c_double
    use, intrinsic :: iso_c_binding     , only:  RK32_DEF => c_float
    use, intrinsic :: iso_c_binding     , only: CK128_DEF => c_long_double_complex
    use, intrinsic :: iso_c_binding     , only:  CK64_DEF => c_double_complex
    use, intrinsic :: iso_c_binding     , only:  CK32_DEF => c_float_complex
    use, intrinsic :: iso_c_binding     , only:    SK_DEF => c_char
    use, intrinsic :: iso_c_binding     , only:    IK_DEF => c_int32_t
   !use, intrinsic :: iso_c_binding     , only:    LK_DEF => c_bool ! c_bool is one byte, which is problematic for some Fortran intrinsics like `execute_command_line()`.
    use, intrinsic :: iso_c_binding     , only:    CK_DEF => c_double_complex
    use, intrinsic :: iso_c_binding     , only:    RK_DEF => c_double
    !>  \endcond excluded

#else

    use, intrinsic :: iso_fortran_env   , only:  IK64_DEF => int64
    use, intrinsic :: iso_fortran_env   , only:  IK32_DEF => int32
    use, intrinsic :: iso_fortran_env   , only:  IK16_DEF => int16
    use, intrinsic :: iso_fortran_env   , only:   IK8_DEF => int8
    use, intrinsic :: iso_fortran_env   , only: RK128_DEF => real128
    use, intrinsic :: iso_fortran_env   , only:  RK64_DEF => real64
    use, intrinsic :: iso_fortran_env   , only:  RK32_DEF => real32
    use, intrinsic :: iso_fortran_env   , only: CK128_DEF => real128
    use, intrinsic :: iso_fortran_env   , only:  CK64_DEF => real64
    use, intrinsic :: iso_fortran_env   , only:  CK32_DEF => real32

    integer     , parameter :: SK_DEF   = kind("a")
    integer     , parameter :: IK_DEF   = kind(1) ! IK64_DEF
    integer     , parameter :: CK_DEF   = kind(1.d0)
    integer     , parameter :: RK_DEF   = kind(1.d0)

#endif

    integer     , parameter :: LK_DEF   = kind(.true.)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Supported kind type parameters of the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    integer , parameter :: SK5 = character_kinds(SK5_ENABLED)
#else
    integer , parameter :: SK5 = -1
#endif
#if SK4_ENABLED
    integer , parameter :: SK4 = character_kinds(SK4_ENABLED)
#else
    integer , parameter :: SK4 = -1
#endif
#if SK3_ENABLED
    integer , parameter :: SK3 = character_kinds(SK3_ENABLED)
#else
    integer , parameter :: SK3 = -1
#endif
#if SK2_ENABLED
    integer , parameter :: SK2 = character_kinds(SK2_ENABLED)
#else
    integer , parameter :: SK2 = -1
#endif
#if SK1_ENABLED
    integer , parameter :: SK1 = character_kinds(SK1_ENABLED)
#else
    integer , parameter :: SK1 = -1
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    integer , parameter :: IK5 = integer_kinds(IK5_ENABLED)
    integer , parameter :: IR5 = range(0_IK5)
#else
    integer , parameter :: IK5 = -1
    integer , parameter :: IR5 = -1
#endif
#if IK4_ENABLED
    integer , parameter :: IK4 = integer_kinds(IK4_ENABLED)
    integer , parameter :: IR4 = range(0_IK4)
#else
    integer , parameter :: IK4 = -1
    integer , parameter :: IR4 = -1
#endif
#if IK3_ENABLED
    integer , parameter :: IK3 = integer_kinds(IK3_ENABLED)
    integer , parameter :: IR3 = range(0_IK3)
#else
    integer , parameter :: IK3 = -1
    integer , parameter :: IR3 = -1
#endif
#if IK2_ENABLED
    integer , parameter :: IK2 = integer_kinds(IK2_ENABLED)
    integer , parameter :: IR2 = range(0_IK2)
#else
    integer , parameter :: IK2 = -1
    integer , parameter :: IR2 = -1
#endif
#if IK1_ENABLED
    integer , parameter :: IK1 = integer_kinds(IK1_ENABLED)
    integer , parameter :: IR1 = range(0_IK1)
#else
    integer , parameter :: IK1 = -1
    integer , parameter :: IR1 = -1
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    integer , parameter :: LK5 = logical_kinds(LK5_ENABLED)
#else
    integer , parameter :: LK5 = -1
#endif
#if LK4_ENABLED
    integer , parameter :: LK4 = logical_kinds(LK4_ENABLED)
#else
    integer , parameter :: LK4 = -1
#endif
#if LK3_ENABLED
    integer , parameter :: LK3 = logical_kinds(LK3_ENABLED)
#else
    integer , parameter :: LK3 = -1
#endif
#if LK2_ENABLED
    integer , parameter :: LK2 = logical_kinds(LK2_ENABLED)
#else
    integer , parameter :: LK2 = -1
#endif
#if LK1_ENABLED
    integer , parameter :: LK1 = logical_kinds(LK1_ENABLED)
#else
    integer , parameter :: LK1 = -1
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    integer , parameter :: CK5 = real_kinds(CK5_ENABLED)
    integer , parameter :: CP5 = precision(0._CK5)
    !>  \cond excluded
#else
    integer , parameter :: CK5 = -1
    integer , parameter :: CP5 = 0
    !>  \endcond excluded
#endif
#if CK4_ENABLED
    integer , parameter :: CK4 = real_kinds(CK4_ENABLED)
    integer , parameter :: CP4 = precision(0._CK4)
    !>  \cond excluded
#else
    integer , parameter :: CK4 = -1
    integer , parameter :: CP4 = 0
    !>  \endcond excluded
#endif
#if CK3_ENABLED
    integer , parameter :: CK3 = real_kinds(CK3_ENABLED)
    integer , parameter :: CP3 = precision(0._CK3)
    !>  \cond excluded
#else
    integer , parameter :: CK3 = -1
    integer , parameter :: CP3 = 0
    !>  \endcond excluded
#endif
#if CK2_ENABLED
    integer , parameter :: CK2 = real_kinds(CK2_ENABLED)
    integer , parameter :: CP2 = precision(0._CK2)
    !>  \cond excluded
#else
    integer , parameter :: CK2 = -1
    integer , parameter :: CP2 = 0
    !>  \endcond excluded
#endif
#if CK1_ENABLED
    integer , parameter :: CK1 = real_kinds(CK1_ENABLED)
    integer , parameter :: CP1 = precision(0._CK1)
    !>  \cond excluded
#else
    integer , parameter :: CK1 = -1
    integer , parameter :: CP1 = 0
    !>  \endcond excluded
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    integer , parameter :: RK5 = real_kinds(RK5_ENABLED)
    integer , parameter :: RP5 = precision(0._RK5)
    !>  \cond excluded
#else
    integer , parameter :: RK5 = -1
    integer , parameter :: RP5 = 0
    !>  \endcond excluded
#endif
#if RK4_ENABLED
    integer , parameter :: RK4 = real_kinds(RK4_ENABLED)
    integer , parameter :: RP4 = precision(0._RK4)
    !>  \cond excluded
#else
    integer , parameter :: RK4 = -1
    integer , parameter :: RP4 = 0
    !>  \endcond excluded
#endif
#if RK3_ENABLED
    integer , parameter :: RK3 = real_kinds(RK3_ENABLED)
    integer , parameter :: RP3 = precision(0._RK3)
    !>  \cond excluded
#else
    integer , parameter :: RK3 = -1
    integer , parameter :: RP3 = 0
    !>  \endcond excluded
#endif
#if RK2_ENABLED
    integer , parameter :: RK2 = real_kinds(RK2_ENABLED)
    integer , parameter :: RP2 = precision(0._RK2)
    !>  \cond excluded
#else
    integer , parameter :: RK2 = -1
    integer , parameter :: RP2 = 0
    !>  \endcond excluded
#endif
#if RK1_ENABLED
    integer , parameter :: RK1 = real_kinds(RK1_ENABLED)
    integer , parameter :: RP1 = precision(0._RK1)
    !>  \cond excluded
#else
    integer , parameter :: RK1 = -1
    integer , parameter :: RP1 = 0
    !>  \endcond excluded
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Default kind type parameters of the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    integer     , parameter ::  SK      =    SK_DEF !<  \public The default `character` kind in the ParaMonte library: `kind("a")`      in Fortran, `c_char`            in \CFI mode.
    integer     , parameter ::  IK      =    IK_DEF !<  \public The default `integer`   kind in the ParaMonte library: `int32`          in Fortran, `c_int32_t`         in \CFI mode.
    integer     , parameter ::  LK      =    LK_DEF !<  \public The default `logical`   kind in the ParaMonte library: `kind(.true.)`   in Fortran, `kind(.true.)`      in \CFI mode.
    integer     , parameter ::  CK      =    CK_DEF !<  \public The default `complex`   kind in the ParaMonte library: `real64`         in Fortran, `c_double_complex`  in \CFI mode.
    integer     , parameter ::  RK      =    RK_DEF !<  \public The default `real`      kind in the ParaMonte library: `real64`         in Fortran, `c_double`          in \CFI mode.

    integer     , parameter ::  IK64    =  IK64_DEF !<  \public The `integer`   kind for a  64-bits container.
    integer     , parameter ::  IK32    =  IK32_DEF !<  \public The `integer`   kind for a  32-bits container.
    integer     , parameter ::  IK16    =  IK16_DEF !<  \public The `integer`   kind for a  16-bits container.
    integer     , parameter ::  IK8     =   IK8_DEF !<  \public The `integer`   kind for an  8-bits container.
    integer     , parameter ::  CK128   = CK128_DEF !<  \public The `complex`   kind for a 128-bits container.
    integer     , parameter ::  CK64    =  CK64_DEF !<  \public The `complex`   kind for a  64-bits container.
    integer     , parameter ::  CK32    =  CK32_DEF !<  \public The `complex`   kind for a  32-bits container.
    integer     , parameter ::  RK128   = RK128_DEF !<  \public The `real`      kind for a 128-bits container.
    integer     , parameter ::  RK64    =  RK64_DEF !<  \public The `real`      kind for a  64-bits container.
    integer     , parameter ::  RK32    =  RK32_DEF !<  \public The `real`      kind for a  32-bits container.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Other popular non-default kind type parameters of the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    integer     , parameter :: SKU  = selected_char_kind("iso_10646")           !<  \public The UNICODE string  kind in Fortran mode.
    integer     , parameter :: SKD  = selected_char_kind("default")             !<  \public The DEFAULT string  kind in Fortran mode.
    integer     , parameter :: SKA  = selected_char_kind("ascii")               !<  \public The ASCII   string  kind in Fortran mode.
    integer     , parameter :: IKS  = kind(1)                                   !<  \public The single-precision integer kind in Fortran mode. On most platforms, this is a 32-bit integer kind.
    integer     , parameter :: IKD  = selected_int_kind(18)                     !<  \public The double precision integer kind in Fortran mode. On most platforms, this is a 64-bit integer kind.
    integer     , parameter :: IKQ  = selected_int_kind(36)                     !<  \public The quadro-precision integer kind in Fortran mode. On most platforms, this is a 128-bit integer kind.
                                                                                !!          There is no guarantee of the existence of this kind (e.g., Intel does not support this).
    integer     , parameter :: RKS  = kind(1.)                                  !<  \public The single-precision real kind in Fortran mode. On most platforms, this is an 32-bit real kind.
    integer     , parameter :: RKD  = kind(1.d0)                                !<  \public The double precision real kind in Fortran mode. On most platforms, this is an 64-bit real kind.
    integer     , parameter :: RKQ  = selected_real_kind(2*precision(1._RKD))   !<  \public The quadro-precision real kind in Fortran mode. On most platforms, this is an 128-bit real kind.
    integer     , parameter :: CKS  = RKS                                       !<  \public The single-precision complex kind in Fortran mode. On most platforms, this is a 32-bit real kind.
    integer     , parameter :: CKD  = RKD                                       !<  \public The double precision complex kind in Fortran mode. On most platforms, this is a 64-bit real kind.
    integer     , parameter :: CKQ  = RKQ                                       !<  \public The quadro-precision complex kind in Fortran mode. On most platforms, this is a 128-bit real kind.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Define the vector of all type kinds supported by the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `integer` vector containing all defined `character` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `character` kinds supported by the processor.<br>
    !>
    !>  \finmain{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: SKALL(*) = pack([SK1, SK2, SK3, SK4, SK5], 0 < [SK1, SK2, SK3, SK4, SK5])

    !>  \brief
    !>  The `integer` vector containing all defined `integer` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `integer` kinds supported by the processor.<br>
    !>
    !>  \finmain{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: IKALL(*) = pack([IK1, IK2, IK3, IK4, IK5], 0 < [IK1, IK2, IK3, IK4, IK5])

    !>  \brief
    !>  The `integer` vector containing all defined `logical` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `logical` kinds supported by the processor.<br>
    !>
    !>  \finmain{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: LKALL(*) = pack([LK1, LK2, LK3, LK4, LK5], 0 < [LK1, LK2, LK3, LK4, LK5])

    !>  \brief
    !>  The `integer` vector containing all defined `complex` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `complex` kinds supported by the processor.<br>
    !>
    !>  \finmain{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: CKALL(*) = pack([CK1, CK2, CK3, CK4, CK5], 0 < [CK1, CK2, CK3, CK4, CK5])

    !>  \brief
    !>  The `integer` vector containing all defined `real` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `real` kinds supported by the processor.<br>
    !>
    !>  \finmain{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: RKALL(*) = pack([RK1, RK2, RK3, RK4, RK5], 0 < [RK1, RK2, RK3, RK4, RK5])

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The highest and the lowest kinds supported by the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     !>  \cond excluded
     integer    , parameter , private   :: IRL = minval([IR1, IR2, IR3, IR4, IR5], mask = 0 < [IR1, IR2, IR3, IR4, IR5])
     integer    , parameter , private   :: RPL = minval([RP1, RP2, RP3, RP4, RP5], mask = 0 < [RP1, RP2, RP3, RP4, RP5])
     integer    , parameter , private   :: IRH = maxval([IR1, IR2, IR3, IR4, IR5])
     integer    , parameter , private   :: RPH = maxval([RP1, RP2, RP3, RP4, RP5])
     !>  \endcond excluded

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest range `integer` kind available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \IKL is the same as the value of \IKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the lowest range `integer` kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-range `integer` kind \IKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-range `integer` kind \IKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-range `integer` kind of the library \IKL, the same does not hold for \IKW when its value is different from \IKL.<br>
    !>
    !>  \finmain{IKL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: IKL  = selected_int_kind(IRL)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest-precision `real` kind available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \RKL is the same as the value of \RKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the lowest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-precision `real` kind \RKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-precision `real` kind \RKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-precision `real` kind of the library \RKL, the same does not hold for \RKW when its value is different from \RKL.<br>
    !>
    !>  \finmain{RKL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: RKL  = selected_real_kind(RPL)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest-precision `complex` kind available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \CKL is the same as the value of \CKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the lowest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-precision `complex` kind \CKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-precision `complex` kind \CKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-precision `complex` kind of the library \CKL, the same does not hold for \CKW when its value is different from \CKL.<br>
    !>
    !>  \finmain{CKL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: CKL  = selected_real_kind(RPL)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest range `integer` kind available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \IKH is the same as the value of \IKB under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the highest range `integer` kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-range `integer` kind \IKH <b>supported by a specific library build</b>  is not necessarily the same as the best-range `integer` kind \IKB <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-range `integer` kind of the library \IKH, the same does not hold for \IKB when its value is different from \IKH.<br>
    !>
    !>  \finmain{IKH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: IKH  = selected_int_kind(IRH)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest-precision `real` kind available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \RKH is the same as the value of \RKB under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the highest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-precision `real` kind \RKH <b>supported by a specific library build</b>  is not necessarily the same as the best-precision `real` kind \RKB <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-precision `real` kind of the library \RKH, the same does not hold for \RKB when its value is different from \RKH.<br>
    !>
    !>  \finmain{RKH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: RKH  = selected_real_kind(RPH)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest-precision `complex` kind available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \CKH is the same as the value of \CKB under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the highest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-precision `complex` kind \CKH <b>supported by a specific library build</b>  is not necessarily the same as the best-precision `complex` kind \CKB <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-precision `complex` kind of the library \CKH, the same does not hold for \CKB when its value is different from \CKH.<br>
    !>
    !>  \finmain{CKH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: CKH  = selected_real_kind(RPH)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The highest and the lowest (the best and the worst) kinds supported by the processor.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>W</b>orst-range `integer` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \IKL is the same as the value of \IKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `integer` kind type parameters that exclude the lowest-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-range `integer` kind \IKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-range `integer` kind \IKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-range `integer` kind of the library \IKL, the same does not hold for \IKW when its value is different from \IKL.<br>
    !>
    !>  \finmain{IKW}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: IKW  = selected_int_kind(1)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>W</b>orst-precision `real` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \RKL is the same as the value of \RKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the lowest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-precision `real` kind \RKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-precision `real` kind \RKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-precision `real` kind of the library \RKL, the same does not hold for \RKW when its value is different from \RKL.<br>
    !>
    !>  \finmain{RKW}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: RKW  = selected_real_kind(1)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>W</b>orst-precision `complex` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \CKL is the same as the value of \CKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the lowest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-precision `complex` kind \CKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-precision `complex` kind \CKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-precision `complex` kind of the library \CKL, the same does not hold for \CKW when its value is different from \CKL.<br>
    !>
    !>  \finmain{CKW}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: CKW  = selected_real_kind(1)

    !>  \cond excluded
    integer     , parameter , private   :: IKB_VEC_RAW(*) = [ selected_int_kind(IRH + 9 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 8 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 4 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 8 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 8 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 7 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 6 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 5 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 4 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 3 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 2 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 1 * range(int(0, IKW))) &
                                                            , selected_int_kind(IRH + 0 * range(int(0, IKW))) &
                                                            ]
    integer     , parameter , private   :: RKB_VEC_RAW(*) = [ selected_real_kind(RPH + 9 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 8 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 4 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 8 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 8 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 7 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 6 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 5 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 4 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 3 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 2 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 1 * precision(real(0, RKW))) &
                                                            , selected_real_kind(RPH + 0 * precision(real(0, RKW))) &
                                                            ]
    integer     , parameter , private   :: IKB_VEC(*) = pack(IKB_VEC_RAW, 0 <= IKB_VEC_RAW)
    integer     , parameter , private   :: RKB_VEC(*) = pack(RKB_VEC_RAW, 0 <= RKB_VEC_RAW)
    integer     , parameter , private   :: CKB_VEC(*) = RKB_VEC
    !>  \endcond excluded

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>B</b>est-range `integer` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \IKH is the same as the value of \IKB under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `integer` kind type parameters that exclude the highest-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-range `integer` kind \IKH <b>supported by a specific library build</b>  is not necessarily the same as the best-range `integer` kind \IKB <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-range `integer` kind of the library \IKH, the same does not hold for \IKB when its value is different from \IKH.<br>
    !>
    !>  The current Fortran standard does not allow automatic selection of the highest-range `integer` kind made available by the processor.<br>
    !>  However, such a kind is essential for defining `integer` constants of highest-range that can be later coerced to `integer` kinds of lower range.<br>
    !>
    !>  \finmain{IKB}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: IKB  = IKB_VEC(1) ! maxval([(selected_int_kind(i), i = 1, 1000)]) ! merge(IKH6, merge(IKH5, merge(IKH4, merge(IKH3, merge(IKH2, IKS, IKH2 > 0), IKH3 > 0), IKH4 > 0), IKH5 > 0), IKH6 > 0) ! only up to `3` times more precise than the double range kind.

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>B</b>est-precision `real` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \RKH is the same as the value of \RKB under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the highest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-precision `real` kind \RKH <b>supported by a specific library build</b>  is not necessarily the same as the best-precision `real` kind \RKB <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-precision `real` kind of the library \RKH, the same does not hold for \RKB when its value is different from \RKH.<br>
    !>
    !>  The current Fortran standard does not allow automatic selection of the highest-precision `real` kind made available by the processor.<br>
    !>  However, such a kind is essential for defining `real` constants of highest-precision that can be later coerced to `real` kinds of lower precision.<br>
    !>
    !>  \finmain{RKB}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: RKB  = RKB_VEC(1) ! maxval([(selected_real_kind(i), i = 1, 1000)]) ! merge(RKH6, merge(RKH5, merge(RKH4, merge(RKH3, merge(RKH2, RKS, RKH2 > 0), RKH3 > 0), RKH4 > 0), RKH5 > 0), RKH6 > 0) ! only up to `3` times more precise than the double precision kind.

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>B</b>est-precision `complex` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \CKH is the same as the value of \CKB under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the highest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-precision `complex` kind \CKH <b>supported by a specific library build</b>  is not necessarily the same as the best-precision `complex` kind \CKB <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-precision `complex` kind of the library \CKH, the same does not hold for \CKB when its value is different from \CKH.<br>
    !>
    !>  The current Fortran standard does not allow automatic selection of the highest-precision `complex` kind made available by the processor.<br>
    !>  However, such a kind is essential for defining `complex` constants of highest-precision that can be later coerced to `complex` kinds of lower precision.<br>
    !>
    !>  \finmain{CKB}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    integer     , parameter             :: CKB  = RKB

    ! Constants for kind precisions.<br>

    !integer     , parameter :: RKP     = precision(1._RK)              !<  \public The precision corresponding to the real kind \RK.<br>
    !integer     , parameter :: RK32P   = precision(1._RK32)            !<  \public The precision corresponding to the real kind \RK1.<br>
    !integer     , parameter :: RK64P   = precision(1._RK64)            !<  \public The precision corresponding to the real kind \RK2.<br>
    !integer     , parameter :: RK128P  = precision(1._RK128)           !<  \public The precision corresponding to the real kind \RK3.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract` derived type for creating objects of class [model_type](@ref pm_kind::model_type)
    !>  that contain the characteristics of the processor representation model used for the requested data object.<br>
    !>
    !>  \details
    !>  The Fortran standard recognizes three data representation models:<br>
    !>
    !>  Integer bit model
    !>  -----------------
    !>
    !>  The Fortran standard (2023) defines intrinsic procedures for manipulating **bits held within integers**.<br>
    !>  They are based on a **model** in which an integer holds \f$s\f$ bits \f$\{w_k : k = 0, \ldots, s - 1\}\f$,
    !>  in a sequence **from right to left**, based on the non-negative value,<br>
    !>  \f{equation}{
    !>      \sum^{s - 1}_{k = 0} w_k \times 2^k ~.
    !>  \f}
    !>  This model is valid only in the context of the standard bit-manipulation intrinsic procedures.<br>
    !>  It is identical to the model for `integer` data object described below when radix \f$r = 2\f$ and \f$w_{s - 1} = 0\f$ (corresponding to the integer-sign bit).<br>
    !>  However, when \f$r /= 2\f$ or when \f$w_{s - 1} = 1\f$ the models do not correspond, and the value expressed as an integer may vary from processor to processor.<br>
    !>
    !>  Integer data model
    !>  ------------------
    !>
    !>  The standard Fortran numeric inquiry and manipulation functions are defined in
    !>  terms of a **model set of integers** for each kind of integer data type implemented.<br>
    !>  For each kind of `integer` type, the model set is,
    !>  \f{equation}{
    !>      i = s \times \sum^{q}_{k = 1} w_k \times r^{k - 1} ~,
    !>  \f}
    !>  where,
    !>  <ol>
    !>      <li>    \f$s\f$ is \f$\pm1\f$, representing the **sign** (returned by the intrinsic `sing()`),<br>
    !>      <li>    \f$r\f$ is an integer exceeding \f$1\f$ (usually \f$2\f$), representing the **radix** (returned by the intrinsic `radix()`),<br>
    !>      <li>    \f$q\f$ is a positive integer, representing the number of **significant (equivalent to `real` mantissa) bits** (returned by the intrinsic `digits()`),<br>
    !>      <li>    \f$w_k\f$ is an integer in the range \f$0\leq w_k < r\f$.<br>
    !>  </ol>
    !>
    !>  Real data model
    !>  ---------------
    !>
    !>  The standard Fortran numeric inquiry and manipulation functions are defined in
    !>  terms of a **model set of integers** for each kind of real data type implemented.<br>
    !>  For each kind of `real` type, the model set is,
    !>  \f{equation}{
    !>      x =
    !>      \begin{cases}
    !>          0, \\
    !>          s \times b^e \times \sum^{p}_{k = 1} f_k \times b^{-k} ~,
    !>      \end{cases}
    !>  \f}
    !>  where,
    !>  <ol>
    !>      <li>    \f$s\f$ is \f$\pm1\f$, representing the **sign** (returned by the intrinsic `sing()`),<br>
    !>      <li>    \f$b\f$ is an integer exceeding \f$1\f$, representing the **radix** (returned by the intrinsic `radix()`),<br>
    !>      <li>    \f$p\f$ is an integer exceeding \f$1\f$, representing the number of **significant (mantissa) bits** (returned by the intrinsic `digits()`),<br>
    !>      <li>    \f$e\f$ is an integer in the range \f$e_{min} \leq e \leq e_{max}\f$, representing the **exponent** of the `real` value (returned by the intrinsic `exponent()`),<br>
    !>      <li>    \f$f_k\f$ is an integer in the range \f$0\leq f_k < b\f$ except that \f$f_1\f$ is nonzero.<br>
    !>  </ol>
    !>  Values of the parameters in these models are chosen for the processor so as best to fit the hardware with the proviso that all model numbers are representable.<br>
    !>  It is likely that there are some machine numbers that lie outside the model.<br>
    !>  For example, many computers represent the integer \f$-r^q\f$, and the IEEE standard for floating-point arithmetic (`ISO/IEC/IEEE 60559:2011`)
    !>  contains reals with \f$f_1 = 0\f$ (called **subnormal numbers**) and register numbers with increased precision and range.<br>
    !>
    !>  \interface{model_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: model_type, modeli_type, modelb_type, modelr_type
    !>      class(model_type), allocatable :: model
    !>
    !>      model = modelr_type(mold)
    !>      model = modelb_type(mold)
    !>      model = modeli_type(mold)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [model_type](@ref pm_kind::model_type)<br>
    !>  [modeln_type](@ref pm_kind::modeln_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>  [modelb_type](@ref pm_kind::modelb_type)<br>
    !>  [modelr_type](@ref pm_kind::modelr_type)<br>
    !>
    !>  \finmain{model_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract :: model_type
        integer(IK) :: kind         !<  \public The scalar `integer` of default kind \IK, representing the processor-dependent `kind` type-parameter of data type of interest.<br>
        integer(IK) :: storage_size !<  \public The scalar `integer` of default kind \IK, whose value is the size, in bits,
                                    !!          that would be taken in memory by a scalar of the same type and kind as the data type of interest.<br>
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract` derived type for creating objects of class [modeln_type](@ref pm_kind::modeln_type)
    !>  that contain the characteristics of the processor representation model used for the requested **numeric** data object.<br>
    !>
    !>  \interface{modeln_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: model_type, modeln_type, modeli_type, modelb_type, modelr_type
    !>      class(modeln_type), allocatable :: modelNum
    !>      class(model_type), allocatable :: model
    !>
    !>      modelNum = modelr_type(mold)
    !>      modelNum = modelb_type(mold)
    !>      modelNum = modeli_type(mold)
    !>
    !>      model = modelr_type(mold)
    !>      model = modelb_type(mold)
    !>      model = modeli_type(mold)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [model_type](@ref pm_kind::model_type)<br>
    !>  [modeln_type](@ref pm_kind::modeln_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>  [modelb_type](@ref pm_kind::modelb_type)<br>
    !>  [modelr_type](@ref pm_kind::modelr_type)<br>
    !>
    !>  \finmain{modeln_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract, extends(model_type) :: modeln_type
        integer(IK) :: digits   !<  \public The scalar `integer` of default kind \IK, whose value is the number of significant bit digits in the model that includes the numeric data type of interest,
                                !!          that is `q` or `p` in the `integer` and `real` models detailed in [model_type](@ref pm_kind::model_type).<br>
        integer(IK) :: radix    !<  \public The scalar `integer` of default kind \IK, whose value is the base in the model that includes the numeric data of interest,
                                !!          that is `r` or `b` in the `integer` and `real` models detailed in [model_type](@ref pm_kind::model_type).<br>
        integer(IK) :: range    !<  \public The scalar `integer` of default kind \IK, holding the equivalent decimal exponent range in the models representing `integer` or `real` data type of interest.<br>
                                !!          The value is `int(log10(huge))` for `integer` data type and `int(min(log10(huge), -log10(tiny)))` for `real` data type,
                                !!          where `huge` and `tiny` are the largest and smallest positive numbers in the corresponding model.<br>
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract` derived type for creating objects of class [modeli_type](@ref pm_kind::modeli_type)
    !>  that contain the characteristics of the processor representation model used for the requested `integer` data object.<br>
    !>
    !>  \details
    !>  See the documentations of [model_type](@ref pm_kind::model_type) for details of model sets for supported data types.<br>
    !>
    !>  \param[in]  mold    :   The input scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                          </ol>
    !>                          whose type and kind will determine the characteristics of the output model object.<br>
    !>                          The specific value of the input `mold` is irrelevant and ignored.<br>
    !>
    !>  \return
    !>  `model`             :   The output scalar or array of type,<br>
    !>                          <ol>
    !>                              <li>    [modeli_type](@ref pm_kind::modeli_type),<br>
    !>                          </ol>
    !>                          containing the characteristics of the processor representation model used for the input data entity.<br>
    !>
    !>  \interface{modeli_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: model_type, modeln_type, modeli_type
    !>      class(modeln_type), allocatable :: modelNum
    !>      class(model_type), allocatable :: model
    !>      type(modeli_type) :: modelInt
    !>
    !>      modelInt = modeli_type(mold)
    !>      modelNum = modeli_type(mold)
    !>      model = modeli_type(mold)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [model_type](@ref pm_kind::model_type)<br>
    !>  [modelb_type](@ref pm_kind::modelb_type)<br>
    !>  [modeln_type](@ref pm_kind::modeln_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>
    !>  \example{modeli_type}
    !>  \include{lineno} example/pm_kind/modeli_type/main.F90
    !>  \compilef{modeli_type}
    !>  \output{modeli_type}
    !>  \include{lineno} example/pm_kind/modeli_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_kind](@ref test_pm_kind)
    !>
    !>  \finmain{modeli_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(modeln_type) :: modeli_type
        integer(IKH) :: huge    !<  \public The scalar `integer` of the highest kind \IKH supported by the ParaMonte library,
                                !!          representing the largest value in the model that includes the `integer` data type of interest.<br>
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract` derived type for creating objects of class [modelb_type](@ref pm_kind::modelb_type)
    !>  that contain the characteristics of the processor representation model used for the requested **numeric** data object.<br>
    !>
    !>  \details
    !>  See the documentations of [model_type](@ref pm_kind::model_type) for details of model sets for supported data types.<br>
    !>
    !>  \param[in]  mold    :   The input scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                          </ol>
    !>                          whose type and kind will determine the characteristics of the output model object.<br>
    !>                          The specific value of the input `mold` is irrelevant and ignored.<br>
    !>
    !>  \return
    !>  `model`             :   The output scalar or array of type,<br>
    !>                          <ol>
    !>                              <li>    [modelb_type](@ref pm_kind::modelb_type),<br>
    !>                          </ol>
    !>                          containing the characteristics of the processor representation model used for the input data entity.<br>
    !>
    !>  \interface{modelb_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: model_type, modeln_type, modeli_type, modelb_type
    !>      class(modelb_type), allocatable :: modelBit
    !>      class(modeli_type), allocatable :: modelInt
    !>      class(modeln_type), allocatable :: modelNum
    !>      class(model_type), allocatable :: model
    !>
    !>      modelBit = modelb_type(mold)
    !>      modelInt = modelb_type(mold)
    !>      modelNum = modelb_type(mold)
    !>      model = modelb_type(mold)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [model_type](@ref pm_kind::model_type)<br>
    !>  [modeln_type](@ref pm_kind::modeln_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>  [modelb_type](@ref pm_kind::modelb_type)<br>
    !>  [modelr_type](@ref pm_kind::modelr_type)<br>
    !>
    !>  \example{modelb_type}
    !>  \include{lineno} example/pm_kind/modelb_type/main.F90
    !>  \compilef{modelb_type}
    !>  \output{modelb_type}
    !>  \include{lineno} example/pm_kind/modelb_type/main.out.F90
    !>
    !>  \finmain{modelb_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(modeli_type) :: modelb_type
        integer(IK) :: bit_size !<  \public The scalar `integer` of default kind \IK, whose value is the number of bits in the model
                                !!          for bits within an integer of the same type parameter the bit (integer) data type of interest.<br>
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the `abstract` derived type for creating objects of class [modelr_type](@ref pm_kind::modelr_type)
    !>  that contain the characteristics of the processor representation model used for the requested `integer` data object.<br>
    !>
    !>  \details
    !>  See the documentations of [model_type](@ref pm_kind::model_type) for details of model sets for supported data types.<br>
    !>
    !>  \param[in]  mold    :   The input scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          whose type and kind will determine the characteristics of the output model object.<br>
    !>                          The specific value of the input `mold` is irrelevant and ignored.<br>
    !>
    !>  \return
    !>  `model`             :   The output scalar or array of type,<br>
    !>                          <ol>
    !>                              <li>    [modelr_type](@ref pm_kind::modelr_type),<br>
    !>                          </ol>
    !>                          containing the characteristics of the processor representation model used for the input data entity.<br>
    !>
    !>  \interface{modelr_type}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: model_type, modeln_type, modelr_type
    !>      type(modeln_type), allocatable :: modelNum
    !>      class(model_type), allocatable :: model
    !>      type(modelr_type) :: modelReal
    !>
    !>      modelReal = modelr_type(mold)
    !>      modelNum = modelr_type(mold)
    !>      model = modelr_type(mold)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [model_type](@ref pm_kind::model_type)<br>
    !>  [modeln_type](@ref pm_kind::modeln_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>  [modelb_type](@ref pm_kind::modelb_type)<br>
    !>  [modelr_type](@ref pm_kind::modelr_type)<br>
    !>
    !>  \example{modelr_type}
    !>  \include{lineno} example/pm_kind/modelr_type/main.F90
    !>  \compilef{modelr_type}
    !>  \output{modelr_type}
    !>  \include{lineno} example/pm_kind/modelr_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_kind](@ref test_pm_kind)
    !>
    !>  \finmain{modelr_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(modeln_type) :: modelr_type
        real(RKH)   :: epsilon          !<  \public The scalar `real` of the highest kind \RKH supported by the ParaMonte library,
                                        !!          whose value is almost negligible compared with the value `1.0` in the model that includes the `real` data type of interest,
                                        !!          that is \f$b^{1-p}\f$ as detailed in the real data model of [model_type](@ref pm_kind::model_type).<br>
                                        !!          \note
                                        !!          EPSILON makes it easy to select a delta for algorithms (such as root locators) that search until the calculation is within delta of an estimate.<br>
                                        !!          If delta is too small (smaller than the decimal resolution of the data type), the algorithm might never halt.<br>
                                        !!          By scaling the value returned by EPSILON to the estimate, you obtain a delta that ensures search termination.<br>
        real(RKH)   :: huge             !<  \public The scalar `real` of the highest kind \RKH supported by the ParaMonte library,
                                        !!          representing the largest value in the model that includes the `real` data type of interest.<br>
        integer(IK) :: maxexponent      !<  \public The scalar `integer` of default kind \IK, representing the maximum exponent in the model that includes the `real` data type of interest.<br>
                                        !!          This corresponds to \f$e_{max}\f$ in the model set for the `real` data type detailed in [model_type](@ref pm_kind::model_type).<br>
        integer(IK) :: minexponent      !<  \public The scalar `integer` of default kind \IK, representing the minimum exponent in the model that includes the `real` data type of interest.<br>
                                        !!          This corresponds to \f$e_{min}\f$ in the model set for the `real` data type detailed in [model_type](@ref pm_kind::model_type).<br>
        integer(IK) :: precision        !<  \public The scalar `integer` of default kind \IK, representing the equivalent decimal precision (number of digits after the decimal point)
                                        !!          in the model representing real numbers with the same type parameter value as the `real` data type of interest.<br>
                                        !!          The value is `int((p - 1) * log10(b)) + k` where `k` is `1` if `b` is an integral power of `10` and `0` otherwise.<br>
                                        !!          The meaning of the `real` data type model parameters are detailed in [model_type](@ref pm_kind::model_type).<br>
        real(RKH)   :: tiny             !<  \public The scalar `real` of the highest kind \RKH supported by the ParaMonte library,
                                        !!          representing the smallest positive number \f$b^{e_min} - 1\f$ in the model that includes the `real` data type of interest.<br>
                                        !!          The meaning of the `real` data type model parameters are detailed in [model_type](@ref pm_kind::model_type).<br>
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    !>  \brief
    !>  Generate and return an object of class [model_type](@ref pm_kind::model_type) containing
    !>  the characteristics of the type and kind of an input scalar or array of arbitrary `rank` of type `integer`.<br>
    !>
    !>  \details
    !>  This generic interface contains the constructors of objects of type,<br>
    !>  <ol>
    !>      <li>    [modeli_type](@ref pm_kind::modeli_type)<br>
    !>  </ol>
    !>
    !>  \param[in]  mold    :   The input scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                          </ol>
    !>                          whose type and kind will determine the characteristics of the output model object.<br>
    !>                          The specific value of the input `mold` is irrelevant and ignored.<br>
    !>
    !>  \return
    !>  `model`             :   The output scalar or array of type,<br>
    !>                          <ol>
    !>                              <li>    [modeli_type](@ref pm_kind::modeli_type),<br>
    !>                          </ol>
    !>                          containing the characteristics of the processor representation model used for the input data entity.<br>
    !>
    !>  \interface{construct_modeli}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: model_type, modeln_type, modeli_type
    !>      class(modeln_type), allocatable :: modelNum
    !>      class(model_type), allocatable :: model
    !>      type(modeli_type) :: modelInt
    !>
    !>      modelInt = modeli_type(mold)
    !>      modelNum = modeli_type(mold)
    !>      model = modeli_type(mold)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [model_type](@ref pm_kind::model_type)<br>
    !>  [modelb_type](@ref pm_kind::modelb_type)<br>
    !>  [modeln_type](@ref pm_kind::modeln_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>
    !>  \example{construct_modeli}
    !>  \include{lineno} example/pm_kind/modeli_type/main.F90
    !>  \compilef{construct_modeli}
    !>  \output{construct_modeli}
    !>  \include{lineno} example/pm_kind/modeli_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_kind](@ref test_pm_kind)
    !>
    !>  \finmain{construct_modeli}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface modeli_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function construct_modeli_IK5(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modeli_IK5
#endif
        integer(IK5), intent(in)    :: mold
        type(modeli_type)           :: model
    end function
#endif

#if IK4_ENABLED
    pure elemental module function construct_modeli_IK4(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modeli_IK4
#endif
        integer(IK4), intent(in)    :: mold
        type(modeli_type)           :: model
    end function
#endif

#if IK3_ENABLED
    pure elemental module function construct_modeli_IK3(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modeli_IK3
#endif
        integer(IK3), intent(in)    :: mold
        type(modeli_type)           :: model
    end function
#endif

#if IK2_ENABLED
    pure elemental module function construct_modeli_IK2(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modeli_IK2
#endif
        integer(IK2), intent(in)    :: mold
        type(modeli_type)           :: model
    end function
#endif

#if IK1_ENABLED
    pure elemental module function construct_modeli_IK1(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modeli_IK1
#endif
        integer(IK1), intent(in)    :: mold
        type(modeli_type)           :: model
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    !>  \brief
    !>  Generate and return an object of type [modelb_type](@ref pm_kind::modelb_type) containing
    !>  the characteristics of the type and kind of an input scalar or array of arbitrary `rank` of type `integer` used as a bit storage.<br>
    !>
    !>  \details
    !>  This generic interface contains the constructors of objects of type,<br>
    !>  <ol>
    !>      <li>    [modelb_type](@ref pm_kind::modelb_type)<br>
    !>  </ol>
    !>
    !>  \param[in]  mold    :   The input scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                          </ol>
    !>                          whose type and kind will determine the characteristics of the output model object.<br>
    !>                          The specific value of the input `mold` is irrelevant and ignored.<br>
    !>
    !>  \return
    !>  `model`             :   The output scalar or array of type,<br>
    !>                          <ol>
    !>                              <li>    [modelb_type](@ref pm_kind::modelb_type),<br>
    !>                          </ol>
    !>                          containing the characteristics of the processor representation model used for the input data entity.<br>
    !>
    !>  \interface{construct_modelb}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: model_type, modeln_type, modelb_type
    !>      class(modeln_type), allocatable :: modelNum
    !>      class(model_type), allocatable :: model
    !>      type(modelb_type) :: modelBit
    !>
    !>      modelBit = modelb_type(mold)
    !>      modelNum = modelb_type(mold)
    !>      model = modelb_type(mold)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [model_type](@ref pm_kind::model_type)<br>
    !>  [modeln_type](@ref pm_kind::modeln_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>  [modelb_type](@ref pm_kind::modelb_type)<br>
    !>  [modelr_type](@ref pm_kind::modelr_type)<br>
    !>
    !>  \example{construct_modelb}
    !>  \include{lineno} example/pm_kind/modelb_type/main.F90
    !>  \compilef{construct_modelb}
    !>  \output{construct_modelb}
    !>  \include{lineno} example/pm_kind/modelb_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_kind](@ref test_pm_kind)
    !>
    !>  \finmain{construct_modelb}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface modelb_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function construct_modelb_IK5(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelb_IK5
#endif
        integer(IK5), intent(in)    :: mold
        type(modelb_type)           :: model
    end function
#endif

#if IK4_ENABLED
    pure elemental module function construct_modelb_IK4(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelb_IK4
#endif
        integer(IK4), intent(in)    :: mold
        type(modelb_type)           :: model
    end function
#endif

#if IK3_ENABLED
    pure elemental module function construct_modelb_IK3(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelb_IK3
#endif
        integer(IK3), intent(in)    :: mold
        type(modelb_type)           :: model
    end function
#endif

#if IK2_ENABLED
    pure elemental module function construct_modelb_IK2(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelb_IK2
#endif
        integer(IK2), intent(in)    :: mold
        type(modelb_type)           :: model
    end function
#endif

#if IK1_ENABLED
    pure elemental module function construct_modelb_IK1(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelb_IK1
#endif
        integer(IK1), intent(in)    :: mold
        type(modelb_type)           :: model
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    !>  \brief
    !>  Generate and return an object of class [model_type](@ref pm_kind::model_type) containing
    !>  the characteristics of the type and kind of an input scalar or array of arbitrary `rank` of type `real`.<br>
    !>
    !>  \details
    !>  This generic interface contains the constructors of objects of type,<br>
    !>  <ol>
    !>      <li>    [modelr_type](@ref pm_kind::modelr_type)<br>
    !>  </ol>
    !>
    !>  \param[in]  mold    :   The input scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          whose type and kind will determine the characteristics of the output model object.<br>
    !>                          The specific value of the input `mold` is irrelevant and ignored.<br>
    !>
    !>  \return
    !>  `model`             :   The output scalar or array of type,<br>
    !>                          <ol>
    !>                              <li>    [modelr_type](@ref pm_kind::modelr_type),<br>
    !>                          </ol>
    !>                          containing the characteristics of the processor representation model used for the input data entity.<br>
    !>
    !>  \interface{construct_modelr}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: model_type, modeln_type, modelr_type
    !>      type(modeln_type), allocatable :: modelNum
    !>      class(model_type), allocatable :: model
    !>      type(modelr_type) :: modelReal
    !>
    !>      modelReal = modelr_type(mold)
    !>      modelNum = modelr_type(mold)
    !>      model = modelr_type(mold)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [model_type](@ref pm_kind::model_type)<br>
    !>  [modeln_type](@ref pm_kind::modeln_type)<br>
    !>  [modeli_type](@ref pm_kind::modeli_type)<br>
    !>  [modelb_type](@ref pm_kind::modelb_type)<br>
    !>  [modelr_type](@ref pm_kind::modelr_type)<br>
    !>
    !>  \example{construct_modelr}
    !>  \include{lineno} example/pm_kind/modelr_type/main.F90
    !>  \compilef{construct_modelr}
    !>  \output{construct_modelr}
    !>  \include{lineno} example/pm_kind/modelr_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_kind](@ref test_pm_kind)
    !>
    !>  \finmain{construct_modelr}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface modelr_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function construct_modelr_RK5(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelr_RK5
#endif
        real(RK5)   , intent(in)    :: mold
        type(modelr_type)           :: model
    end function
#endif

#if RK4_ENABLED
    pure elemental module function construct_modelr_RK4(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelr_RK4
#endif
        real(RK4)   , intent(in)    :: mold
        type(modelr_type)           :: model
    end function
#endif

#if RK3_ENABLED
    pure elemental module function construct_modelr_RK3(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelr_RK3
#endif
        real(RK3)   , intent(in)    :: mold
        type(modelr_type)           :: model
    end function
#endif

#if RK2_ENABLED
    pure elemental module function construct_modelr_RK2(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelr_RK2
#endif
        real(RK2)   , intent(in)    :: mold
        type(modelr_type)           :: model
    end function
#endif

#if RK1_ENABLED
    pure elemental module function construct_modelr_RK1(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct_modelr_RK1
#endif
        real(RK1)   , intent(in)    :: mold
        type(modelr_type)           :: model
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_kind ! LCOV_EXCL_LINE