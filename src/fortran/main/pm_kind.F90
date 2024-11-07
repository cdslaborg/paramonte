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
!>  \details
!>
!>  Default Kind Type Parameters
!>  ----------------------------
!>
!>  The following list contains the default kind type parameters of interest to the end users of the ParaMonte library.<br>
!>  Note that the default kind type parameters may not necessarily correspond to the default kind type parameters used by the processor.<br>
!>
!>  <ol>
!>      <li>    The constant [SK](@ref pm_kind::SK) represents the default `character` kind type parameter used within the entire ParaMonte library.
!>      <li>    The constant [IK](@ref pm_kind::IK) represents the default `integer` kind type parameter used within the entire ParaMonte library.
!>      <li>    The constant [LK](@ref pm_kind::LK) represents the default `logical` kind type parameter used within the entire ParaMonte library.
!>      <li>    The constant [CK](@ref pm_kind::CK) represents the default `complex` kind type parameter used within the entire ParaMonte library.
!>      <li>    The constant [RK](@ref pm_kind::RK) represents the default `real` kind type parameter used within the entire ParaMonte library.
!>  </ol>
!>
!>  Single-Range/Precision Kind Type Parameters
!>  -------------------------------------------
!>
!>  <ol>
!>      <li>    The constant [IKS](@ref pm_kind::IKS) represents the single-range `integer` kind type parameter used within the entire ParaMonte library.<br>
!>              This corresponds to the default `integer` kind type parameter of the processor.<br>
!>      <li>    The constant [CKS](@ref pm_kind::CKS) represents the single-precision `complex` kind type parameter used within the entire ParaMonte library.<br>
!>              This corresponds to the default `real` kind type parameter of the processor.<br>
!>      <li>    The constant [RKS](@ref pm_kind::RKS) represents the single-precision `real` kind type parameter used within the entire ParaMonte library.<br>
!>              This corresponds to the default `real` kind type parameter of the processor.<br>
!>  </ol>
!>  All of the above kind type parameters are guaranteed by the standard to exist and be supported.<br>
!>
!>  Double-Range/Precision Kind Type Parameters
!>  -------------------------------------------
!>
!>  <ol>
!>      <li>    The constant [IKD](@ref pm_kind::IKD) represents the double-range `integer` kind type parameter used within the entire ParaMonte library.<br>
!>              This corresponds to the default `integer` kind type parameter of the processor.<br>
!>      <li>    The constant [CKD](@ref pm_kind::CKD) represents the double-precision `complex` kind type parameter used within the entire ParaMonte library.<br>
!>              This corresponds to the default `real` kind type parameter of the processor.<br>
!>      <li>    The constant [RKD](@ref pm_kind::RKD) represents the double-precision `real` kind type parameter used within the entire ParaMonte library.<br>
!>              This corresponds to the default `real` kind type parameter of the processor.<br>
!>  </ol>
!>  All of the above kind type parameters are guaranteed by the standard to exist and be supported.<br>
!>
!>  Quadru-Range/Precision Kind Type Parameters
!>  -------------------------------------------
!>
!>  <ol>
!>      <li>    The constant [IKQ](@ref pm_kind::IKQ) represents the quadru-range `integer` kind type parameter used within the entire ParaMonte library.<br>
!>              This corresponds to the default `integer` kind type parameter of the processor.<br>
!>      <li>    The constant [CKQ](@ref pm_kind::CKQ) represents the quadru-precision `complex` kind type parameter used within the entire ParaMonte library.<br>
!>              This corresponds to the default `real` kind type parameter of the processor.<br>
!>      <li>    The constant [RKQ](@ref pm_kind::RKQ) represents the quadru-precision `real` kind type parameter used within the entire ParaMonte library.<br>
!>              This corresponds to the default `real` kind type parameter of the processor.<br>
!>  </ol>
!>  There is no guarantee of the existence of the above kind type parameters on all processors and compilers.<br>
!>  If a kind type parameter does not exist, the processor assigns a negative value to it.<br>
!>
!>  Lowest-Range Kind Type Parameters
!>  ---------------------------------
!>
!>  This module specifies the kind type parameters that yield the lowest range intrinsic `integer`, `complex`,
!>  and `real` types that are supported by the specific ParaMonte library build as in the following:
!>
!>  <ol>
!>      <li>    The constant [IKL](@ref pm_kind::IKL) represents the lowest-range `integer` kind type parameter supported by the ParaMonte library.
!>      <li>    The constant [CKLR](@ref pm_kind::CKLR) represents the lowest-range `complex` kind type parameter supported by the ParaMonte library.
!>      <li>    The constant [RKLR](@ref pm_kind::RKLR) represents the lowest-range `real` kind type parameter supported by the ParaMonte library.
!>  </ol>
!>
!>  Highest-Range Kind Type Parameters
!>  ----------------------------------
!>
!>  This module specifies the kind type parameters that yield the highest range intrinsic `integer`, `complex`,
!>  and `real` types that are supported by the specific ParaMonte library build as in the following:
!>
!>  <ol>
!>      <li>    The constant [IKH](@ref pm_kind::IKH) represents the highest-range `integer` kind type parameter supported by the ParaMonte library.
!>      <li>    The constant [CKHR](@ref pm_kind::CKHR) represents the highest-range `complex` kind type parameter supported by the ParaMonte library.
!>      <li>    The constant [RKHR](@ref pm_kind::RKHR) represents the highest-range `real` kind type parameter supported by the ParaMonte library.
!>  </ol>
!>
!>  Lowest-Precision Kind Type Parameters
!>  -------------------------------------
!>
!>  This module specifies the kind type parameters that yield the lowest precision intrinsic `complex` and
!>  and `real` types that are supported by the specific ParaMonte library build as in the following:
!>
!>  <ol>
!>      <li>    The constant [CKL](@ref pm_kind::CKL) represents the lowest-precision `complex` kind type parameter supported by the ParaMonte library.
!>      <li>    The constant [RKL](@ref pm_kind::RKL) represents the lowest-precision `real` kind type parameter supported by the ParaMonte library.
!>  </ol>
!>
!>  Highest-Precision Kind Type Parameters
!>  --------------------------------------
!>
!>  This module specifies the kind type parameters that yield the highest precision intrinsic `complex` and
!>  and `real` types that are supported by the specific ParaMonte library build as in the following:
!>
!>  <ol>
!>      <li>    The constant [CKH](@ref pm_kind::CKH) represents the highest-precision `complex` kind type parameter supported by the ParaMonte library.
!>      <li>    The constant [RKH](@ref pm_kind::RKH) represents the highest-precision `real` kind type parameter supported by the ParaMonte library.
!>  </ol>
!>
!>  Worst-Range Kind Type Parameters
!>  --------------------------------
!>
!>  This module specifies the kind type parameters that yield the worst range intrinsic `integer`, `complex`,
!>  and `real` types that supported by the processor (whether or not supported by the specific ParaMonte library build) as in the following:
!>
!>  <ol>
!>      <li>    The constant [IKW](@ref pm_kind::IKW) represents the worst-range `integer` kind type parameter supported by the processor.
!>      <li>    The constant [CKWR](@ref pm_kind::CKWR) represents the worst-range `complex` kind type parameter supported by the processor.
!>      <li>    The constant [RKWR](@ref pm_kind::RKWR) represents the worst-range `real` kind type parameter supported by the processor.
!>  </ol>
!>
!>  Best-Range Kind Type Parameters
!>  -------------------------------
!>
!>  This module specifies the kind type parameters that yield the highest range intrinsic `integer`, `complex`,
!>  and `real` types that supported by the processor (whether or not supported by the specific ParaMonte library build) as in the following:
!>
!>  <ol>
!>      <li>    The constant [IKB](@ref pm_kind::IKB) represents the highest-range `integer` kind type parameter supported by the processor.
!>      <li>    The constant [CKBR](@ref pm_kind::CKBR) represents the highest-range `complex` kind type parameter supported by the processor.
!>      <li>    The constant [RKBR](@ref pm_kind::RKBR) represents the highest-range `real` kind type parameter supported by the processor.
!>  </ol>
!>
!>  Worst-Precision Kind Type Parameters
!>  ------------------------------------
!>
!>  This module specifies the kind type parameters that yield the worst precision intrinsic `complex` and
!>  and `real` types that supported by the processor (whether or not supported by the specific ParaMonte library build) as in the following:
!>
!>  <ol>
!>      <li>    The constant [CKW](@ref pm_kind::CKW) represents the worst-precision `complex` kind type parameter supported by the processor.
!>      <li>    The constant [RKW](@ref pm_kind::RKW) represents the worst-precision `real` kind type parameter supported by the processor.
!>  </ol>
!>
!>  Best-Precision Kind Type Parameters
!>  -----------------------------------
!>
!>  This module specifies the kind type parameters that yield the highest precision intrinsic `complex` and
!>  and `real` types that supported by the processor (whether or not supported by the specific ParaMonte library build) as in the following:
!>
!>  <ol>
!>      <li>    The constant [CKB](@ref pm_kind::CKB) represents the highest-precision `complex` kind type parameter supported by the processor.
!>      <li>    The constant [RKB](@ref pm_kind::RKB) represents the highest-precision `real` kind type parameter supported by the processor.
!>  </ol>
!>
!>  Developer Guidelines
!>  ====================
!>
!>  The current implementation of the ParaMonte library can recognize up to five kind type parameters
!>  for each of the five intrinsic types in the latest Fortran programming language standard 2023:<br>
!>  <ol>
!>      <li>    The `character` kinds are prefixed by `SK` standing for **string kind** type:
!>              <ol>
!>                  <li>    [SK5](@ref pm_kind::SK5)
!>                  <li>    [SK4](@ref pm_kind::SK4)
!>                  <li>    [SK3](@ref pm_kind::SK3)
!>                  <li>    [SK2](@ref pm_kind::SK2)
!>                  <li>    [SK1](@ref pm_kind::SK1)
!>              </ol>
!>      <li>    The `integer` kinds are prefixed by `IK` standing for **integer kind** type:
!>              <ol>
!>                  <li>    [IK5](@ref pm_kind::IK5)
!>                  <li>    [IK4](@ref pm_kind::IK4)
!>                  <li>    [IK3](@ref pm_kind::IK3)
!>                  <li>    [IK2](@ref pm_kind::IK2)
!>                  <li>    [IK1](@ref pm_kind::IK1)
!>              </ol>
!>      <li>    The `logical` kinds are prefixed by `LK` standing for **logical kind** type:
!>              <ol>
!>                  <li>    [LK5](@ref pm_kind::LK5)
!>                  <li>    [LK4](@ref pm_kind::LK4)
!>                  <li>    [LK3](@ref pm_kind::LK3)
!>                  <li>    [LK2](@ref pm_kind::LK2)
!>                  <li>    [LK1](@ref pm_kind::LK1)
!>              </ol>
!>      <li>    The `complex` kinds are prefixed by `CK` standing for **complex kind** type:
!>              <ol>
!>                  <li>    [CK5](@ref pm_kind::CK5)
!>                  <li>    [CK4](@ref pm_kind::CK4)
!>                  <li>    [CK3](@ref pm_kind::CK3)
!>                  <li>    [CK2](@ref pm_kind::CK2)
!>                  <li>    [CK1](@ref pm_kind::CK1)
!>              </ol>
!>      <li>    The `real` kinds are prefixed by `CK` standing for **real kind** type:
!>              <ol>
!>                  <li>    [RK5](@ref pm_kind::RK5)
!>                  <li>    [RK4](@ref pm_kind::RK4)
!>                  <li>    [RK3](@ref pm_kind::RK3)
!>                  <li>    [RK2](@ref pm_kind::RK2)
!>                  <li>    [RK1](@ref pm_kind::RK1)
!>              </ol>
!>  </ol>
!>
!>  A few remarks are in order:
!>  <ol>
!>      <li>    Any unsupported or non-existing kind type parameter is given a negative value.<br>
!>      <li>    Not all five kinds per intrinsic type may be available on a given platform or compiler or desired for a particular library build.<br>
!>      <li>    The availability of the above kind type parameters can be also controlled via the
!>              [optional library build flags](https://github.com/cdslaborg/paramonte/blob/main/install.config.md#TIER-4-ParaMonte-library-build-configuration-flags).<br>
!>      <li>    The above kind type parameters are considered internal library implementations.<br>
!>      <li>    An end user must never use these kind type parameters from this module.<br>
!>  </ol>
!>
!>  Precision Parameters
!>  --------------------
!>
!>  The precisions of the corresponding `real` and `complex` kind type parameters are represented by the letter `R`.
!>  <ol>
!>      <li>    `complex` precisions are prefixed by `CP` standing for **complex precision**:
!>              <ol>
!>                  <li>    [CP5](@ref pm_kind::CP5)
!>                  <li>    [CP4](@ref pm_kind::CP4)
!>                  <li>    [CP3](@ref pm_kind::CP3)
!>                  <li>    [CP2](@ref pm_kind::CP2)
!>                  <li>    [CP1](@ref pm_kind::CP1)
!>              </ol>
!>      <li>    `real` kinds are prefixed by `CK` standing for **real precision**:
!>              <ol>
!>                  <li>    [RP5](@ref pm_kind::RP5)
!>                  <li>    [RP4](@ref pm_kind::RP4)
!>                  <li>    [RP3](@ref pm_kind::RP3)
!>                  <li>    [RP2](@ref pm_kind::RP2)
!>                  <li>    [RP1](@ref pm_kind::RP1)
!>              </ol>
!>  </ol>
!>  The precision of non-existing kind type parameters are set to zero.<br>
!>
!>  Intrinsic Type Representation Models
!>  ------------------------------------
!>
!>  The following representation model derived types also exist in this module:<br>
!>  <ol>
!>      <li>    [modelb_type](@ref pm_kind::modelb_type) can hold information about binary (bit) containers.<br>
!>      <li>    [modeli_type](@ref pm_kind::modeli_type) can hold information about `integer` containers of various kinds.<br>
!>      <li>    [modelr_type](@ref pm_kind::modelr_type) can hold information about `complex` and `real` containers of various kinds.<br>
!>  </ol>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

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
    integer , parameter :: CR5 = range(0._CK5)
    !>  \cond excluded
#else
    integer , parameter :: CK5 = -1
    integer , parameter :: CP5 = 0
    integer , parameter :: CR5 = 0
    !>  \endcond excluded
#endif
#if CK4_ENABLED
    integer , parameter :: CK4 = real_kinds(CK4_ENABLED)
    integer , parameter :: CP4 = precision(0._CK4)
    integer , parameter :: CR4 = range(0._CK4)
    !>  \cond excluded
#else
    integer , parameter :: CK4 = -1
    integer , parameter :: CP4 = 0
    integer , parameter :: CR4 = 0
    !>  \endcond excluded
#endif
#if CK3_ENABLED
    integer , parameter :: CK3 = real_kinds(CK3_ENABLED)
    integer , parameter :: CP3 = precision(0._CK3)
    integer , parameter :: CR3 = range(0._CK3)
    !>  \cond excluded
#else
    integer , parameter :: CK3 = -1
    integer , parameter :: CP3 = 0
    integer , parameter :: CR3 = 0
    !>  \endcond excluded
#endif
#if CK2_ENABLED
    integer , parameter :: CK2 = real_kinds(CK2_ENABLED)
    integer , parameter :: CP2 = precision(0._CK2)
    integer , parameter :: CR2 = range(0._CK2)
    !>  \cond excluded
#else
    integer , parameter :: CK2 = -1
    integer , parameter :: CP2 = 0
    integer , parameter :: CR2 = 0
    !>  \endcond excluded
#endif
#if CK1_ENABLED
    integer , parameter :: CK1 = real_kinds(CK1_ENABLED)
    integer , parameter :: CP1 = precision(0._CK1)
    integer , parameter :: CR1 = range(0._CK1)
    !>  \cond excluded
#else
    integer , parameter :: CK1 = -1
    integer , parameter :: CP1 = 0
    integer , parameter :: CR1 = 0
    !>  \endcond excluded
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    integer , parameter :: RK5 = real_kinds(RK5_ENABLED)
    integer , parameter :: RP5 = precision(0._RK5)
    integer , parameter :: RR5 = range(0._RK5)
    !>  \cond excluded
#else
    integer , parameter :: RK5 = -1
    integer , parameter :: RP5 = 0
    integer , parameter :: RR5 = 0
    !>  \endcond excluded
#endif
#if RK4_ENABLED
    integer , parameter :: RK4 = real_kinds(RK4_ENABLED)
    integer , parameter :: RP4 = precision(0._RK4)
    integer , parameter :: RR4 = range(0._RK4)
    !>  \cond excluded
#else
    integer , parameter :: RK4 = -1
    integer , parameter :: RP4 = 0
    integer , parameter :: RR4 = 0
    !>  \endcond excluded
#endif
#if RK3_ENABLED
    integer , parameter :: RK3 = real_kinds(RK3_ENABLED)
    integer , parameter :: RP3 = precision(0._RK3)
    integer , parameter :: RR3 = range(0._RK3)
    !>  \cond excluded
#else
    integer , parameter :: RK3 = -1
    integer , parameter :: RP3 = 0
    integer , parameter :: RR3 = 0
    !>  \endcond excluded
#endif
#if RK2_ENABLED
    integer , parameter :: RK2 = real_kinds(RK2_ENABLED)
    integer , parameter :: RP2 = precision(0._RK2)
    integer , parameter :: RR2 = range(0._RK2)
    !>  \cond excluded
#else
    integer , parameter :: RK2 = -1
    integer , parameter :: RP2 = 0
    integer , parameter :: RR2 = 0
    !>  \endcond excluded
#endif
#if RK1_ENABLED
    integer , parameter :: RK1 = real_kinds(RK1_ENABLED)
    integer , parameter :: RP1 = precision(0._RK1)
    integer , parameter :: RR1 = range(0._RK1)
    !>  \cond excluded
#else
    integer , parameter :: RK1 = -1
    integer , parameter :: RP1 = 0
    integer , parameter :: RR1 = 0
    !>  \endcond excluded
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Default kind type parameters of the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    integer     , parameter ::  SK      =    SK_DEF !<  \public The default `character` kind in the ParaMonte library: `kind("a")`      in Fortran, `c_char`            in \cfi mode.
    integer     , parameter ::  IK      =    IK_DEF !<  \public The default `integer`   kind in the ParaMonte library: `int32`          in Fortran, `c_int32_t`         in \cfi mode.
    integer     , parameter ::  LK      =    LK_DEF !<  \public The default `logical`   kind in the ParaMonte library: `kind(.true.)`   in Fortran, `kind(.true.)`      in \cfi mode.
    integer     , parameter ::  CK      =    CK_DEF !<  \public The default `complex`   kind in the ParaMonte library: `real64`         in Fortran, `c_double_complex`  in \cfi mode.
    integer     , parameter ::  RK      =    RK_DEF !<  \public The default `real`      kind in the ParaMonte library: `real64`         in Fortran, `c_double`          in \cfi mode.

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
    integer     , parameter :: IKQ  = selected_int_kind(36)                     !<  \public The quadru-precision integer kind in Fortran mode. On most platforms, this is a 128-bit integer kind.
                                                                                !!          There is no guarantee of the existence of this kind (e.g., Intel does not support this).
    integer     , parameter :: RKS  = kind(1.)                                  !<  \public The single-precision real kind in Fortran mode. On most platforms, this is an 32-bit real kind.
    integer     , parameter :: RKD  = kind(1.d0)                                !<  \public The double precision real kind in Fortran mode. On most platforms, this is an 64-bit real kind.
    integer     , parameter :: RKQ  = selected_real_kind(2*precision(1._RKD))   !<  \public The quadru-precision real kind in Fortran mode. On most platforms, this is an 128-bit real kind.
    integer     , parameter :: CKS  = RKS                                       !<  \public The single-precision complex kind in Fortran mode. On most platforms, this is a 32-bit real kind.
    integer     , parameter :: CKD  = RKD                                       !<  \public The double precision complex kind in Fortran mode. On most platforms, this is a 64-bit real kind.
    integer     , parameter :: CKQ  = RKQ                                       !<  \public The quadru-precision complex kind in Fortran mode. On most platforms, this is a 128-bit real kind.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Define the vector of all type kinds supported by the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The `integer` vector containing all defined `character` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `character` kinds supported by the processor.<br>
    !>
    !>  \final{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: SKALL(*) = pack([SK1, SK2, SK3, SK4, SK5], 0 < [SK1, SK2, SK3, SK4, SK5])

    !>  \brief
    !>  The `integer` vector containing all defined `integer` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `integer` kinds supported by the processor.<br>
    !>
    !>  \final{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: IKALL(*) = pack([IK1, IK2, IK3, IK4, IK5], 0 < [IK1, IK2, IK3, IK4, IK5])

    !>  \brief
    !>  The `integer` vector containing all defined `logical` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `logical` kinds supported by the processor.<br>
    !>
    !>  \final{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: LKALL(*) = pack([LK1, LK2, LK3, LK4, LK5], 0 < [LK1, LK2, LK3, LK4, LK5])

    !>  \brief
    !>  The `integer` vector containing all defined `complex` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `complex` kinds supported by the processor.<br>
    !>
    !>  \final{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: CKALL(*) = pack([CK1, CK2, CK3, CK4, CK5], 0 < [CK1, CK2, CK3, CK4, CK5])

    !>  \brief
    !>  The `integer` vector containing all defined `real` kinds within the ParaMonte library.<br>
    !>
    !>  \details
    !>  Note that this vector may only be a subset of the `real` kinds supported by the processor.<br>
    !>
    !>  \final{SKALL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: RKALL(*) = pack([RK1, RK2, RK3, RK4, RK5], 0 < [RK1, RK2, RK3, RK4, RK5])

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The highest and the lowest rnages and precisions for the kind type parameters supported by the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest range among all `integer` kinds supported by the specific library build.<br>
    !>
    !>  \final{IRL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: IRL = minval([IR1, IR2, IR3, IR4, IR5], mask = 0 < [IR1, IR2, IR3, IR4, IR5]) ! Integer Range Lowest

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest range among all `integer` kinds supported by the specific library build.<br>
    !>
    !>  \final{IRH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: IRH = maxval([IR1, IR2, IR3, IR4, IR5], mask = 0 < [IR1, IR2, IR3, IR4, IR5]) ! Integer Range Highest

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest precision among all `complex` kinds supported by the specific library build.<br>
    !>
    !>  \final{CPL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: CPL = minval([CP1, CP2, CP3, CP4, CP5], mask = 0 < [CP1, CP2, CP3, CP4, CP5]) ! Complex Precision Lowest

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest precision among all `complex` kinds supported by the specific library build.<br>
    !>
    !>  \final{CPH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: CPH = maxval([CP1, CP2, CP3, CP4, CP5], mask = 0 < [CP1, CP2, CP3, CP4, CP5]) ! Complex Precision Highest

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest decimal exponent range among all `complex` kinds supported by the specific library build.<br>
    !>
    !>  \final{CRL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: CRL = minval([CR1, CR2, CR3, CR4, CR5], mask = 0 < [CR1, CR2, CR3, CR4, CR5]) ! Complex Range Lowest

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest decimal exponent range among all `complex` kinds supported by the specific library build.<br>
    !>
    !>  \final{CRH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: CRH = maxval([CR1, CR2, CR3, CR4, CR5], mask = 0 < [CR1, CR2, CR3, CR4, CR5]) ! Complex Range Lowest

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest precision among all `real` kinds supported by the specific library build.<br>
    !>
    !>  \final{RPL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: RPL = minval([RP1, RP2, RP3, RP4, RP5], mask = 0 < [RP1, RP2, RP3, RP4, RP5]) ! Real Precision Lowest

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest precision among all `real` kinds supported by the specific library build.<br>
    !>
    !>  \final{RPH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: RPH = maxval([RP1, RP2, RP3, RP4, RP5], mask = 0 < [RP1, RP2, RP3, RP4, RP5]) ! Real Precision Highest

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest decimal exponent range among all `real` kinds supported by the specific library build.<br>
    !>
    !>  \final{RRL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: RRL = minval([RR1, RR2, RR3, RR4, RR5], mask = 0 < [RR1, RR2, RR3, RR4, RR5]) ! Real Range Lowest

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest decimal exponent range among all `real` kinds supported by the specific library build.<br>
    !>
    !>  \final{RRH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
     integer    , parameter             :: RRH = maxval([RR1, RR2, RR3, RR4, RR5], mask = 0 < [RR1, RR2, RR3, RR4, RR5]) ! Real Range Lowest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The lowest precision/range kinds supported by the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest range `integer` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \IKL is the same as the value of \IKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the lowest range `integer` kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-range `integer` kind \IKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-range `integer` kind \IKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-range `integer` kind of the library \IKL, the same does not hold for \IKW when its value is different from \IKL.<br>
    !>
    !>  \final{IKL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: IKL  = selected_int_kind(IRL)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest-precision `complex` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \CKL is the same as the value of \CKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the lowest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-precision `complex` kind \CKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-precision `complex` kind \CKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-precision `complex` kind of the library \CKL, the same does not hold for \CKW when its value is different from \CKL.<br>
    !>
    !>  \final{CKL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: CKL  = selected_real_kind(CPL)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>l</b>owest-precision `real` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \RKL is the same as the value of \RKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the lowest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-precision `real` kind \RKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-precision `real` kind \RKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-precision `real` kind of the library \RKL, the same does not hold for \RKW when its value is different from \RKL.<br>
    !>
    !>  \final{RKL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: RKL  = selected_real_kind(RPL)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest-decimal-exponent-range `complex` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \CKLR is the same as the value of \CKWR under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the highest-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-range `complex` kind \CKLR <b>supported by a specific library build</b>  is not necessarily the same as the best-range `complex` kind \CKWR <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-range `complex` kind of the library \CKLR, the same does not hold for \CKWR when its value is different from \CKLR.<br>
    !>
    !>  \final{CKL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: CKLR = selected_real_kind(r = CRL)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest-decimal-exponent-range `real` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \RKLR is the same as the value of \RKWR under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the highest-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-range `real` kind \RKLR <b>supported by a specific library build</b>  is not necessarily the same as the best-range `real` kind \RKWR <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-range `real` kind of the library \RKLR, the same does not hold for \RKWR when its value is different from \RKLR.<br>
    !>
    !>  \final{RKL}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: RKLR = selected_real_kind(r = RRL)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The highest precision/range kinds supported by the ParaMonte library.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest range `integer` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \IKH is the same as the value of \IKB under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the highest range `integer` kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-range `integer` kind \IKH <b>supported by a specific library build</b>  is not necessarily the same as the best-range `integer` kind \IKB <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-range `integer` kind of the library \IKH, the same does not hold for \IKB when its value is different from \IKH.<br>
    !>
    !>  \final{IKH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: IKH  = selected_int_kind(IRH)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest-precision `complex` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \CKH is the same as the value of \CKB under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the highest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-precision `complex` kind \CKH <b>supported by a specific library build</b>  is not necessarily the same as the best-precision `complex` kind \CKB <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-precision `complex` kind of the library \CKH, the same does not hold for \CKB when its value is different from \CKH.<br>
    !>
    !>  \final{CKH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: CKH  = selected_real_kind(CPH)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest-precision `real` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \RKH is the same as the value of \RKB under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the highest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-precision `real` kind \RKH <b>supported by a specific library build</b>  is not necessarily the same as the best-precision `real` kind \RKB <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-precision `real` kind of the library \RKH, the same does not hold for \RKB when its value is different from \RKH.<br>
    !>
    !>  \final{RKH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: RKH  = selected_real_kind(RPH)

    integer     , parameter , private   :: CKHR_VEC_RAW(*)= [ selected_real_kind(merge(huge(1), CP5, CP5 < 1), r = CRH) &
                                                            , selected_real_kind(merge(huge(1), CP4, CP4 < 1), r = CRH) &
                                                            , selected_real_kind(merge(huge(1), CP3, CP3 < 1), r = CRH) &
                                                            , selected_real_kind(merge(huge(1), CP2, CP2 < 1), r = CRH) &
                                                            , selected_real_kind(merge(huge(1), CP1, CP1 < 1), r = CRH) ]
    integer     , parameter , private   :: CKHR_VEC(*) = pack(CKHR_VEC_RAW, 0 <= CKHR_VEC_RAW)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest-decimal-exponent-range `complex` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \CKHR is the same as the value of \CKBR under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the highest-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-range `complex` kind \CKHR <b>supported by a specific library build</b>  is not necessarily the same as the best-range `complex` kind \CKBR <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-range `complex` kind of the library \CKHR, the same does not hold for \CKBR when its value is different from \CKHR.<br>
    !>
    !>  \final{CKH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: CKHR = CKHR_VEC(1)

    integer     , parameter , private   :: RKHR_VEC_RAW(*)= [ selected_real_kind(merge(huge(1), RP5, RP5 < 1), r = RRH) &
                                                            , selected_real_kind(merge(huge(1), RP4, RP4 < 1), r = RRH) &
                                                            , selected_real_kind(merge(huge(1), RP3, RP3 < 1), r = RRH) &
                                                            , selected_real_kind(merge(huge(1), RP2, RP2 < 1), r = RRH) &
                                                            , selected_real_kind(merge(huge(1), RP1, RP1 < 1), r = RRH) ]
    integer     , parameter , private   :: RKHR_VEC(*) = pack(RKHR_VEC_RAW, 0 <= RKHR_VEC_RAW)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>h</b>ighest-decimal-exponent-range `real` kind type parameter available in the specific library build.<br>
    !>
    !>  \details
    !>  Although the value of \RKHR is the same as the value of \RKBR under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the highest-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-range `real` kind \RKHR <b>supported by a specific library build</b>  is not necessarily the same as the best-range `real` kind \RKBR <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-range `real` kind of the library \RKHR, the same does not hold for \RKBR when its value is different from \RKHR.<br>
    !>
    !>  \final{RKH}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: RKHR = RKHR_VEC(1)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The worst range/precision kind type parameters supported by the processor.
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
    !>  \final{IKW}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: IKW  = selected_int_kind(1)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>W</b>orst-precision `complex` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \CKL is the same as the value of \CKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the lowest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-precision `complex` kind \CKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-precision `complex` kind \CKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-precision `complex` kind of the library \CKL, the same does not hold for \CKW when its value is different from \CKL.<br>
    !>
    !>  \final{CKW}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: CKW  = selected_real_kind(1)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>W</b>orst-precision `real` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \RKL is the same as the value of \RKW under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the lowest-precision kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-precision `real` kind \RKL <b>supported by a specific library build</b>  is not necessarily the same as the worst-precision `real` kind \RKW <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-precision `real` kind of the library \RKL, the same does not hold for \RKW when its value is different from \RKL.<br>
    !>
    !>  \final{RKW}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: RKW  = selected_real_kind(1)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>W</b>orst-decimal-exponent-range `complex` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \CKLR is the same as the value of \CKWR under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the lowest-decimal-exponent-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-decimal-exponent-range `complex` kind \CKLR <b>supported by a specific library build</b>  is not necessarily the same as the worst-decimal-exponent-range `complex` kind \CKWR <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-decimal-exponent-range `complex` kind of the library \CKLR, the same does not hold for \CKWR when its value is different from \CKLR.<br>
    !>
    !>  \final{CKWR}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: CKWR = selected_real_kind(r = 1)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>W</b>orst-decimal-exponent-range `real` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \RKLR is the same as the value of \RKWR under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the lowest-decimal-exponent-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the lowest-decimal-exponent-range `real` kind \RKLR <b>supported by a specific library build</b>  is not necessarily the same as the worst-decimal-exponent-range `real` kind \RKWR <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the lowest-decimal-exponent-range `real` kind of the library \RKLR, the same does not hold for \RKWR when its value is different from \RKLR.<br>
    !>
    !>  \final{RKWR}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: RKWR = selected_real_kind(r = 1)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The best range/precision kind type parameters supported by the processor.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The vector of `integer` constants of intrinsic default kind, representing the vector of `50` possible `integer`
    !>  ranges of the kinds provided in `integer_kinds` of the intrinsic Fortran module `iso_fortran_env`.<br>
    !>
    !>  \details
    !>  This nightmare is necessary to identify the highest-range `integer` kind type parameter supported by the processor.
    !>
    !>  \warning
    !>  The current implementation assumes the `integer_kinds` vector of the intrinsic Fortran module `iso_fortran_env` has a maximum length of `50`.<br>
    !>  While this assumption will likely hold for many more decades to come, it is bound to fail in the distant future.<br>
    !>
    !>  \final{integer_kinds_range}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter , private   :: integer_kinds_range(*) = [ range(int(0, kind = integer_kinds(min(size(integer_kinds),  1)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds),  2)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds),  3)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds),  4)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds),  5)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds),  6)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds),  7)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds),  8)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds),  9)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 10)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 11)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 12)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 13)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 14)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 15)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 16)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 17)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 18)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 19)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 20)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 21)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 22)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 23)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 24)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 25)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 26)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 27)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 28)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 29)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 30)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 31)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 32)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 33)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 34)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 35)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 36)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 37)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 38)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 39)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 40)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 40)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 41)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 42)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 43)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 44)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 45)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 46)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 47)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 48)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 49)))) &
                                                                    , range(int(0, kind = integer_kinds(min(size(integer_kinds), 50)))) &
                                                                    ]

    !>  \brief
    !>  The vector of `integer` constants of intrinsic default kind, representing the vector of `50` possible `real`
    !>  precisions of the kinds provided in `real_kinds` of the intrinsic Fortran module `iso_fortran_env`.<br>
    !>
    !>  \details
    !>  This nightmare is necessary to identify the highest-range/precision `real` kind type parameter supported by the processor.
    !>
    !>  \warning
    !>  The current implementation assumes the `real_kinds` vector of the intrinsic Fortran module `iso_fortran_env` has a maximum length of `50`.<br>
    !>  While this assumption will likely hold for many more decades to come, it is bound to fail in the distant future.<br>
    !>
    !>  \final{real_kinds_range}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter , private   :: real_kinds_range(*) =    [ range(real(0, kind = real_kinds(min(size(real_kinds),  1)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds),  2)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds),  3)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds),  4)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds),  5)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds),  6)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds),  7)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds),  8)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds),  9)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 10)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 11)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 12)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 13)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 14)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 15)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 16)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 17)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 18)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 19)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 20)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 21)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 22)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 23)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 24)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 25)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 26)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 27)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 28)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 29)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 30)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 31)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 32)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 33)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 34)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 35)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 36)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 37)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 38)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 39)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 40)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 40)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 41)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 42)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 43)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 44)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 45)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 46)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 47)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 48)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 49)))) &
                                                                    , range(real(0, kind = real_kinds(min(size(real_kinds), 50)))) &
                                                                    ]

    !>  \brief
    !>  The vector of `integer` constants of intrinsic default kind, representing the vector of `50` possible `real`
    !>  highest-range precisions of the kinds provided in `real_kinds` of the intrinsic Fortran module `iso_fortran_env`.<br>
    !>
    !>  \details
    !>  This nightmare is necessary to identify the highest-range highest-precision `real` kind type parameter supported by the processor.
    !>
    !>  \warning
    !>  The current implementation assumes the `real_kinds` vector of the intrinsic Fortran module `iso_fortran_env` has a maximum length of `50`.<br>
    !>  While this assumption will likely hold for many more decades to come, it is bound to fail in the distant future.<br>
    !>
    !>  \final{real_kinds_precision}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter , private   :: real_kinds_precision(*)= [ precision(real(0, kind = real_kinds(min(size(real_kinds),  1)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds),  2)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds),  3)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds),  4)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds),  5)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds),  6)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds),  7)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds),  8)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds),  9)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 10)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 11)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 12)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 13)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 14)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 15)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 16)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 17)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 18)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 19)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 20)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 21)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 22)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 23)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 24)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 25)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 26)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 27)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 28)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 29)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 30)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 31)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 32)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 33)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 34)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 35)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 36)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 37)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 38)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 39)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 40)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 40)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 41)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 42)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 43)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 44)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 45)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 46)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 47)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 48)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 49)))) &
                                                                    , precision(real(0, kind = real_kinds(min(size(real_kinds), 50)))) &
                                                                    ]

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the
    !>  highest-decimal-exponent-range of `real` types made available by the processor.
    !>
    !>  \final{real_kinds_range_max}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter , private   :: real_kinds_range_max = maxval(real_kinds_range, dim = 1)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing
    !>  the highest-precision of `real` types made available by the processor.
    !>
    !>  \final{real_kinds_precision_max}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter , private   :: real_kinds_precision_max = maxval(real_kinds_precision, dim = 1)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing
    !>  the lowest-precision of `real` types made available by the processor.
    !>
    !>  \final{real_kinds_precision_min}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter , private   :: real_kinds_precision_min = minval(real_kinds_precision, dim = 1)

    !>  \brief
    !>  The scalar `real` constant of intrinsic default kind, representing
    !>  the step size between the highest and lowest-precision of `real` types made available by the processor.
    !>
    !>  \details
    !>  This constant is internally used within the module to identify the highest-range highest-precision `real` kind type parameter.
    !>
    !>  \final{real_kinds_precision_hop}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    real        , parameter , private   :: real_kinds_precision_hop = real(real_kinds_precision_max - real_kinds_precision_min) / size(real_kinds_precision)

    !>  \brief
    !>  The vector of `integer` constants of intrinsic default kind, representing the `50` possible `real`
    !>  highest-range precisions of the kinds provided in `real_kinds` of the intrinsic Fortran module `iso_fortran_env`.<br>
    !>
    !>  \details
    !>  This nightmare is necessary to identify the highest-range highest-precision `real` kind type parameter supported by the processor.
    !>
    !>  \warning
    !>  The current implementation assumes the `real_kinds` vector of the intrinsic Fortran module `iso_fortran_env` has a maximum length of `50`.<br>
    !>  While this assumption will likely hold for many more decades to come, it is bound to fail in the distant future.<br>
    !>
    !>  \final{real_kinds_prmax_kind}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter , private   :: real_kinds_prmax_kind(*) =   [ selected_real_kind(p = real_kinds_precision_max, r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  0)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  1)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  2)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  3)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  4)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  5)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  6)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  7)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  8)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop *  9)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 10)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 11)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 12)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 13)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 14)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 15)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 16)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 17)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 18)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 19)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 20)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 21)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 22)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 23)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 24)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 25)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 26)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 27)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 28)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 29)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 30)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 31)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 32)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 33)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 34)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 35)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 36)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 37)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 38)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 39)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 40)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 41)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 42)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 43)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 44)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 45)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 46)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 47)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 48)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 49)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = max(real_kinds_precision_min, nint(real_kinds_precision_max - real_kinds_precision_hop * 50)), r = real_kinds_range_max) &
                                                                        , selected_real_kind(p = real_kinds_precision_min, r = real_kinds_range_max) &
                                                                        ]

    !>  \brief
    !>  The vector of `integer` constants of intrinsic default kind, containing all maximum-range `real` kind type parameters
    !>  corresponding to various precisions provided in `real_kinds` of the intrinsic Fortran module `iso_fortran_env`.<br>
    !>
    !>  \details
    !>  This nightmare is necessary to identify the highest-range highest-precision `real` kind type parameter supported by the processor.
    !>
    !>  \final{real_kinds_prmax_kind_avail}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter , private   :: real_kinds_prmax_kind_avail(*) = pack(real_kinds_prmax_kind, 0 <= real_kinds_prmax_kind)

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
    !>  \final{IKB}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    ! maxval([(selected_int_kind(i), i = 1, 1000)]) ! merge(IKH6, merge(IKH5, merge(IKH4, merge(IKH3, merge(IKH2, IKS, IKH2 > 0), IKH3 > 0), IKH4 > 0), IKH5 > 0), IKH6 > 0) ! only up to `3` times more precise than the double range kind.
    integer     , parameter             :: IKB = selected_int_kind(maxval(integer_kinds_range, dim = 1))

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
    !>  \final{CKB}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: CKB = selected_real_kind(maxval(real_kinds_precision, dim = 1))

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
    !>  \final{RKB}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    ! maxval([(selected_real_kind(i), i = 1, 1000)]) ! merge(RKH6, merge(RKH5, merge(RKH4, merge(RKH3, merge(RKH2, RKS, RKH2 > 0), RKH3 > 0), RKH4 > 0), RKH5 > 0), RKH6 > 0) ! only up to `3` times more precise than the double precision kind.
    integer     , parameter             :: RKB = selected_real_kind(maxval(real_kinds_precision, dim = 1))

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>B</b>est-decimal-exponent-range `complex` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \CKHR is the same as the value of \CKBR under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `complex` kind type parameters that exclude the highest-decimal-exponent-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-decimal-exponent-range `complex` kind \CKHR <b>supported by a specific library build</b>  is not necessarily the same as the best-decimal-exponent-range `complex` kind \CKBR <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-decimal-exponent-range `complex` kind of the library \CKHR, the same does not hold for \CKBR when its value is different from \CKHR.<br>
    !>
    !>  The current Fortran standard does not allow automatic selection of the highest-decimal-exponent-range `complex` kind made available by the processor.<br>
    !>  However, such a kind is essential for defining `complex` constants of highest-decimal-exponent-range that can be later coerced to `complex` kinds of lower decimal-exponent-range.<br>
    !>
    !>  \final{CKBR}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: CKBR = real_kinds_prmax_kind_avail(1)

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind, representing the <b>B</b>est-decimal-exponent-range `real` kind **supported by the processor**.<br>
    !>
    !>  \details
    !>  Although the value of \RKHR is the same as the value of \RKBR under normal (default) library builds, the two are not necessarily the same.<br>
    !>  This situation occurs when the library is built for `real` kind type parameters that exclude the highest-decimal-exponent-range kind <b>supported by the processor</b>.<br>
    !>  In other words, the highest-decimal-exponent-range `real` kind \RKHR <b>supported by a specific library build</b>  is not necessarily the same as the best-decimal-exponent-range `real` kind \RKBR <b>supported by the processor</b>.<br>
    !>  While all relevant routines of the library are guaranteed to support the highest-decimal-exponent-range `real` kind of the library \RKHR, the same does not hold for \RKBR when its value is different from \RKHR.<br>
    !>
    !>  The current Fortran standard does not allow automatic selection of the highest-decimal-exponent-range `real` kind made available by the processor.<br>
    !>  However, such a kind is essential for defining `real` constants of highest-decimal-exponent-range that can be later coerced to `real` kinds of lower decimal-exponent-range.<br>
    !>
    !>  \final{RKBR}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    integer     , parameter             :: RKBR = real_kinds_prmax_kind_avail(1)

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
    !>  \final{model_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
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
    !>  \final{modeln_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
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
    !>  \final{modeli_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
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
    !>  \final{modelb_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
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
    !>  \final{modelr_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(modeln_type) :: modelr_type
        real(RKH)   :: epsilon          !<  \public The scalar `real` of the highest kind \RKH supported by the ParaMonte library,
                                        !!          whose value is almost negligible compared with the value `1.0` in the model that includes the `real` data type of interest,
                                        !!          that is \f$b^{1-p}\f$ as detailed in the real data model of [model_type](@ref pm_kind::model_type).<br>
                                        !!          \note
                                        !!          EPSILON makes it easy to select a delta for algorithms (such as root locators) that search until the calculation is within delta of an estimate.<br>
                                        !!          If delta is too small (smaller than the decimal resolution of the data type), the algorithm might never halt.<br>
                                        !!          By scaling the value returned by EPSILON to the estimate, you obtain a delta that ensures search termination.<br>
        real(RKHR)  :: huge             !<  \public The scalar `real` of the highest kind \RKHR supported by the ParaMonte library,
                                        !!          representing the largest value in the model that includes the `real` data type of interest.<br>
        integer(IK) :: maxexponent      !<  \public The scalar `integer` of default kind \IK, representing the maximum exponent in the model that includes the `real` data type of interest.<br>
                                        !!          This corresponds to \f$e_{max}\f$ in the model set for the `real` data type detailed in [model_type](@ref pm_kind::model_type).<br>
        integer(IK) :: minexponent      !<  \public The scalar `integer` of default kind \IK, representing the minimum exponent in the model that includes the `real` data type of interest.<br>
                                        !!          This corresponds to \f$e_{min}\f$ in the model set for the `real` data type detailed in [model_type](@ref pm_kind::model_type).<br>
        integer(IK) :: precision        !<  \public The scalar `integer` of default kind \IK, representing the equivalent decimal precision (number of digits after the decimal point)
                                        !!          in the model representing real numbers with the same type parameter value as the `real` data type of interest.<br>
                                        !!          The value is `int((p - 1) * log10(b)) + k` where `k` is `1` if `b` is an integral power of `10` and `0` otherwise.<br>
                                        !!          The meaning of the `real` data type model parameters are detailed in [model_type](@ref pm_kind::model_type).<br>
        real(RKHR)  :: tiny             !<  \public The scalar `real` of the highest kind \RKHR supported by the ParaMonte library,
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
    !>  \interface{modeli_typer}
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
    !>  \example{modeli_typer}
    !>  \include{lineno} example/pm_kind/modeli_type/main.F90
    !>  \compilef{modeli_typer}
    !>  \output{modeli_typer}
    !>  \include{lineno} example/pm_kind/modeli_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_kind](@ref test_pm_kind)
    !>
    !>  \final{modeli_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface modeli_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function modeli_typer_IK5(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modeli_typer_IK5
#endif
        integer(IK5), intent(in)    :: mold
        type(modeli_type)           :: model
    end function
#endif

#if IK4_ENABLED
    pure elemental module function modeli_typer_IK4(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modeli_typer_IK4
#endif
        integer(IK4), intent(in)    :: mold
        type(modeli_type)           :: model
    end function
#endif

#if IK3_ENABLED
    pure elemental module function modeli_typer_IK3(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modeli_typer_IK3
#endif
        integer(IK3), intent(in)    :: mold
        type(modeli_type)           :: model
    end function
#endif

#if IK2_ENABLED
    pure elemental module function modeli_typer_IK2(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modeli_typer_IK2
#endif
        integer(IK2), intent(in)    :: mold
        type(modeli_type)           :: model
    end function
#endif

#if IK1_ENABLED
    pure elemental module function modeli_typer_IK1(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modeli_typer_IK1
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
    !>  \interface{modelb_typer}
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
    !>  \example{modelb_typer}
    !>  \include{lineno} example/pm_kind/modelb_type/main.F90
    !>  \compilef{modelb_typer}
    !>  \output{modelb_typer}
    !>  \include{lineno} example/pm_kind/modelb_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_kind](@ref test_pm_kind)
    !>
    !>  \final{modelb_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface modelb_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function modelb_typer_IK5(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelb_typer_IK5
#endif
        integer(IK5), intent(in)    :: mold
        type(modelb_type)           :: model
    end function
#endif

#if IK4_ENABLED
    pure elemental module function modelb_typer_IK4(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelb_typer_IK4
#endif
        integer(IK4), intent(in)    :: mold
        type(modelb_type)           :: model
    end function
#endif

#if IK3_ENABLED
    pure elemental module function modelb_typer_IK3(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelb_typer_IK3
#endif
        integer(IK3), intent(in)    :: mold
        type(modelb_type)           :: model
    end function
#endif

#if IK2_ENABLED
    pure elemental module function modelb_typer_IK2(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelb_typer_IK2
#endif
        integer(IK2), intent(in)    :: mold
        type(modelb_type)           :: model
    end function
#endif

#if IK1_ENABLED
    pure elemental module function modelb_typer_IK1(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelb_typer_IK1
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
    !>  \interface{modelr_typer}
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
    !>  \example{modelr_typer}
    !>  \include{lineno} example/pm_kind/modelr_type/main.F90
    !>  \compilef{modelr_typer}
    !>  \output{modelr_typer}
    !>  \include{lineno} example/pm_kind/modelr_type/main.out.F90
    !>
    !>  \test
    !>  [test_pm_kind](@ref test_pm_kind)
    !>
    !>  \final{modelr_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface modelr_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function modelr_typer_RK5(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelr_typer_RK5
#endif
        real(RK5)   , intent(in)    :: mold
        type(modelr_type)           :: model
    end function
#endif

#if RK4_ENABLED
    pure elemental module function modelr_typer_RK4(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelr_typer_RK4
#endif
        real(RK4)   , intent(in)    :: mold
        type(modelr_type)           :: model
    end function
#endif

#if RK3_ENABLED
    pure elemental module function modelr_typer_RK3(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelr_typer_RK3
#endif
        real(RK3)   , intent(in)    :: mold
        type(modelr_type)           :: model
    end function
#endif

#if RK2_ENABLED
    pure elemental module function modelr_typer_RK2(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelr_typer_RK2
#endif
        real(RK2)   , intent(in)    :: mold
        type(modelr_type)           :: model
    end function
#endif

#if RK1_ENABLED
    pure elemental module function modelr_typer_RK1(mold) result(model)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: modelr_typer_RK1
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