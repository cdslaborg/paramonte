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
!>  This module contains procedures and generic interfaces for convenient allocation and filling of arrays of
!>  arbitrary intrinsic types (i.e., `character`, `integer`, `logical`, `complex`, `real`), kinds, and non-zero ranks (up to `3`).<br>
!>
!>  \details
!>  The functionalities of the procedures of this module are equivalent and expand
!>  the functionalities of the convenience MATLAB functions `zeros()` and `ones()`.<br>
!>
!>  \note
!>  Use the intrinsic Fortran function `repeat()` to generate and initialize scalars of type `character`.<br>
!>
!>  \see
!>  [pm_arrayInit](@ref pm_arrayInit)<br>
!>  [pm_matrixInit](@ref pm_matrixInit)<br>
!>
!>  \test
!>  [test_pm_arrayFill](@ref test_pm_arrayFill)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayFill

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayFill"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an array of the specified rank and shape of arbitrary intrinsic
    !>  type and kind with all its elements initialized to the user specified value.<br>
    !>
    !>  \param[in]      val     :   The input scalar of the same type and kind as the output `array`,
    !>                              containing the initialization value for the elements of the array.<br>
    !>  \param[in]      s1      :   The input scalar of type `integer` of default kind \IK, representing the size of the output `array` along its first dimension.<br>
    !>  \param[in]      s2      :   The input scalar of type `integer` of default kind \IK, representing the size of the output `array` array along its second dimension.<br>
    !>                              (**optional**. It must be present if `s3` is present.)
    !>  \param[in]      s3      :   The input scalar of type `integer` of default kind \IK, representing the size of the output `array` array along its third dimension.<br>
    !>                              (**optional**. It can be present only if `s2` is present.)
    !>
    !>  \return
    !>  `array`                 :   The output scalar of,<br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL of the specified length type parameter `size` or,<br>
    !>                              </ul>
    !>                              the output array of rank `1`, `2`, or `3` of shape `(s1)` or `(s1, s2)` or `(s1, s2, s3)`, of either,<br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or<br>
    !>                                  <li>    type `integer` of kind \IKALL, or<br>
    !>                                  <li>    type `logical` of kind \LKALL, or<br>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL.<br>
    !>                              </ul>
    !>                              On output, the `array` is allocated and initialized to the requested input value.<br>
    !>
    !>  \interface{getFilled}
    !>  \code{.F90}
    !>
    !>      use pm_arrayFill, only: getFilled
    !>
    !>      array = getFilled(val, s1) ! array(1:s1)
    !>      array = getFilled(val, s1, s2) ! array(1:s1, 1:s2)
    !>      array = getFilled(val, s1, s2, s3) ! array(1:s1, 1:s2, 1:s3)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= s1` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= s2` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= s3` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  This generic interface is merely a convenient functional around the Fortran `allocate` statement.<br>
    !>  There is no reason to use this generic interface in performance-critical applications.<br>
    !>
    !>  \see
    !>  [getMatInit](@ref pm_matrixInit::getMatInit)<br>
    !>  [getCoreHalo](@ref pm_arrayInit::getCoreHalo)<br>
    !>  [setCoreHalo](@ref pm_arrayInit::setCoreHalo)<br>
    !>  [setResized](@ref pm_arrayResize::setResized)<br>
    !>  [setRebound](@ref pm_arrayRebind::setRebound)<br>
    !>  [setRebilled](@ref pm_arrayRebill::setRebilled)<br>
    !>  [setRefilled](@ref pm_arrayRefill::setRefilled)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>  [getPadded](@ref pm_arrayPad::getPadded)<br>
    !>  [setPadded](@ref pm_arrayPad::setPadded)<br>
    !>
    !>  \example{getFilled}
    !>  \include{lineno} example/pm_arrayFill/getFilled/main.F90
    !>  \compilef{getFilled}
    !>  \output{getFilled}
    !>  \include{lineno} example/pm_arrayFill/getFilled/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayFill](@ref test_pm_arrayFill)
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be extended to `array` arguments of higher ranks.<br>
    !>
    !>  \final{getFilled}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getFilled

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getFilled_D1_SK5(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(in)                        :: s1
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1)
    end function
#endif

#if SK4_ENABLED
    PURE module function getFilled_D1_SK4(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(in)                        :: s1
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1)
    end function
#endif

#if SK3_ENABLED
    PURE module function getFilled_D1_SK3(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(in)                        :: s1
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1)
    end function
#endif

#if SK2_ENABLED
    PURE module function getFilled_D1_SK2(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(in)                        :: s1
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1)
    end function
#endif

#if SK1_ENABLED
    PURE module function getFilled_D1_SK1(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(in)                        :: s1
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getFilled_D1_IK5(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(in)                        :: s1
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1)
    end function
#endif

#if IK4_ENABLED
    PURE module function getFilled_D1_IK4(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(in)                        :: s1
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1)
    end function
#endif

#if IK3_ENABLED
    PURE module function getFilled_D1_IK3(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(in)                        :: s1
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1)
    end function
#endif

#if IK2_ENABLED
    PURE module function getFilled_D1_IK2(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(in)                        :: s1
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1)
    end function
#endif

#if IK1_ENABLED
    PURE module function getFilled_D1_IK1(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(in)                        :: s1
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getFilled_D1_LK5(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(in)                        :: s1
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1)
    end function
#endif

#if LK4_ENABLED
    PURE module function getFilled_D1_LK4(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(in)                        :: s1
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1)
    end function
#endif

#if LK3_ENABLED
    PURE module function getFilled_D1_LK3(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(in)                        :: s1
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1)
    end function
#endif

#if LK2_ENABLED
    PURE module function getFilled_D1_LK2(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(in)                        :: s1
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1)
    end function
#endif

#if LK1_ENABLED
    PURE module function getFilled_D1_LK1(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(in)                        :: s1
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getFilled_D1_CK5(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(in)                        :: s1
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1)
    end function
#endif

#if CK4_ENABLED
    PURE module function getFilled_D1_CK4(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(in)                        :: s1
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1)
    end function
#endif

#if CK3_ENABLED
    PURE module function getFilled_D1_CK3(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(in)                        :: s1
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1)
    end function
#endif

#if CK2_ENABLED
    PURE module function getFilled_D1_CK2(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(in)                        :: s1
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1)
    end function
#endif

#if CK1_ENABLED
    PURE module function getFilled_D1_CK1(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(in)                        :: s1
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getFilled_D1_RK5(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(in)                        :: s1
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1)
    end function
#endif

#if RK4_ENABLED
    PURE module function getFilled_D1_RK4(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(in)                        :: s1
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1)
    end function
#endif

#if RK3_ENABLED
    PURE module function getFilled_D1_RK3(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(in)                        :: s1
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1)
    end function
#endif

#if RK2_ENABLED
    PURE module function getFilled_D1_RK2(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(in)                        :: s1
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1)
    end function
#endif

#if RK1_ENABLED
    PURE module function getFilled_D1_RK1(val, s1) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(in)                        :: s1
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getFilled_D2_SK5(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(in)                        :: s1, s2
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2)
    end function
#endif

#if SK4_ENABLED
    PURE module function getFilled_D2_SK4(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(in)                        :: s1, s2
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2)
    end function
#endif

#if SK3_ENABLED
    PURE module function getFilled_D2_SK3(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(in)                        :: s1, s2
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2)
    end function
#endif

#if SK2_ENABLED
    PURE module function getFilled_D2_SK2(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(in)                        :: s1, s2
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2)
    end function
#endif

#if SK1_ENABLED
    PURE module function getFilled_D2_SK1(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(in)                        :: s1, s2
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getFilled_D2_IK5(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(in)                        :: s1, s2
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2)
    end function
#endif

#if IK4_ENABLED
    PURE module function getFilled_D2_IK4(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(in)                        :: s1, s2
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2)
    end function
#endif

#if IK3_ENABLED
    PURE module function getFilled_D2_IK3(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(in)                        :: s1, s2
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2)
    end function
#endif

#if IK2_ENABLED
    PURE module function getFilled_D2_IK2(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(in)                        :: s1, s2
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2)
    end function
#endif

#if IK1_ENABLED
    PURE module function getFilled_D2_IK1(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(in)                        :: s1, s2
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getFilled_D2_LK5(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(in)                        :: s1, s2
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2)
    end function
#endif

#if LK4_ENABLED
    PURE module function getFilled_D2_LK4(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(in)                        :: s1, s2
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2)
    end function
#endif

#if LK3_ENABLED
    PURE module function getFilled_D2_LK3(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(in)                        :: s1, s2
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2)
    end function
#endif

#if LK2_ENABLED
    PURE module function getFilled_D2_LK2(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(in)                        :: s1, s2
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2)
    end function
#endif

#if LK1_ENABLED
    PURE module function getFilled_D2_LK1(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(in)                        :: s1, s2
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getFilled_D2_CK5(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(in)                        :: s1, s2
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2)
    end function
#endif

#if CK4_ENABLED
    PURE module function getFilled_D2_CK4(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(in)                        :: s1, s2
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2)
    end function
#endif

#if CK3_ENABLED
    PURE module function getFilled_D2_CK3(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(in)                        :: s1, s2
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2)
    end function
#endif

#if CK2_ENABLED
    PURE module function getFilled_D2_CK2(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(in)                        :: s1, s2
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2)
    end function
#endif

#if CK1_ENABLED
    PURE module function getFilled_D2_CK1(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(in)                        :: s1, s2
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getFilled_D2_RK5(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(in)                        :: s1, s2
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2)
    end function
#endif

#if RK4_ENABLED
    PURE module function getFilled_D2_RK4(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(in)                        :: s1, s2
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2)
    end function
#endif

#if RK3_ENABLED
    PURE module function getFilled_D2_RK3(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(in)                        :: s1, s2
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2)
    end function
#endif

#if RK2_ENABLED
    PURE module function getFilled_D2_RK2(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(in)                        :: s1, s2
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2)
    end function
#endif

#if RK1_ENABLED
    PURE module function getFilled_D2_RK1(val, s1, s2) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D2_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(in)                        :: s1, s2
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getFilled_D3_SK5(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2, s3)
    end function
#endif

#if SK4_ENABLED
    PURE module function getFilled_D3_SK4(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2, s3)
    end function
#endif

#if SK3_ENABLED
    PURE module function getFilled_D3_SK3(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2, s3)
    end function
#endif

#if SK2_ENABLED
    PURE module function getFilled_D3_SK2(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2, s3)
    end function
#endif

#if SK1_ENABLED
    PURE module function getFilled_D3_SK1(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        character(*,SKG)            , intent(in)                        :: val
        character(len(val,IK),SKG)                                      :: array(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getFilled_D3_IK5(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if IK4_ENABLED
    PURE module function getFilled_D3_IK4(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if IK3_ENABLED
    PURE module function getFilled_D3_IK3(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if IK2_ENABLED
    PURE module function getFilled_D3_IK2(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if IK1_ENABLED
    PURE module function getFilled_D3_IK1(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        integer(IKG)                , intent(in)                        :: val
        integer(IKG)                                                    :: array(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getFilled_D3_LK5(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if LK4_ENABLED
    PURE module function getFilled_D3_LK4(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if LK3_ENABLED
    PURE module function getFilled_D3_LK3(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if LK2_ENABLED
    PURE module function getFilled_D3_LK2(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if LK1_ENABLED
    PURE module function getFilled_D3_LK1(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        logical(LKG)                , intent(in)                        :: val
        logical(LKG)                                                    :: array(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getFilled_D3_CK5(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if CK4_ENABLED
    PURE module function getFilled_D3_CK4(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if CK3_ENABLED
    PURE module function getFilled_D3_CK3(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if CK2_ENABLED
    PURE module function getFilled_D3_CK2(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2, s3)
    end function
#endif

#if CK1_ENABLED
    PURE module function getFilled_D3_CK1(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        complex(CKG)                , intent(in)                        :: val
        complex(CKG)                                                    :: array(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getFilled_D3_RK5(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2, s3)
    end function
#endif

#if RK4_ENABLED
    PURE module function getFilled_D3_RK4(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2, s3)
    end function
#endif

#if RK3_ENABLED
    PURE module function getFilled_D3_RK3(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2, s3)
    end function
#endif

#if RK2_ENABLED
    PURE module function getFilled_D3_RK2(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2, s3)
    end function
#endif

#if RK1_ENABLED
    PURE module function getFilled_D3_RK1(val, s1, s2, s3) result(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFilled_D3_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(in)                        :: s1, s2, s3
        real(RKG)                   , intent(in)                        :: val
        real(RKG)                                                       :: array(s1, s2, s3)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayFill ! LCOV_EXCL_LINE