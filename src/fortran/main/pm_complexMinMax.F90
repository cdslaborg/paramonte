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
!>  This module contains procedures and generic interfaces for computing element-wise minimum/maximum value/location
!>  of the real and imaginary components of scalars and arrays of arbitrary ranks of type `complex` of arbitrary kinds.
!>
!>  \remark
!>  The primary purpose of the procedures in this module is to provide a convenient
!>  element-wise evaluation of the minimum of the real and imaginary components of complex numbers.<br>
!>  Such cases frequently occur in various library testing scenarios.<br>
!>
!>  \test
!>  [test_pm_complexMinMax](@ref test_pm_complexMinMax)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \cond excluded
    ! Note: pointer rank remapping does not work for non-contiguous `array`.
#if 1
#define __CONTIGUOUS
#else
#define __CONTIGUOUS, contiguous
#endif
!>  \endcond excluded

module pm_complexMinMax

    use pm_kind, only: SK, IK
    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_complexMinMax"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the component-wise minimum value of (both real and imaginary parts of) the input `complex`.
    !>
    !>  \param[in]  a1  :   The input scalar or array of arbitrary rank of type `complex` of kind \CKALL.
    !>  \param[in]  a2  :   The input scalar or array of arbitrary rank of the same type and kind as `a1`.
    !>
    !>  \return
    !>  `val`           :   The output scalar or array of the same rank and shape as `a1` whose real and imaginary parts
    !>                      are set to the minimum values of the corresponding real and imaginary parts of the input `a1` and `a2`.
    !>
    !>  \interface{min}
    !>  \code{.F90}
    !>
    !>      use pm_complexMinMax, only: min
    !>
    !>      val = min(a1, a2) ! = cmplx(min(a1%re, a2%re), min(a1%im, a2%im), kind(a1))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [min](@ref pm_complexMinMax::min)<br>
    !>  [max](@ref pm_complexMinMax::max)<br>
    !>  [minval](@ref pm_complexMinMax::minval)<br>
    !>  [maxval](@ref pm_complexMinMax::maxval)<br>
    !>  [minloc](@ref pm_complexMinMax::minloc)<br>
    !>  [maxloc](@ref pm_complexMinMax::maxloc)<br>
    !>  [pm_complexCompareAll](@ref pm_complexCompareAll)<br>
    !>  [pm_complexCompareAny](@ref pm_complexCompareAny)<br>
    !>  [pm_complexCompareLex](@ref pm_complexCompareLex)<br>
    !>  [pm_arrayMinMax](@ref pm_arrayMinMax)<br>
    !>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
    !>
    !>  \example{min}
    !>  \include{lineno} example/pm_complexMinMax/min/main.F90
    !>  \compilef{min}
    !>  \output{min}
    !>  \include{lineno} example/pm_complexMinMax/min/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexMinMax](@ref test_pm_complexMinMax)
    !>
    !>  \final{min}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface min

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function min_D0_CK5(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: min_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

#if CK4_ENABLED
    pure elemental module function min_D0_CK4(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: min_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

#if CK3_ENABLED
    pure elemental module function min_D0_CK3(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: min_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

#if CK2_ENABLED
    pure elemental module function min_D0_CK2(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: min_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

#if CK1_ENABLED
    pure elemental module function min_D0_CK1(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: min_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the component-wise maximum value of (both real and imaginary parts of) the input `complex`.
    !>
    !>  \param[in]  a1  :   The input scalar or array of arbitrary rank of type `complex` of kind \CKALL.
    !>  \param[in]  a2  :   The input scalar or array of arbitrary rank of the same type and kind as `a1`.
    !>
    !>  \return
    !>  `val`           :   The output scalar or array of the same rank and shape as `a1` whose real and imaginary parts
    !>                      are set to the maximum values of the corresponding real and imaginary parts of the input `a1` and `a2`.
    !>
    !>  \interface{max}
    !>  \code{.F90}
    !>
    !>      use pm_complexMinMax, only: max
    !>
    !>      val = max(a1, a2) ! = cmplx(max(a1%re, a2%re), max(a1%im, a2%im), kind(a1))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [min](@ref pm_complexMinMax::min)<br>
    !>  [max](@ref pm_complexMinMax::max)<br>
    !>  [minval](@ref pm_complexMinMax::minval)<br>
    !>  [maxval](@ref pm_complexMinMax::maxval)<br>
    !>  [minloc](@ref pm_complexMinMax::minloc)<br>
    !>  [maxloc](@ref pm_complexMinMax::maxloc)<br>
    !>  [pm_complexCompareAll](@ref pm_complexCompareAll)<br>
    !>  [pm_complexCompareAny](@ref pm_complexCompareAny)<br>
    !>  [pm_complexCompareLex](@ref pm_complexCompareLex)<br>
    !>  [pm_arrayMinMax](@ref pm_arrayMinMax)<br>
    !>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
    !>
    !>  \example{max}
    !>  \include{lineno} example/pm_complexMinMax/max/main.F90
    !>  \compilef{max}
    !>  \output{max}
    !>  \include{lineno} example/pm_complexMinMax/max/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexMinMax](@ref test_pm_complexMinMax)
    !>
    !>  \final{max}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface max

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function max_D0_CK5(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: max_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

#if CK4_ENABLED
    pure elemental module function max_D0_CK4(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: max_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

#if CK3_ENABLED
    pure elemental module function max_D0_CK3(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: max_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

#if CK2_ENABLED
    pure elemental module function max_D0_CK2(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: max_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

#if CK1_ENABLED
    pure elemental module function max_D0_CK1(a1, a2) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: max_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)                        :: a1, a2
        complex(CKG)                                    :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the component-wise minimum value of (both real and imaginary parts of) the input `complex`.
    !>
    !>  \details
    !>  This generic interfaces extends the behavior of the intrinsic Fortran `minval(a, dim = dim)` for arguments of `real` type to arguments of type `complex`.<br>
    !>
    !>  \param[in]  array   :   The input array of rank `(1..2)` of type `complex` of kind \CKALL.
    !>  \param[in]  dim     :   The input scalar of type `integer` of default kind \IK of value `1 <= dim .and. dim <= rank(array)`,
    !>                          representing the dimension along which the minimum must be computed.<br>
    !>                          (**optional**. If missing, the output value is a scalar containing the minimum of the entire input `array`.)
    !>
    !>  \return
    !>  `val`               :   The output object of the same type and kind as the input `array` of rank `rank(array) - 1` whose real and
    !>                          imaginary parts are the minimum values of the corresponding real and imaginary parts of the input `array`.<br>
    !>                          The returned values are `+huge(real(0, kind(array))` if the input `array` has zero size along the dimension of interest.<br>
    !>
    !>  \interface{minval}
    !>  \code{.F90}
    !>
    !>      use pm_complexMinMax, only: minval
    !>
    !>      val = minval(array) ! = cmplx(minval(array%re), minval(array%im), kind(array))
    !>      val = minval(array, dim) ! = cmplx(minval(array%re, dim), minval(array%im, dim), kind(array))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [min](@ref pm_complexMinMax::min)<br>
    !>  [max](@ref pm_complexMinMax::max)<br>
    !>  [minval](@ref pm_complexMinMax::minval)<br>
    !>  [maxval](@ref pm_complexMinMax::maxval)<br>
    !>  [minloc](@ref pm_complexMinMax::minloc)<br>
    !>  [maxloc](@ref pm_complexMinMax::maxloc)<br>
    !>  [pm_complexCompareAll](@ref pm_complexCompareAll)<br>
    !>  [pm_complexCompareAny](@ref pm_complexCompareAny)<br>
    !>  [pm_complexCompareLex](@ref pm_complexCompareLex)<br>
    !>  [pm_arrayMinMax](@ref pm_arrayMinMax)<br>
    !>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
    !>
    !>  \example{minval}
    !>  \include{lineno} example/pm_complexMinMax/minval/main.F90
    !>  \compilef{minval}
    !>  \output{minval}
    !>  \include{lineno} example/pm_complexMinMax/minval/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexMinMax](@ref test_pm_complexMinMax)
    !>
    !>  \final{minval}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface minval

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function minvalALL_D1_CK5(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK4_ENABLED
    pure module function minvalALL_D1_CK4(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK3_ENABLED
    pure module function minvalALL_D1_CK3(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK2_ENABLED
    pure module function minvalALL_D1_CK2(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK1_ENABLED
    pure module function minvalALL_D1_CK1(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function minvalALL_D2_CK5(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK4_ENABLED
    module function minvalALL_D2_CK4(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK3_ENABLED
    module function minvalALL_D2_CK3(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK2_ENABLED
    module function minvalALL_D2_CK2(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK1_ENABLED
    module function minvalALL_D2_CK1(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalALL_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function minvalDIM_D1_CK5(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK4_ENABLED
    PURE module function minvalDIM_D1_CK4(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK3_ENABLED
    PURE module function minvalDIM_D1_CK3(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK2_ENABLED
    PURE module function minvalDIM_D1_CK2(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK1_ENABLED
    PURE module function minvalDIM_D1_CK1(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function minvalDIM_D2_CK5(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function minvalDIM_D2_CK4(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function minvalDIM_D2_CK3(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function minvalDIM_D2_CK2(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function minvalDIM_D2_CK1(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minvalDIM_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the component-wise maximum value of (both real and imaginary parts of) the input `complex`.
    !>
    !>  \details
    !>  This generic interfaces extends the behavior of the intrinsic Fortran `maxval(a, dim = dim)` for arguments of `real` type to arguments of type `complex`.<br>
    !>
    !>  \param[in]  array   :   The input array of rank `(1..2)` of type `complex` of kind \CKALL.
    !>  \param[in]  dim     :   The input scalar of type `integer` of default kind \IK of value `1 <= dim <= rank(array)`,
    !>                          representing the dimension along which the maximum must be computed.<br>
    !>                          (**optional**. If missing, the output value is a scalar containing the maximum of the entire input `array`.)
    !>
    !>  \return
    !>  `val`               :   The output object of the same type and kind as the input `array` of rank `rank(array) - 1` whose real and
    !>                          imaginary parts are the maximum values of the corresponding real and imaginary parts of the input `array`.<br>
    !>                          The returned values are `-huge(real(0, kind(array))` if the input `array` has zero size along the dimension of interest.<br>
    !>
    !>  \interface{maxval}
    !>  \code{.F90}
    !>
    !>      use pm_complexMinMax, only: maxval
    !>
    !>      val = maxval(array) ! = cmplx(maxval(array%re), maxval(array%im), kind(array))
    !>      val = maxval(array, dim) ! = cmplx(maxval(array%re, dim), maxval(array%im, dim), kind(array))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [min](@ref pm_complexMinMax::min)<br>
    !>  [max](@ref pm_complexMinMax::max)<br>
    !>  [minval](@ref pm_complexMinMax::minval)<br>
    !>  [maxval](@ref pm_complexMinMax::maxval)<br>
    !>  [minloc](@ref pm_complexMinMax::minloc)<br>
    !>  [maxloc](@ref pm_complexMinMax::maxloc)<br>
    !>  [pm_complexCompareAll](@ref pm_complexCompareAll)<br>
    !>  [pm_complexCompareAny](@ref pm_complexCompareAny)<br>
    !>  [pm_complexCompareLex](@ref pm_complexCompareLex)<br>
    !>  [pm_arrayMinMax](@ref pm_arrayMinMax)<br>
    !>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
    !>
    !>  \example{maxval}
    !>  \include{lineno} example/pm_complexMinMax/maxval/main.F90
    !>  \compilef{maxval}
    !>  \output{maxval}
    !>  \include{lineno} example/pm_complexMinMax/maxval/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexMinMax](@ref test_pm_complexMinMax)
    !>
    !>  \final{maxval}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface maxval

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function maxvalALL_D1_CK5(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK4_ENABLED
    pure module function maxvalALL_D1_CK4(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK3_ENABLED
    pure module function maxvalALL_D1_CK3(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK2_ENABLED
    pure module function maxvalALL_D1_CK2(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK1_ENABLED
    pure module function maxvalALL_D1_CK1(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function maxvalALL_D2_CK5(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK4_ENABLED
    module function maxvalALL_D2_CK4(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK3_ENABLED
    module function maxvalALL_D2_CK3(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK2_ENABLED
    module function maxvalALL_D2_CK2(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK1_ENABLED
    module function maxvalALL_D2_CK1(array) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalALL_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        complex(CKG)                                        :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function maxvalDIM_D1_CK5(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK4_ENABLED
    PURE module function maxvalDIM_D1_CK4(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK3_ENABLED
    PURE module function maxvalDIM_D1_CK3(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK2_ENABLED
    PURE module function maxvalDIM_D1_CK2(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

#if CK1_ENABLED
    PURE module function maxvalDIM_D1_CK1(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        complex(CKG)                                        :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function maxvalDIM_D2_CK5(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function maxvalDIM_D2_CK4(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function maxvalDIM_D2_CK3(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function maxvalDIM_D2_CK2(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function maxvalDIM_D2_CK1(array, dim) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxvalDIM_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        complex(CKG)                                        :: val(size(array, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the component-wise location of the minimum value of (both real and imaginary parts of) the input `complex`.
    !>
    !>  \details
    !>  This generic interfaces extends the behavior of the intrinsic Fortran `minloc(a, dim = dim)` for arguments of `real` type to arguments of type `complex`.<br>
    !>
    !>  \param[in]  array   :   The input array of rank `(1..2)` of type `complex` of kind \CKALL.
    !>  \param[in]  dim     :   The input scalar of type `integer` of default kind \IK of value `1 <= dim <= rank(array)`,
    !>                          representing the dimension along which the location of the minimum must be computed.<br>
    !>                          (**optional**. If missing, the output value is a scalar containing the location of the minimum of the entire input `array`.)
    !>
    !>  \return
    !>  `val`               :   The output,
    !>                          <ol>
    !>                              <li>    vector of size `2`, if `array` is of rank `1` or if the input argument `dim` is missing,
    !>                              <li>    array of shape `(1:2, size(array, 3 - dim))`, if `array` is of rank `2` and the input argument `dim` is present,
    !>                          </ol>
    !>                          of type `integer` of default kind \IK.<br>
    !>                          <ol>
    !>                              <li>    The elements in the first  row of `loc` contain the locations of the first occurrences of the minimum values of the corresponding real      parts of the input `array`.<br>
    !>                              <li>    The elements in the second row of `loc` contain the locations of the first occurrences of the minimum values of the corresponding imaginary parts of the input `array`.<br>
    !>                          </ol>
    !>                          The returned locations are `0` if the input array has a zero size.<br>
    !>
    !>  \interface{minloc}
    !>  \code{.F90}
    !>
    !>      use pm_complexMinMax, only: minloc
    !>
    !>      loc(1:2) = minloc(array) ! = [minloc(array%re), minloc(array%im)]
    !>      loc(1:2) = minloc(array(:), dim) ! = [minloc(array%re, dim), minloc(array%im, dim)]
    !>      loc(1:2, 1:size(array, 3 - dim)) = minloc(array(:,:), dim)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [min](@ref pm_complexMinMax::min)<br>
    !>  [max](@ref pm_complexMinMax::max)<br>
    !>  [minval](@ref pm_complexMinMax::minval)<br>
    !>  [maxval](@ref pm_complexMinMax::maxval)<br>
    !>  [minloc](@ref pm_complexMinMax::minloc)<br>
    !>  [maxloc](@ref pm_complexMinMax::maxloc)<br>
    !>  [pm_complexCompareAll](@ref pm_complexCompareAll)<br>
    !>  [pm_complexCompareAny](@ref pm_complexCompareAny)<br>
    !>  [pm_complexCompareLex](@ref pm_complexCompareLex)<br>
    !>  [pm_arrayMinMax](@ref pm_arrayMinMax)<br>
    !>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
    !>
    !>  \example{minloc}
    !>  \include{lineno} example/pm_complexMinMax/minloc/main.F90
    !>  \compilef{minloc}
    !>  \output{minloc}
    !>  \include{lineno} example/pm_complexMinMax/minloc/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexMinMax](@ref test_pm_complexMinMax)
    !>
    !>  \todo
    !>  \pmed
    !>  This generic interface can be extended with optional `mask` and `back` arguments to match those of the intrinsic `minloc()`.<br>
    !>
    !>  \final{minloc}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface minloc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function minlocALL_D1_CK5(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK4_ENABLED
    pure module function minlocALL_D1_CK4(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK3_ENABLED
    pure module function minlocALL_D1_CK3(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK2_ENABLED
    pure module function minlocALL_D1_CK2(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK1_ENABLED
    pure module function minlocALL_D1_CK1(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function minlocALL_D2_CK5(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK4_ENABLED
    module function minlocALL_D2_CK4(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK3_ENABLED
    module function minlocALL_D2_CK3(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK2_ENABLED
    module function minlocALL_D2_CK2(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK1_ENABLED
    module function minlocALL_D2_CK1(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocALL_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function minlocDIM_D1_CK5(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK4_ENABLED
    PURE module function minlocDIM_D1_CK4(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK3_ENABLED
    PURE module function minlocDIM_D1_CK3(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK2_ENABLED
    PURE module function minlocDIM_D1_CK2(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK1_ENABLED
    PURE module function minlocDIM_D1_CK1(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function minlocDIM_D2_CK5(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function minlocDIM_D2_CK4(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function minlocDIM_D2_CK3(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function minlocDIM_D2_CK2(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function minlocDIM_D2_CK1(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: minlocDIM_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the component-wise location of the maximum value of (both real and imaginary parts of) the input `complex`.
    !>
    !>  \details
    !>  This generic interfaces extends the behavior of the intrinsic Fortran `maxloc(a, dim = dim)` for arguments of `real` type to arguments of type `complex`.<br>
    !>
    !>  \param[in]  array   :   The input array of rank `(1..2)` of type `complex` of kind \CKALL.
    !>  \param[in]  dim     :   The input scalar of type `integer` of default kind \IK of value `1 <= dim <= rank(array)`,
    !>                          representing the dimension along which the location of the maximum must be computed.<br>
    !>                          (**optional**. If missing, the output value is a scalar containing the location of the maximum of the entire input `array`.)
    !>
    !>  \return
    !>  `val`               :   The output,
    !>                          <ol>
    !>                              <li>    vector of size `2`, if `array` is of rank `1` or if the input argument `dim` is missing,
    !>                              <li>    array of shape `(1:2, size(array, 3 - dim))`, if `array` is of rank `2` and the input argument `dim` is present,
    !>                          </ol>
    !>                          of type `integer` of default kind \IK.<br>
    !>                          <ol>
    !>                              <li>    The elements in the first  row of `loc` contain the locations of the first occurrences of the maximum values of the corresponding real      parts of the input `array`.<br>
    !>                              <li>    The elements in the second row of `loc` contain the locations of the first occurrences of the maximum values of the corresponding imaginary parts of the input `array`.<br>
    !>                          </ol>
    !>                          The returned locations are `0` if the input array has a zero size.<br>
    !>
    !>  \interface{maxloc}
    !>  \code{.F90}
    !>
    !>      use pm_complexMinMax, only: maxloc
    !>
    !>      loc(1:2) = maxloc(array) ! = [maxloc(array%re), maxloc(array%im)]
    !>      loc(1:2) = maxloc(array(:), dim) ! = [maxloc(array%re, dim), maxloc(array%im, dim)]
    !>      loc(1:2, 1:size(array, 3 - dim)) = maxloc(array(:,:), dim)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [min](@ref pm_complexMinMax::min)<br>
    !>  [max](@ref pm_complexMinMax::max)<br>
    !>  [minval](@ref pm_complexMinMax::minval)<br>
    !>  [maxval](@ref pm_complexMinMax::maxval)<br>
    !>  [minloc](@ref pm_complexMinMax::minloc)<br>
    !>  [maxloc](@ref pm_complexMinMax::maxloc)<br>
    !>  [pm_complexCompareAll](@ref pm_complexCompareAll)<br>
    !>  [pm_complexCompareAny](@ref pm_complexCompareAny)<br>
    !>  [pm_complexCompareLex](@ref pm_complexCompareLex)<br>
    !>  [pm_arrayMinMax](@ref pm_arrayMinMax)<br>
    !>  [pm_mathMinMax](@ref pm_mathMinMax)<br>
    !>
    !>  \example{maxloc}
    !>  \include{lineno} example/pm_complexMinMax/maxloc/main.F90
    !>  \compilef{maxloc}
    !>  \output{maxloc}
    !>  \include{lineno} example/pm_complexMinMax/maxloc/main.out.F90
    !>
    !>  \test
    !>  [test_pm_complexMinMax](@ref test_pm_complexMinMax)
    !>
    !>  \todo
    !>  \pmed
    !>  This generic interface can be extended with optional `mask` and `back` arguments to match those of the intrinsic `maxloc()`.<br>
    !>
    !>  \final{maxloc}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface maxloc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function maxlocALL_D1_CK5(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK4_ENABLED
    pure module function maxlocALL_D1_CK4(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK3_ENABLED
    pure module function maxlocALL_D1_CK3(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK2_ENABLED
    pure module function maxlocALL_D1_CK2(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK1_ENABLED
    pure module function maxlocALL_D1_CK1(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function maxlocALL_D2_CK5(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK4_ENABLED
    module function maxlocALL_D2_CK4(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK3_ENABLED
    module function maxlocALL_D2_CK3(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK2_ENABLED
    module function maxlocALL_D2_CK2(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK1_ENABLED
    module function maxlocALL_D2_CK1(array) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocALL_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous, target    :: array(:,:)
        integer(IK)                                         :: loc(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function maxlocDIM_D1_CK5(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK4_ENABLED
    PURE module function maxlocDIM_D1_CK4(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK3_ENABLED
    PURE module function maxlocDIM_D1_CK3(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK2_ENABLED
    PURE module function maxlocDIM_D1_CK2(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

#if CK1_ENABLED
    PURE module function maxlocDIM_D1_CK1(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:)
        integer(IK)                                         :: loc(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function maxlocDIM_D2_CK5(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D2_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function maxlocDIM_D2_CK4(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D2_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function maxlocDIM_D2_CK3(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D2_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function maxlocDIM_D2_CK2(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D2_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function maxlocDIM_D2_CK1(array, dim) result(loc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: maxlocDIM_D2_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                            :: dim
        complex(CKG), intent(in)    __CONTIGUOUS            :: array(:,:)
        integer(IK)                                         :: loc(2, size(array, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_complexMinMax ! LCOV_EXCL_LINE