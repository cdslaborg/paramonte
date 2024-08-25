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
!>  This module contains procedures and generic interfaces for computing `sum(x)`
!>  accurately when `x` is a long array of wildly varying `real` or `complex` values.<br>
!>
!>  \details
!>  Some of the algorithms of this module are inspired by the work of
!>  Blanchard, Higham, and Mary, 2020, A Class of Fast and Accurate Summation Algorithms.<br>
!>
!>  \see
!>  [pm_mathCumSum](@ref pm_mathCumSum)<br>
!>
!>  \benchmarks
!>
!>  \benchmark{getSum, The effects of `method` on runtime efficiency}
!>  The following program compares the runtime performance of [getSum](@ref pm_mathSum::getSum)
!>  algorithms with the default Fortran `sum()` intrinsic function.<br>
!>
!>  \include{lineno} benchmark/pm_mathSum/getSum/main.F90
!>  \compilefb{getSum}
!>  \postprocb{getSum}
!>  \include{lineno} benchmark/pm_mathSum/getSum/main.py
!>  \visb{getSum}
!>  \image html benchmark/pm_mathSum/getSum/benchmark.getSum.runtime.png width=1000
!>  \image html benchmark/pm_mathSum/getSum/benchmark.getSum.runtime.ratio.png width=1000
!>  \moralb{getSum}
!>      -#  Among all summation algorithms, [fablocked_type](@ref pm_mathSum::fablocked_type)
!>          appears to offer the most accurate result while also being even faster than the
!>          default Fortran `sum()` and all other implemented summation algorithms for array sizes `> 100`.<br>
!>
!>  \benchmark{getDot, The effects of `method` on runtime efficiency}
!>  The following program compares the runtime performance of [getDot](@ref pm_mathSum::getDot)
!>  algorithms with the default Fortran `dot_product()` intrinsic function.<br>
!>
!>  \include{lineno} benchmark/pm_mathSum/getDot/main.F90
!>  \compilefb{getDot}
!>  \postprocb{getDot}
!>  \include{lineno} benchmark/pm_mathSum/getDot/main.py
!>  \visb{getDot}
!>  \image html benchmark/pm_mathSum/getDot/benchmark.getDot.runtime.png width=1000
!>  \image html benchmark/pm_mathSum/getDot/benchmark.getDot.runtime.ratio.png width=1000
!>  \moralb{getDot}
!>      -#  Among all summation algorithms, [fablocked_type](@ref pm_mathSum::fablocked_type)
!>          appears to offer the most accurate result while also being even faster than the
!>          default Fortran `dot_product()` and all other implemented summation algorithms for array sizes `> 100`.<br>
!>
!>  \final
!>
!>  \test
!>  [test_pm_mathSum](@ref test_pm_mathSum)
!>
!>  \author
!>  \AmirShahmoradi, August 8, 2024, 8:45 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathSum

    use pm_kind, only: SK, IK
    use pm_control, only: iteration_type, iteration
    use pm_control, only: recursion_type, recursion

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathSum"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type, abstract :: sum_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  the **Fast Accurate Blocked summation** method of Blanchard, Higham, and Mary, 2020
    !>  within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [fablocked](@ref pm_mathSum::fablocked)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [fablocked](@ref pm_mathSum::fablocked)<br>
    !>  [kahanbabu](@ref pm_mathSum::kahanbabu)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [fablocked_type](@ref pm_mathSum::fablocked_type)<br>
    !>  [kahanbabu_type](@ref pm_mathSum::kahanbabu_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>
    !>  \final{fablocked_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, August 8, 2024, 8:45 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
    type, extends(sum_type) :: fablocked_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a scalar `parameter` object of type [fablocked_type](@ref pm_mathSum::fablocked_type).<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [fablocked](@ref pm_mathSum::fablocked)<br>
    !>  [kahanbabu](@ref pm_mathSum::kahanbabu)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [fablocked_type](@ref pm_mathSum::fablocked_type)<br>
    !>  [kahanbabu_type](@ref pm_mathSum::kahanbabu_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>
    !>  \final{fablocked}
    !>
    !>  \author
    !>  \AmirShahmoradi, August 8, 2024, 8:45 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
    type(fablocked_type), parameter :: fablocked = fablocked_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: fablocked
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  the **Naive Blocked summation** method of Blanchard, Higham, and Mary, 2020
    !>  within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [nablocked](@ref pm_mathSum::nablocked)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [fablocked](@ref pm_mathSum::fablocked)<br>
    !>  [nablocked](@ref pm_mathSum::nablocked)<br>
    !>  [kahanbabu](@ref pm_mathSum::kahanbabu)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [fablocked_type](@ref pm_mathSum::fablocked_type)<br>
    !>  [nablocked_type](@ref pm_mathSum::nablocked_type)<br>
    !>  [kahanbabu_type](@ref pm_mathSum::kahanbabu_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>
    !>  \final{nablocked_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, August 8, 2024, 8:45 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
    type, extends(sum_type) :: nablocked_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a scalar `parameter` object of type [nablocked_type](@ref pm_mathSum::nablocked_type).<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [fablocked](@ref pm_mathSum::fablocked)<br>
    !>  [kahanbabu](@ref pm_mathSum::kahanbabu)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [fablocked_type](@ref pm_mathSum::fablocked_type)<br>
    !>  [kahanbabu_type](@ref pm_mathSum::kahanbabu_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>
    !>  \final{nablocked}
    !>
    !>  \author
    !>  \AmirShahmoradi, August 8, 2024, 8:45 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
    type(nablocked_type), parameter :: nablocked = nablocked_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: nablocked
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to request
    !>  the **Kahan-Babuska summation** method within the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [kahanbabu](@ref pm_mathSum::kahanbabu)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [fablocked](@ref pm_mathSum::fablocked)<br>
    !>  [kahanbabu](@ref pm_mathSum::kahanbabu)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [kahanbabu_type](@ref pm_mathSum::kahanbabu_type)<br>
    !>  [fablocked_type](@ref pm_mathSum::fablocked_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>
    !>  \final{kahanbabu_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, August 8, 2024, 8:45 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
    type, extends(sum_type) :: kahanbabu_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a scalar `parameter` object of type [kahanbabu_type](@ref pm_mathSum::kahanbabu_type).<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [fablocked](@ref pm_mathSum::fablocked)<br>
    !>  [kahanbabu](@ref pm_mathSum::kahanbabu)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [fablocked_type](@ref pm_mathSum::fablocked_type)<br>
    !>  [kahanbabu_type](@ref pm_mathSum::kahanbabu_type)<br>
    !>  [iteration_type](@ref pm_control::iteration_type)<br>
    !>  [recursion_type](@ref pm_control::recursion_type)<br>
    !>
    !>  \final{kahanbabu}
    !>
    !>  \author
    !>  \AmirShahmoradi, August 8, 2024, 8:45 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
    type(kahanbabu_type), parameter :: kahanbabu = kahanbabu_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: kahanbabu
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the expression `sum(x)` accurately (almost without roundoff error).<br>
    !>
    !>  \param[in]  x           :   The input vector of arbitrary size, of<br>
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              representing the `x` value whose sum of elements is to be returned.<br>
    !>  \param[in]  method      :   The input scalar object that can be,
    !>                              <ol>
    !>                                  <li>    the constant [iteration](@ref pm_control::iteration) or equivalently,
    !>                                          an object of type [iteration_type](@ref pm_control::iteration_type).<br>
    !>                                          Specifying this option forces the use of normal iterative approach to summing all array elements.<br>
    !>                                          This option is equivalent to the default implementations of the Fortran intrinsic `sum()`.<br>
    !>                                          This approach is the fastest serial method among all, but also generally the least accurate.<br>
    !>                                  <li>    the constant [recursion](@ref pm_control::recursion) or equivalently,
    !>                                          an object of type [recursion_type](@ref pm_control::recursion_type).<br>
    !>                                          Specifying this option forces the use of recursive pairwise approach to summing all array elements.<br>
    !>                                  <li>    the constant [kahanbabu](@ref pm_mathSum::kahanbabu) or equivalently,
    !>                                          an object of type [kahanbabu_type](@ref pm_mathSum::kahanbabu_type).<br>
    !>                                          Specifying this option forces the use of the Kahan-Babuska compensated approach to summing all array elements.<br>
    !>                                          This algorithm, while accurate, can be up to 2-4 times more expensive than the iterative approach discussed above.
    !>                                  <li>    the constant [fablocked](@ref pm_mathSum::fablocked) or equivalently,
    !>                                          an object of type [fablocked_type](@ref pm_mathSum::fablocked_type).<br>
    !>                                          Specifying this option forces the use of the Fast Accurate Blocked approach to summing all array elements.<br>
    !>                                  <li>    the constant [nablocked](@ref pm_mathSum::nablocked) or equivalently,
    !>                                          an object of type [nablocked_type](@ref pm_mathSum::nablocked_type).<br>
    !>                                          Specifying this option forces the use of the Naive Blocked approach to summing all array elements.<br>
    !>                              </ol>
    !>                              The presence of this argument is merely for compile-time resolution of the procedures of this generic interface.<br>
    !>                              (**optional**, default = [fablocked](@ref pm_mathSum::fablocked).)
    !>
    !>  \return
    !>  `sumres`                :   The output scalar of the same type and kind
    !>                              containing the result of summing all elements of the input `x`.<br>
    !>
    !>  \interface{getSum}
    !>  \code{.F90}
    !>
    !>      use pm_mathSum, only: getSum, iteration, recursion, kahanbabu, fablocked, nablocked
    !>
    !>      sumres = getSum(x)
    !>      sumres = getSum(x, method)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The returned value is `0` is the size of the condition `size(x) == 0` holds.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
    !>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
    !>  [getCumSum](@ref pm_mathCumSum::getCumSum)<br>
    !>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
    !>  [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp)<br>
    !>  [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)<br>
    !>
    !>  \example{getSum}
    !>  \include{lineno} example/pm_mathSum/getSum/main.F90
    !>  \compilef{getSum}
    !>  \output{getSum}
    !>  \include{lineno} example/pm_mathSum/getSum/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{getSum, The effects of `method` on runtime efficiency}
    !>  The following program compares the runtime performance of [getSum](@ref pm_mathSum::getSum)
    !>  algorithms with the default Fortran `sum()` intrinsic function.<br>
    !>
    !>  \include{lineno} benchmark/pm_mathSum/getSum/main.F90
    !>  \compilefb{getSum}
    !>  \postprocb{getSum}
    !>  \include{lineno} benchmark/pm_mathSum/getSum/main.py
    !>  \visb{getSum}
    !>  \image html benchmark/pm_mathSum/getSum/benchmark.getSum.runtime.png width=1000
    !>  \image html benchmark/pm_mathSum/getSum/benchmark.getSum.runtime.ratio.png width=1000
    !>  \moralb{getSum}
    !>      -#  Among all summation algorithms, [fablocked_type](@ref pm_mathSum::fablocked_type)
    !>          appears to offer the most accurate result while also being even faster than the
    !>          default Fortran `sum()` and all other implemented summation algorithms for array sizes `> 100`.<br>
    !>
    !>  \test
    !>  [test_pm_mathSum](@ref test_pm_mathSum)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, August 8, 2024, 10:23 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>

    ! Default

    interface getSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSumDef_CK5(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
    end function
#endif

#if CK4_ENABLED
    PURE module function getSumDef_CK4(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
    end function
#endif

#if CK3_ENABLED
    PURE module function getSumDef_CK3(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
    end function
#endif

#if CK2_ENABLED
    PURE module function getSumDef_CK2(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
    end function
#endif

#if CK1_ENABLED
    PURE module function getSumDef_CK1(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSumDef_RK5(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
    end function
#endif

#if RK4_ENABLED
    PURE module function getSumDef_RK4(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
    end function
#endif

#if RK3_ENABLED
    PURE module function getSumDef_RK3(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
    end function
#endif

#if RK2_ENABLED
    PURE module function getSumDef_RK2(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
    end function
#endif

#if RK1_ENABLED
    PURE module function getSumDef_RK1(x) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumDef_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Iteration

    interface getSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSumIte_CK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE module function getSumIte_CK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE module function getSumIte_CK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE module function getSumIte_CK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE module function getSumIte_CK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSumIte_RK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getSumIte_RK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getSumIte_RK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getSumIte_RK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getSumIte_RK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumIte_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Recursion

    interface getSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module function getSumRec_CK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE recursive module function getSumRec_CK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE recursive module function getSumRec_CK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE recursive module function getSumRec_CK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE recursive module function getSumRec_CK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module function getSumRec_RK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE recursive module function getSumRec_RK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE recursive module function getSumRec_RK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE recursive module function getSumRec_RK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE recursive module function getSumRec_RK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumRec_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! fablocked

    interface getSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSumFAB_CK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE module function getSumFAB_CK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE module function getSumFAB_CK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE module function getSumFAB_CK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE module function getSumFAB_CK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSumFAB_RK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getSumFAB_RK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getSumFAB_RK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getSumFAB_RK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getSumFAB_RK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumFAB_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! nablocked

    interface getSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSumNAB_CK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE module function getSumNAB_CK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE module function getSumNAB_CK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE module function getSumNAB_CK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE module function getSumNAB_CK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSumNAB_RK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getSumNAB_RK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getSumNAB_RK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getSumNAB_RK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getSumNAB_RK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumNAB_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! kahanbabu

    interface getSum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSumKAB_CK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE module function getSumKAB_CK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE module function getSumKAB_CK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE module function getSumKAB_CK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE module function getSumKAB_CK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:)
        complex(TKG)                                        :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSumKAB_RK5(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getSumKAB_RK4(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getSumKAB_RK3(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getSumKAB_RK2(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getSumKAB_RK1(x, method) result(sumres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSumKAB_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:)
        real(TKG)                                           :: sumres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the expression `dot_product(x, y)` accurately (almost without roundoff error).<br>
    !>
    !>  \param[in]  x           :   The input `contiguous` vector of arbitrary size, of<br>
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              representing the `x` value whose dot product with `y` is to be returned.<br>
    !>  \param[in]  y           :   The input `contiguous` vector of the same size, type, and kind as the input `x`, of<br>
    !>                              representing the `y` values whose dot product with `x` is to be returned.<br>
    !>  \param[in]  method      :   The input scalar object that can be,
    !>                              <ol>
    !>                                  <li>    the constant [iteration](@ref pm_control::iteration) or equivalently,
    !>                                          an object of type [iteration_type](@ref pm_control::iteration_type).<br>
    !>                                          Specifying this option forces the use of normal iterative approach to summing all array elements.<br>
    !>                                          This option is equivalent to the default implementations of the Fortran intrinsic `sum()`.<br>
    !>                                          This approach is the fastest serial method among all, but also generally the least accurate.<br>
    !>                                  <li>    the constant [recursion](@ref pm_control::recursion) or equivalently,
    !>                                          an object of type [recursion_type](@ref pm_control::recursion_type).<br>
    !>                                          Specifying this option forces the use of recursive pairwise approach to summing all array elements.<br>
    !>                                  <li>    the constant [kahanbabu](@ref pm_mathSum::kahanbabu) or equivalently,
    !>                                          an object of type [kahanbabu_type](@ref pm_mathSum::kahanbabu_type).<br>
    !>                                          Specifying this option forces the use of the Kahan-Babuska compensated approach to summing all array elements.<br>
    !>                                          This algorithm, while accurate, can be up to 2-4 times more expensive than the iterative approach discussed above.
    !>                                  <li>    the constant [fablocked](@ref pm_mathSum::fablocked) or equivalently,
    !>                                          an object of type [fablocked_type](@ref pm_mathSum::fablocked_type).<br>
    !>                                          Specifying this option forces the use of the Fast Accurate Blocked approach to summing all array elements.<br>
    !>                                  <li>    the constant [nablocked](@ref pm_mathSum::nablocked) or equivalently,
    !>                                          an object of type [nablocked_type](@ref pm_mathSum::nablocked_type).<br>
    !>                                          Specifying this option forces the use of the Naive Blocked approach to summing all array elements.<br>
    !>                              </ol>
    !>                              The presence of this argument is merely for compile-time resolution of the procedures of this generic interface.<br>
    !>                              (**optional**, default = [fablocked](@ref pm_mathSum::fablocked).)
    !>
    !>  \return
    !>  `dotres`                :   The output scalar of the same type and kind
    !>                              containing the result of the dot product of the input `x` and `y` vectors.<br>
    !>
    !>  \interface{getSum}
    !>  \code{.F90}
    !>
    !>      use pm_mathSum, only: getDot, iteration, recursion, kahanbabu, fablocked, nablocked
    !>
    !>      dotres = getDot(x, y)
    !>      dotres = getDot(x, y, method)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The returned value is `0` is the size of the condition `size(x) == 0` holds.<br>
    !>  The condition `size(x) == size(y)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
    !>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
    !>  [getCumSum](@ref pm_mathCumSum::getCumSum)<br>
    !>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
    !>  [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp)<br>
    !>  [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)<br>
    !>
    !>  \example{getSum}
    !>  \include{lineno} example/pm_mathSum/getDot/main.F90
    !>  \compilef{getSum}
    !>  \output{getSum}
    !>  \include{lineno} example/pm_mathSum/getDot/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{getDot, The effects of `method` on runtime efficiency}
    !>  The following program compares the runtime performance of [getDot](@ref pm_mathSum::getDot)
    !>  algorithms with the default Fortran `dot_product()` intrinsic function.<br>
    !>
    !>  \include{lineno} benchmark/pm_mathSum/getDot/main.F90
    !>  \compilefb{getDot}
    !>  \postprocb{getDot}
    !>  \include{lineno} benchmark/pm_mathSum/getDot/main.py
    !>  \visb{getDot}
    !>  \image html benchmark/pm_mathSum/getDot/benchmark.getDot.runtime.png width=1000
    !>  \image html benchmark/pm_mathSum/getDot/benchmark.getDot.runtime.ratio.png width=1000
    !>  \moralb{getDot}
    !>      -#  Among all summation algorithms, [fablocked_type](@ref pm_mathSum::fablocked_type)
    !>          appears to offer the most accurate result while also being even faster than the
    !>          default Fortran `dot_product()` and all other implemented summation algorithms for array sizes `> 100`.<br>
    !>
    !>  \test
    !>  [test_pm_mathSum](@ref test_pm_mathSum)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, August 8, 2024, 10:23 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>

    ! Default

    interface getDot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDotDef_CK5(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
    end function
#endif

#if CK4_ENABLED
    PURE module function getDotDef_CK4(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
    end function
#endif

#if CK3_ENABLED
    PURE module function getDotDef_CK3(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
    end function
#endif

#if CK2_ENABLED
    PURE module function getDotDef_CK2(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
    end function
#endif

#if CK1_ENABLED
    PURE module function getDotDef_CK1(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDotDef_RK5(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
    end function
#endif

#if RK4_ENABLED
    PURE module function getDotDef_RK4(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
    end function
#endif

#if RK3_ENABLED
    PURE module function getDotDef_RK3(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
    end function
#endif

#if RK2_ENABLED
    PURE module function getDotDef_RK2(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
    end function
#endif

#if RK1_ENABLED
    PURE module function getDotDef_RK1(x, y) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotDef_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Iteration

    interface getDot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDotIte_CK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE module function getDotIte_CK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE module function getDotIte_CK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE module function getDotIte_CK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE module function getDotIte_CK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDotIte_RK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDotIte_RK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDotIte_RK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDotIte_RK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDotIte_RK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotIte_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(iteration_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Recursion

    interface getDot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module function getDotRec_CK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE recursive module function getDotRec_CK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE recursive module function getDotRec_CK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE recursive module function getDotRec_CK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE recursive module function getDotRec_CK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module function getDotRec_RK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE recursive module function getDotRec_RK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE recursive module function getDotRec_RK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE recursive module function getDotRec_RK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE recursive module function getDotRec_RK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotRec_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(recursion_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! fablocked

    interface getDot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDotFAB_CK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE module function getDotFAB_CK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE module function getDotFAB_CK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE module function getDotFAB_CK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE module function getDotFAB_CK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDotFAB_RK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDotFAB_RK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDotFAB_RK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDotFAB_RK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDotFAB_RK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotFAB_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(fablocked_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! nablocked

    interface getDot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDotNAB_CK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE module function getDotNAB_CK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE module function getDotNAB_CK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE module function getDotNAB_CK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE module function getDotNAB_CK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDotNAB_RK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDotNAB_RK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDotNAB_RK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDotNAB_RK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDotNAB_RK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotNAB_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(nablocked_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! kahanbabu

    interface getDot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDotKAB_CK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    PURE module function getDotKAB_CK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    PURE module function getDotKAB_CK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    PURE module function getDotKAB_CK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    PURE module function getDotKAB_CK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)        , intent(in)    , contiguous    :: x(:), y(:)
        complex(TKG)                                        :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDotKAB_RK5(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getDotKAB_RK4(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getDotKAB_RK3(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getDotKAB_RK2(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getDotKAB_RK1(x, y, method) result(dotres)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDotKAB_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)           , intent(in)    , contiguous    :: x(:), y(:)
        real(TKG)                                           :: dotres
        type(kahanbabu_type), intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathSum ! LCOV_EXCL_LINE