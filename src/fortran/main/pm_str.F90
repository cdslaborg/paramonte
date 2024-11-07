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
!>  This module contains classes and procedures for various string manipulations and inquiries.
!>
!>  \devnote
!>  This module and its generic interfaces could have been written for strings of all kinds.
!>  However, there is no clear usage for the extended functionality of wrapping non-string arrays.
!>  Consequently, it is implemented specifically for strings.
!>
!>  \devnote
!>  Do **not** change the double back-ticks in <pre>``"(g0,:,',')"``</pre> to single back-ticks in any documentations
!>  in this module. Doxygen version 1.9 has difficultly parsing the code sections with single back-tick when the code contains
!>  advanced features of modern Fortran `g0` edit descriptor.
!>
!>  \see
!>  [pm_val2str](@ref pm_val2str)<br>
!>  [pm_strANSI](@ref pm_strANSI)<br>
!>  [pm_strASCII](@ref pm_strASCII)<br>
!>  [pm_arrayPad](@ref pm_arrayPad)<br>
!>  [pm_arrayFind](@ref pm_arrayFind)<br>
!>  [pm_arraySplit](@ref pm_arraySplit)<br>
!>  [pm_arrayCenter](@ref pm_arrayCenter)<br>
!>
!>  \test
!>  [test_pm_str](@ref test_pm_str)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_str

    use pm_kind, only: SK, IK, LK
    use pm_container, only: css_type
#if PDT_ENABLED
    use pm_container, only: css_pdt
#endif

!#if SK5_ENABLED
!    use pm_arraySplit, only: setSplit_D0_D0_SK5, setSplitCustom_D0_D0_SK5, setSplitInstance_D0_D0_SK5, setSplitCustomInstance_D0_D0_SK5
!#endif
!
!#if SK4_ENABLED
!    use pm_arraySplit, only: setSplit_D0_D0_SK4, setSplitCustom_D0_D0_SK4, setSplitInstance_D0_D0_SK4, setSplitCustomInstance_D0_D0_SK4
!#endif
!
!#if SK3_ENABLED
!    use pm_arraySplit, only: setSplit_D0_D0_SK3, setSplitCustom_D0_D0_SK3, setSplitInstance_D0_D0_SK3, setSplitCustomInstance_D0_D0_SK3
!#endif
!
!#if SK2_ENABLED
!    use pm_arraySplit, only: setSplit_D0_D0_SK2, setSplitCustom_D0_D0_SK2, setSplitInstance_D0_D0_SK2, setSplitCustomInstance_D0_D0_SK2
!#endif
!
!#if SK1_ENABLED
!    use pm_arraySplit, only: setSplit_D0_D0_SK1, setSplitCustom_D0_D0_SK1, setSplitInstance_D0_D0_SK1, setSplitCustomInstance_D0_D0_SK1
!#endif

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_str"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    character(*, SK), parameter :: SPACE = SK_" "                               !<  \public The space `character` of default kind \SK of length `1`.
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: SPACE
#endif

    character(*, SK), parameter :: NLC = new_line(SK_"a")                       !<  \public The newline `character` of default kind \SK as returned by `new_line(SK_"a")`.
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: NLC
#endif

    character(*, SK), parameter :: UNDEFINED = SK_"UNDEFINED"                   !<  \public The scalar `character` parameter of default kind \SK containing `"UNDEFINED"`.
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: UNDEFINED
#endif

    ! Intel compilers cannot apply `minval()` or `maxval()` to empty character array constants as of 2022.
    ! Until this internal compiler bug is resolved, the following constants are commented out.
   !character(*, SK), parameter :: COL_SEQ_BEG = maxval([character(1,SK)::],1)  !<  \public The scalar `character` parameter of default kind \SK containing the first character in the processor collating sequence.
   !character(*, SK), parameter :: COL_SEQ_END = minval([character(1,SK)::],1)  !<  \public The scalar `character` parameter of default kind \SK containing the final character in the processor collating sequence.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_str_alleq
    !>  Generate and return `.true.` if the input scalar strings are equal (non-lexically), otherwise, return `.false.`.
    !>
    !>  \brief
    !>  The procedures of this generic interface return `.true.` if and only if,
    !>  <ol>
    !>      <li>    the two input scalar strings have the same value and length type parameter, or
    !>      <li>    the length of at least one of the input strings is `1` and
    !>              each character of the other string equals the input string
    !>              with the length of `1`.
    !>  </ol>
    !>  The functionality of this generic interface is equivalent to the Fortran intrinsic `all()` for array entities.
    !>
    !>  \param[in]  str1    :   The input scalar `character` of kind \SKALL of arbitrary length type parameter.
    !>  \param[in]  str2    :   The input scalar `character` of the same kind as `str1` of arbitrary length type parameter.
    !>
    !>  \return
    !>  `allEqual`          :   The output scalar logical of default kind \LK that is `.true.` if and only if,
    !>                          <ol>
    !>                              <li>    the two input scalar strings have the same value and length type parameter, or
    !>                              <li>    the length of at least one of the input strings is `1` and
    !>                                      each character of the other string equals the input string
    !>                                      with the length of `1`.
    !>                          </ol>
    !>
    !>  \interface{alleq}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_str, only: operator(.alleq.)
    !>      logical(LK) :: allEqual
    !>
    !>      allEqual = str1 .alleq. str2
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [pm_arrayCompareLex](@ref pm_arrayCompareLex)<br>
    !>
    !>  \example{alleq}
    !>  \include{lineno} example/pm_str/alleq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_str/alleq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface operator(.alleq.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function alleq_SK5(str1, str2) result(allEqual)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: alleq_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str1, str2
        logical(LK)                                 :: allEqual
    end function
#endif

#if SK4_ENABLED
    pure elemental module function alleq_SK4(str1, str2) result(allEqual)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: alleq_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str1, str2
        logical(LK)                                 :: allEqual
    end function
#endif

#if SK3_ENABLED
    pure elemental module function alleq_SK3(str1, str2) result(allEqual)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: alleq_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str1, str2
        logical(LK)                                 :: allEqual
    end function
#endif

#if SK2_ENABLED
    pure elemental module function alleq_SK2(str1, str2) result(allEqual)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: alleq_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str1, str2
        logical(LK)                                 :: allEqual
    end function
#endif

#if SK1_ENABLED
    pure elemental module function alleq_SK1(str1, str2) result(allEqual)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: alleq_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str1, str2
        logical(LK)                                 :: allEqual
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the location of the first occurrence of the minimum of the input string(s).
    !>
    !>  \param[in]  str         :   The input scalar (or array of the same rank and shape as other array-like input arguments) of
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary `len` type parameter,
    !>                                  <li>    type string container [css_pdt](@ref pm_container::css_pdt) of kind \SKALL,
    !>                                  <li>    type string container [css_type](@ref pm_container::css_type) of default kind \SK,
    !>                              </ol>
    !>                              containing the string whose end must be checked for equivalence with the input `suffix`.
    !>  `endedWith`             :   The output scalar or array of the same shape as the input `str` of type `integer` of default
    !>  \param[in]  suffix      :   The input scalar of type `logical` of default kind \LK.<br>
    !>                              If present and `.true.`, then the location of the last element of minimum value in `str` will be output.<br>
    !>                              (**optional**, default = `.false.`)
    !>
    !>  \return
    !>  `endedWith`             :   The output scalar (or array of the same rank and shape as any array-like input argument) of type `logical` of default kind \LK
    !>                              that is `.true.` if and only if the input `str` ends with the specified input `suffix`.<br>
    !>
    !>  \interface{isEndedWith}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: isEndedWith
    !>
    !>      endedWith = isEndedWith(str, suffix)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This generic interface always yields `.true.` for an input empty `suffix`.<br>
    !>
    !>  \remark
    !>  Individual characters are compared using the Fortran intrinsic relational operators for objects of type `character`.<br>
    !>
    !>  \pure
    !>
    !>  \example{isEndedWith}
    !>  \include{lineno} example/pm_str/isEndedWith/main.F90
    !>  \compilef{isEndedWith}
    !>  \output{isEndedWith}
    !>  \include{lineno} example/pm_str/isEndedWith/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.11.1 Build 20231117_000000}
    !>  \desc
    !>  The following example fails on Windows and Linux with \ifort{2021.11.1 Build 20231117_000000}.
    !>  \code{.F90}
    !>      call disp%show("isEndedWith(css_type('ParaMonte'), [css_type('ParaMonte'), css_type('Monte  '), css_type('Monte  ', trimmed = .true._LK), css_type('monte')])")
    !>      call disp%show( isEndedWith(css_type('ParaMonte'), [css_type('ParaMonte'), css_type('Monte  '), css_type('Monte  ', trimmed = .true._LK), css_type('monte')]) )
    !>  \endcode
    !>  with the following message:
    !>  \code{.sh}
    !>      forrtl: severe (408): fort: (8): Attempt to fetch from allocatable variable VAL when it is not allocated.
    !>  \endcode
    !>  Attempts to reproduce this error in an isolated simple case were unsuccessful.<br>
    !>  The code successfully compiles and runs using \gfortran.
    !>  \remedy
    !>  For now, the above example is commented out in the example script.<br>
    !>
    !>  \final{isEndedWith}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isEndedWith

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isEndedWith_SK5(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)    , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isEndedWith_SK4(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)    , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isEndedWith_SK3(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)    , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isEndedWith_SK2(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)    , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isEndedWith_SK1(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)    , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure elemental module function isEndedWith_PSSK5(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        type(css_pdt(SKG))  , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isEndedWith_PSSK4(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        type(css_pdt(SKG))  , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isEndedWith_PSSK3(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        type(css_pdt(SKG))  , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isEndedWith_PSSK2(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        type(css_pdt(SKG))  , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isEndedWith_PSSK1(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        type(css_pdt(SKG))  , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental module function isEndedWith_BSSK(str, suffix) result(endedWith)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isEndedWith_BSSK
#endif
        use pm_kind, only: SKG => SK
        type(css_type)      , intent(in)            :: str, suffix
        logical(LK)                                 :: endedWith
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the location of the first occurrence of the minimum of the input string(s).
    !>
    !>  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of arbitrary `len` type parameter.
    !>  \param[in]  mask    :   The input vector of type `logical` of default kind \LK of the same length as the `len` type parameter of the input `str`.<br>
    !>                          If present, then only those characters of `str` corresponding to a `.true.` in `mask` will be considered to find `minLoc`.<br>
    !>                          When `mask` is present, `str` must be a scalar (the procedures are non-elemental when `mask` is present).<br>
    !>                          (**optional**, default = `[(.true., i = 1, len(str))]`)
    !>  \param[in]  back    :   The input scalar of type `logical` of default kind \LK.<br>
    !>                          If present and `.true.`, then the location of the last element of minimum value in `str` will be output.<br>
    !>                          (**optional**, default = `.false.`)
    !>
    !>  \return
    !>  `minLoc`            :   The output scalar or array of the same shape as the input `str` of type `integer` of default
    !>                          kind \IK containing the location of the first character of minimum value in the input `str`.<br>
    !>                          If `len(str) == 0` on input, then `minLoc = 1_IK` on output, similar to the behavior of the
    !>                          Fortran intrinsic functions `minloc` and `maxloc` in the Intel Fortran compiler.<br>
    !>
    !>  \interface{getMinLoc}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: getMinLoc
    !>
    !>      minLoc = getMinLoc(str, back = back)
    !>      minLoc = getMinLoc(str, mask, back = back)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input `str` must have non-zero `len` type-parameter.<br>
    !>  The length of the input `mask` vector must equal the length of the input `str`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  Individual characters are compared using the Fortran intrinsic relational operators for objects of type `character`.<br>
    !>
    !>  \note
    !>  The procedures under this generic interface are intended to (partially) overload the
    !>  functionality of the Fortran intrinsic function `minloc` to also accept input scalar string arguments.<br>
    !>
    !>  \see
    !>  [getMinLoc](@ref pm_str::getMinLoc)<br>
    !>  [getMaxLoc](@ref pm_str::getMaxLoc)<br>
    !>  [getMinVal](@ref pm_str::getMinVal)<br>
    !>  [getMaxVal](@ref pm_str::getMaxVal)<br>
    !>
    !>  \example{getMinLoc}
    !>  \include{lineno} example/pm_str/getMinLoc/main.F90
    !>  \compilef{getMaxLoc}
    !>  \output{getMaxLoc}
    !>  \include{lineno} example/pm_str/getMinLoc/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \final{getMaxLoc}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getMinLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module function getMinLocNoMask_SK5(str, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocNoMask_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

#if SK4_ENABLED
    PURE elemental module function getMinLocNoMask_SK4(str, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocNoMask_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

#if SK3_ENABLED
    PURE elemental module function getMinLocNoMask_SK3(str, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocNoMask_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

#if SK2_ENABLED
    PURE elemental module function getMinLocNoMask_SK2(str, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocNoMask_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

#if SK1_ENABLED
    PURE elemental module function getMinLocNoMask_SK1(str, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocNoMask_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMinLocMasked_SK5(str, mask, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocMasked_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

#if SK4_ENABLED
    PURE module function getMinLocMasked_SK4(str, mask, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocMasked_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

#if SK3_ENABLED
    PURE module function getMinLocMasked_SK3(str, mask, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocMasked_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

#if SK2_ENABLED
    PURE module function getMinLocMasked_SK2(str, mask, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocMasked_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

#if SK1_ENABLED
    PURE module function getMinLocMasked_SK1(str, mask, back) result(minLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinLocMasked_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: minLoc
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the location of the first occurrence of the maximum of the input string(s).<br>
    !>
    !>  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of arbitrary `len` type parameter.<br>
    !>  \param[in]  mask    :   The input vector of type `logical` of default kind \LK of the same length as the `len` type parameter of the input `str`.<br>
    !>                          If present, then only those characters of `str` corresponding to a `.true.` in `mask` will be considered to find `maxLoc`.<br>
    !>                          When `mask` is present, `str` must be a scalar (the procedures are non-elemental when `mask` is present).<br>
    !>                          (**optional**, default = `[(.true., i = 1, len(str))]`)
    !>  \param[in]  back    :   The input scalar of type `logical` of default kind \LK.<br>
    !>                          If present and `.true.`, then the location of the last element of maximum value in `str` will be output.<br>
    !>                          (**optional**, default = `.false.`)
    !>
    !>  \return
    !>  `maxLoc`            :   The output scalar or array of the same shape as the input `str` of type `integer` of default
    !>                          kind \IK containing the location of the first character of maximum value in the input `str`.<br>
    !>                          If `len(str) == 0` on input, then `maxLoc = 1_IK` on output, similar to the behavior of the
    !>                          Fortran intrinsic functions `minloc` and `maxloc` in the Intel Fortran compiler.<br>
    !>
    !>  \interface{getMaxLoc}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: getMaxLoc
    !>
    !>      maxLoc = getMaxLoc(str, back = back)
    !>      maxLoc = getMaxLoc(str, mask, back = back)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input `str` must have non-zero `len` type-parameter.<br>
    !>  The length of the input `mask` vector must equal the length of the input `str`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  Individual characters are compared using the Fortran intrinsic relational operators for objects of type `character`.<br>
    !>
    !>  \note
    !>  The procedures under this generic interface are intended to (partially) overload the
    !>  functionality of the Fortran intrinsic function `maxloc` to also accept input scalar string arguments.<br>
    !>
    !>  \see
    !>  [getMinLoc](@ref pm_str::getMinLoc)<br>
    !>  [getMaxLoc](@ref pm_str::getMaxLoc)<br>
    !>  [getMinVal](@ref pm_str::getMinVal)<br>
    !>  [getMaxVal](@ref pm_str::getMaxVal)<br>
    !>
    !>  \example{getMaxLoc}
    !>  \include{lineno} example/pm_str/getMaxLoc/main.F90
    !>  \compilef{getMaxLoc}
    !>  \output{getMaxLoc}
    !>  \include{lineno} example/pm_str/getMaxLoc/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \final{getMaxLoc}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getMaxLoc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module function getMaxLocNoMask_SK5(str, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocNoMask_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

#if SK4_ENABLED
    PURE elemental module function getMaxLocNoMask_SK4(str, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocNoMask_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

#if SK3_ENABLED
    PURE elemental module function getMaxLocNoMask_SK3(str, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocNoMask_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

#if SK2_ENABLED
    PURE elemental module function getMaxLocNoMask_SK2(str, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocNoMask_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

#if SK1_ENABLED
    PURE elemental module function getMaxLocNoMask_SK1(str, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocNoMask_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMaxLocMasked_SK5(str, mask, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocMasked_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

#if SK4_ENABLED
    PURE module function getMaxLocMasked_SK4(str, mask, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocMasked_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

#if SK3_ENABLED
    PURE module function getMaxLocMasked_SK3(str, mask, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocMasked_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

#if SK2_ENABLED
    PURE module function getMaxLocMasked_SK2(str, mask, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocMasked_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

#if SK1_ENABLED
    PURE module function getMaxLocMasked_SK1(str, mask, back) result(maxLoc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxLocMasked_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        logical(LK)     , intent(in), optional      :: back
        integer(IK)                                 :: maxLoc
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return minimum character value in the input string(s).<br>
    !>
    !>  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of arbitrary `len` type parameter.<br>
    !>  \param[in]  mask    :   The input vector of type `logical` of default kind \LK of the same length as the `len` type parameter of the input `str`.<br>
    !>                          If present, then only those characters of `str` corresponding to a `.true.` in `mask` will be considered to find `minVal`.<br>
    !>                          When `mask` is present, `str` must be a scalar (the procedures are non-elemental when `mask` is present).<br>
    !>                          (**optional**, default = `[(.true., i = 1, len(str))]`)
    !>
    !>  \return
    !>  `minVal`            :   The output scalar or array of the same type, kind, and rank as the input `str` but with
    !>                          a `len` type-parameter of `1` containing the character of minimum value in the input `str`.<br>
    !>
    !>  \interface{getMinVal}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: getMinVal
    !>      character(1, kind(str)) :: minVal
    !>
    !>      minVal = getMinVal(str)
    !>      minVal = getMinVal(str, mask)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input `str` must have non-zero `len` type-parameter.<br>
    !>  The length of the input `mask` vector must equal the length of the input `str`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  Individual characters are compared using the Fortran intrinsic relational operators for objects of type `character`.<br>
    !>
    !>  \note
    !>  The procedures under this generic interface are intended to (partially) overload the
    !>  functionality of the Fortran intrinsic function `minval` to also accept input scalar string arguments.<br>
    !>
    !>  \see
    !>  [getMinLoc](@ref pm_str::getMinLoc)<br>
    !>  [getMaxLoc](@ref pm_str::getMaxLoc)<br>
    !>  [getMinVal](@ref pm_str::getMinVal)<br>
    !>  [getMaxVal](@ref pm_str::getMaxVal)<br>
    !>
    !>  \example{getMinVal}
    !>  \include{lineno} example/pm_str/getMinVal/main.F90
    !>  \compilef{getMinVal}
    !>  \output{getMinVal}
    !>  \include{lineno} example/pm_str/getMinVal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \todo
    !>  \phigh
    !>  The restriction on the non-zero length of the input string can be lifted to make the behavior compatible with the corresponding intrinsic routine.<br>
    !>
    !>  \final{getMinVal}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getMinVal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module function getMinValNoMask_SK5(str) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValNoMask_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: minVal
    end function
#endif

#if SK4_ENABLED
    PURE elemental module function getMinValNoMask_SK4(str) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValNoMask_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: minVal
    end function
#endif

#if SK3_ENABLED
    PURE elemental module function getMinValNoMask_SK3(str) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValNoMask_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: minVal
    end function
#endif

#if SK2_ENABLED
    PURE elemental module function getMinValNoMask_SK2(str) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValNoMask_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: minVal
    end function
#endif

#if SK1_ENABLED
    PURE elemental module function getMinValNoMask_SK1(str) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValNoMask_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: minVal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMinValMasked_SK5(str, mask) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValMasked_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: minVal
    end function
#endif

#if SK4_ENABLED
    PURE module function getMinValMasked_SK4(str, mask) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValMasked_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: minVal
    end function
#endif

#if SK3_ENABLED
    PURE module function getMinValMasked_SK3(str, mask) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValMasked_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: minVal
    end function
#endif

#if SK2_ENABLED
    PURE module function getMinValMasked_SK2(str, mask) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValMasked_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: minVal
    end function
#endif

#if SK1_ENABLED
    PURE module function getMinValMasked_SK1(str, mask) result(minVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinValMasked_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: minVal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return maximum character value in the input string(s).<br>
    !>
    !>  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of arbitrary `len` type parameter.<br>
    !>  \param[in]  mask    :   The input vector of type `logical` of default kind \LK of the same length as the `len` type parameter of the input `str`.<br>
    !>                          If present, then only those characters of `str` corresponding to a `.true.` in `mask` will be considered to find `maxVal`.<br>
    !>                          When `mask` is present, `str` must be a scalar (the procedures are non-elemental when `mask` is present).<br>
    !>                          (**optional**, default = `[(.true., i = 1, len(str))]`)
    !>
    !>  \return
    !>  `maxVal`            :   The output scalar or array of the same type, kind, and rank as the input `str` but with
    !>                          a `len` type-parameter of `1` containing the character of maximum value in the input `str`.<br>
    !>
    !>  \interface{getMaxVal}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: getMaxVal
    !>      character(1, kind(str)) :: maxVal
    !>
    !>      maxVal = getMaxVal(str)
    !>      maxVal = getMaxVal(str, mask)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input `str` must have non-zero `len` type-parameter.<br>
    !>  The length of the input `mask` vector must equal the length of the input `str`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  Individual characters are compared using the Fortran intrinsic relational operators for objects of type `character`.<br>
    !>
    !>  \note
    !>  The procedures under this generic interface are intended to (partially) overload the
    !>  functionality of the Fortran intrinsic function `maxval` to also accept input scalar string arguments.<br>
    !>
    !>  \see
    !>  [getMinLoc](@ref pm_str::getMinLoc)<br>
    !>  [getMaxLoc](@ref pm_str::getMaxLoc)<br>
    !>  [getMinVal](@ref pm_str::getMinVal)<br>
    !>  [getMaxVal](@ref pm_str::getMaxVal)<br>
    !>
    !>  \example{getMaxVal}
    !>  \include{lineno} example/pm_str/getMaxVal/main.F90
    !>  \compilef{getMaxVal}
    !>  \output{getMaxVal}
    !>  \include{lineno} example/pm_str/getMaxVal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \todo
    !>  \phigh
    !>  The restriction on the non-zero length of the input string can be lifted to make the behavior compatible with the corresponding intrinsic routine.<br>
    !>
    !>  \final{getMaxVal}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getMaxVal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE elemental module function getMaxValNoMask_SK5(str) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValNoMask_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: maxVal
    end function
#endif

#if SK4_ENABLED
    PURE elemental module function getMaxValNoMask_SK4(str) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValNoMask_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: maxVal
    end function
#endif

#if SK3_ENABLED
    PURE elemental module function getMaxValNoMask_SK3(str) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValNoMask_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: maxVal
    end function
#endif

#if SK2_ENABLED
    PURE elemental module function getMaxValNoMask_SK2(str) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValNoMask_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: maxVal
    end function
#endif

#if SK1_ENABLED
    PURE elemental module function getMaxValNoMask_SK1(str) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValNoMask_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: maxVal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMaxValMasked_SK5(str, mask) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValMasked_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: maxVal
    end function
#endif

#if SK4_ENABLED
    PURE module function getMaxValMasked_SK4(str, mask) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValMasked_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: maxVal
    end function
#endif

#if SK3_ENABLED
    PURE module function getMaxValMasked_SK3(str, mask) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValMasked_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: maxVal
    end function
#endif

#if SK2_ENABLED
    PURE module function getMaxValMasked_SK2(str, mask) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValMasked_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: maxVal
    end function
#endif

#if SK1_ENABLED
    PURE module function getMaxValMasked_SK1(str, mask) result(maxVal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxValMasked_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str
        logical(LK)     , intent(in), contiguous    :: mask(:)
        character(1,SKG)                            :: maxVal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a vector of single-characters each element of which corresponds to one character in the input scalar string.<br>
    !>
    !>  \param[in]  str :   The input scalar `character` of arbitrary `len` type parameter of kind \SKALL.<br>
    !>
    !>  \return
    !>  `charVec`       :   The output vector of the same type and kind as the input `str` of `len` type-parameter `1`,
    !>                      composed of single characters corresponding to each character in the input scalar `str`.<br>
    !>
    !>  \interface{getCharVec}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: getCharVec
    !>      character(1, kind(str)) :: charVec(len(str))
    !>
    !>      charVec(1:len(str)) = getCharVec(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [getCharSeq](@ref pm_str::getCharSeq)<br>
    !>  [getCharVec](@ref pm_str::getCharVec)<br>
    !>
    !>  \example{getCharVec}
    !>  \include{lineno} example/pm_str/getCharVec/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_str/getCharVec/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getCharVec

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getCharVec_SK5(str) result(charVec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharVec_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: charVec(len(str, IK))
    end function
#endif

#if SK4_ENABLED
    pure module function getCharVec_SK4(str) result(charVec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharVec_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: charVec(len(str, IK))
    end function
#endif

#if SK3_ENABLED
    pure module function getCharVec_SK3(str) result(charVec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharVec_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: charVec(len(str, IK))
    end function
#endif

#if SK2_ENABLED
    pure module function getCharVec_SK2(str) result(charVec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharVec_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: charVec(len(str, IK))
    end function
#endif

#if SK1_ENABLED
    pure module function getCharVec_SK1(str) result(charVec)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharVec_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                :: str
        character(1,SKG)                            :: charVec(len(str, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar string resulting from the concatenation of
    !>  the individual elements (single-characters) of the input vector of characters.<br>
    !>
    !>  \param[in]  charVec :   The input `contiguous` vector of rank `1` of type `character` of default \SKALL
    !>                          of `len` type parameter `1`, each element of which contains a single character.<br>
    !>
    !>  \return
    !>  `string`            :   The output scalar of the same type and kind as `charVec` but with `len` type-parameter of `size(charVec)`
    !>                          whose content results from the concatenation of the single characters in the input `charVec`.<br>
    !>
    !>  \interface{getCharSeq}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: getCharSeq
    !>      character(size(charVec), kind(charVec)) :: str
    !>
    !>      str = getCharSeq(charVec(:))
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [getCharSeq](@ref pm_str::getCharSeq)<br>
    !>  [getCharVec](@ref pm_str::getCharVec)<br>
    !>
    !>  \example{getCharSeq}
    !>  \include{lineno} example/pm_str/getCharSeq/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_str/getCharSeq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getCharSeq

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getCharSeq_SK5(charVec) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharSeq_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG), intent(in), contiguous    :: charVec(:)
        character(size(charVec,1,IK),SKG)           :: str
    end function
#endif

#if SK4_ENABLED
    pure module function getCharSeq_SK4(charVec) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharSeq_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG), intent(in), contiguous    :: charVec(:)
        character(size(charVec,1,IK),SKG)           :: str
    end function
#endif

#if SK3_ENABLED
    pure module function getCharSeq_SK3(charVec) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharSeq_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG), intent(in), contiguous    :: charVec(:)
        character(size(charVec,1,IK),SKG)           :: str
    end function
#endif

#if SK2_ENABLED
    pure module function getCharSeq_SK2(charVec) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharSeq_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG), intent(in), contiguous    :: charVec(:)
        character(size(charVec,1,IK),SKG)           :: str
    end function
#endif

#if SK1_ENABLED
    pure module function getCharSeq_SK1(charVec) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharSeq_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG), intent(in), contiguous    :: charVec(:)
        character(size(charVec,1,IK),SKG)           :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a vector of single-characters each element of which corresponds to one character in the input scalar string.<br>
    !>
    !>  \param[in]  str  :   The input scalar `character` of arbitrary `len` type parameter of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strt`           :   The output vector of the same type and kind as the input `str` of `len` type-parameter `1`,
    !>                          composed of single characters corresponding to each character in the input scalar `str`.<br>
    !>
    !>  \interface{getTrimmedTZ}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: getTrimmedTZ
    !>      character(1, kind(str)) :: strt(len(str))
    !>
    !>      strt(1:len(str)) = getTrimmedTZ(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [pm_arrayStrip](@ref pm_arrayStrip)<br>
    !>
    !>  \example{getTrimmedTZ}
    !>  \include{lineno} example/pm_str/getTrimmedTZ/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_str/getTrimmedTZ/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getTrimmedTZ

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getTrimmedTZ_SK5(str) result(strt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTrimmedTZ_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG), intent(in)                    :: str
        character(:,SKG)                , allocatable   :: strt
    end function
#endif

#if SK4_ENABLED
    pure module function getTrimmedTZ_SK4(str) result(strt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTrimmedTZ_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG), intent(in)                    :: str
        character(:,SKG)                , allocatable   :: strt
    end function
#endif

#if SK3_ENABLED
    pure module function getTrimmedTZ_SK3(str) result(strt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTrimmedTZ_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG), intent(in)                    :: str
        character(:,SKG)                , allocatable   :: strt
    end function
#endif

#if SK2_ENABLED
    pure module function getTrimmedTZ_SK2(str) result(strt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTrimmedTZ_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG), intent(in)                    :: str
        character(:,SKG)                , allocatable   :: strt
    end function
#endif

#if SK1_ENABLED
    pure module function getTrimmedTZ_SK1(str) result(strt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getTrimmedTZ_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG), intent(in)                    :: str
        character(:,SKG)                , allocatable   :: strt
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the input scalar string wrapped within the specified width while
    !>  optionally prefixing each line and respecting the initial indentation in all subsequent lines.<br>
    !>
    !>  \param[in]  str         :   The input scalar of type `character` of kind \SKALL that is to be wrapped.<br>
    !>  \param[in]  prefix      :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter,
    !>                              whose contents will be prefixed to each new line (including the first line) before an indentation occurs.<br>
    !>                              Note that `prefix` does **not** count toward filling the specified `width` or `maxwidth` of each wrapped line.<br>
    !>                              (**optional**, default = `""`)
    !>  \param[in]  indent      :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter, whose contents
    !>                              are the pattern (a set of characters) to search for at the beginning of each paragraph in the input `str`.<br>
    !>                              The repetitions of the `indent` at the beginning of each paragraph (including the beginning of `str`) collectively
    !>                              construct the indentation to be applied to all wrapped lines in the specific paragraph within `str`.<br>
    !>                              For example, the following choice of input arguments,
    !>                              \code{.F90}
    !>                                  strWrapped = getStrWrapped(str = "~+~+This is a test.", indent = "~+")
    !>                              \endcode
    !>                              will lead to prefixing each wrapped line in the output with `~+~+`.<br>
    !>                              If `indent` is as an empty string, then no search for indentation character will be made and no line will be indented.<br>
    !>                              (**optional**, default = `" "`. If `indent = ""`, then no indentation of the lines will occur.)
    !>  \param[in]  break       :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter,
    !>                              whose contents are the set of characters at which breaking the input `str` is allowed.<br>
    !>                              Setting `break = ""` will lead to wrapping the input `str` at any arbitrary character
    !>                              as soon as the length of a given line reaches the specified `width`.<br>
    !>                              (**optional**, default = `" "`)
    !>  \param[in]  newline     :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter,
    !>                              whose contents represent the dividing point between the wrapped lines in the input `str`.<br>
    !>                              For example, if the string contains the C-style newline character \f$\ms{\n}\f$, then set \f$\ms{newline = '\\n'}\f$.<br>
    !>                              Note that `newline` does **not** count toward the `width` or `maxwidth` of the wrapped line.<br>
    !>                              The presence of two adjacent `newline` subsrtings signals the beginning of a new paragraph.<br>
    !>                              In such a case, the indentation of the new paragraph will be computed afresh by couting the
    !>                              number of repetitions of the input `indent` pattern at the beginning of the new paragraph.<br>
    !>                              (**optional**, default = `new_line(SKG_"a")` where `SKG` refers to the kind of `str`.)
    !>  \param[in]  linefeed    :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter,
    !>                              whose contents represent the linefeed to use for dividing `str` into separate wrapped lines.<br>
    !>                              Each line that is to be wrapped will end with a `linefeed`.<br>
    !>                              Additionally, all instances of `newline` in the input `str` will be replaced with `linefeed` in the wrapped `str`.<br>
    !>                              Note that `linefeed` does **not** count toward the `width` or `maxwidth` of the wrapped line.<br>
    !>                              This argument is extremely useful for fast conversion of C-style newline instances to linefeed characters while wrapping `str`.<br>
    !>                              See below for example usages of this input argument.<br>
    !>                              (**optional**, default = the input argument `newline`)
    !>  \param[in]  width       :   The input scalar of type `integer` of default kind \IK, representing the length of
    !>                              each line in the output wrapping, **with the possibility of word overflows**.<br>
    !>                              If a word overflows the specified input line width, it is allowed to continue
    !>                              in the same line until the next `break` character is encountered.<br>
    !>                              Note that the lengths of the input arguments `prefix`, `newline`, and `linefeed` do **not** count toward the `width` or `maxwidth` of the wrapped line.<br>
    !>                              (**optional**, default = `132_IK`)
    !>  \param[in]  maxwidth    :   The input scalar of type `integer` of default kind \IK,
    !>                              representing the maximum length of each wrapped line, beyond which no character can appear (except `newline`).<br>
    !>                              If a word overflows the specified input `maxwidth` line limit, it will be broken and continued on the next line.<br>
    !>                              If `break = ""`, then the value of `maxwidth` becomes irrelevant because line wrapping will strictly occur at `width`.<br>
    !>                              (**optional**, default = `huge(0_IK) / 2_IK`)
    !>
    !>  \return
    !>  `strWrapped`            :   The output `allocatable` scalar `character` of the same kind as the input `str`,
    !>                              containing the input `str` whose contents is divided into lines of maximum width `maxwidth`
    !>                              by inserting the input `newline` at the appropriate positions specified via `break`.<br>
    !>
    !>  \interface{getStrWrapped}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: getStrWrapped
    !>      character(:,kind(str)), allocatable :: strWrapped
    !>
    !>      strWrapped = getStrWrapped(str, prefix = prefix, indent = indent, break = break, newline = newline, linefeed = linefeed, width = width, maxwidth = maxwidth)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `len(indent) < width <= maxwidth` must hold at all times.<br>
    !>  Otherwise, appropriate adjustments to the values of `width` and `maxWidth` will be made to infinite loops.<br>
    !>
    !>  \impure
    !>
    !>  \remark
    !>  This procedure exists primarily to facilitate implementation of other procedures in other modules of the library, e.g., [pm_err](@ref pm_err).<br>
    !>
    !>  \remark
    !>  It is sometimes desirable to create a list of wrapped lines instead of a single string of all wrapped lines.<br>
    !>  Generating a list of lines can be much more expensive than a single string. But if needed, it can be done by calling
    !>  [setSplit](@ref pm_arraySplit::setSplit) and setting the argument `sep` to the `newline` argument passed
    !>  to the procedures under this generic interface.<br>
    !>
    !>  \see
    !>  [display_type](@ref pm_io::display_type)<br>
    !>  [getCentered](@ref pm_arrayCenter::getCentered)<br>
    !>  [setCentered](@ref pm_arrayCenter::setCentered)<br>
    !>  [isFailedGetShellWidth](@ref pm_sysShell::isFailedGetShellWidth)<br>
    !>  [isFailedGetShellHeight](@ref pm_sysShell::isFailedGetShellHeight)<br>
    !>
    !>  \example{getStrWrapped}
    !>  \include{lineno} example/pm_str/getStrWrapped/main.F90
    !>  \compilef{getStrWrapped}
    !>  \output{getStrWrapped}
    !>  \include{lineno} example/pm_str/getStrWrapped/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getStrWrapped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getStrWrapped_SK5(str, prefix, indent, break, newline, linefeed, width, maxwidth) result(strWrapped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrapped_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: str
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline, linefeed
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        character(:,SKG)                            , allocatable   :: strWrapped
    end function
#endif

#if SK4_ENABLED
    impure module function getStrWrapped_SK4(str, prefix, indent, break, newline, linefeed, width, maxwidth) result(strWrapped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrapped_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: str
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline, linefeed
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        character(:,SKG)                            , allocatable   :: strWrapped
    end function
#endif

#if SK3_ENABLED
    impure module function getStrWrapped_SK3(str, prefix, indent, break, newline, linefeed, width, maxwidth) result(strWrapped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrapped_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: str
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline, linefeed
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        character(:,SKG)                            , allocatable   :: strWrapped
    end function
#endif

#if SK2_ENABLED
    impure module function getStrWrapped_SK2(str, prefix, indent, break, newline, linefeed, width, maxwidth) result(strWrapped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrapped_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: str
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline, linefeed
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        character(:,SKG)                            , allocatable   :: strWrapped
    end function
#endif

#if SK1_ENABLED
    impure module function getStrWrapped_SK1(str, prefix, indent, break, newline, linefeed, width, maxwidth) result(strWrapped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrapped_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: str
        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline, linefeed
        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
        character(:,SKG)                            , allocatable   :: strWrapped
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the number of non-overlapping repetitions of the input `pattern` in the input `str`
    !>  until the first occurrence of a mismatch between `pattern` and the relevant segment of `str`.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface are meant to be used primarily within the wrapping procedures in the parent module.<br>
    !>  The procedures under this generic interface count the number (`countPattern`) of non-overlapping repetitions of the input `pattern` at
    !>  the beginning of the input `str`. The returned result the length of the indentation at the beginning of `str`, **in units of `len(pattern)`**.<br>
    !>
    !>  \param[in]  str         :   The input scalar of type `character` of kind \SKALL whose initial indentation is to be found.<br>
    !>  \param[in]  pattern     :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter.<br>
    !>                              It represents the building block of the potential indentation (at the beginning of the input `str`)
    !>                              whose length is to be found.<br>
    !>
    !>  \return
    !>  `lenIndent`             :   The output scalar `integer` of default kind \IK, representing the length of the indentation
    !>                              identified at the beginning of the input `str`, **in units of `len(pattern)`**.<br>
    !>
    !>  \interface{getLenIndent}
    !>  \code{.F90}
    !>
    !>      use pm_str, only: getLenIndent
    !>      integer(IK) :: lenIndent
    !>
    !>      lenIndent = getLenIndent(str, pattern)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \warning
    !>  If the input `str` or `pattern` has a zero length or `len(str) < len(pattern)`, the returned result is `0`.<br>
    !>
    !>  \see
    !>  [getStrWrapped](@ref pm_str::getStrWrapped)<br>
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getBin](@ref pm_arraySearch::getBin)<br>
    !>
    !>  \example{getLenIndent}
    !>  \include{lineno} example/pm_str/getLenIndent/main.F90
    !>  \compilef{getLenIndent}
    !>  \output{getLenIndent}
    !>  \include{lineno} example/pm_str/getLenIndent/main.out.F90
    !>
    !>  \test
    !>  [test_pm_str](@ref test_pm_str)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLenIndent

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getLenIndent_SK5(str, pattern) result(lenIndent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLenIndent_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: str, pattern
        integer(IK)                                                 :: lenIndent
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getLenIndent_SK4(str, pattern) result(lenIndent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLenIndent_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: str, pattern
        integer(IK)                                                 :: lenIndent
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getLenIndent_SK3(str, pattern) result(lenIndent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLenIndent_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: str, pattern
        integer(IK)                                                 :: lenIndent
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getLenIndent_SK2(str, pattern) result(lenIndent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLenIndent_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: str, pattern
        integer(IK)                                                 :: lenIndent
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getLenIndent_SK1(str, pattern) result(lenIndent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLenIndent_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: str, pattern
        integer(IK)                                                 :: lenIndent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!   >  \brief
!   >  Wrap the input scalar string to within the specified width while
!   >  optionally respecting the initial indentation in all subsequent lines.<br>
!   >
!   >  \param[in]  str         :   The input scalar of type `character` of kind \SKALL that is to be wrapped.<br>
!   >  \param[in]  prefix      :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter,
!   >                              whose contents will be prefixed to each new line (including the first line) before an indentation occurs.<br>
!   >                              Note that `prefix` does **not** count toward filling the specified `width` or `maxwidth` of each wrapped line.<br>
!   >                              (**optional**, default = `SKG_""`)
!   >  \param[in]  indent      :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter,
!   >                              whose contents are the set of characters to search for at the beginning of the input `str` representing the indentation.<br>
!   >                              If specified as a non-empty input argument, then all wrappings of `str` will be indented with as many instances of `indent`
!   >                              that are found at the beginning of the input `str`.<br>
!   >                              For example, the following choice of input arguments,
!   >                              \code{.F90}
!   >                                  List = getStrWrappedList(str = "~+~+This is a test.", indent = "~+")
!   >                              \endcode
!   >                              will lead to prefixing each wrapped line in the output `List` with `~+~+`.<br>
!   >                              (**optional**, default = `" "`. If `indent = ""`, then no indentation of the lines will occur.)
!   >  \param[in]  break       :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter,
!   >                              whose contents are the set of characters at which breaking the input `str` is allowed.<br>
!   >                              Setting `break = ""` will lead to wrapping the input `str` at any arbitrary character
!   >                              as soon as the length of a given line reaches the specified `width`.<br>
!   >                              (**optional**, default = " ")
!   >  \param[in]  newline    :   The input scalar of the same type and kind as the input `str` of arbitrary length type parameter.<br>
!   >                              Upon encountering a pattern that matches `newline`, the pattern will be removed and a new line will be created
!   >                              starting with `prefix//indent` followed by the rest of the `str` that appears immediately after `newline`.<br>
!   >                              (**optional**, default = `new_line(SKG_" ")`)
!   >  \param[in]  width       :   The input scalar of type `integer` of default kind \IK, representing the length of each line in the
!   >                              output wrapping, **with the possibility of word overflows**.<br>
!   !>                             If a word overflows the specified input line width,
!   >                              it is allowed to continue in the same line until the next `break` character is encountered.<br>
!   >                              (**optional**, default = `132_IK`)
!   >  \param[in]  maxwidth    :   The input scalar of type `integer` of default kind \IK,
!   >                              representing the maximum length of each line in `List`, beyond which no character can appear in any line.<br>
!   >                              If a word overflows the specified input `maxwidth` line limit, it will be broken and continued on the next line.<br>
!   >                              If `break = ""`, then the value of `maxwidth` becomes irrelevant because line wrapping will strictly occur at `width`.<br>
!   >                              (**optional**, default = `huge(0_IK) / 2_IK`)
!   >
!   >  \return
!   >  `List`                  :   The output `allocatable` array rank `1` of type string-container [css_pdt](@ pm_container::css_pdt)
!   >                              whose string contents are of the same kind as the kind of the input `str`, each element of which
!   >                              contains a segment of the input `str` wrapped within the specified width.<br>
!   >
!   >  \interface{getStrWrappedList}
!   >  \code{.F90}
!   >
!   >      use pm_str, only: getStrWrappedList
!   >      use pm_container, only: css_pdt
!   >
!   >      type(css_pdt(kind(str))). allocatable :: List(:)
!   >
!   >      List = getStrWrappedList(str)
!   >
!   >  \endcode
!   >
!   >  \warning
!   >  An input `str` satisfying `len(str,IK) >= maxwidth` can potentially cause arithmetic overflow in this algorithm.<br>
!   >  When \IK refers to the 32-bit `integer` kind, this issue becomes relevant only for strings longer than \f$1^9\f$ characters.<br>
!   >
!   >  \warning
!   >  The condition `len(indent) < width <= maxwidth` must hold at all times.<br>
!   >  \vericons
!   >
!   >  \warnpure
!   >
!   >  \remark
!   >  This procedure exists primarily to facilitate implementation of other procedures in other modules of the library, e.g., [pm_err](@ref pm_err).<br>
!   >  Other algorithms in this library listed below can provide more general solutions to the problem solved by this procedure.<br>
!   >
!   >  \see
!   >  [getLocNonSpace](@ref pm_str::getLocNonSpace)<br>
!   >  [setLoc](@ref pm_arrayFind::setLoc)<br>
!   >  [getLoc](@ref pm_arrayFind::getLoc)<br>
!   >  [getBin](@ref pm_arraySearch::getBin)<br>
!   >
!   >  \example{getStrWrappedList}
!   >  \include{lineno} example/pm_str/getStrWrappedList/main.F90
!   >  \compile{getStrWrappedList}
!   >  \output{getStrWrappedList}
!   >  \include{lineno} example/pm_str/getStrWrappedList/main.out.F90
!   >
!   >  \test
!   >  [test_pm_str](@ref test_pm_str)
!   >
!   >  \final
!   >
!   >  \author
!   >  Amir Shahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!    interface getStrWrappedList
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    impure module function getStrWrappedList_SK5(str, prefix, indent, break, newline, width, maxwidth) result(List)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrappedList_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        character(*,SKG)            , intent(in)                    :: str
!        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
!        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
!        type(css_pdt(SKG))                       , allocatable   :: List(:)
!    end function
!#endif
!
!#if SK4_ENABLED
!    impure module function getStrWrappedList_SK4(str, prefix, indent, break, newline, width, maxwidth) result(List)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrappedList_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        character(*,SKG)            , intent(in)                    :: str
!        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
!        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
!        type(css_pdt(SKG))                       , allocatable   :: List(:)
!    end function
!#endif
!
!#if SK3_ENABLED
!    impure module function getStrWrappedList_SK3(str, prefix, indent, break, newline, width, maxwidth) result(List)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrappedList_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        character(*,SKG)            , intent(in)                    :: str
!        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
!        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
!        type(css_pdt(SKG))                       , allocatable   :: List(:)
!    end function
!#endif
!
!#if SK2_ENABLED
!    impure module function getStrWrappedList_SK2(str, prefix, indent, break, newline, width, maxwidth) result(List)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrappedList_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        character(*,SKG)            , intent(in)                    :: str
!        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
!        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
!        type(css_pdt(SKG))                       , allocatable   :: List(:)
!    end function
!#endif
!
!#if SK1_ENABLED
!    impure module function getStrWrappedList_SK1(str, prefix, indent, break, newline, width, maxwidth) result(List)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getStrWrappedList_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        character(*,SKG)            , intent(in)                    :: str
!        character(*,SKG)            , intent(in)    , optional      :: prefix, indent, break, newline
!        integer(IK)                 , intent(in)    , optional      :: width, maxwidth
!        type(css_pdt(SKG))                       , allocatable   :: List(:)
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!contains
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !>  Generate and return a vector of single-characters each element of which corresponds to one character in the input scalar string.<br>
!    !>
!    !>  \param[in]  string  :   The input scalar `character` of arbitrary `len` type parameter of default kind \SK.<br>
!    !>
!    !>  \return
!    !>  `charVec`           :   The output vector of the same type and kind as the input `string`, composed
!    !>                          of single characters corresponding to each character in the input scalar string.<br>
!    !>
!    !>  \interface{getCharVec}
!    !>  \code{.F90}
!    !>
!    !>      use pm_str, only: getCharVec
!    !>      character(*, SK) :: string
!    !>      character(1, SK) :: charVec(len(string))
!    !>
!    !>      charVec = getCharVec(string)
!    !>
!    !>  \endcode
!    !>
!    !>  \pure
!    !>
!    !>  \see
!    !>  [getCharSeq](@ref pm_str::getCharSeq)<br>
!    !>
!    !>  \example{getCharVec}
!    !>  \include{lineno} example/pm_str/getCharVec/main.F90
!    !>  \compilef
!    !>  \output
!    !>  \include{lineno} example/pm_str/getCharVec/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_str](@ref test_pm_str)
!    !>
!    !>  \final
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!    pure function getCharVec(str) result(charVec)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getCharVec
!#endif
!        use pm_kind, only: SK, IK
!        character(*, SK), intent(in)    :: str
!        character(1, SK)                :: charVec(len(str, kind=IK))
!        integer(IK)                     :: i
!        do concurrent(i = 1_IK : len(str, kind = IK))
!            charVec(i) = str(i:i)
!        end do
!    end function
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !>  Generate and return a scalar string resulting from the concatenation of
!    !>  the individual elements (single-characters) of the input vector of characters.<br>
!    !>
!    !>  \param[in]  charVec :   The input `contiguous` vector of rank `1` of type `character` of default kind \SK
!    !>                          of `len` type parameter `1`, each element of which contains a single character.<br>
!    !>
!    !>  \return
!    !>  `string`            :   The output string scalar of the same kind as `charVec` whose content
!    !>                          is the concatenation of the single characters in the input `charVec`.<br>
!    !>
!    !>  \interface{getCharSeq}
!    !>  \code{.F90}
!    !>
!    !>      use pm_str, only: getCharSeq
!    !>      character(1, SK) :: charVec(*)
!    !>      character(size(charVec), SK) :: string
!    !>
!    !>      string = getCharSeq(charVec)
!    !>
!    !>  \endcode
!    !>
!    !>  \pure
!    !>
!    !>  \see
!    !>  [getCharVec](@ref pm_str::getCharVec)<br>
!    !>
!    !>  \example{getCharSeq}
!    !>  \include{lineno} example/pm_str/getCharSeq/main.F90
!    !>  \compilef
!    !>  \output
!    !>  \include{lineno} example/pm_str/getCharSeq/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_str](@ref test_pm_str)
!    !>
!    !>  \final
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!    pure function getCharSeq(charVec) result(string)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getCharSeq
!#endif
!        use pm_kind, only: SK, IK
!        character(1, SK), intent(in), contiguous    :: charVec(:)
!        character(size(charVec, kind = IK), SK)     :: string
!        integer(IK)                                 :: i
!        do concurrent(i = 1_IK : size(charVec, kind = IK))
!            string(i:i) = charVec(i)
!        end do
!    end function
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !>  Pad the input assumed-size character vector `str` with the input assumed-size
!    !>  character vector `fill` for a length of `lenPadded` and return the resulting new string.<br>
!    !>
!    !>  \param[in]  str         :   The input assumed-size character vector to be padded.<br>
!    !>  \param[in]  fill        :   The input assumed-size character vector to be used as the fill for padding.<br>
!    !>  \param[in]  lenPadded   :   The length of the resulting final allocatable character vector.<br>
!    !>
!    !>  \return
!    !>  `strPadded`             :   The output string padded with `fill`.<br>
!    !>
!    !>  \remark
!    !>  Note that `fill` is an assumed-size character vector of any length. However, if the full length of fill does
!    !>  not fit at the end of the padded output `strPadded`, the `fill` will be cut at the end of the output `strPadded`.<br>
!    pure function getPadded(str, fill, lenPadded) result(strPadded)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getPadded
!#endif
!        use pm_kind, only: SK, IK
!        implicit none
!        character(*, SK), intent(in)    :: str
!        character(*, SK), intent(in)    :: fill
!        integer(IK)     , intent(in)    :: lenPadded
!        character(lenPadded)            :: strPadded
!        character(:, SK), allocatable   :: pad
!        integer(IK)                     :: lenStr, lenFill, fillCount, diff ! LCOV_EXCL_LINE
!        lenStr = len(str, kind = IK)
!        if (lenStr >= lenPadded) then
!            strPadded = str
!            return
!        end if
!        lenFill = len(fill, kind = IK)
!        diff = lenPadded - lenStr
!        fillCount = diff / lenFill + 1_IK
!        pad = repeat(fill,fillCount)
!        strPadded = str // pad(1:diff)
!    end function getPadded
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_str ! LCOV_EXCL_LINE