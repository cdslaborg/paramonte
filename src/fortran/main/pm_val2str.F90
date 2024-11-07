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
!>  This module contains the generic procedures for converting values of different types and kinds to Fortran strings.
!>
!>  \details
!>  The functionalities of the generic interfaces of this module are similar to and go beyond the Python `join()` string method.<br>
!>  To replicate the Python `join()` functionality, simply call the generic interfaces with `format` input argument.<br>
!>  For example, to create a concatenated `, ` separated string from an integer vector, one can try,<br>
!>  \code{.F90}
!>      use pm_val2str, only: getStr
!>      character(:), allocatable :: str
!>      str = getStr([1, 2, 3, 4])
!>  \endcode
!>  or to separate the items via an arbitrary separator, e.g., `_paramone_`, try,
!>  \code{.F90}
!>      use pm_val2str, only: getStr
!>      character(:), allocatable :: str
!>      str = getStr([1, 2, 3, 4], format = "(*(g0,:,'_paramone_'))")
!>  \endcode
!>
!>  \devnote
!>  Do **not** change the double back-ticks in <pre>``"(g0,:,',')"``</pre> to single back-ticks in any documentations in this module.<br>
!>  Doxygen version 1.9 has difficultly parsing the code sections with single back-tick
!>  when the code contains advanced features of modern Fortran `g0` edit descriptor.<br>
!>
!>  \see
!>  [pm_val2str](@ref pm_val2str)<br>
!>  [pm_val2int](@ref pm_val2int)<br>
!>  [pm_val2logical](@ref pm_val2logical)<br>
!>  [pm_val2complex](@ref pm_val2complex)<br>
!>  [pm_val2real](@ref pm_val2real)<br>
!>
!>  \test
!>  [test_pm_val2str](@ref test_pm_val2str)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_val2str

    use pm_kind, only: SK, IK, LK

    implicit none

    public

    character(*, SK), parameter :: MODULE_NAME = "@pm_val2str"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the conversion of the input value to an output Fortran string,
    !>  possibly with the requested format and sign symbol, if provided.
    !>
    !>  \param[in]  val     :   The input scalar or `contiguous` array of rank 0, 1, or 2 of either <br>
    !>                          <ol>
    !>                              <li>    type [css_type](@ref pm_container::css_type) (derived type container of scalar string of default kind \SK) or<br>
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized derived type container of scalar string of kind \SKALL) or<br>
    !>                              <li>    type `character` of kind \SKALL or<br>
    !>                              <li>    type `integer` of kind \IKALL or<br>
    !>                              <li>    type `logical` of kind \LKALL or<br>
    !>                              <li>    type `complex` of kind \CKALL or<br>
    !>                              <li>    type `real` of kind \RKALL.<br>
    !>                          </ol>
    !>                          containing the value to be converted to string.
    !>  \param[in]  format  :   The input scalar of type `character` of default kind \SK representing the Fortran
    !>                          IO format to be used when writing the value to the output allocatable character <br>
    !>                          (**optional**, the default value depends on the type of the input `val`:<br>
    !>                          <ol>
    !>                              <li>    If `val` is of type `complex` and the `sign` argument is missing, then <code>format = "(*('(',g0,',',g0,')',:,','))"</code>.<br>
    !>                              <li>    If `val` is of type `complex` and `sign = .true.`, then <code>format = "(*('(',sp,g0,',',g0,')',:,','))"</code>.<br>
    !>                              <li>    If `val` is of type `integer` or `real` and the `sign` argument is missing, then <code>format = "(*(g0,:,','))"</code>.<br>
    !>                              <li>    If `val` is of type `integer` or `real` and the `sign` argument is missing, then <code>format = "(*(sp,g0,:,','))"</code>.<br>
    !>                              <li>    If `val` is neither `integer` nor `real` nor `complex`, then <code>format = "(*(g0,:,','))"</code>. However,<br>
    !>                                      <ol>
    !>                                          <li>    If `val` is of type `character`, then the blank characters are trimmed form the right side of each value.<br>
    !>                                          <li>    If `val` is of type `logical`, then the values `TRUE` and `FALSE` are output, corresponding to `.true.` and `.false.`.<br>
    !>                                      </ol>
    !>                          </ol>
    !>                          In all cases, the default separator is a comma `,`.)
    !>  \param[in]  length  :   The input scalar of type `integer` of default kind \IK representing the length of the output string.<br>
    !>                          If provided, the output string will be preallocated with the requested input length. Then, the value will be written to it.<br>
    !>                          If missing, the blank characters on both ends of the output string will be trimmed.<br>
    !>                          (**optional**, default = automatically determined such that all leading and trailing blank characters are removed)
    !>  \param[in]  signed  :   The input scalar of type `logical` of default kind \LK. If `.true.`, the output, if it is a number, will be signed (`+` or `-`).<br>
    !>                          Otherwise, the output numbers will be unsigned. The value of `signed` is ignored if the optional input argument `format` is present.<br>
    !>                          (**optional**, default = `.false.`)
    !>
    !>  \return
    !>  `str`               :   The output `allocatable` scalar of type `character` of default kind \SK containing the input value as a string with the specified format.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_val2str, only: getStr
    !>
    !>      str = getStr(val, format = format, length = length, signed = signed)
    !>      str = getStr(val(:), format = format, length = length, signed = signed)
    !>      str = getStr(val(:,:), format = format, length = length, signed = signed)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \warning
    !>  If the input argument `length` is present, its value must be large enough to contain the full input value.
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setStr](@ref pm_val2str::setStr)<br>
    !>  [getInt](@ref pm_val2int::getInt)<br>
    !>  [setInt](@ref pm_val2int::setInt)<br>
    !>  [getLogical](@ref pm_val2logical::getLogical)<br>
    !>  [setLogical](@ref pm_val2logical::setLogical)<br>
    !>  [getComplex](@ref pm_val2complex::getComplex)<br>
    !>  [setComplex](@ref pm_val2complex::setComplex)<br>
    !>  [getReal](@ref pm_val2real::getReal)<br>
    !>  [setReal](@ref pm_val2real::setReal)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_val2str/getStr/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_val2str/getStr/main.out.F90
    !>
    !>  \test
    !>  [test_pm_val2str](@ref test_pm_val2str)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-11}
    !>  \desc
    !>  The \gfortran cannot handle IO of string arrays of zero-length of non-zero size in `getStr()` in assertion check blocks.<br>
    !>  The issue happens to be with deferred-length strings `character(:,SKG), allocatable :: Array(:)`.<br>
    !>  The \ifort can handle this properly.<br>
    !>  \remedy
    !>  Avoid deferred-length strings where \gfortran cannot compile.
    !>
    !>  \bug
    !>
    !>  \todo
    !>  \plow This generic interface can be expanded to arrays of higher dimensions than two.
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getStr

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStr_D0_SK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_SK5_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK5
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK4_ENABLED
    PURE module function getStr_D0_SK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_SK4_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK4
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK3_ENABLED
    PURE module function getStr_D0_SK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_SK3_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK3
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK2_ENABLED
    PURE module function getStr_D0_SK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_SK2_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK2
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK1_ENABLED
    PURE module function getStr_D0_SK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_SK1_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK1
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getStr_D0_IK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_IK5_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK5
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK4_ENABLED
    PURE module function getStr_D0_IK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_IK4_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK4
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK3_ENABLED
    PURE module function getStr_D0_IK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_IK3_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK3
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK2_ENABLED
    PURE module function getStr_D0_IK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_IK2_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK2
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK1_ENABLED
    PURE module function getStr_D0_IK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_IK1_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK1
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getStr_D0_LK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_LK5_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK5
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK4_ENABLED
    PURE module function getStr_D0_LK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_LK4_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK4
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK3_ENABLED
    PURE module function getStr_D0_LK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_LK3_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK3
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK2_ENABLED
    PURE module function getStr_D0_LK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_LK2_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK2
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK1_ENABLED
    PURE module function getStr_D0_LK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_LK1_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK1
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getStr_D0_CK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_CK5_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK5
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK4_ENABLED
    PURE module function getStr_D0_CK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_CK4_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK4
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK3_ENABLED
    PURE module function getStr_D0_CK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_CK3_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK3
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK2_ENABLED
    PURE module function getStr_D0_CK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_CK2_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK2
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK1_ENABLED
    PURE module function getStr_D0_CK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_CK1_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK1
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getStr_D0_RK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_RK5_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK5
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK4_ENABLED
    PURE module function getStr_D0_RK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_RK4_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK4
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK3_ENABLED
    PURE module function getStr_D0_RK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_RK3_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK3
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK2_ENABLED
    PURE module function getStr_D0_RK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_RK2_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK2
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK1_ENABLED
    PURE module function getStr_D0_RK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_RK1_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK1
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getStr_D0_PSSK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_PSSK5_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK5
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK4_ENABLED
    PURE module function getStr_D0_PSSK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_PSSK4_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK4
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK3_ENABLED
    PURE module function getStr_D0_PSSK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_PSSK3_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK3
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK2_ENABLED
    PURE module function getStr_D0_PSSK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_PSSK2_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK2
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK1_ENABLED
    PURE module function getStr_D0_PSSK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_PSSK1_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK1
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getStr_D0_BSSK_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D0_BSSK_SK
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK, SKG => SK1
        type(css_type)          , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStr_D1_SK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_SK5_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK4_ENABLED
    PURE module function getStr_D1_SK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_SK4_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK3_ENABLED
    PURE module function getStr_D1_SK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_SK3_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK2_ENABLED
    PURE module function getStr_D1_SK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_SK2_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK1_ENABLED
    PURE module function getStr_D1_SK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_SK1_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getStr_D1_IK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_IK5_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK4_ENABLED
    PURE module function getStr_D1_IK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_IK4_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK3_ENABLED
    PURE module function getStr_D1_IK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_IK3_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK2_ENABLED
    PURE module function getStr_D1_IK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_IK2_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK1_ENABLED
    PURE module function getStr_D1_IK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_IK1_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getStr_D1_LK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_LK5_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK4_ENABLED
    PURE module function getStr_D1_LK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_LK4_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK3_ENABLED
    PURE module function getStr_D1_LK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_LK3_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK2_ENABLED
    PURE module function getStr_D1_LK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_LK2_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK1_ENABLED
    PURE module function getStr_D1_LK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_LK1_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getStr_D1_CK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_CK5_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK4_ENABLED
    PURE module function getStr_D1_CK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_CK4_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK3_ENABLED
    PURE module function getStr_D1_CK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_CK3_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK2_ENABLED
    PURE module function getStr_D1_CK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_CK2_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK1_ENABLED
    PURE module function getStr_D1_CK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_CK1_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getStr_D1_RK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_RK5_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK4_ENABLED
    PURE module function getStr_D1_RK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_RK4_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK3_ENABLED
    PURE module function getStr_D1_RK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_RK3_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK2_ENABLED
    PURE module function getStr_D1_RK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_RK2_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK1_ENABLED
    PURE module function getStr_D1_RK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_RK1_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getStr_D1_PSSK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_PSSK5_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK4_ENABLED
    PURE module function getStr_D1_PSSK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_PSSK4_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK3_ENABLED
    PURE module function getStr_D1_PSSK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_PSSK3_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK2_ENABLED
    PURE module function getStr_D1_PSSK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_PSSK2_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK1_ENABLED
    PURE module function getStr_D1_PSSK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_PSSK1_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getStr_D1_BSSK_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D1_BSSK_SK
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getStr_D2_SK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_SK5_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK4_ENABLED
    PURE module function getStr_D2_SK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_SK4_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK3_ENABLED
    PURE module function getStr_D2_SK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_SK3_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK2_ENABLED
    PURE module function getStr_D2_SK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_SK2_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK1_ENABLED
    PURE module function getStr_D2_SK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_SK1_SK
#endif
        use pm_kind, only: SKO => SK, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getStr_D2_IK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_IK5_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK4_ENABLED
    PURE module function getStr_D2_IK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_IK4_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK3_ENABLED
    PURE module function getStr_D2_IK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_IK3_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK2_ENABLED
    PURE module function getStr_D2_IK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_IK2_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if IK1_ENABLED
    PURE module function getStr_D2_IK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_IK1_SK
#endif
        use pm_kind, only: SKO => SK, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getStr_D2_LK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_LK5_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK4_ENABLED
    PURE module function getStr_D2_LK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_LK4_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK3_ENABLED
    PURE module function getStr_D2_LK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_LK3_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK2_ENABLED
    PURE module function getStr_D2_LK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_LK2_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if LK1_ENABLED
    PURE module function getStr_D2_LK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_LK1_SK
#endif
        use pm_kind, only: SKO => SK, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getStr_D2_CK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_CK5_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK4_ENABLED
    PURE module function getStr_D2_CK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_CK4_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK3_ENABLED
    PURE module function getStr_D2_CK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_CK3_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK2_ENABLED
    PURE module function getStr_D2_CK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_CK2_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if CK1_ENABLED
    PURE module function getStr_D2_CK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_CK1_SK
#endif
        use pm_kind, only: SKO => SK, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getStr_D2_RK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_RK5_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK4_ENABLED
    PURE module function getStr_D2_RK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_RK4_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK3_ENABLED
    PURE module function getStr_D2_RK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_RK3_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK2_ENABLED
    PURE module function getStr_D2_RK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_RK2_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if RK1_ENABLED
    PURE module function getStr_D2_RK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_RK1_SK
#endif
        use pm_kind, only: SKO => SK, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getStr_D2_PSSK5_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_PSSK5_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK4_ENABLED
    PURE module function getStr_D2_PSSK4_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_PSSK4_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK3_ENABLED
    PURE module function getStr_D2_PSSK3_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_PSSK3_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK2_ENABLED
    PURE module function getStr_D2_PSSK2_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_PSSK2_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#if SK1_ENABLED
    PURE module function getStr_D2_PSSK1_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_PSSK1_SK
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getStr_D2_BSSK_SK(val, format, length, signed) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStr_D2_BSSK_SK
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        integer(IK)             , intent(in)    , optional      :: length
        logical(LK)             , intent(in)    , optional      :: signed
        character(:,SKO)        , allocatable                   :: str
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the conversion of the input value to an output Fortran string,
    !>  possibly with the requested format and sign symbol, if provided.
    !>
    !>  \return
    !>  \param[out] str     :   The output scalar of type `character` of kind \SKALL containing the input value as a string with the specified format.<br>
    !>  \param[out] length  :   The output scalar of type `integer` of default kind \IK representing the length of the record written to the output string (i.e., `str(1:length)` is the record).<br>
    !>  \param[in]  val     :   The input scalar or `contiguous` array of rank 0, 1, or 2 of either <br>
    !>                          <ol>
    !>                              <li>    type [css_type](@ref pm_container::css_type) (derived type container of scalar string of default kind \SK) or<br>
    !>                              <li>    type [css_pdt](@ref pm_container::css_pdt) (parameterized derived type container of scalar string of kind \SKALL) or<br>
    !>                              <li>    type `character` of kind \SKALL or<br>
    !>                              <li>    type `integer` of kind \IKALL or<br>
    !>                              <li>    type `logical` of kind \LKALL or<br>
    !>                              <li>    type `complex` of kind \CKALL or<br>
    !>                              <li>    type `real` of kind \RKALL.<br>
    !>                          </ol>
    !>                          containing the value to be converted to string.
    !>  \param[in]  format  :   The input scalar of type `character` of default kind \SK representing the Fortran
    !>                          IO format to be used when writing the value to the output allocatable character <br>
    !>                          (**optional**, the default value depends on the type of the input `val`: <br>
    !>                          <ol>
    !>                              <li>    If `val` is of type `complex` and the `sign` argument is missing, then <code>format = "(*('(',g0,', ',g0,')',:,', '))"</code>.<br>
    !>                              <li>    If `val` is of type `complex` and `sign = .true.`, then <code>format = "(*('(',sp,g0,', ',g0,')',:,', '))"</code>.<br>
    !>                              <li>    If `val` is of type `integer` or `real` and the `sign` argument is missing, then <code>format = "(*(g0,:,', '))"</code>.<br>
    !>                              <li>    If `val` is of type `integer` or `real` and the `sign` argument is missing, then <code>format = "(*(sp,g0,:,', '))"</code>.<br>
    !>                              <li>    If `val` is neither `integer` nor `real` nor `complex`, then <code>format = "(*(g0,:,','))"</code>. However,<br>
    !>                                      <ol>
    !>                                          <li>    If `val` is of type `character`, then the blank characters are trimmed form the right side of each value.<br>
    !>                                          <li>    If `val` is of type `logical`, then the values `TRUE` and `FALSE` are output, corresponding to `.true.` and `.false.`.<br>
    !>                                      </ol>
    !>                          </ol>
    !>                          In all cases, the default separator is a comma followed by white space character `, `.)
    !>  \param[in]  signed  :   The input scalar of type `logical` of default kind \LK. If `.true.`, the output, if it is a number, will be signed (`+` or `-`).<br>
    !>                          Otherwise, the output numbers will be unsigned. The value of `signed` is ignored if the optional input argument `format` is present.<br>
    !>                          (**optional**, default = `.false.`)
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_val2str, only: setStr
    !>
    !>      call setStr(str, val, format = format, length = length, signed = signed)
    !>      call setStr(str, val(:), format = format, length = length, signed = signed)
    !>      call setStr(str, val(:,:), format = format, length = length, signed = signed)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \warning
    !>  The length of the input argument `str` must be large enough to contain the full input value.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [setStr](@ref pm_val2str::setStr)<br>
    !>  [getInt](@ref pm_val2int::getInt)<br>
    !>  [setInt](@ref pm_val2int::setInt)<br>
    !>  [getLogical](@ref pm_val2logical::getLogical)<br>
    !>  [setLogical](@ref pm_val2logical::setLogical)<br>
    !>  [getComplex](@ref pm_val2complex::getComplex)<br>
    !>  [setComplex](@ref pm_val2complex::setComplex)<br>
    !>  [getReal](@ref pm_val2real::getReal)<br>
    !>  [setReal](@ref pm_val2real::setReal)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_val2str/setStr/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_val2str/setStr/main.out.F90
    !>
    !>  \test
    !>  [test_pm_val2str](@ref test_pm_val2str)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-11}
    !>  \desc
    !>  The \gfortran cannot handle IO of string arrays of zero-length of non-zero size in `getStr()` in assertion check blocks.<br>
    !>  The issue happens to be with deferred-length strings `character(:,SKG), allocatable :: Array(:)`.<br>
    !>  The \ifort can handle this properly.<br>
    !>  \remedy
    !>  Avoid deferred-length strings where \gfortran cannot compile.
    !>
    !>  \todo
    !>  \plow This generic interface can be expanded to arrays of higher dimensions than two.
    !>
    !>  \todo
    !>  \pvhigh
    !>  The default length of fields must be generically defined.<br>
    !>  Currently, the maximum length of the fields is set to `127`.<br>
    !>  While this number is more than enough precisions of much higher orders available by compilers,
    !>  it is bound to fail in the distant future.<br>
    !>  Therefore, it must be generically set based on the `precision` of the fields.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! SKO => SK5

#if SK5_ENABLED
    interface setStr

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D0_SK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK5_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK5
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_SK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK4_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK4
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_SK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK3_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK3
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_SK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK2_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK2
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_SK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK1_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK1
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D0_IK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK5_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK5
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D0_IK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK4_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK4
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D0_IK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK3_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK3
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D0_IK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK2_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK2
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D0_IK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK1_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK1
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D0_LK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK5_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK5
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D0_LK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK4_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK4
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D0_LK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK3_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK3
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D0_LK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK2_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK2
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D0_LK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK1_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK1
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D0_CK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK5_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK5
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D0_CK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK4_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK4
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D0_CK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK3_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK3
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D0_CK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK2_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK2
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D0_CK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK1_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK1
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D0_RK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK5_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK5
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D0_RK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK4_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK4
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D0_RK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK3_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK3
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D0_RK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK2_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK2
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D0_RK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK1_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK1
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D0_PSSK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK5_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK5
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_PSSK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK4_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK4
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_PSSK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK3_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK3
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_PSSK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK2_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK2
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_PSSK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK1_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK1
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D0_BSSK_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_BSSK_SK5
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK5, SKG => SK1
        type(css_type)          , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D1_SK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK5_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_SK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK4_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_SK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK3_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_SK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK2_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_SK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK1_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D1_IK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK5_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D1_IK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK4_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D1_IK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK3_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D1_IK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK2_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D1_IK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK1_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D1_LK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK5_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D1_LK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK4_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D1_LK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK3_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D1_LK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK2_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D1_LK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK1_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D1_CK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK5_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D1_CK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK4_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D1_CK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK3_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D1_CK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK2_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D1_CK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK1_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D1_RK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK5_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D1_RK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK4_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D1_RK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK3_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D1_RK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK2_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D1_RK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK1_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D1_PSSK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK5_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_PSSK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK4_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_PSSK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK3_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_PSSK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK2_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_PSSK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK1_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D1_BSSK_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_BSSK_SK5
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK5, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D2_SK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK5_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_SK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK4_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_SK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK3_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_SK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK2_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_SK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK1_SK5
#endif
        use pm_kind, only: SKO => SK5, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D2_IK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK5_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D2_IK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK4_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D2_IK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK3_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D2_IK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK2_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D2_IK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK1_SK5
#endif
        use pm_kind, only: SKO => SK5, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D2_LK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK5_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D2_LK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK4_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D2_LK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK3_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D2_LK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK2_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D2_LK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK1_SK5
#endif
        use pm_kind, only: SKO => SK5, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D2_CK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK5_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D2_CK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK4_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D2_CK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK3_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D2_CK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK2_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D2_CK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK1_SK5
#endif
        use pm_kind, only: SKO => SK5, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D2_RK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK5_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D2_RK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK4_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D2_RK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK3_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D2_RK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK2_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D2_RK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK1_SK5
#endif
        use pm_kind, only: SKO => SK5, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D2_PSSK5_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK5_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_PSSK4_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK4_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_PSSK3_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK3_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_PSSK2_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK2_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_PSSK1_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK1_SK5
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK5, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D2_BSSK_SK5(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_BSSK_SK5
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK5, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
#endif

    ! SKO => SK4

#if SK4_ENABLED
    interface setStr

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D0_SK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK5_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK5
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_SK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK4_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK4
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_SK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK3_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK3
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_SK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK2_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK2
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_SK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK1_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK1
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D0_IK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK5_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK5
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D0_IK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK4_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK4
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D0_IK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK3_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK3
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D0_IK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK2_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK2
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D0_IK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK1_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK1
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D0_LK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK5_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK5
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D0_LK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK4_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK4
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D0_LK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK3_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK3
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D0_LK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK2_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK2
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D0_LK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK1_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK1
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D0_CK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK5_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK5
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D0_CK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK4_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK4
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D0_CK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK3_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK3
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D0_CK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK2_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK2
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D0_CK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK1_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK1
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D0_RK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK5_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK5
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D0_RK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK4_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK4
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D0_RK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK3_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK3
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D0_RK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK2_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK2
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D0_RK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK1_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK1
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D0_PSSK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK5_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK5
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_PSSK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK4_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK4
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_PSSK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK3_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK3
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_PSSK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK2_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK2
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_PSSK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK1_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK1
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D0_BSSK_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_BSSK_SK4
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK4, SKG => SK1
        type(css_type)          , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D1_SK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK5_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_SK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK4_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_SK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK3_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_SK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK2_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_SK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK1_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D1_IK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK5_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D1_IK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK4_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D1_IK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK3_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D1_IK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK2_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D1_IK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK1_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D1_LK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK5_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D1_LK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK4_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D1_LK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK3_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D1_LK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK2_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D1_LK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK1_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D1_CK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK5_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D1_CK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK4_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D1_CK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK3_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D1_CK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK2_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D1_CK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK1_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D1_RK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK5_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D1_RK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK4_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D1_RK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK3_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D1_RK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK2_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D1_RK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK1_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D1_PSSK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK5_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_PSSK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK4_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_PSSK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK3_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_PSSK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK2_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_PSSK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK1_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D1_BSSK_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_BSSK_SK4
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK4, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D2_SK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK5_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_SK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK4_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_SK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK3_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_SK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK2_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_SK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK1_SK4
#endif
        use pm_kind, only: SKO => SK4, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D2_IK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK5_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D2_IK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK4_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D2_IK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK3_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D2_IK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK2_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D2_IK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK1_SK4
#endif
        use pm_kind, only: SKO => SK4, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D2_LK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK5_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D2_LK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK4_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D2_LK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK3_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D2_LK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK2_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D2_LK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK1_SK4
#endif
        use pm_kind, only: SKO => SK4, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D2_CK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK5_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D2_CK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK4_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D2_CK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK3_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D2_CK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK2_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D2_CK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK1_SK4
#endif
        use pm_kind, only: SKO => SK4, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D2_RK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK5_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D2_RK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK4_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D2_RK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK3_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D2_RK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK2_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D2_RK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK1_SK4
#endif
        use pm_kind, only: SKO => SK4, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D2_PSSK5_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK5_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_PSSK4_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK4_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_PSSK3_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK3_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_PSSK2_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK2_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_PSSK1_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK1_SK4
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK4, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D2_BSSK_SK4(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_BSSK_SK4
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK4, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
#endif

    ! SKO => SK3

#if SK3_ENABLED
    interface setStr

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D0_SK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK5_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK5
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_SK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK4_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK4
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_SK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK3_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK3
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_SK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK2_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK2
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_SK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK1_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK1
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D0_IK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK5_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK5
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D0_IK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK4_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK4
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D0_IK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK3_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK3
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D0_IK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK2_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK2
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D0_IK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK1_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK1
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D0_LK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK5_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK5
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D0_LK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK4_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK4
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D0_LK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK3_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK3
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D0_LK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK2_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK2
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D0_LK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK1_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK1
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D0_CK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK5_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK5
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D0_CK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK4_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK4
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D0_CK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK3_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK3
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D0_CK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK2_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK2
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D0_CK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK1_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK1
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D0_RK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK5_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK5
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D0_RK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK4_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK4
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D0_RK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK3_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK3
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D0_RK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK2_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK2
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D0_RK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK1_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK1
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D0_PSSK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK5_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK5
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_PSSK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK4_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK4
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_PSSK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK3_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK3
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_PSSK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK2_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK2
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_PSSK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK1_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK1
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D0_BSSK_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_BSSK_SK3
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK3, SKG => SK1
        type(css_type)          , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D1_SK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK5_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_SK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK4_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_SK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK3_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_SK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK2_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_SK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK1_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D1_IK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK5_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D1_IK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK4_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D1_IK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK3_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D1_IK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK2_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D1_IK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK1_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D1_LK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK5_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D1_LK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK4_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D1_LK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK3_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D1_LK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK2_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D1_LK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK1_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D1_CK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK5_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D1_CK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK4_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D1_CK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK3_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D1_CK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK2_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D1_CK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK1_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D1_RK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK5_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D1_RK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK4_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D1_RK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK3_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D1_RK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK2_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D1_RK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK1_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D1_PSSK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK5_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_PSSK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK4_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_PSSK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK3_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_PSSK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK2_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_PSSK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK1_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D1_BSSK_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_BSSK_SK3
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK3, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D2_SK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK5_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_SK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK4_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_SK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK3_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_SK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK2_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_SK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK1_SK3
#endif
        use pm_kind, only: SKO => SK3, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D2_IK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK5_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D2_IK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK4_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D2_IK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK3_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D2_IK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK2_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D2_IK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK1_SK3
#endif
        use pm_kind, only: SKO => SK3, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D2_LK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK5_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D2_LK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK4_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D2_LK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK3_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D2_LK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK2_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D2_LK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK1_SK3
#endif
        use pm_kind, only: SKO => SK3, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D2_CK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK5_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D2_CK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK4_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D2_CK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK3_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D2_CK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK2_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D2_CK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK1_SK3
#endif
        use pm_kind, only: SKO => SK3, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D2_RK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK5_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D2_RK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK4_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D2_RK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK3_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D2_RK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK2_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D2_RK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK1_SK3
#endif
        use pm_kind, only: SKO => SK3, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D2_PSSK5_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK5_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_PSSK4_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK4_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_PSSK3_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK3_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_PSSK2_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK2_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_PSSK1_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK1_SK3
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK3, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D2_BSSK_SK3(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_BSSK_SK3
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK3, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
#endif

    ! SKO => SK2

#if SK2_ENABLED
    interface setStr

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D0_SK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK5_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK5
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_SK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK4_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK4
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_SK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK3_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK3
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_SK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK2_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK2
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_SK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK1_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK1
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D0_IK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK5_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK5
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D0_IK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK4_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK4
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D0_IK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK3_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK3
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D0_IK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK2_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK2
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D0_IK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK1_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK1
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D0_LK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK5_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK5
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D0_LK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK4_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK4
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D0_LK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK3_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK3
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D0_LK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK2_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK2
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D0_LK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK1_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK1
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D0_CK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK5_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK5
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D0_CK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK4_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK4
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D0_CK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK3_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK3
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D0_CK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK2_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK2
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D0_CK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK1_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK1
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D0_RK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK5_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK5
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D0_RK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK4_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK4
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D0_RK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK3_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK3
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D0_RK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK2_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK2
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D0_RK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK1_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK1
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D0_PSSK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK5_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK5
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_PSSK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK4_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK4
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_PSSK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK3_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK3
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_PSSK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK2_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK2
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_PSSK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK1_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK1
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D0_BSSK_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_BSSK_SK2
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK2, SKG => SK1
        type(css_type)          , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D1_SK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK5_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_SK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK4_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_SK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK3_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_SK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK2_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_SK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK1_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D1_IK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK5_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D1_IK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK4_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D1_IK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK3_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D1_IK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK2_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D1_IK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK1_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D1_LK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK5_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D1_LK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK4_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D1_LK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK3_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D1_LK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK2_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D1_LK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK1_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D1_CK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK5_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D1_CK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK4_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D1_CK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK3_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D1_CK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK2_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D1_CK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK1_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D1_RK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK5_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D1_RK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK4_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D1_RK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK3_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D1_RK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK2_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D1_RK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK1_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D1_PSSK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK5_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_PSSK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK4_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_PSSK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK3_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_PSSK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK2_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_PSSK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK1_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D1_BSSK_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_BSSK_SK2
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK2, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D2_SK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK5_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_SK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK4_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_SK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK3_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_SK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK2_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_SK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK1_SK2
#endif
        use pm_kind, only: SKO => SK2, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D2_IK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK5_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D2_IK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK4_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D2_IK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK3_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D2_IK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK2_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D2_IK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK1_SK2
#endif
        use pm_kind, only: SKO => SK2, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D2_LK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK5_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D2_LK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK4_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D2_LK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK3_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D2_LK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK2_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D2_LK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK1_SK2
#endif
        use pm_kind, only: SKO => SK2, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D2_CK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK5_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D2_CK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK4_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D2_CK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK3_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D2_CK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK2_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D2_CK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK1_SK2
#endif
        use pm_kind, only: SKO => SK2, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D2_RK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK5_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D2_RK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK4_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D2_RK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK3_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D2_RK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK2_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D2_RK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK1_SK2
#endif
        use pm_kind, only: SKO => SK2, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D2_PSSK5_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK5_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_PSSK4_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK4_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_PSSK3_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK3_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_PSSK2_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK2_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_PSSK1_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK1_SK2
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK2, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D2_BSSK_SK2(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_BSSK_SK2
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK2, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
#endif

    ! SKO => SK1

#if SK1_ENABLED
    interface setStr

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D0_SK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK5_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK5
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_SK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK4_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK4
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_SK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK3_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK3
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_SK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK2_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK2
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_SK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_SK1_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK1
        character(*,SKG)        , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D0_IK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK5_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK5
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D0_IK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK4_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK4
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D0_IK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK3_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK3
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D0_IK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK2_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK2
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D0_IK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_IK1_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK1
        integer(IKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D0_LK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK5_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK5
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D0_LK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK4_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK4
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D0_LK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK3_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK3
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D0_LK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK2_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK2
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D0_LK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_LK1_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK1
        logical(LKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D0_CK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK5_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK5
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D0_CK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK4_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK4
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D0_CK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK3_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK3
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D0_CK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK2_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK2
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D0_CK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_CK1_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK1
        complex(CKG)            , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D0_RK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK5_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK5
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D0_RK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK4_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK4
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D0_RK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK3_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK3
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D0_RK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK2_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK2
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D0_RK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_RK1_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK1
        real(RKG)               , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D0_PSSK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK5_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK5
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D0_PSSK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK4_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK4
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D0_PSSK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK3_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK3
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D0_PSSK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK2_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK2
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D0_PSSK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_PSSK1_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK1
        type(css_pdt(SKG))      , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D0_BSSK_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D0_BSSK_SK1
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK1, SKG => SK1
        type(css_type)          , intent(in)                    :: val
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D1_SK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK5_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_SK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK4_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_SK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK3_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_SK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK2_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_SK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_SK1_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D1_IK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK5_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D1_IK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK4_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D1_IK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK3_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D1_IK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK2_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D1_IK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_IK1_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D1_LK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK5_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D1_LK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK4_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D1_LK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK3_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D1_LK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK2_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D1_LK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_LK1_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D1_CK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK5_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D1_CK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK4_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D1_CK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK3_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D1_CK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK2_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D1_CK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_CK1_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D1_RK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK5_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D1_RK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK4_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D1_RK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK3_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D1_RK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK2_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D1_RK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_RK1_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D1_PSSK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK5_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D1_PSSK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK4_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D1_PSSK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK3_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D1_PSSK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK2_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D1_PSSK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_PSSK1_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D1_BSSK_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D1_BSSK_SK1
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK1, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setStr_D2_SK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK5_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK5
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_SK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK4_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK4
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_SK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK3_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK3
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_SK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK2_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK2
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_SK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_SK1_SK1
#endif
        use pm_kind, only: SKO => SK1, SKG => SK1
        character(*,SKG), target, intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setStr_D2_IK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK5_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setStr_D2_IK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK4_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setStr_D2_IK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK3_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setStr_D2_IK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK2_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setStr_D2_IK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_IK1_SK1
#endif
        use pm_kind, only: SKO => SK1, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setStr_D2_LK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK5_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setStr_D2_LK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK4_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setStr_D2_LK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK3_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setStr_D2_LK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK2_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setStr_D2_LK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_LK1_SK1
#endif
        use pm_kind, only: SKO => SK1, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setStr_D2_CK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK5_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setStr_D2_CK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK4_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setStr_D2_CK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK3_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setStr_D2_CK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK2_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setStr_D2_CK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_CK1_SK1
#endif
        use pm_kind, only: SKO => SK1, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setStr_D2_RK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK5_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setStr_D2_RK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK4_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setStr_D2_RK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK3_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setStr_D2_RK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK2_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setStr_D2_RK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_RK1_SK1
#endif
        use pm_kind, only: SKO => SK1, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setStr_D2_PSSK5_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK5_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK5
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setStr_D2_PSSK4_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK4_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK4
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setStr_D2_PSSK3_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK3_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK3
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setStr_D2_PSSK2_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK2_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK2
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setStr_D2_PSSK1_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_PSSK1_SK1
#endif
        use pm_container, only: css_pdt
        use pm_kind, only: SKO => SK1, SKG => SK1
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setStr_D2_BSSK_SK1(str, length, val, format, signed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStr_D2_BSSK_SK1
#endif
        use pm_container, only: css_type
        use pm_kind, only: SKO => SK1, SKG => SK1
        type(css_type)          , intent(in)    , contiguous    :: val(:,:)
        character(*, SK)        , intent(in)    , optional      :: format
        logical(LK)             , intent(in)    , optional      :: signed
        integer(IK)             , intent(out)                   :: length
        character(*,SKO)        , intent(out)                   :: str
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_val2str ! LCOV_EXCL_LINE