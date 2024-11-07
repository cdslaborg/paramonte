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
!>  This module contains procedures and types for facilitating the conversion of values
!>  of different types (e.g., intrinsic Fortran string and logical) to complex values of different kinds.
!>
!>  \devnote
!>  Do **not** change the double back-ticks in <pre>``"(g0,:,',')"``</pre> to single back-ticks in any documentations in this module.<br>
!>  Doxygen version 1.9 has difficultly parsing the code sections with single back-tick when the code contains
!>  advanced features of modern Fortran `g0` edit descriptor.<br>
!>
!>  \see
!>  [pm_val2str](@ref pm_val2str)<br>
!>  [pm_val2int](@ref pm_val2int)<br>
!>  [pm_val2logical](@ref pm_val2logical)<br>
!>  [pm_val2complex](@ref pm_val2complex)<br>
!>  [pm_val2real](@ref pm_val2real)<br>
!>
!>  \test
!>  [test_pm_val2complex](@ref test_pm_val2complex)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_val2complex

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_val2complex"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the conversion of the input value to a `complex` value of default kind \CK.
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of either,<br>
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL, or<br>
    !>                              <li>    type `logical` of kind \LKALL, or<br>
    !>                          </ol>
    !>                          whose value will be converted to an output of type `complex` of default kind \CK.
    !>
    !>  \return
    !>  `conversion`        :   The scalar or array of the same shape as the input `val` of the type `complex`
    !>                          of default kind \CK, containing the conversion of the input value `val`.<br>
    !>                          <ol>
    !>                              <li>    If the input `val` is a `logical` that evaluates to `.true.`, then `conversion = (1., 1.)`, otherwise `conversion = (0., 0.)`.
    !>                              <li>    If the input `val` is a `character`, then `conversion` is the result of the List-directed I/O action `read(val,*) conversion`.
    !>                          </ol>
    !>
    !>  \interface{getComplex}
    !>  \code{.F90}
    !>
    !>      use pm_val2complex, only: getComplex
    !>
    !>      conversion = getComplex(val)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [getInt](@ref pm_val2int::getInt)<br>
    !>  [setInt](@ref pm_val2int::setInt)<br>
    !>  [getReal](@ref pm_val2real::getReal)<br>
    !>  [setReal](@ref pm_val2real::setReal)<br>
    !>  [getComplex](@ref pm_val2complex::getComplex)<br>
    !>  [setComplex](@ref pm_val2complex::setComplex)<br>
    !>  [getLogical](@ref pm_val2logical::getLogical)<br>
    !>  [setLogical](@ref pm_val2logical::setLogical)<br>
    !>
    !>  \example{getComplex}
    !>  \include{lineno} example/pm_val2complex/getComplex/main.F90
    !>  \compilef{getComplex}
    !>  \output{getComplex}
    !>  \include{lineno} example/pm_val2complex/getComplex/main.out.F90
    !>
    !>  \test
    !>  [test_pm_val2complex](@ref test_pm_val2complex)
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface can be extended to support conversion to `complex` of arbitrary kind
    !>  via an optional input `mold` argument whose type and kind dictates that of the output.<br>
    !>
    !>  \final{getComplex}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getComplex

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function getComplexDef_LK5(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_LK5
#endif
        use pm_kind, only: CKG => CK, LKG => LK5
        logical(LKG)                , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

#if LK4_ENABLED
    pure elemental module function getComplexDef_LK4(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_LK4
#endif
        use pm_kind, only: CKG => CK, LKG => LK4
        logical(LKG)                , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

#if LK3_ENABLED
    pure elemental module function getComplexDef_LK3(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_LK3
#endif
        use pm_kind, only: CKG => CK, LKG => LK3
        logical(LKG)                , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

#if LK2_ENABLED
    pure elemental module function getComplexDef_LK2(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_LK2
#endif
        use pm_kind, only: CKG => CK, LKG => LK2
        logical(LKG)                , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

#if LK1_ENABLED
    pure elemental module function getComplexDef_LK1(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_LK1
#endif
        use pm_kind, only: CKG => CK, LKG => LK1
        logical(LKG)                , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getComplexDef_SK5(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_SK5
#endif
        use pm_kind, only: CKG => CK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getComplexDef_SK4(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_SK4
#endif
        use pm_kind, only: CKG => CK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getComplexDef_SK3(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_SK3
#endif
        use pm_kind, only: CKG => CK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getComplexDef_SK2(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_SK2
#endif
        use pm_kind, only: CKG => CK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getComplexDef_SK1(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplexDef_SK1
#endif
        use pm_kind, only: CKG => CK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: val
        complex(CKG)                                                :: conversion
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the conversion of the input value to a `complex` value of arbitrary kind.
    !>
    !>  \param[out] conversion  :   The output scalar or array of the same shape as other array-like arguments,
    !>                              of type `complex` of kind \CKALL, containing the conversion of the input value `val`.<br>
    !>                              <ol>
    !>                                  <li>    If the input `val` is a `logical` that evaluates to `.true.`, then `conversion = (1., 1.)`, otherwise `conversion = (0., 0.)`.
    !>                                  <li>    If the input `val` is a `character`, then `conversion` is the result of the List-directed I/O action `read(val,*) conversion`.
    !>                              </ol>
    !>  \param[in]      val     :   The input scalar or array of the same shape as other array-like arguments, of either<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL, or<br>
    !>                                  <li>    type `logical` of kind \LKALL,<br>
    !>                              </ol>
    !>                              whose value will be converted to an output of type `complex`.
    !>  \param[out]     iostat  :   The output `integer` of default kind \IK, of the same rank as the input `val` containing
    !>                              the Fortran IO status error code returned by the `read()` statement of the Fortran standard.<br>
    !>                              On return, it is zero <b>if and only if</b> no error occurs during the execution of the Fortran `read()` statement.<br>
    !>                              Refer to the Fortran standard `read/write` statements for the meaning of different non-zero output values for `iostat`.<br>
    !>                              (**optional**, it can be present only if `val` is of type `character`. if missing, the program will halt upon an IO error.)
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_val2complex, only: setComplex
    !>
    !>      call setComplex(conversion, val)
    !>      call setComplex(conversion, val, iostat)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  When the input `val` is an array, the output `iostat`, if present, must also be a vector of the same shape and size as `val`.
    !>
    !>  \pure
    !>
    !>  \remark
    !>  The purity of the procedures under this generic interface break when the output argument `iostat` is present.
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getStr](@ref pm_val2str::getStr)<br>
    !>  [getInt](@ref pm_val2int::getInt)<br>
    !>  [setInt](@ref pm_val2int::setInt)<br>
    !>  [getReal](@ref pm_val2real::getReal)<br>
    !>  [setReal](@ref pm_val2real::setReal)<br>
    !>  [getComplex](@ref pm_val2complex::getComplex)<br>
    !>  [setComplex](@ref pm_val2complex::setComplex)<br>
    !>  [getLogical](@ref pm_val2logical::getLogical)<br>
    !>  [setLogical](@ref pm_val2logical::setLogical)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_val2complex/setComplex/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_val2complex/setComplex/main.out.F90
    !>
    !>  \test
    !>  [test_pm_val2complex](@ref test_pm_val2complex)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setComplex

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED && LK5_ENABLED
    pure elemental module subroutine setComplexDef_CK5_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_LK5
#endif
        use pm_kind, only: CKG => CK5, LKG => LK5
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK5_ENABLED && LK4_ENABLED
    pure elemental module subroutine setComplexDef_CK5_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_LK4
#endif
        use pm_kind, only: CKG => CK5, LKG => LK4
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK5_ENABLED && LK3_ENABLED
    pure elemental module subroutine setComplexDef_CK5_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_LK3
#endif
        use pm_kind, only: CKG => CK5, LKG => LK3
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK5_ENABLED && LK2_ENABLED
    pure elemental module subroutine setComplexDef_CK5_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_LK2
#endif
        use pm_kind, only: CKG => CK5, LKG => LK2
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK5_ENABLED && LK1_ENABLED
    pure elemental module subroutine setComplexDef_CK5_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_LK1
#endif
        use pm_kind, only: CKG => CK5, LKG => LK1
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexDef_CK5_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_SK5
#endif
        use pm_kind, only: CKG => CK5, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK5_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexDef_CK5_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_SK4
#endif
        use pm_kind, only: CKG => CK5, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK5_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexDef_CK5_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_SK3
#endif
        use pm_kind, only: CKG => CK5, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK5_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexDef_CK5_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_SK2
#endif
        use pm_kind, only: CKG => CK5, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK5_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexDef_CK5_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK5_SK1
#endif
        use pm_kind, only: CKG => CK5, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK4_ENABLED && LK5_ENABLED
    pure elemental module subroutine setComplexDef_CK4_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_LK5
#endif
        use pm_kind, only: CKG => CK4, LKG => LK5
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK4_ENABLED && LK4_ENABLED
    pure elemental module subroutine setComplexDef_CK4_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_LK4
#endif
        use pm_kind, only: CKG => CK4, LKG => LK4
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK4_ENABLED && LK3_ENABLED
    pure elemental module subroutine setComplexDef_CK4_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_LK3
#endif
        use pm_kind, only: CKG => CK4, LKG => LK3
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK4_ENABLED && LK2_ENABLED
    pure elemental module subroutine setComplexDef_CK4_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_LK2
#endif
        use pm_kind, only: CKG => CK4, LKG => LK2
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK4_ENABLED && LK1_ENABLED
    pure elemental module subroutine setComplexDef_CK4_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_LK1
#endif
        use pm_kind, only: CKG => CK4, LKG => LK1
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK4_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexDef_CK4_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_SK5
#endif
        use pm_kind, only: CKG => CK4, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK4_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexDef_CK4_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_SK4
#endif
        use pm_kind, only: CKG => CK4, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK4_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexDef_CK4_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_SK3
#endif
        use pm_kind, only: CKG => CK4, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK4_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexDef_CK4_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_SK2
#endif
        use pm_kind, only: CKG => CK4, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK4_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexDef_CK4_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK4_SK1
#endif
        use pm_kind, only: CKG => CK4, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED && LK5_ENABLED
    pure elemental module subroutine setComplexDef_CK3_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_LK5
#endif
        use pm_kind, only: CKG => CK3, LKG => LK5
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK3_ENABLED && LK4_ENABLED
    pure elemental module subroutine setComplexDef_CK3_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_LK4
#endif
        use pm_kind, only: CKG => CK3, LKG => LK4
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK3_ENABLED && LK3_ENABLED
    pure elemental module subroutine setComplexDef_CK3_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_LK3
#endif
        use pm_kind, only: CKG => CK3, LKG => LK3
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK3_ENABLED && LK2_ENABLED
    pure elemental module subroutine setComplexDef_CK3_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_LK2
#endif
        use pm_kind, only: CKG => CK3, LKG => LK2
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK3_ENABLED && LK1_ENABLED
    pure elemental module subroutine setComplexDef_CK3_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_LK1
#endif
        use pm_kind, only: CKG => CK3, LKG => LK1
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexDef_CK3_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_SK5
#endif
        use pm_kind, only: CKG => CK3, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK3_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexDef_CK3_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_SK4
#endif
        use pm_kind, only: CKG => CK3, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK3_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexDef_CK3_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_SK3
#endif
        use pm_kind, only: CKG => CK3, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK3_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexDef_CK3_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_SK2
#endif
        use pm_kind, only: CKG => CK3, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK3_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexDef_CK3_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK3_SK1
#endif
        use pm_kind, only: CKG => CK3, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK2_ENABLED && LK5_ENABLED
    pure elemental module subroutine setComplexDef_CK2_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_LK5
#endif
        use pm_kind, only: CKG => CK2, LKG => LK5
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK2_ENABLED && LK4_ENABLED
    pure elemental module subroutine setComplexDef_CK2_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_LK4
#endif
        use pm_kind, only: CKG => CK2, LKG => LK4
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK2_ENABLED && LK3_ENABLED
    pure elemental module subroutine setComplexDef_CK2_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_LK3
#endif
        use pm_kind, only: CKG => CK2, LKG => LK3
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK2_ENABLED && LK2_ENABLED
    pure elemental module subroutine setComplexDef_CK2_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_LK2
#endif
        use pm_kind, only: CKG => CK2, LKG => LK2
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK2_ENABLED && LK1_ENABLED
    pure elemental module subroutine setComplexDef_CK2_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_LK1
#endif
        use pm_kind, only: CKG => CK2, LKG => LK1
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK2_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexDef_CK2_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_SK5
#endif
        use pm_kind, only: CKG => CK2, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK2_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexDef_CK2_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_SK4
#endif
        use pm_kind, only: CKG => CK2, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK2_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexDef_CK2_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_SK3
#endif
        use pm_kind, only: CKG => CK2, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK2_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexDef_CK2_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_SK2
#endif
        use pm_kind, only: CKG => CK2, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK2_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexDef_CK2_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK2_SK1
#endif
        use pm_kind, only: CKG => CK2, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK1_ENABLED && LK5_ENABLED
    pure elemental module subroutine setComplexDef_CK1_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_LK5
#endif
        use pm_kind, only: CKG => CK1, LKG => LK5
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK1_ENABLED && LK4_ENABLED
    pure elemental module subroutine setComplexDef_CK1_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_LK4
#endif
        use pm_kind, only: CKG => CK1, LKG => LK4
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK1_ENABLED && LK3_ENABLED
    pure elemental module subroutine setComplexDef_CK1_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_LK3
#endif
        use pm_kind, only: CKG => CK1, LKG => LK3
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK1_ENABLED && LK2_ENABLED
    pure elemental module subroutine setComplexDef_CK1_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_LK2
#endif
        use pm_kind, only: CKG => CK1, LKG => LK2
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if CK1_ENABLED && LK1_ENABLED
    pure elemental module subroutine setComplexDef_CK1_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_LK1
#endif
        use pm_kind, only: CKG => CK1, LKG => LK1
        complex(CKG)                , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK1_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexDef_CK1_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_SK5
#endif
        use pm_kind, only: CKG => CK1, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK1_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexDef_CK1_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_SK4
#endif
        use pm_kind, only: CKG => CK1, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK1_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexDef_CK1_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_SK3
#endif
        use pm_kind, only: CKG => CK1, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK1_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexDef_CK1_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_SK2
#endif
        use pm_kind, only: CKG => CK1, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if CK1_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexDef_CK1_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexDef_CK1_SK1
#endif
        use pm_kind, only: CKG => CK1, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexIO_CK5_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK5_SK5
#endif
        use pm_kind, only: CKG => CK5, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK5_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexIO_CK5_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK5_SK4
#endif
        use pm_kind, only: CKG => CK5, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK5_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexIO_CK5_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK5_SK3
#endif
        use pm_kind, only: CKG => CK5, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK5_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexIO_CK5_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK5_SK2
#endif
        use pm_kind, only: CKG => CK5, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK5_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexIO_CK5_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK5_SK1
#endif
        use pm_kind, only: CKG => CK5, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK4_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexIO_CK4_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK4_SK5
#endif
        use pm_kind, only: CKG => CK4, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK4_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexIO_CK4_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK4_SK4
#endif
        use pm_kind, only: CKG => CK4, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK4_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexIO_CK4_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK4_SK3
#endif
        use pm_kind, only: CKG => CK4, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK4_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexIO_CK4_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK4_SK2
#endif
        use pm_kind, only: CKG => CK4, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK4_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexIO_CK4_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK4_SK1
#endif
        use pm_kind, only: CKG => CK4, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK3_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexIO_CK3_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK3_SK5
#endif
        use pm_kind, only: CKG => CK3, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK3_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexIO_CK3_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK3_SK4
#endif
        use pm_kind, only: CKG => CK3, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK3_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexIO_CK3_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK3_SK3
#endif
        use pm_kind, only: CKG => CK3, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK3_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexIO_CK3_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK3_SK2
#endif
        use pm_kind, only: CKG => CK3, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK3_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexIO_CK3_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK3_SK1
#endif
        use pm_kind, only: CKG => CK3, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK2_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexIO_CK2_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK2_SK5
#endif
        use pm_kind, only: CKG => CK2, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK2_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexIO_CK2_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK2_SK4
#endif
        use pm_kind, only: CKG => CK2, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK2_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexIO_CK2_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK2_SK3
#endif
        use pm_kind, only: CKG => CK2, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK2_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexIO_CK2_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK2_SK2
#endif
        use pm_kind, only: CKG => CK2, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK2_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexIO_CK2_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK2_SK1
#endif
        use pm_kind, only: CKG => CK2, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK1_ENABLED && SK5_ENABLED
    pure elemental module subroutine setComplexIO_CK1_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK1_SK5
#endif
        use pm_kind, only: CKG => CK1, SKG => SK5
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK1_ENABLED && SK4_ENABLED
    pure elemental module subroutine setComplexIO_CK1_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK1_SK4
#endif
        use pm_kind, only: CKG => CK1, SKG => SK4
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK1_ENABLED && SK3_ENABLED
    pure elemental module subroutine setComplexIO_CK1_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK1_SK3
#endif
        use pm_kind, only: CKG => CK1, SKG => SK3
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK1_ENABLED && SK2_ENABLED
    pure elemental module subroutine setComplexIO_CK1_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK1_SK2
#endif
        use pm_kind, only: CKG => CK1, SKG => SK2
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if CK1_ENABLED && SK1_ENABLED
    pure elemental module subroutine setComplexIO_CK1_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setComplexIO_CK1_SK1
#endif
        use pm_kind, only: CKG => CK1, SKG => SK1
        complex(CKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_val2complex ! LCOV_EXCL_LINE