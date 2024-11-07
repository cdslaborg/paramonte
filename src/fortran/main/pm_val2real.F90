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
!>  of different types (e.g., intrinsic Fortran string and logical) to real values of different kinds.
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
!>  [test_pm_val2real](@ref test_pm_val2real)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_val2real

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_val2real"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the conversion of the input value to a `real` value of default kind \RK.
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of either<br>
    !>                          &nbsp; type `character` of kind \SKALL, or<br>
    !>                          &nbsp; type `logical` of kind \LKALL, or<br>
    !>                          &nbsp; whose value will be converted to an output of type `real` of default kind \RK.
    !>
    !>  \return
    !>  `conversion`        :   The scalar or array of the same shape as the input `val` of the type `real`
    !>                          of default kind \RK, containing the conversion of the input value `val`.<br>
    !>                              -#  If the input `val` is a `logical` that evaluates to `.true.`, then `conversion = 1._RK`, otherwise `conversion = 0._RK`.
    !>                              -#  If the input `val` is a `character`, then `conversion` is the result of the List-directed I/O action `read(val,*) conversion`.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_val2real, only: getReal
    !>
    !>      conversion = getReal(val)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
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
    !>  \include{lineno} example/pm_val2real/getReal/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_val2real/getReal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_val2real](@ref test_pm_val2real)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getReal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure elemental module function getRealDef_LK5(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_LK5
#endif
        use pm_kind, only: RKG => RK, LKG => LK5
        logical(LKG)                , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

#if LK4_ENABLED
    pure elemental module function getRealDef_LK4(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_LK4
#endif
        use pm_kind, only: RKG => RK, LKG => LK4
        logical(LKG)                , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

#if LK3_ENABLED
    pure elemental module function getRealDef_LK3(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_LK3
#endif
        use pm_kind, only: RKG => RK, LKG => LK3
        logical(LKG)                , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

#if LK2_ENABLED
    pure elemental module function getRealDef_LK2(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_LK2
#endif
        use pm_kind, only: RKG => RK, LKG => LK2
        logical(LKG)                , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

#if LK1_ENABLED
    pure elemental module function getRealDef_LK1(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_LK1
#endif
        use pm_kind, only: RKG => RK, LKG => LK1
        logical(LKG)                , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getRealDef_SK5(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_SK5
#endif
        use pm_kind, only: RKG => RK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getRealDef_SK4(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_SK4
#endif
        use pm_kind, only: RKG => RK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getRealDef_SK3(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_SK3
#endif
        use pm_kind, only: RKG => RK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getRealDef_SK2(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_SK2
#endif
        use pm_kind, only: RKG => RK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getRealDef_SK1(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRealDef_SK1
#endif
        use pm_kind, only: RKG => RK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: val
        real(RKG)                                                   :: conversion
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the conversion of the input value to a `real` value of arbitrary kind.
    !>
    !>  \param[out] conversion  :   The output scalar or array of the same shape as other array-like arguments,
    !>                              of type `real` of kind \RKALL, containing the conversion of the input value `val`.<br>
    !>                              <ol>
    !>                                  <li>    If the input `val` is a `logical` that evaluates to `.true.`, then `conversion = 1.`, otherwise `conversion = 0.`.
    !>                                  <li>    If the input `val` is a `character`, then `conversion` is the result of the List-directed I/O action `read(val,*) conversion`.
    !>                              </ol>
    !>  \param[in]      val     :   The input scalar or array of the same shape as other array-like arguments, of either<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL, or<br>
    !>                                  <li>    type `logical` of kind \LKALL, or<br>
    !>                              </ol>
    !>                              whose value will be converted to an output of type `real`.
    !>  \param[out]     iostat  :   The output `integer` of default kind \IK, of the same rank as the input `val` containing
    !>                              the Fortran IO status error code returned by the `read()` statement of the Fortran standard.<br>
    !>                              On return, it is zero <b>if and only if</b> no error occurs during the execution of the Fortran `read()` statement.<br>
    !>                              Refer to the Fortran standard `read/write` statements for the meaning of different non-zero output values for `iostat`.<br>
    !>                              (**optional**, it can be present only if `val` is of type `character`. if missing, the program will halt upon an IO error.)
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_val2real, only: setReal
    !>
    !>      call setReal(conversion, val)
    !>      call setReal(conversion, val, iostat)
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
    !>  \include{lineno} example/pm_val2real/setReal/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_val2real/setReal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_val2real](@ref test_pm_val2real)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setReal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && LK5_ENABLED
    pure elemental module subroutine setRealDef_RK5_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_LK5
#endif
        use pm_kind, only: RKG => RK5, LKG => LK5
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK5_ENABLED && LK4_ENABLED
    pure elemental module subroutine setRealDef_RK5_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_LK4
#endif
        use pm_kind, only: RKG => RK5, LKG => LK4
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK5_ENABLED && LK3_ENABLED
    pure elemental module subroutine setRealDef_RK5_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_LK3
#endif
        use pm_kind, only: RKG => RK5, LKG => LK3
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK5_ENABLED && LK2_ENABLED
    pure elemental module subroutine setRealDef_RK5_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_LK2
#endif
        use pm_kind, only: RKG => RK5, LKG => LK2
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK5_ENABLED && LK1_ENABLED
    pure elemental module subroutine setRealDef_RK5_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_LK1
#endif
        use pm_kind, only: RKG => RK5, LKG => LK1
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealDef_RK5_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_SK5
#endif
        use pm_kind, only: RKG => RK5, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK5_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealDef_RK5_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_SK4
#endif
        use pm_kind, only: RKG => RK5, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK5_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealDef_RK5_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_SK3
#endif
        use pm_kind, only: RKG => RK5, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK5_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealDef_RK5_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_SK2
#endif
        use pm_kind, only: RKG => RK5, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK5_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealDef_RK5_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK5_SK1
#endif
        use pm_kind, only: RKG => RK5, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK4_ENABLED && LK5_ENABLED
    pure elemental module subroutine setRealDef_RK4_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_LK5
#endif
        use pm_kind, only: RKG => RK4, LKG => LK5
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK4_ENABLED && LK4_ENABLED
    pure elemental module subroutine setRealDef_RK4_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_LK4
#endif
        use pm_kind, only: RKG => RK4, LKG => LK4
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK4_ENABLED && LK3_ENABLED
    pure elemental module subroutine setRealDef_RK4_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_LK3
#endif
        use pm_kind, only: RKG => RK4, LKG => LK3
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK4_ENABLED && LK2_ENABLED
    pure elemental module subroutine setRealDef_RK4_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_LK2
#endif
        use pm_kind, only: RKG => RK4, LKG => LK2
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK4_ENABLED && LK1_ENABLED
    pure elemental module subroutine setRealDef_RK4_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_LK1
#endif
        use pm_kind, only: RKG => RK4, LKG => LK1
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK4_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealDef_RK4_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_SK5
#endif
        use pm_kind, only: RKG => RK4, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK4_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealDef_RK4_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_SK4
#endif
        use pm_kind, only: RKG => RK4, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK4_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealDef_RK4_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_SK3
#endif
        use pm_kind, only: RKG => RK4, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK4_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealDef_RK4_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_SK2
#endif
        use pm_kind, only: RKG => RK4, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK4_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealDef_RK4_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK4_SK1
#endif
        use pm_kind, only: RKG => RK4, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED && LK5_ENABLED
    pure elemental module subroutine setRealDef_RK3_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_LK5
#endif
        use pm_kind, only: RKG => RK3, LKG => LK5
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK3_ENABLED && LK4_ENABLED
    pure elemental module subroutine setRealDef_RK3_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_LK4
#endif
        use pm_kind, only: RKG => RK3, LKG => LK4
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK3_ENABLED && LK3_ENABLED
    pure elemental module subroutine setRealDef_RK3_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_LK3
#endif
        use pm_kind, only: RKG => RK3, LKG => LK3
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK3_ENABLED && LK2_ENABLED
    pure elemental module subroutine setRealDef_RK3_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_LK2
#endif
        use pm_kind, only: RKG => RK3, LKG => LK2
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK3_ENABLED && LK1_ENABLED
    pure elemental module subroutine setRealDef_RK3_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_LK1
#endif
        use pm_kind, only: RKG => RK3, LKG => LK1
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealDef_RK3_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_SK5
#endif
        use pm_kind, only: RKG => RK3, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK3_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealDef_RK3_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_SK4
#endif
        use pm_kind, only: RKG => RK3, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK3_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealDef_RK3_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_SK3
#endif
        use pm_kind, only: RKG => RK3, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK3_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealDef_RK3_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_SK2
#endif
        use pm_kind, only: RKG => RK3, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK3_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealDef_RK3_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK3_SK1
#endif
        use pm_kind, only: RKG => RK3, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK2_ENABLED && LK5_ENABLED
    pure elemental module subroutine setRealDef_RK2_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_LK5
#endif
        use pm_kind, only: RKG => RK2, LKG => LK5
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK2_ENABLED && LK4_ENABLED
    pure elemental module subroutine setRealDef_RK2_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_LK4
#endif
        use pm_kind, only: RKG => RK2, LKG => LK4
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK2_ENABLED && LK3_ENABLED
    pure elemental module subroutine setRealDef_RK2_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_LK3
#endif
        use pm_kind, only: RKG => RK2, LKG => LK3
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK2_ENABLED && LK2_ENABLED
    pure elemental module subroutine setRealDef_RK2_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_LK2
#endif
        use pm_kind, only: RKG => RK2, LKG => LK2
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK2_ENABLED && LK1_ENABLED
    pure elemental module subroutine setRealDef_RK2_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_LK1
#endif
        use pm_kind, only: RKG => RK2, LKG => LK1
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK2_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealDef_RK2_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_SK5
#endif
        use pm_kind, only: RKG => RK2, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK2_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealDef_RK2_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_SK4
#endif
        use pm_kind, only: RKG => RK2, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK2_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealDef_RK2_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_SK3
#endif
        use pm_kind, only: RKG => RK2, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK2_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealDef_RK2_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_SK2
#endif
        use pm_kind, only: RKG => RK2, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK2_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealDef_RK2_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK2_SK1
#endif
        use pm_kind, only: RKG => RK2, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK1_ENABLED && LK5_ENABLED
    pure elemental module subroutine setRealDef_RK1_LK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_LK5
#endif
        use pm_kind, only: RKG => RK1, LKG => LK5
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK1_ENABLED && LK4_ENABLED
    pure elemental module subroutine setRealDef_RK1_LK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_LK4
#endif
        use pm_kind, only: RKG => RK1, LKG => LK4
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK1_ENABLED && LK3_ENABLED
    pure elemental module subroutine setRealDef_RK1_LK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_LK3
#endif
        use pm_kind, only: RKG => RK1, LKG => LK3
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK1_ENABLED && LK2_ENABLED
    pure elemental module subroutine setRealDef_RK1_LK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_LK2
#endif
        use pm_kind, only: RKG => RK1, LKG => LK2
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

#if RK1_ENABLED && LK1_ENABLED
    pure elemental module subroutine setRealDef_RK1_LK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_LK1
#endif
        use pm_kind, only: RKG => RK1, LKG => LK1
        real(RKG)                   , intent(out)                   :: conversion
        logical(LKG)                , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK1_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealDef_RK1_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_SK5
#endif
        use pm_kind, only: RKG => RK1, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK1_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealDef_RK1_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_SK4
#endif
        use pm_kind, only: RKG => RK1, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK1_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealDef_RK1_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_SK3
#endif
        use pm_kind, only: RKG => RK1, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK1_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealDef_RK1_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_SK2
#endif
        use pm_kind, only: RKG => RK1, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if RK1_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealDef_RK1_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealDef_RK1_SK1
#endif
        use pm_kind, only: RKG => RK1, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
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

#if RK5_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealErr_RK5_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK5_SK5
#endif
        use pm_kind, only: RKG => RK5, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK5_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealErr_RK5_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK5_SK4
#endif
        use pm_kind, only: RKG => RK5, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK5_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealErr_RK5_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK5_SK3
#endif
        use pm_kind, only: RKG => RK5, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK5_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealErr_RK5_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK5_SK2
#endif
        use pm_kind, only: RKG => RK5, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK5_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealErr_RK5_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK5_SK1
#endif
        use pm_kind, only: RKG => RK5, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK4_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealErr_RK4_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK4_SK5
#endif
        use pm_kind, only: RKG => RK4, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK4_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealErr_RK4_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK4_SK4
#endif
        use pm_kind, only: RKG => RK4, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK4_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealErr_RK4_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK4_SK3
#endif
        use pm_kind, only: RKG => RK4, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK4_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealErr_RK4_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK4_SK2
#endif
        use pm_kind, only: RKG => RK4, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK4_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealErr_RK4_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK4_SK1
#endif
        use pm_kind, only: RKG => RK4, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealErr_RK3_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK3_SK5
#endif
        use pm_kind, only: RKG => RK3, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK3_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealErr_RK3_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK3_SK4
#endif
        use pm_kind, only: RKG => RK3, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK3_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealErr_RK3_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK3_SK3
#endif
        use pm_kind, only: RKG => RK3, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK3_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealErr_RK3_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK3_SK2
#endif
        use pm_kind, only: RKG => RK3, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK3_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealErr_RK3_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK3_SK1
#endif
        use pm_kind, only: RKG => RK3, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK2_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealErr_RK2_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK2_SK5
#endif
        use pm_kind, only: RKG => RK2, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK2_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealErr_RK2_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK2_SK4
#endif
        use pm_kind, only: RKG => RK2, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK2_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealErr_RK2_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK2_SK3
#endif
        use pm_kind, only: RKG => RK2, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK2_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealErr_RK2_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK2_SK2
#endif
        use pm_kind, only: RKG => RK2, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK2_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealErr_RK2_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK2_SK1
#endif
        use pm_kind, only: RKG => RK2, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK1_ENABLED && SK5_ENABLED
    pure elemental module subroutine setRealErr_RK1_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK1_SK5
#endif
        use pm_kind, only: RKG => RK1, SKG => SK5
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK1_ENABLED && SK4_ENABLED
    pure elemental module subroutine setRealErr_RK1_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK1_SK4
#endif
        use pm_kind, only: RKG => RK1, SKG => SK4
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK1_ENABLED && SK3_ENABLED
    pure elemental module subroutine setRealErr_RK1_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK1_SK3
#endif
        use pm_kind, only: RKG => RK1, SKG => SK3
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK1_ENABLED && SK2_ENABLED
    pure elemental module subroutine setRealErr_RK1_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK1_SK2
#endif
        use pm_kind, only: RKG => RK1, SKG => SK2
        real(RKG)                   , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if RK1_ENABLED && SK1_ENABLED
    pure elemental module subroutine setRealErr_RK1_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRealErr_RK1_SK1
#endif
        use pm_kind, only: RKG => RK1, SKG => SK1
        real(RKG)                   , intent(out)                   :: conversion
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

end module pm_val2real ! LCOV_EXCL_LINE