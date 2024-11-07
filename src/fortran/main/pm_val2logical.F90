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
!>  of different types (e.g., intrinsic Fortran strings) to logical values of different kinds.
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
!>  [test_pm_val2logical](@ref test_pm_val2logical)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_val2logical

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_val2logical"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the conversion of the input value to a `logical` value of default kind \LK.
    !>
    !>  \details
    !>  See the documentation of [setLogical](@ref pm_val2logical::setLogical) for methodological details.<br>
    !>
    !>  \param[in]  val     :   The input scalar or array of arbitrary rank of<br>
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL,<br>
    !>                          </ol>
    !>                          whose value will be converted to an output of type `logical` of default kind \LK.
    !>
    !>  \return
    !>  `conversion`        :   The scalar or array of the same shape as the input `val` of the type `logical`
    !>                          of default kind \LK, containing the conversion of the input value `val`.<br>
    !>
    !>  \interface{getLogical}
    !>  \code{.F90}
    !>
    !>      use pm_val2logical, only: getLogical
    !>
    !>      conversion = getLogical(val)
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
    !>  \example{getLogical}
    !>  \include{lineno} example/pm_val2logical/getLogical/main.F90
    !>  \compilef{getLogical}
    !>  \output{getLogical}
    !>  \include{lineno} example/pm_val2logical/getLogical/main.out.F90
    !>
    !>  \test
    !>  [test_pm_val2logical](@ref test_pm_val2logical)
    !>
    !>  \final{getLogical}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLogical

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getLogicalDef_SK5(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogicalDef_SK5
#endif
        use pm_kind, only: LKG => LK, SKG => SK5
        character(*,SKG)            , intent(in)                    :: val
        logical(LKG)                                                :: conversion
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getLogicalDef_SK4(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogicalDef_SK4
#endif
        use pm_kind, only: LKG => LK, SKG => SK4
        character(*,SKG)            , intent(in)                    :: val
        logical(LKG)                                                :: conversion
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getLogicalDef_SK3(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogicalDef_SK3
#endif
        use pm_kind, only: LKG => LK, SKG => SK3
        character(*,SKG)            , intent(in)                    :: val
        logical(LKG)                                                :: conversion
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getLogicalDef_SK2(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogicalDef_SK2
#endif
        use pm_kind, only: LKG => LK, SKG => SK2
        character(*,SKG)            , intent(in)                    :: val
        logical(LKG)                                                :: conversion
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getLogicalDef_SK1(val) result(conversion)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogicalDef_SK1
#endif
        use pm_kind, only: LKG => LK, SKG => SK1
        character(*,SKG)            , intent(in)                    :: val
        logical(LKG)                                                :: conversion
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the conversion of the input value to a `logical` value of arbitrary kind.
    !>
    !>  \param[out] conversion  :   The output scalar or array of the same shape as other array-like arguments,
    !>                              of type `logical` of kind \LKALL, containing the conversion of the input value `val`.<br>
    !>                              <ol>
    !>                                  <li>    If the input `val` is a `character`, then `conversion` is the result of the List-directed I/O action `read(val,*) conversion`.
    !>                              </ol>
    !>  \param[in]      val     :   The input scalar or array of the same shape as other array-like arguments, of<br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL,<br>
    !>                              </ol>
    !>                              whose value will be converted to an output of type `logical`.
    !>  \param[out]     iostat  :   The output `integer` of default kind \IK, of the same rank as the input `val` containing
    !>                              the Fortran IO status error code returned by the `read()` statement of the Fortran standard.<br>
    !>                              On return, it is zero <b>if and only if</b> no error occurs during the execution of the Fortran `read()` statement.<br>
    !>                              Refer to the Fortran standard `read/write` statements for the meaning of different non-zero output values for `iostat`.<br>
    !>                              (**optional**, it can be present only if `val` is of type `character`. if missing, the program will halt upon an IO error.)
    !>
    !>  \interface{setLogical}
    !>  \code{.F90}
    !>
    !>      use pm_val2logical, only: setLogical
    !>
    !>      call setLogical(conversion, val)
    !>      call setLogical(conversion, val, iostat)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  When the input `val` is an array, the output `iostat`,
    !>  if present, must also be a vector of the same shape and size as `val`.<br>
    !>
    !>  \pure
    !>
    !>  \remark
    !>  The purity of the procedures under this generic interface break when the output argument `iostat` is present.<br>
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
    !>  \example{setLogical}
    !>  \include{lineno} example/pm_val2logical/setLogical/main.F90
    !>  \compilef{setLogical}
    !>  \output{setLogical}
    !>  \include{lineno} example/pm_val2logical/setLogical/main.out.F90
    !>
    !>  \test
    !>  [test_pm_val2logical](@ref test_pm_val2logical)
    !>
    !>  \final{setLogical}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setLogical

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalDef_LK5_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK5_SK5
#endif
        use pm_kind, only: LKG => LK5, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK5_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalDef_LK5_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK5_SK4
#endif
        use pm_kind, only: LKG => LK5, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK5_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalDef_LK5_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK5_SK3
#endif
        use pm_kind, only: LKG => LK5, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK5_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalDef_LK5_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK5_SK2
#endif
        use pm_kind, only: LKG => LK5, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK5_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalDef_LK5_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK5_SK1
#endif
        use pm_kind, only: LKG => LK5, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK4_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalDef_LK4_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK4_SK5
#endif
        use pm_kind, only: LKG => LK4, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK4_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalDef_LK4_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK4_SK4
#endif
        use pm_kind, only: LKG => LK4, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK4_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalDef_LK4_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK4_SK3
#endif
        use pm_kind, only: LKG => LK4, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK4_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalDef_LK4_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK4_SK2
#endif
        use pm_kind, only: LKG => LK4, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK4_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalDef_LK4_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK4_SK1
#endif
        use pm_kind, only: LKG => LK4, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK3_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalDef_LK3_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK3_SK5
#endif
        use pm_kind, only: LKG => LK3, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK3_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalDef_LK3_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK3_SK4
#endif
        use pm_kind, only: LKG => LK3, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK3_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalDef_LK3_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK3_SK3
#endif
        use pm_kind, only: LKG => LK3, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK3_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalDef_LK3_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK3_SK2
#endif
        use pm_kind, only: LKG => LK3, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK3_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalDef_LK3_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK3_SK1
#endif
        use pm_kind, only: LKG => LK3, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK2_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalDef_LK2_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK2_SK5
#endif
        use pm_kind, only: LKG => LK2, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK2_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalDef_LK2_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK2_SK4
#endif
        use pm_kind, only: LKG => LK2, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK2_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalDef_LK2_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK2_SK3
#endif
        use pm_kind, only: LKG => LK2, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK2_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalDef_LK2_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK2_SK2
#endif
        use pm_kind, only: LKG => LK2, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK2_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalDef_LK2_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK2_SK1
#endif
        use pm_kind, only: LKG => LK2, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK1_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalDef_LK1_SK5(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK1_SK5
#endif
        use pm_kind, only: LKG => LK1, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK1_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalDef_LK1_SK4(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK1_SK4
#endif
        use pm_kind, only: LKG => LK1, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK1_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalDef_LK1_SK3(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK1_SK3
#endif
        use pm_kind, only: LKG => LK1, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK1_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalDef_LK1_SK2(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK1_SK2
#endif
        use pm_kind, only: LKG => LK1, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
    end subroutine
#endif

#if LK1_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalDef_LK1_SK1(conversion, val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalDef_LK1_SK1
#endif
        use pm_kind, only: LKG => LK1, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
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

#if LK5_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalErr_LK5_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK5_SK5
#endif
        use pm_kind, only: LKG => LK5, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK5_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalErr_LK5_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK5_SK4
#endif
        use pm_kind, only: LKG => LK5, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK5_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalErr_LK5_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK5_SK3
#endif
        use pm_kind, only: LKG => LK5, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK5_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalErr_LK5_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK5_SK2
#endif
        use pm_kind, only: LKG => LK5, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK5_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalErr_LK5_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK5_SK1
#endif
        use pm_kind, only: LKG => LK5, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK4_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalErr_LK4_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK4_SK5
#endif
        use pm_kind, only: LKG => LK4, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK4_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalErr_LK4_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK4_SK4
#endif
        use pm_kind, only: LKG => LK4, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK4_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalErr_LK4_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK4_SK3
#endif
        use pm_kind, only: LKG => LK4, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK4_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalErr_LK4_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK4_SK2
#endif
        use pm_kind, only: LKG => LK4, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK4_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalErr_LK4_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK4_SK1
#endif
        use pm_kind, only: LKG => LK4, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK3_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalErr_LK3_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK3_SK5
#endif
        use pm_kind, only: LKG => LK3, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK3_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalErr_LK3_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK3_SK4
#endif
        use pm_kind, only: LKG => LK3, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK3_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalErr_LK3_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK3_SK3
#endif
        use pm_kind, only: LKG => LK3, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK3_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalErr_LK3_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK3_SK2
#endif
        use pm_kind, only: LKG => LK3, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK3_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalErr_LK3_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK3_SK1
#endif
        use pm_kind, only: LKG => LK3, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK2_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalErr_LK2_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK2_SK5
#endif
        use pm_kind, only: LKG => LK2, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK2_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalErr_LK2_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK2_SK4
#endif
        use pm_kind, only: LKG => LK2, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK2_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalErr_LK2_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK2_SK3
#endif
        use pm_kind, only: LKG => LK2, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK2_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalErr_LK2_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK2_SK2
#endif
        use pm_kind, only: LKG => LK2, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK2_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalErr_LK2_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK2_SK1
#endif
        use pm_kind, only: LKG => LK2, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK1_ENABLED && SK5_ENABLED
    pure elemental module subroutine setLogicalErr_LK1_SK5(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK1_SK5
#endif
        use pm_kind, only: LKG => LK1, SKG => SK5
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK1_ENABLED && SK4_ENABLED
    pure elemental module subroutine setLogicalErr_LK1_SK4(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK1_SK4
#endif
        use pm_kind, only: LKG => LK1, SKG => SK4
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK1_ENABLED && SK3_ENABLED
    pure elemental module subroutine setLogicalErr_LK1_SK3(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK1_SK3
#endif
        use pm_kind, only: LKG => LK1, SKG => SK3
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK1_ENABLED && SK2_ENABLED
    pure elemental module subroutine setLogicalErr_LK1_SK2(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK1_SK2
#endif
        use pm_kind, only: LKG => LK1, SKG => SK2
        logical(LKG)                , intent(out)                   :: conversion
        character(*,SKG)            , intent(in)                    :: val
        integer(IK)                 , intent(out)                   :: iostat
    end subroutine
#endif

#if LK1_ENABLED && SK1_ENABLED
    pure elemental module subroutine setLogicalErr_LK1_SK1(conversion, val, iostat)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setLogicalErr_LK1_SK1
#endif
        use pm_kind, only: LKG => LK1, SKG => SK1
        logical(LKG)                , intent(out)                   :: conversion
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

end module pm_val2logical ! LCOV_EXCL_LINE