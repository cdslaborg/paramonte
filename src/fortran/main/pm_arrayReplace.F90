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
!>  This module contains procedures and generic interfaces for replacing patterns within arrays of various types.<br>
!>
!>  \benchmarks
!>
!>  \benchmark{getReplaced_vs_setReplaced, The runtime performance of [getReplaced](@ref pm_arrayReplace::getReplaced) vs. [setReplaced](@ref pm_arrayReplace::setReplaced)}
!>  \include{lineno} benchmark/pm_arrayReplace/getReplaced_vs_setReplaced/main.F90
!>  \compilefb{getReplaced_vs_setReplaced}
!>  \postprocb{getReplaced_vs_setReplaced}
!>  \include{lineno} benchmark/pm_arrayReplace/getReplaced_vs_setReplaced/main.py
!>  \visb{getReplaced_vs_setReplaced}
!>  \image html benchmark/pm_arrayReplace/getReplaced_vs_setReplaced/benchmark.getReplaced_vs_setReplaced.runtime.png width=1000
!>  \image html benchmark/pm_arrayReplace/getReplaced_vs_setReplaced/benchmark.getReplaced_vs_setReplaced.runtime.ratio.png width=1000
!>  \moralb{getReplaced_vs_setReplaced}
!>      -#  The procedures under the generic interface [getReplaced](@ref pm_arrayReplace::getReplaced) are functions while
!>          the procedures under the generic interface [setReplaced](@ref pm_arrayReplace::setReplaced) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs slightly less efficiently than the subroutine interface.<br>
!>          Note that this benchmark does not even include the cost of repeated reallcations, that is, the allocation of `Replaced` happen only once in all tests.<br>
!>      -#  Furthermore, the recursive `getReplaced()` implementation with recursive allocations appears to be 3-33 times
!>          slower than the subroutine implementation, depending on the size of the array within which the pattern is to be replaced.<br>
!>      -#  Note that this benchmark considers the worst-case scenario where all elements of the input `array` match the
!>          input `pattern` and must be therefore, replaced.<br>
!>
!>  \see
!>  [pm_arrayInsert](@ref pm_arrayInsert)<br>
!>  [pm_arrayRemove](@ref pm_arrayRemove)<br>
!>
!>  \test
!>  [test_pm_arrayReplace](@ref test_pm_arrayReplace)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayReplace

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayReplace"

!>  \cond excluded
!   \bug
!   The following bypasses the bug reported below that creates a conflict between Intel and gfortran.
#if     __INTEL_COMPILER
#define LEN_ARRAY :
#else
#define LEN_ARRAY len(array)
#endif
!>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   !abstract interface
   !function iseq_proc_CK3(object1, object2) result(equivalent)
   !    use pm_kind, only: CK => CK3
   !    complex(CK) , intent(in)    :: Object1(:), Object2(:)
   !    logical(LK)                 :: equivalent
   !end function
   !end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
#if 0
    !   gfortran 11 cannot run examples that contain procedure dummy arguments with explicit interface,
    !   yielding the following error: Fortran runtime error: array bound mismatch for dimension 1 of array 'segment' (0/7021782950660276225)
    !   Intel ifort 2021.6 can successfully compile and run the examples.
    abstract interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    function iseq_D0_D0_SK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG)        , intent(in)                    :: segment
        character(*,SKG)        , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if SK4_ENABLED
    function iseq_D0_D0_SK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG)        , intent(in)                    :: segment
        character(*,SKG)        , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if SK3_ENABLED
    function iseq_D0_D0_SK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG)        , intent(in)                    :: segment
        character(*,SKG)        , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if SK2_ENABLED
    function iseq_D0_D0_SK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG)        , intent(in)                    :: segment
        character(*,SKG)        , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if SK1_ENABLED
    function iseq_D0_D0_SK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG)        , intent(in)                    :: segment
        character(*,SKG)        , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    function iseq_D0_D0_IK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)            , intent(in)                    :: segment
        integer(IKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if IK4_ENABLED
    function iseq_D0_D0_IK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)            , intent(in)                    :: segment
        integer(IKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if IK3_ENABLED
    function iseq_D0_D0_IK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)            , intent(in)                    :: segment
        integer(IKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if IK2_ENABLED
    function iseq_D0_D0_IK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)            , intent(in)                    :: segment
        integer(IKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if IK1_ENABLED
    function iseq_D0_D0_IK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)            , intent(in)                    :: segment
        integer(IKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    function iseq_D0_D0_LK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)            , intent(in)                    :: segment
        logical(LKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if LK4_ENABLED
    function iseq_D0_D0_LK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)            , intent(in)                    :: segment
        logical(LKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if LK3_ENABLED
    function iseq_D0_D0_LK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)            , intent(in)                    :: segment
        logical(LKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if LK2_ENABLED
    function iseq_D0_D0_LK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)            , intent(in)                    :: segment
        logical(LKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if LK1_ENABLED
    function iseq_D0_D0_LK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)            , intent(in)                    :: segment
        logical(LKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    function iseq_D0_D0_CK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)            , intent(in)                    :: segment
        complex(CKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if CK4_ENABLED
    function iseq_D0_D0_CK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)            , intent(in)                    :: segment
        complex(CKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if CK3_ENABLED
    function iseq_D0_D0_CK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)            , intent(in)                    :: segment
        complex(CKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if CK2_ENABLED
    function iseq_D0_D0_CK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)            , intent(in)                    :: segment
        complex(CKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if CK1_ENABLED
    function iseq_D0_D0_CK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)            , intent(in)                    :: segment
        complex(CKG)            , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    function iseq_D0_D0_RK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)               , intent(in)                    :: segment
        real(RKG)               , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if RK4_ENABLED
    function iseq_D0_D0_RK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)               , intent(in)                    :: segment
        real(RKG)               , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if RK3_ENABLED
    function iseq_D0_D0_RK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)               , intent(in)                    :: segment
        real(RKG)               , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if RK2_ENABLED
    function iseq_D0_D0_RK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)               , intent(in)                    :: segment
        real(RKG)               , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

#if RK1_ENABLED
    function iseq_D0_D0_RK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D0_D0_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)               , intent(in)                    :: segment
        real(RKG)               , intent(in)                    :: pattern
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    function iseq_D1_D1_SK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: segment(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if SK4_ENABLED
    function iseq_D1_D1_SK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: segment(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if SK3_ENABLED
    function iseq_D1_D1_SK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: segment(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if SK2_ENABLED
    function iseq_D1_D1_SK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: segment(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if SK1_ENABLED
    function iseq_D1_D1_SK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: segment(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    function iseq_D1_D1_IK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: segment(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if IK4_ENABLED
    function iseq_D1_D1_IK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: segment(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if IK3_ENABLED
    function iseq_D1_D1_IK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: segment(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if IK2_ENABLED
    function iseq_D1_D1_IK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: segment(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if IK1_ENABLED
    function iseq_D1_D1_IK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: segment(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    function iseq_D1_D1_LK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: segment(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if LK4_ENABLED
    function iseq_D1_D1_LK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: segment(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if LK3_ENABLED
    function iseq_D1_D1_LK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: segment(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if LK2_ENABLED
    function iseq_D1_D1_LK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: segment(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if LK1_ENABLED
    function iseq_D1_D1_LK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: segment(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    function iseq_D1_D1_CK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: segment(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if CK4_ENABLED
    function iseq_D1_D1_CK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: segment(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if CK3_ENABLED
    function iseq_D1_D1_CK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: segment(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if CK2_ENABLED
    function iseq_D1_D1_CK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: segment(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if CK1_ENABLED
    function iseq_D1_D1_CK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: segment(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    function iseq_D1_D1_RK5(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: segment(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if RK4_ENABLED
    function iseq_D1_D1_RK4(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: segment(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if RK3_ENABLED
    function iseq_D1_D1_RK3(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: segment(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if RK2_ENABLED
    function iseq_D1_D1_RK2(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: segment(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

#if RK1_ENABLED
    function iseq_D1_D1_RK1(segment, pattern) result(equivalent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: iseq_D1_D1_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: segment(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        logical(LK)                                             :: equivalent
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
#endif
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an `arrayNew` of the same type and kind as the input `array`,
    !>  in which the requested instances of the input `pattern` have been replaced with the input `replacement`.
    !>
    !>  \details
    !>  If an input vector of `instance` is specified, representing the specific instances of pattern to change,
    !>  then only those specific instances will be changed.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or<br>
    !>                                  <li>    type `logical` of kind \LKALL, or<br>
    !>                                  <li>    type `integer` of kind \IKALL, or<br>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              or
    !>                              <ul>
    !>                                  <li>    a scalar `character` of kind \SKALL,<br>
    !>                              </ul>
    !>                              within which the specific instances of the input `pattern` must be replaced.<br>
    !>  \param[in]  pattern     :   The input `contiguous` array of rank `1` of the same type and kind as the input `array`,
    !>                              containing the pattern whose instances will have to be replaced in the input `array`.<br>
    !>  \param[in]  replacement :   The input `contiguous` array of rank `1` of the same type and kind as the input `array`,
    !>                              containing the replacement that will replace in the instances of `pattern` in the input `array`.<br>
    !>  \param      iseq        :   The `external` user-specified function that takes either two input assumed-length `character` arguments
    !>                              (if the input `array` is also an assumed-length `character`) or two array-valued **explicit-shape**
    !>                              arguments of the same type and kind as the input `array`.<br>
    !>                              It must return a scalar of type `logical` of default kind \LK that is `.true.` if all elements of
    !>                              the two input arguments are equivalent (e.g., equal) according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The the input `array` is array-valued, then the last argument to `iseq` is the length of the input `pattern`.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is array-valued,<br>
    !>                              \code{.F90}
    !>                                  function iseq(segment, pattern, lenPattern) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      integer(IK) , intent(in)    :: lenPattern
    !>                                      TYPE(KIND)  , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  character(*, SK), intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                  integer(IK)     , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                  logical(LK)     , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                  complex(CK)     , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                  real(RK)        , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                              \endcode
    !>                              where the kinds `SKG`, `IKG`, `LKG`, `CKG`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),<br>
    !>                              \code{.F90}
    !>                                  function iseq(segment, pattern) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, LK, CK, RK
    !>                                  character(*, SK), intent(in)    :: segment, pattern
    !>                                  integer(IK)     , intent(in)    :: segment, pattern
    !>                                  logical(LK)     , intent(in)    :: segment, pattern
    !>                                  complex(CK)     , intent(in)    :: segment, pattern
    !>                                  real(RK)        , intent(in)    :: segment, pattern
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where an equivalence test other than exact identity is needed,
    !>                              for example, when the array segments should match the input `pattern` only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`)
    !>  \param[in]  instance    :   The input `contiguous` array of rank `1` of type `integer` of default kind \IK,
    !>                              containing the instances of the input `pattern` in the input `array` that should be replaced with the input `replacement`.<br>
    !>                              Any element of `instance` that points to an out-of-scope instance of `pattern` in the input `array` will be ignored.<br>
    !>                              Any element of `instance` that is negatively valued will be counted from end of the input `array`.<br>
    !>                              For example, `instance = [2,-1]` requests replacing the second instance of `pattern` in `array` from the beginning and
    !>                              replacing the first instance of `pattern` starting from the end of `array`.<br>
    !>                              (**optional**, the default value corresponds to replacing all instances of `pattern` with `replacement` in `array`)
    !>  \param[in]  sorted      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all in ascending-order.<br>
    !>                              This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from
    !>                              the beginning of the input `array`.<br>Setting `sorted = .true.` will lead to faster runtime of the procedure.<br>
    !>                              However, the onus will be strictly on the user to ensure all elements of `instance` are in ascending-order.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `sorted = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>  \param[in]  unique      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all unique.<br>
    !>                              This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from
    !>                              the beginning of the input `array`.<br>Setting `unique = .true.` will lead to faster runtime of the procedure.<br>
    !>                              However, the onus will be strictly on the user to ensure all elements of `instance` are unique.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `unique = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>
    !>  \return
    !>  `arrayNew`              :   The output `allocatable` array of the same type and kind as the input `array` in which all
    !>                              requested instances of the input `pattern` have been replaced with the input `replacement`.
    !>
    !>  \interface{getReplaced}
    !>  \code{.F90}
    !>
    !>      use pm_arrayReplace, only: getReplaced
    !>
    !>      arrayNew = getReplaced(array, pattern, replacement)
    !>      arrayNew = getReplaced(array, pattern, replacement, iseq)
    !>      arrayNew = getReplaced(array, pattern, replacement, instance, sorted = sorted, unique = unique)
    !>      arrayNew = getReplaced(array, pattern, replacement, iseq, instance, sorted = sorted, unique = unique)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.<br>
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.<br>
    !>
    !>  \remark
    !>  The functions under this generic interface are slightly slower than the [setReplaced](@ref pm_arrayReplace::setReplaced) subroutine implementations.<br>
    !>  See [pm_arrayReplace](@ref pm_arrayReplace) for the relevant benchmarks.<br>
    !>
    !>  \see
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setRemoved](@ref pm_arrayRemove::setRemoved)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{getReplaced}
    !>  \include{lineno} example/pm_arrayReplace/getReplaced/main.F90
    !>  \compilef{getReplaced}
    !>  \output{getReplaced}
    !>  \include{lineno} example/pm_arrayReplace/getReplaced/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayReplace](@ref test_pm_arrayReplace)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.2.0}, \gfortran{10, 11}
    !>  \desc
    !>  The \ifort{2021.2.0} has a bug for the following interface definition,<br>
    !>  \code{.F90}
    !>      character(len(array))   , allocatable   :: arrayNew(:)
    !>  \endcode
    !>  leading to an **internal compiler error**.<br>
    !>  \remedy
    !>  For now, the remedy seems to be to redefine the interface as,
    !>  \code{.F90}
    !>      character(:, SK), allocatable   :: arrayNew(:)
    !>  \endcode
    !>  and changing the allocation method accordingly in the implementation to,
    !>  \code{.F90}
    !>      allocate(character(len(array, kind = IK)) :: arrayNew(lenArrayNew))
    !>  \endcode
    !>  However, this introduces `internal compiler error: Segmentation fault` with gfortran versions 10 and 11.
    !>  Here is a code snippet to regenerate the bug in Intel ifort (uncomment the commented line to reproduce the gfortran bug),
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(arrayNew)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array),SK)    , allocatable   :: arrayNew(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>                 !character(:, SK)            , allocatable   :: arrayNew(:) ! catastrophic internal error with gfortran 10.3. Fine with ifort 2021.2
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(arrayNew, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: arrayNew(:)
    !>          arrayNew = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>  It turns out that both gfortran and Intel do not tolerate the separation of interface from implementation in the above code snippet.
    !>  If one duplicates the interface in the implementation submodule, then both compilers compile and run the code with no errors.
    !>  This is the remedy that is currently used in this [getReplaced](@ref pm_arrayReplace::getReplaced) generic interface
    !>  (interface duplication where the bug exists). Here is a workaround example for the bug in the above code snippet,
    !>  \code{.F90}
    !>
    !>      module pm_explicitLenResult
    !>          implicit none
    !>          interface
    !>              pure module function bug(array) result(arrayNew)
    !>                  character(*, SK), intent(in), contiguous    :: array(:)
    !>                  character(len(array))   , allocatable   :: arrayNew(:) ! catastrophic internal error with ifort 2021.2. Fine with gfortran 10.3
    !>              end function
    !>          end interface
    !>      end module pm_explicitLenResult
    !>
    !>      submodule (pm_explicitLenResult) routines
    !>          implicit none
    !>      contains
    !>          module procedure bug
    !>             allocate(arrayNew, source = array)
    !>          end procedure
    !>      end submodule routines
    !>
    !>      program main
    !>          use pm_explicitLenResult, only: bug
    !>          character(2) :: array(3) = ["AA", "BB", "CC"]
    !>          character(2), allocatable :: arrayNew(:)
    !>          arrayNew = bug(array)
    !>      end program main
    !>
    !>  \endcode
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10, 11}
    !>  \desc
    !>  The \gfortran cannot run examples that use procedure dummy arguments `iseq()` with explicit interface,
    !>  yielding the following error:<br>
    !>  \code{.sh}
    !>      Fortran runtime error: array bound mismatch for dimension 1 of array 'segment' (0/7021782950660276225)
    !>  \endcode
    !>  The \ifort{2021.6} can successfully compile and run the examples.<br>
    !>  \remedy
    !>  For now, all `iseq()` procedure dummy argument interfaces are kept implicit.<br>
    !>
    !>  \todo
    !>  The internal compiler error with \ifort and \gfortran has to be fixed in the future versions.<br>
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \todo
    !>  \phigh 
    !>  This generic interface can be extended to scalar input `pattern` and `replacement` arguments.<br>
    !>  Such an extension will likely lead to runtime performance gain for cases where `pattern` and `replacement` are arrays of length `1`.<br>
    !>  See [pm_arraySplit](@ref pm_arraySplit) for a relevant benchmark about this matter.<br>
    !>
    !>  \todo
    !>  \phigh 
    !>  A benchmark comparing the performance of [setReplaced](@ref pm_arrayReplace::setReplaced) with and without `sorted, unique`
    !>  optional input arguments should be added.<br>
    !>
    !>  \final{getReplaced}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getReplaced

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComDefIns_D0_D0_D0_SK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D0_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComDefIns_D0_D0_D0_SK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D0_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComDefIns_D0_D0_D0_SK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D0_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComDefIns_D0_D0_D0_SK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D0_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(:,SKG)                        , allocatable   :: arrayNew
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComDefIns_D0_D0_D0_SK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D0_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_SK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_SK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_SK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_SK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_SK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_IK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_IK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_IK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_IK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_IK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_LK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_LK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_LK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_LK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_LK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_CK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_CK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_CK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_CK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_CK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_RK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_RK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_RK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_RK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D0_RK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_SK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_SK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_SK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_SK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_SK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_IK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_IK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_IK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_IK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_IK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_LK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_LK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_LK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_LK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_LK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_CK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_CK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_CK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_CK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_CK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_RK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_RK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_RK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_RK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D0_D1_RK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D0_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_SK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_SK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_SK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_SK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_SK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_IK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_IK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_IK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_IK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_IK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_LK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_LK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_LK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_LK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_LK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_CK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_CK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_CK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_CK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_CK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_RK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_RK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_RK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_RK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D0_RK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_SK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_SK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_SK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_SK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_SK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_IK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_IK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_IK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_IK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_IK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_LK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_LK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_LK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_LK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_LK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_CK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_CK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_CK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_CK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_CK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_RK5(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_RK4(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_RK3(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_RK2(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    PURE module function getReplacedDefComDefIns_D1_D1_D1_RK1(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComDefIns_D1_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComDefIns_D0_D0_D0_SK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D0_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComDefIns_D0_D0_D0_SK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D0_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComDefIns_D0_D0_D0_SK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D0_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComDefIns_D0_D0_D0_SK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D0_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: arrayNew
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComDefIns_D0_D0_D0_SK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D0_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_SK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_SK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_SK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_SK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_SK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_IK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_IK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_IK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_IK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_IK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_LK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_LK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_LK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_LK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_LK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_CK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_CK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_CK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_CK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_CK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_RK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_RK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_RK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_RK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D0_RK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_SK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_SK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_SK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_SK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_SK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_IK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_IK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_IK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_IK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_IK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_LK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_LK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_LK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_LK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_LK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_CK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_CK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_CK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_CK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_CK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_RK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_RK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_RK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_RK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    module function getReplacedCusComDefIns_D1_D0_D1_RK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D0_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_SK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_SK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_SK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_SK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_SK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_IK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_IK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_IK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_IK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_IK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_LK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_LK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_LK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_LK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_LK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_CK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_CK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_CK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_CK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_CK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_RK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_RK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_RK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_RK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D0_RK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_SK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_SK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_SK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_SK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_SK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_IK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_IK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_IK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_IK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_IK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_LK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_LK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_LK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_LK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_LK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_CK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_CK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_CK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_CK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_CK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_RK5(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_RK4(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_RK3(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_RK2(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    module function getReplacedCusComDefIns_D1_D1_D1_RK1(array, pattern, replacement, iseq) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComDefIns_D1_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComCusIns_D0_D0_D0_SK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D0_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComCusIns_D0_D0_D0_SK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D0_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComCusIns_D0_D0_D0_SK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D0_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComCusIns_D0_D0_D0_SK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D0_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComCusIns_D0_D0_D0_SK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D0_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_SK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_SK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_SK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_SK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_SK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_IK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_IK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_IK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_IK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_IK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_LK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_LK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_LK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_LK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_LK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_CK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_CK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_CK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_CK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_CK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_RK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_RK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_RK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_RK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D0_RK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_SK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_SK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_SK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_SK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_SK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_IK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_IK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_IK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_IK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_IK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_LK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_LK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_LK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_LK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_LK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_CK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_CK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_CK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_CK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_CK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_RK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_RK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_RK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_RK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D0_D1_RK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D0_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_SK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_SK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_SK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_SK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_SK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_IK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_IK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_IK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_IK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_IK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_LK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_LK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_LK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_LK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_LK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_CK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_CK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_CK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_CK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_CK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_RK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_RK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_RK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_RK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D0_RK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_SK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_SK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_SK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_SK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_SK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_IK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_IK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_IK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_IK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_IK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_LK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_LK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_LK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_LK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_LK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_CK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_CK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_CK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_CK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_CK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_RK5(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_RK4(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_RK3(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_RK2(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    PURE module function getReplacedDefComCusIns_D1_D1_D1_RK1(array, pattern, replacement, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedDefComCusIns_D1_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComCusIns_D0_D0_D0_SK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D0_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComCusIns_D0_D0_D0_SK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D0_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComCusIns_D0_D0_D0_SK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D0_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComCusIns_D0_D0_D0_SK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D0_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComCusIns_D0_D0_D0_SK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D0_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(:,SKG)                        , allocatable   :: arrayNew
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_SK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_SK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_SK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_SK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_SK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_IK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_IK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_IK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_IK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_IK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_LK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_LK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_LK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_LK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_LK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_CK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_CK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_CK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_CK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_CK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_RK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_RK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_RK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_RK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D0_RK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_SK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_SK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_SK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_SK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_SK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_IK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_IK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_IK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_IK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_IK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_LK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_LK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_LK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_LK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_LK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_CK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_CK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_CK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_CK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_CK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_RK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_RK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_RK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_RK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    module function getReplacedCusComCusIns_D1_D0_D1_RK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D0_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_SK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_SK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_SK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_SK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_SK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_IK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_IK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_IK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_IK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_IK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_LK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_LK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_LK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_LK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_LK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_CK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_CK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_CK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_CK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_CK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_RK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_RK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_RK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_RK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D0_RK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_SK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_SK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_SK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

#if SK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_SK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function

#endif

#if SK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_SK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        character(LEN_ARRAY,SKG)                , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_IK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_IK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_IK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if IK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_IK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if IK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_IK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        integer(IKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_LK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_LK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_LK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if LK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_LK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if LK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_LK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        logical(LKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_CK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_CK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_CK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

#if CK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_CK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function

#endif

#if CK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_CK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        complex(CKG)                            , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_RK5(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK4_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_RK4(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK3_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_RK3(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

#if RK2_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_RK2(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function

#endif

#if RK1_ENABLED
    module function getReplacedCusComCusIns_D1_D1_D1_RK1(array, pattern, replacement, iseq, instance, sorted, unique) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplacedCusComCusIns_D1_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
        real(RKG)                               , allocatable   :: arrayNew(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Replace the requested instances of the input `pattern` with the input `replacement` in the allocatable input/output `array`.
    !>
    !>  \details
    !>  If an input vector of `instance` is specified, representing the specific instances of pattern to change,
    !>  then only those specific instances will be changed.
    !>
    !>  \param[inout]   array       :   The input/output `allocatable` array of rank `1` of either <br>
    !>                                  <ul>
    !>                                      <li>    type `character` of kind \SKALL,
    !>                                      <li>    type `logical` of kind \LKALL, or
    !>                                      <li>    type `integer` of kind \IKALL, or
    !>                                      <li>    type `complex` of kind \CKALL, or
    !>                                      <li>    type `real` of kind \RKALL, or
    !>                                  </ul>
    !>                                  scalar `character` of kind \SKALL,<br>
    !>                                  within which the specific instances of the input `pattern` must be replaced.<br>
    !>                                  On output, the array will be possibly (if needed) reallocated to its new value<br>
    !>                                  with the requested instances of `pattern` replaced with `replacement.`<br>
    !>  \param[in]      pattern     :   The input `contiguous` array of rank `1` of the same type and kind as the input `array`,
    !>                                  containing the pattern whose instances will have to be replaced in the input `array`.<br>
    !>  \param[in]      replacement :   The input `contiguous` array of rank `1` of the same type and kind as the input `array`,
    !>                                  containing the replacement that will replace in the instances of `pattern` in the input `array`.<br>
    !>  \param          iseq        :   The `external` user-specified function that takes either two input assumed-length `character` arguments
    !>                                  (if the input `array` is also an assumed-length `character`) or two array-valued **explicit-shape**
    !>                                  arguments of the same type and kind as the input `array`.<br>
    !>                                  It must return a scalar of type `logical` of default kind \LK that is `.true.` if all elements of
    !>                                  the two input arguments are equivalent (e.g., equal) according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                                  The the input `array` is array-valued, then the last argument to `iseq` is the length of the input `pattern`.<br>
    !>                                  The following illustrates the generic interface of `iseq` where `pattern` is array-valued,
    !>                                  \code{.F90}
    !>                                      function iseq(segment, pattern, lenPattern) result(equivalent)
    !>                                          use pm_kind, only: IK, LK
    !>                                          integer(IK) , intent(in)    :: lenPattern
    !>                                          TYPE(KIND)  , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                          logical(LK)                 :: equivalent
    !>                                      end function
    !>                                  \endcode
    !>                                  where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                                  \code{.F90}
    !>                                      use pm_kind, only: SK, IK, LK, CK, RK
    !>                                      character(*, SK), intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                      integer(IK)     , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                      logical(LK)     , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                      complex(CK)     , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                      real(RK)        , intent(in)    :: segment(lenPattern), pattern(lenPattern)
    !>                                  \endcode
    !>                                  where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                                  The following illustrates the generic interface of `iseq` where `pattern` is scalar-valued (**including Fortran scalar strings**),
    !>                                  \code{.F90}
    !>                                      function iseq(segment, pattern) result(equivalent)
    !>                                          use pm_kind, only: LK
    !>                                          TYPE(KIND)  , intent(in)    :: segment, pattern
    !>                                          logical(LK)                 :: equivalent
    !>                                      end function
    !>                                  \endcode
    !>                                  where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                                  \code{.F90}
    !>                                      character(*, SK), intent(in)    :: segment, pattern
    !>                                      integer(IK)     , intent(in)    :: segment, pattern
    !>                                      logical(LK)     , intent(in)    :: segment, pattern
    !>                                      complex(CK)     , intent(in)    :: segment, pattern
    !>                                      real(RK)        , intent(in)    :: segment, pattern
    !>                                  \endcode
    !>                                  where the kinds `SK`, `IK`, `LK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                                  This user-defined equivalence check is extremely useful where an equivalence test other than exact identity is needed,
    !>                                  for example, when the array segments should match the input `pattern` only within a given threshold or,
    !>                                  when the case-sensitivity in character comparisons do not matter.<br>
    !>                                  In such cases, user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                                  (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`)
    !>  \param[in]      instance    :   The input `contiguous` array of rank `1` of type `integer` of default kind \IK,
    !>                                  containing the instances of the input `pattern` in the input `array` that should be replaced with the input `replacement`.<br>
    !>                                  Any element of `instance` that points to an out-of-scope instance of `pattern` in the input `array` will be ignored.<br>
    !>                                  Any element of `instance` that is negatively valued will be counted from end of the input `array`.<br>
    !>                                  For example, `instance = [2,-1]` requests replacing the second instance of `pattern` in `array` from the beginning and
    !>                                  replacing the first instance of `pattern` starting from the end of `array`.<br>
    !>                                  (**optional**, the default value corresponds to replacing all instances of `pattern` with `replacement` in `array`)
    !>  \param[in]      sorted      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all in ascending-order.<br>
    !>                                  This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from
    !>                                  the beginning of the input `array`.<br>Setting `sorted = .true.` will lead to faster runtime of the procedure.<br>
    !>                                  However, the onus will be strictly on the user to ensure all elements of `instance` are in ascending-order.<br>
    !>                                  This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                                  Therefore, set `sorted = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                                  (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>  \param[in]      unique      :   The input `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all unique.<br>
    !>                                  This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from
    !>                                  the beginning of the input `array`.<br>Setting `unique = .true.` will lead to faster runtime of the procedure.<br>
    !>                                  However, the onus will be strictly on the user to ensure all elements of `instance` are unique.<br>
    !>                                  This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                                  Therefore, set `unique = .true.` **only if** you can guarantee the validity of the condition.<br>
    !>                                  (**optional**, default = `.false.`. It can be present as input argument **only if** the input argument `instance` is present.)
    !>
    !>  \interface{setReplaced}
    !>  \code{.F90}
    !>
    !>      use pm_arrayReplace, only: setReplaced
    !>
    !>      call setReplaced(array, pattern, replacement)
    !>      call setReplaced(array, pattern, replacement, iseq)
    !>      call setReplaced(array, pattern, replacement, instance, sorted = sorted, unique = unique)
    !>      call setReplaced(array, pattern, replacement, iseq, instance, sorted = sorted, unique = unique)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.<br>
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.<br>
    !>
    !>  \remark
    !>  The logic behind making the argument `array` an `allocatable` with `intent(inout)` is to make the procedures potentially more efficient with a well-defined behavior.<br>
    !>  The alternative to the current interface is to pass an extra `intent(out), allocatable :: arrayNew(:)` argument along with an `intent(in), contiguous :: array(:)`.<br>
    !>  This, however, leads to unnecessary copies of `array` to `arrayNew` when no replacement practically has occurred.<br>
    !>  Furthermore, the size of the output `arrayNew` must be set a priori by the user which is only determined within the algorithm.<br>
    !>
    !>  \remark
    !>  The subroutines under this interface are slightly faster than the [getReplaced](@ref pm_arrayReplace::getReplaced) function implementations.<br>
    !>  See [pm_arrayReplace](@ref pm_arrayReplace) for the relevant benchmarks.
    !>
    !>  \note
    !>  The procedures under this generic interface perform replacements within the input `array` in-place.<br>
    !>  To return the result in a new array, use [getReplaced](@ref pm_arrayReplace::getReplaced).<br>
    !>  See [pm_arrayReplace](@ref pm_arrayReplace) for the relevant performance benchmarks.<br>
    !>
    !>  \note
    !>  Upon return, the `allocatable, intent(inout) :: array` argument is guaranteed to have the same lower bound as before,
    !>  although its upper bound is potentially different.
    !>
    !>  \see
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setRemoved](@ref pm_arrayRemove::setRemoved)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{setReplaced}
    !>  \include{lineno} example/pm_arrayReplace/setReplaced/main.F90
    !>  \compilef{setReplaced}
    !>  \output{setReplaced}
    !>  \include{lineno} example/pm_arrayReplace/setReplaced/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayReplace](@ref test_pm_arrayReplace)
    !>
    !>  \todo
    !>  \plow This logic behind using `intent(inout), allocatable :: array` argument as opposed to returning the result
    !>  in a new array is that the size of the output array is not known a priori, but is determined within the algorithm.<br>
    !>  If the use of `allocatable` dummy arguments is to be truly avoided, a new interface can be added where the output is
    !>  returned in a `contiguous` `ArrayReplaced` while the input array remains intact.<br>
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \todo
    !>  \phigh This generic interface can be extended to scalar non-string input `pattern` and `replacement` arguments.<br>
    !>  Such an extension will likely lead to runtime performance gain for cases where `pattern` and `replacement` are arrays of length `1`.
    !>  See [pm_arraySplit](@ref pm_arraySplit) for an example relevant benchmark about this matter.<br>
    !>
    !>  \todo
    !>  \phigh A benchmark comparing the performance of [setReplaced](@ref pm_arrayReplace::setReplaced) with and without `sorted, unique`
    !>  optional input arguments should be added.
    !>
    !>  \final{setReplaced}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setReplaced

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D0_D0_D0_SK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D0_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D0_D0_D0_SK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D0_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D0_D0_D0_SK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D0_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D0_D0_D0_SK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D0_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D0_D0_D0_SK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D0_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_SK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_SK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_SK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_SK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_SK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_IK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_IK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_IK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_IK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine

#endif

#if IK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_IK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_LK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_LK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_LK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_LK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine

#endif

#if LK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_LK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_CK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_CK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_CK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_CK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine

#endif

#if CK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_CK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_RK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_RK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_RK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_RK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
    end subroutine

#endif

#if RK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D0_RK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_SK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_SK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_SK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_SK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_SK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_IK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_IK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_IK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_IK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if IK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_IK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_LK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_LK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_LK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_LK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if LK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_LK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_CK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_CK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_CK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_CK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if CK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_CK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_RK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_RK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_RK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_RK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if RK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D0_D1_RK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D0_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_SK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_SK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_SK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_SK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_SK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_IK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_IK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_IK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_IK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine

#endif

#if IK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_IK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_LK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_LK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_LK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_LK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine

#endif

#if LK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_LK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_CK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_CK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_CK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_CK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine

#endif

#if CK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_CK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_RK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_RK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_RK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_RK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
    end subroutine

#endif

#if RK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D0_RK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_SK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_SK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_SK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_SK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_SK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_IK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_IK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_IK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_IK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if IK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_IK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_LK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_LK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_LK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_LK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if LK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_LK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_CK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_CK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_CK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_CK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if CK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_CK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_RK5(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_RK4(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_RK3(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_RK2(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
    end subroutine

#endif

#if RK1_ENABLED
    PURE module subroutine setReplacedDefComDefIns_D1_D1_D1_RK1(array, pattern, replacement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComDefIns_D1_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComDefIns_D0_D0_D0_SK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D0_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComDefIns_D0_D0_D0_SK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D0_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComDefIns_D0_D0_D0_SK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D0_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComDefIns_D0_D0_D0_SK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D0_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComDefIns_D0_D0_D0_SK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D0_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_SK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_SK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_SK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_SK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_SK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_IK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_IK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_IK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_IK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if IK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_IK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_LK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_LK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_LK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_LK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if LK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_LK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_CK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_CK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_CK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_CK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if CK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_CK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_RK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_RK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_RK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_RK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if RK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D0_RK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_SK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_SK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_SK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_SK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_SK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_IK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_IK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_IK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_IK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if IK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_IK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_LK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_LK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_LK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_LK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if LK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_LK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_CK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_CK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_CK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_CK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if CK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_CK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_RK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_RK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_RK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_RK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if RK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D0_D1_RK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D0_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_SK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_SK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_SK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_SK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_SK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_IK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_IK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_IK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_IK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if IK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_IK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_LK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_LK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_LK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_LK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if LK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_LK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_CK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_CK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_CK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_CK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if CK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_CK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_RK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_RK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_RK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_RK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if RK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D0_RK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_SK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_SK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_SK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_SK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_SK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_IK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_IK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_IK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_IK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if IK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_IK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_LK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_LK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_LK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_LK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if LK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_LK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_CK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_CK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_CK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_CK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if CK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_CK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_RK5(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_RK4(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_RK3(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_RK2(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
    end subroutine

#endif

#if RK1_ENABLED
    module subroutine setReplacedCusComDefIns_D1_D1_D1_RK1(array, pattern, replacement, iseq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComDefIns_D1_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D0_D0_D0_SK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D0_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D0_D0_D0_SK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D0_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D0_D0_D0_SK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D0_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D0_D0_D0_SK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D0_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D0_D0_D0_SK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D0_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_SK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_SK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_SK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_SK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_SK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_IK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_IK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_IK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_IK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if IK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_IK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_LK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_LK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_LK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_LK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if LK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_LK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_CK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_CK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_CK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_CK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if CK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_CK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_RK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_RK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_RK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_RK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if RK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D0_RK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_SK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_SK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_SK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_SK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_SK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_IK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_IK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_IK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_IK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if IK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_IK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_LK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_LK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_LK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_LK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if LK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_LK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_CK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_CK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_CK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_CK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if CK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_CK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_RK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_RK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_RK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_RK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if RK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D0_D1_RK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D0_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_SK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_SK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_SK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_SK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_SK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_IK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_IK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_IK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_IK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if IK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_IK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_LK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_LK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_LK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_LK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if LK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_LK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_CK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_CK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_CK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_CK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if CK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_CK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_RK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_RK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_RK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_RK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if RK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D0_RK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_SK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_SK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_SK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_SK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_SK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_IK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_IK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_IK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_IK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if IK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_IK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_LK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_LK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_LK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_LK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if LK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_LK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_CK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_CK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_CK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_CK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if CK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_CK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_RK5(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_RK4(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_RK3(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_RK2(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if RK1_ENABLED
    PURE module subroutine setReplacedDefComCusIns_D1_D1_D1_RK1(array, pattern, replacement, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedDefComCusIns_D1_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComCusIns_D0_D0_D0_SK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D0_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComCusIns_D0_D0_D0_SK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D0_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComCusIns_D0_D0_D0_SK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D0_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComCusIns_D0_D0_D0_SK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D0_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComCusIns_D0_D0_D0_SK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D0_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)        , intent(inout) , allocatable   :: array
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_SK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_SK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_SK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_SK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_SK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_IK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_IK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_IK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_IK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if IK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_IK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_LK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_LK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_LK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_LK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if LK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_LK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_CK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_CK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_CK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_CK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if CK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_CK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_RK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_RK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_RK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_RK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if RK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D0_RK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D0_D0_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_SK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_SK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_SK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_SK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_SK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)                    :: pattern
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_IK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_IK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_IK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_IK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if IK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_IK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)                    :: pattern
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_LK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_LK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_LK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_LK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if LK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_LK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)                    :: pattern
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_CK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_CK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_CK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_CK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if CK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_CK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)                    :: pattern
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_RK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_RK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_RK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_RK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if RK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D0_D1_RK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D0_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)                    :: pattern
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D0_D0_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_SK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_SK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_SK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_SK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_SK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_IK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_IK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_IK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_IK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if IK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_IK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_LK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_LK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_LK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_LK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if LK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_LK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_CK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_CK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_CK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_CK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if CK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_CK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_RK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_RK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_RK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_RK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if RK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D0_RK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)                    :: replacement
       !procedure(iseq_D1_D1_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_SK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_SK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_SK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_SK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if SK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_SK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , allocatable   :: array(:)
        character(*,SKG)        , intent(in)    , contiguous    :: pattern(:)
        character(*,SKG)        , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_SK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_IK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_IK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_IK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_IK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if IK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_IK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , allocatable   :: array(:)
        integer(IKG)            , intent(in)    , contiguous    :: pattern(:)
        integer(IKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_IK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_LK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_LK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_LK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_LK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if LK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_LK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , allocatable   :: array(:)
        logical(LKG)            , intent(in)    , contiguous    :: pattern(:)
        logical(LKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_LK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_CK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_CK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_CK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_CK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if CK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_CK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , allocatable   :: array(:)
        complex(CKG)            , intent(in)    , contiguous    :: pattern(:)
        complex(CKG)            , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_CK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_RK5(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK5)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_RK4(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK4)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_RK3(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK3)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_RK2(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK2)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
    end subroutine

#endif

#if RK1_ENABLED
    module subroutine setReplacedCusComCusIns_D1_D1_D1_RK1(array, pattern, replacement, iseq, instance, sorted, unique)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setReplacedCusComCusIns_D1_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , allocatable   :: array(:)
        real(RKG)               , intent(in)    , contiguous    :: pattern(:)
        real(RKG)               , intent(in)    , contiguous    :: replacement(:)
       !procedure(iseq_D1_D1_RK1)                               :: iseq
        procedure(logical(LK))                                  :: iseq
        integer(IK)             , intent(in)    , contiguous    :: instance(:)
        logical(LK)             , intent(in)    , optional      :: sorted
        logical(LK)             , intent(in)    , optional      :: unique
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \legacy
    !>  Generate and return a copy of the input `array` where all instances of the input `pattern` in it are replaced with the input `replacement`.<br>
    !>
    !>  \param[in]  array       :   The input scalar `character` of default kind \SK.<br>
    !>  \param[in]  pattern     :   The input scalar the same type and kind as the input `array` representing the pattern to be searched within `array`.<br>
    !>  \param[in]  replacement :   The input scalar the same type and kind as the input `array` representing the replacement for instances of `pattern` in `array`.<br>
    !>
    !>  \return
    !>  `arrayNew`              :   The output copy of the input `array` where all instances of the input `pattern` in it are replaced with the input `replacement`.<br>
    !>
    !>  \warning
    !>  This procedure is solely provided for benchmarking purposes to compare the performance difference between
    !>  a naive recursive allocation method of replacement and the optimized method implemented under the generic
    !>  interfaces in [getReplaced](@ref pm_arrayReplace::getReplaced) and [setReplaced](@ref pm_arrayReplace::setReplaced).<br>
    !>
    !>  \pure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>
    !>  [test_pm_arrayReplace](@ref test_pm_arrayReplace)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    pure recursive function getReplaced_recurs_alloc(array, pattern, replacement) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getReplaced_recurs_alloc
#endif
        use pm_kind, only: IK
        implicit none
        character(*, SK), intent(in)  :: array, pattern, replacement
        character(:, SK), allocatable :: arrayNew
        integer(IK)                   :: i, lenArray, lenPattern
        lenArray = len(array)
        lenPattern = len(pattern)
        if (lenArray == 0_IK .or. lenPattern == 0_IK) then
            arrayNew = ""
            return
        elseif (lenArray < lenPattern) then
            arrayNew = array
            return
        end if
        i = 1_IK
        do
            if (array(i:i+lenPattern-1) == pattern) then
                arrayNew = array(1:i-1) // replacement // getReplaced_recurs_alloc(array(i+lenPattern:lenArray), pattern, replacement)
                exit
            end if
            if (i+lenPattern > lenArray) then
                arrayNew = array
                exit
            end if
            i = i + 1_IK
            cycle
        end do
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayReplace ! LCOV_EXCL_LINE