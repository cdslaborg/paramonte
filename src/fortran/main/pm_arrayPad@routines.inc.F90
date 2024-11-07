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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_arrayPad](@ref pm_arrayPad).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>>>>>>>>>>>>>>>>>>>>>>>>
!!>  \bug
!!>  gfortran as of version 10.3 cannot handle regular allocation for assumed-length allocatable character types and returns the following error:
!!>  Fortran runtime error: Integer overflow when calculating the amount of memory to allocate
!!>  The following preprocessor condition bypasses gfortran's bug.
!#if     setPaddedAsisSB_D0_SK_ENABLED || setPaddedMargSB_D0_SK_ENABLED || getPaddedAsisSB_D0_SK_ENABLED || getPaddedMargSB_D0_SK_ENABLED
!#define ALLOCATE_NEW_WITH_STAT allocate(character(lenarrayPaddedMinusOne+1_IK,SK) :: arrayPadded, stat = stat)
!#define ALLOCATE_NEW allocate(character(lenarrayPaddedMinusOne+1_IK,SK) :: arrayPadded)
!#elif   setPaddedAsisSB_D1_SK_ENABLED || setPaddedMargSB_D1_SK_ENABLED || getPaddedAsisSB_D1_SK_ENABLED || getPaddedMargSB_D1_SK_ENABLED
!#define ALLOCATE_NEW_WITH_STAT allocate(character(len(array)) :: arrayPadded(lb : lb + lenarrayPaddedMinusOne), stat = stat)
!#define ALLOCATE_NEW allocate(character(len(array)) :: arrayPadded(lb : lb + lenarrayPaddedMinusOne))
!#else
!#define ALLOCATE_NEW_WITH_STAT allocate(arrayPadded(lb : lb + lenarrayPaddedMinusOne), stat = stat)
!#define ALLOCATE_NEW allocate(arrayPadded(lb : lb + lenarrayPaddedMinusOne))
!#endif
!<<<<<<<<<<<<<<<<<<<<<<<

        integer(IK) :: i, lenArray, lenarrayPaddedMinusOne
#if     Asis_ENABLED && SB_ENABLED
        integer(IK), parameter :: lmsize = 0_IK, rmsize = 0_IK
#elif 	Asis_ENABLED && SL_ENABLED
        integer(IK), parameter :: lmsize = 0_IK
#elif 	Asis_ENABLED && SR_ENABLED
        integer(IK), parameter :: rmsize = 0_IK
#elif   !Marg_ENABLED
#error  "Unrecognized interface."
#endif
#if     setPadded_ENABLED
        integer(IK) :: stat
#if     SK_ENABLED && D0_ENABLED
        character(:,SKG), allocatable :: arrayPadded
#elif   SK_ENABLED && D1_ENABLED
        character(len(array,IK),SKG), allocatable :: arrayPadded(:)
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG), allocatable :: arrayPadded(:)
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG), allocatable :: arrayPadded(:)
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG), allocatable :: arrayPadded(:)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG), allocatable :: arrayPadded(:)
#else
#error  "Unrecognized interface."
#endif
#elif   !getPadded_ENABLED
#error  "Unrecognized interface."
#endif
        ! Set the array bounds.
#if     SK_ENABLED && D0_ENABLED
#define GET_INDEX(i) i:i
        integer(IK) , parameter :: lb = 1_IK
        lenArray = len(array, kind = IK)
#elif   D1_ENABLED
#define GET_INDEX(i) i
#if     getPadded_ENABLED
        integer(IK) , parameter :: lb = 1_IK
#elif   setPadded_ENABLED
        integer(IK) :: lb
        lb = lbound(array, dim = 1, kind = IK)
#endif
        lenArray = size(array, kind = IK)
#else
#error  "Unrecognized interface."
#endif
        ! Verify the validity of the input.
#if     SB_ENABLED || SL_ENABLED
        CHECK_ASSERTION(__LINE__, lpsize >= 0_IK, SK_"The condition `lpsize >= 0_IK` must hold. lpsize = "//getStr(lpsize))
#endif
#if     SB_ENABLED || SR_ENABLED
        CHECK_ASSERTION(__LINE__, rpsize >= 0_IK, SK_"The condition `rpsize >= 0_IK` must hold. rpsize = "//getStr(rpsize))
#endif
#if     Marg_ENABLED && (SB_ENABLED || SL_ENABLED)
        CHECK_ASSERTION(__LINE__, lmsize >= 0_IK, SK_"The condition `lmsize >= 0_IK` must hold. lmsize = "//getStr(lmsize))
#endif
#if     Marg_ENABLED && (SB_ENABLED || SR_ENABLED)
        CHECK_ASSERTION(__LINE__, rmsize >= 0_IK, SK_"The condition `rmsize >= 0_IK` must hold. rmsize = "//getStr(rmsize))
#endif
        ! Set the length of padded array.
#if     SB_ENABLED
        lenarrayPaddedMinusOne = lenArray + lpsize + rpsize + lmsize + rmsize - 1_IK
#elif   SL_ENABLED
        lenarrayPaddedMinusOne = lenArray + lpsize + lmsize - 1_IK
#elif   SR_ENABLED
        lenarrayPaddedMinusOne = lenArray + rpsize + rmsize - 1_IK
#else
#error  "Unrecognized interface."
#endif
        ! Allocate the new array for the subroutine interface.
#if     setPadded_ENABLED
        if (present(failed)) then
#if         SK_ENABLED && D0_ENABLED
            allocate(character(lenarrayPaddedMinusOne + 1_IK, SKG) :: arrayPadded, stat = stat)
#elif       SK_ENABLED && D1_ENABLED
            allocate(character(len(array,IK),SKG) :: arrayPadded(lb : lb + lenarrayPaddedMinusOne), stat = stat)
#else
            allocate(arrayPadded(lb : lb + lenarrayPaddedMinusOne), stat = stat)
#endif
            failed = logical(stat > 0_IK, LK)
            if (failed) return ! LCOV_EXCL_LINE
        else
#if         SK_ENABLED && D0_ENABLED
            allocate(character(lenarrayPaddedMinusOne + 1_IK, SKG) :: arrayPadded)
#elif       SK_ENABLED && D1_ENABLED
            allocate(character(len(array,IK),SKG) :: arrayPadded(lb : lb + lenarrayPaddedMinusOne))
#else
            allocate(arrayPadded(lb : lb + lenarrayPaddedMinusOne))
#endif
        end if
#endif

        ! Fill the left margin, if any.

#if     Marg_ENABLED && (SB_ENABLED || SL_ENABLED)
        if (present(lmfill)) then
            do concurrent(i = lb : lb + lmsize - 1_IK)
                arrayPadded(GET_INDEX(i)) = lmfill
            end do
        end if
#endif

        ! Pad the array contents in the new array.

#if     SB_ENABLED || SL_ENABLED
        do concurrent(i = lb + lmsize : lb + lmsize + lpsize - 1_IK)
            arrayPadded(GET_INDEX(i)) = lpfill
        end do
        arrayPadded(lb + lmsize + lpsize : lb + lmsize + lpsize + lenArray - 1_IK) = array
#endif
#if     SB_ENABLED
        do concurrent(i = lb + lmsize + lpsize + lenArray : lb + lenarrayPaddedMinusOne - rmsize)
            arrayPadded(GET_INDEX(i)) = rpfill
        end do
#elif   SR_ENABLED
        arrayPadded(lb : lb + lenArray - 1_IK) = array
        do concurrent(i = lb + lenArray : lb + lenArrayPaddedMinusOne - rmsize)
            arrayPadded(GET_INDEX(i)) = rpfill
        end do
#endif
        ! Fill the right margin, if any.

#if     Marg_ENABLED && (SB_ENABLED || SR_ENABLED)
        if (present(rmfill)) then
            do concurrent(i = lb + lenarrayPaddedMinusOne - rmsize + 1_IK : lb + lenarrayPaddedMinusOne)
                arrayPadded(GET_INDEX(i)) = rmfill
            end do
        end if
#endif

#if     setPadded_ENABLED
        call move_alloc(arrayPadded, array)
#endif

#undef  GET_INDEX
