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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_arrayCenter](@ref pm_arrayCenter).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \bug
!>  gfortran as of version 10.3 cannot handle regular allocation for assumed-length allocatable character types and returns the following error:
!>  Fortran runtime error: Integer overflow when calculating the amount of memory to allocate
!>  The following preprocessor condition bypasses gfortran's bug.
!#if     getCenteredAsis_D0_SK_ENABLED || getCenteredMarg_D0_SK_ENABLED || setCenteredAsis_D0_SK_ENABLED || setCenteredMarg_D0_SK_ENABLED
!#define ALLOCATE_NEW_WITH_STAT allocate(character(size + lmsize + rmsize, SKG) :: arrayCentered, stat = stat)
!#define ALLOCATE_NEW allocate(character(size + lmsize + rmsize, SKG) :: arrayCentered)
!#elif   getCenteredAsis_D1_SK_ENABLED || getCenteredMarg_D1_SK_ENABLED || setCenteredAsis_D1_SK_ENABLED || setCenteredMarg_D1_SK_ENABLED
!#define ALLOCATE_NEW_WITH_STAT allocate(character(len(array),SKG) :: arrayCentered(lb : lb + size + lmsize + rmsize - 1_IK), stat = stat)
!#define ALLOCATE_NEW allocate(character(len(array),SKG) :: arrayCentered(lb : lb + size + lmsize + rmsize - 1_IK))
!#else
!#define ALLOCATE_NEW_WITH_STAT allocate(arrayCentered(lb : lb + size + lmsize + rmsize - 1_IK), stat = stat)
!#define ALLOCATE_NEW allocate(arrayCentered(lb : lb + size + lmsize + rmsize - 1_IK))
!#endif
        integer(IK) :: i, lenArray, cbeg, cend ! content beginning location in the new array.
        ! Set the default margins.
#if     Asis_ENABLED
        integer(IK), parameter :: lmsize = 0_IK, rmsize = 0_IK
#elif   !Marg_ENABLED
#error  "Unrecognized interface."
#endif
        ! Set the array bounds.
#if     D0_ENABLED && SK_ENABLED
#define GET_SIZE(ARR) len(ARR, kind = IK)
#define GET_INDEX(i) i:i
#elif   D1_ENABLED
#define GET_SIZE(ARR) ubound(ARR, 1, IK)
#define GET_INDEX(i) i
#else
#error  "Unrecognized interface."
#endif
        integer(IK) , parameter :: lb = 1_IK
#if     setCentered_ENABLED
        integer(IK) :: size ! size of arrayCentered
        size = GET_SIZE(arrayCentered) - lmsize - rmsize ! fpp
#endif
        lenArray = GET_SIZE(array) ! fpp
        ! Verify the validity of the input.
        CHECK_ASSERTION(__LINE__, size >= 0_IK, SK_"@getCentered()/setCentered(): The condition `size >= 0_IK` must hold. size = "//getStr(size)) ! fpp
#if     Marg_ENABLED
        CHECK_ASSERTION(__LINE__, lmsize >= 0_IK, SK_"@getCentered()/setCentered(): The condition `lmsize >= 0_IK` must hold. lmsize = "//getStr(lmsize)) ! fpp
        CHECK_ASSERTION(__LINE__, rmsize >= 0_IK, SK_"@getCentered()/setCentered(): The condition `rmsize >= 0_IK` must hold. rmsize = "//getStr(rmsize)) ! fpp
#endif
#if     setCentered_ENABLED
        CHECK_ASSERTION(__LINE__, size >= 0_IK, \
        SK_"@setCentered(): `GET_SIZE(arrayCentered) - lmsize - rmsize >= 0` must hold. GET_SIZE(arrayCentered), lmsize, rmsize = "\
        //getStr([GET_SIZE(arrayCentered), lmsize, rmsize])) ! fpp
#endif
#if     setCentered_ENABLED && D1_ENABLED && SK_ENABLED
        CHECK_ASSERTION(__LINE__, len(array,IK) == len(arrayCentered,IK), \
        SK_"@setCentered(): `len(array) == len(arrayCentered)` must hold. len(array), len(arrayCentered) = "//\
        getStr([len(array,IK), len(arrayCentered,IK)])) ! fpp
#endif
        ! Fill the left margin, if any.
#if     Marg_ENABLED
        if (present(lmfill)) then
            do concurrent(i = 1 : lb + lmsize - 1_IK)
                arrayCentered(GET_INDEX(i)) = lmfill
            end do
        end if
#endif
        ! Center the array contents in the new array.
        if (size > lenArray) then
            cbeg = lb + lmsize + (size - lenArray) / 2_IK
            cend = cbeg + lenArray - 1_IK
            if (present(fill)) then
                do concurrent(i = lb + lmsize : cbeg - 1_IK)
                    arrayCentered(GET_INDEX(i)) = fill
                end do
                arrayCentered(cbeg:cend) = array
                do concurrent(i = cend + 1_IK : lb + lmsize + size - 1_IK)
                    arrayCentered(GET_INDEX(i)) = fill
                end do
            else
                arrayCentered(cbeg:cend) = array
            end if
        else
            cbeg = lb + (lenArray - size) / 2_IK
            cend = cbeg + size - 1_IK
            arrayCentered(lb + lmsize : lb + lmsize + size - 1_IK) = array(cbeg:cend)
        end if
        ! Fill the right margin, if any.
#if     Marg_ENABLED
        if (present(rmfill)) then
            do concurrent(i = lb + lmsize + size : lb + lmsize + size + rmsize - 1_IK)
                arrayCentered(GET_INDEX(i)) = rmfill
            end do
        end if
#endif

#undef  GET_INDEX
#undef  GET_SIZE
