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
!>  This file contains the implementation details of the routines under the generic interfaces in [pm_arrayInsert](@ref pm_arrayInsert).
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Determine assumed-length scalar character vs. array input arguments.
#if     D0_D0_ENABLED && SK_ENABLED
#define GET_SIZE len
#define ALL
#elif   D1_D0_ENABLED || D1_D1_ENABLED
#define GET_SIZE size
#define ALL all
#else
#error  "Unrecognized interface."
#endif
        ! Define logical vs. normal equivalence operators
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#elif   SK_ENABLED || IK_ENABLED || CK_ENABLED || RK_ENABLED
#define IS_EQUAL ==
#else
#error  "Unrecognized interface."
#endif
        integer(IK)                 :: lenArray, lenArrayNew, lenIndex
        integer(IK)                 :: indexNew(size(index))
        integer(IK)                 :: start, stop, i
        logical(LK)                 :: negative
        logical(LK)                 :: unsorted
        ! Scalar vs. vector `insertion` argument.
#if     D1_D0_ENABLED
        integer(IK) , parameter     :: lenInsertion = 1_IK
#elif   D0_D0_ENABLED   || D1_D1_ENABLED
        integer(IK)                 :: lenInsertion
        lenInsertion = GET_SIZE(insertion, kind = IK) ! \warning GET_SIZE is a preprocessor macro.
        if (lenInsertion == 0_IK) then
            arrayNew = array
            return
        end if
#else
#error  "Unrecognized interface."
#endif
        lenIndex = size(index, kind = IK)
        if (lenIndex > 0_IK) then
            lenArray = GET_SIZE(array, kind = IK) ! \warning GET_SIZE is a preprocessor macro.
            lenArrayNew = GET_SIZE(arrayNew, kind = IK) ! \warning GET_SIZE is a preprocessor macro.
            CHECK_ASSERTION(__LINE__, All(abs(index) <= lenArray + 1_IK), \
            SK_": The condition `All(abs(index) <= lenArray + 1_IK)` must hold. lenArray, index = "\
            //getStr([lenArray,index])) ! fpp
#if         SK_ENABLED && (D1_D0_ENABLED || D1_D1_ENABLED)
            CHECK_ASSERTION(__LINE__, len(array, IK) == len(insertion, IK), \
            SK_": The condition `len(array) == len(insertion)` must hold. len(array), len(insertion) = "\
            //getStr([len(array, IK), len(insertion, IK)])) ! fpp
#endif
#if         SK_ENABLED && (D1_D0_ENABLED || D1_D1_ENABLED) && setInserted_ENABLED
            CHECK_ASSERTION(__LINE__, len(array, IK) == len(arrayNew, IK), \
            SK_": The condition `len(array) == len(arrayNew)` must hold. len(array) == len(arrayNew) = "\
            //getStr([len(array), len(arrayNew)])) ! fpp
#endif
#if         (D0_D0_ENABLED || D1_D1_ENABLED) && setInserted_ENABLED
            CHECK_ASSERTION(__LINE__, lenArrayNew == lenArray + lenIndex * lenInsertion, \
            SK_": The condition `lenArrayNew /= lenArray + lenIndex * lenInsertion` must hold. lenArrayNew, lenArray, lenIndex, lenInsertion = "\
            //getStr([lenArrayNew, lenArray, lenIndex, lenInsertion])) ! fpp
#endif
#if         D1_D0_ENABLED || D1_D1_ENABLED
            CHECK_ASSERTION(__LINE__, lenArrayNew == lenArray + lenIndex * lenInsertion, \
            SK_": The condition lenArrayNew == lenArray + lenIndex * lenInsertion must hold. lenArrayNew, lenArray, lenIndex, lenInsertion = "\
            //getStr([lenArrayNew, lenArray, lenIndex, lenInsertion])) ! fpp
#endif
            ! Set the positivity of index.
            if (present(positive)) then
                negative = .not. positive
            else
                negative = .true._LK
            end if
            if (negative) then
                do i = 1_IK, lenIndex
                    if (index(i) > 0_IK) then
                        indexNew(i) = index(i)
                    else
                        indexNew(i) = lenArray + index(i) + 1_IK
                    end if
                end do
            else
                indexNew(1_IK:lenIndex) = index
            end if
            ! Sort the index if needed.
            if (present(sorted)) then
                unsorted = .not. sorted
            else
                unsorted = .true._LK
            end if
            if (unsorted) call setSorted(indexNew)
            ! Insert the `insertion` into the proper locations.
            if (indexNew(1) > 1_IK) arrayNew(1_IK:indexNew(1)-1_IK) = array(1_IK:indexNew(1)-1_IK)
            start = indexNew(1)
            do i = 2_IK, lenIndex
                stop = start + lenInsertion - 1_IK
                arrayNew(start:stop) = insertion
                start = stop + 1_IK
                stop = start + indexNew(i) - indexNew(i-1_IK) - 1_IK
                arrayNew(start:stop) = array(indexNew(i-1_IK):indexNew(i)-1_IK)
                start = stop + 1_IK
            end do
            stop = start + lenInsertion - 1_IK
            arrayNew(start:stop) = insertion
            start = stop + 1_IK
            arrayNew(start:lenArrayNew) = array(indexNew(lenIndex):lenArray)
        else ! index is empty.
            arrayNew = array
            return
        end if

#undef  INDEXNEW
#undef  IS_EQUAL
#undef  GET_SIZE
#undef  ALL