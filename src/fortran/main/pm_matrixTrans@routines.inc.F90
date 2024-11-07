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
!>  This file contains procedure implementations of [pm_matrixTrans](@ref pm_matrixTrans).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#if     setMatTrans_ENABLED
        !%%%%%%%%%%%%%%%%%%

        ! Set transposition type.
#if     Herm_ENABLED
#define OPERATION, operation
#define GET_CONJG(X) conjg(X)
#elif   Symm_ENABLED
#define OPERATION
#define GET_CONJG(X) X
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: imid ! middle row/column of source.
        integer(IK) :: nrow, ncol ! number of row and column of source.
#if     Fix_ENABLED
#define BSIZE_ARG
        integer(IK) , parameter :: BSIZE = 32_IK
#elif   Arb_ENABLED
#define BSIZE_ARG , bsize
        CHECK_ASSERTION(__LINE__, 0_IK < bsize, SK_"@setMatTrans(): The condition `0 < bsize` must hold. bsize = "//getStr(bsize))
#else
#error  "Unrecognized interface."
#endif
        ! Set in-place transposition rule.
#if     Old_ENABLED
#define destin source
        CHECK_ASSERTION(__LINE__, size(source, 1, IK) == size(source, 2, IK), \
        SK_"@setMatTrans(): The condition `size(source, 1) == size(source, 2)` must hold. shape(source) = "//getStr(shape(source, IK)))
#elif   New_ENABLED
        CHECK_ASSERTION(__LINE__, size(source, 1, IK) == size(destin, 2, IK) .and. size(source, 2, IK) == size(destin, 1, IK), \
        SK_"@setMatTrans(): The condition `size(source, 1) == size(destin, 2) .and. size(source, 2) == size(destin,1)` must hold. shape(source), shape(destin) = "//\
        getStr([shape(source, IK), shape(destin, IK)]))
#else
#error  "Unrecognized interface."
#endif
        nrow = size(source, 1, IK)
        ncol = size(source, 2, IK)
        if (nrow <= BSIZE .and. ncol <= BSIZE) then
            destin = transpose(GET_CONJG(source))
        elseif (nrow < ncol) then
            imid = ncol / 2_IK
            call setMatTrans(source(1 : nrow, 1 : imid), destin(1 : imid, 1 : nrow)OPERATION BSIZE_ARG)
            call setMatTrans(source(1 : nrow, imid + 1 : ncol), destin(imid + 1 : ncol, 1 : nrow)OPERATION BSIZE_ARG)
        else
            imid = nrow / 2_IK
            call setMatTrans(source(1 : imid, 1 : ncol), destin(1 : ncol, 1 : imid)OPERATION BSIZE_ARG)
            call setMatTrans(source(imid + 1 : nrow, 1 : ncol), destin(1 : ncol, imid + 1 : nrow)OPERATION BSIZE_ARG)
        end if
#undef  OPERATION
#undef  GET_CONJG
#undef  BSIZE_ARG
#undef  destin

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif