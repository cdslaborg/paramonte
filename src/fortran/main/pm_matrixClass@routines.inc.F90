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
!>  This include file contains procedure implementations of [pm_matrixClass](@ref pm_matrixClass).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the logical equality check.
#if     LK_ENABLED
#define ISEQ .eqv.
#else
#define ISEQ ==
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     isMatClass_ENABLED && Herm_ENABLED || Symm_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     Herm_ENABLED && CK_ENABLED
#define GET_CONJG(X) conjg(X)
#elif   Symm_ENABLED || Herm_ENABLED
#define GET_CONJG(X) X
#else
#error  "Unrecorgnized interface."
#endif
        integer(IK) :: irow, icol
        itis = logical(size(mat, 1, IK) == size(mat, 2, IK), LK) ! check square shape.
        if (itis) then
            loopOverCol: do icol = 1, size(mat, 2, IK)
                loopOverRow: do irow = 1, size(mat, 1, IK)
                    if (mat(irow, icol) ISEQ GET_CONJG(mat(icol, irow))) cycle loopOverRow
                    itis = .false._LK
                    return
                end do loopOverRow
            end do loopOverCol
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isMatClass_ENABLED && PosDef_ENABLED && Ful_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        itis = logical(isMatClass(mat, hermitian), LK)
        if (itis) itis = isMatClass(mat, posdefmat, uppDia, rdpack)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isMatClass_ENABLED && PosDef_ENABLED && (Upp_ENABLED || Low_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
#if     CK_ENABLED
        complex(CKG) :: chol(size(mat, 1, IK), size(mat, 2, IK))
#elif   RK_ENABLED
        real(RKG) :: chol(size(mat, 1, IK), size(mat, 2, IK))
#else
#error  "Unrecorgnized interface."
#endif
#if     RDP_ENABLED
        itis = logical(size(mat, 1, IK) == size(mat, 2, IK), LK) ! check square shape.
#elif   RFP_ENABLED
        itis = isMatPack(pack, shape(mat, IK))
#else
#error  "Unrecorgnized interface."
#endif
        if (itis) then
            !call setMatCopy(chol, rdpack, mat, pack, subset)
            !call setMatChol(chol, subset, info, recursion)
            call setMatChol(mat, subset, info, chol, nothing)
        end if
        itis = logical(info == 0_IK, LK)
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecorgnized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_CONJG
#undef  ISEQ