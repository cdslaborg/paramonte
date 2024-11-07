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
!>  This module contains implementations of the tests of the procedures under the generic interfaces [pm_arrayFill](@ref pm_arrayFill).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif
#if     SK_ENABLED
        character(2,SKG), allocatable   :: vector(:), matrix(:,:), cuboid(:,:,:)
        character(2,SKG), parameter     :: fill = SKG_"**"
#elif   LK_ENABLED
        logical(LKG)    , allocatable   :: vector(:), matrix(:,:), cuboid(:,:,:)
        logical(LKG)    , parameter     :: fill = .false._LKG
#elif   IK_ENABLED
        integer(IKG)    , allocatable   :: vector(:), matrix(:,:), cuboid(:,:,:)
        integer(IKG)    , parameter     :: fill = huge(1_IKG)
#elif   CK_ENABLED
        complex(CKG)    , allocatable   :: vector(:), matrix(:,:), cuboid(:,:,:)
        complex(CKG)    , parameter     :: fill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
#elif   RK_ENABLED
        real(RKG)       , allocatable   :: vector(:), matrix(:,:), cuboid(:,:,:)
        real(RKG)       , parameter     :: fill = huge(0._RKG)
#else
#error  "Unrecognized interface."
#endif
        type(display_type) :: disp
        character(6) :: objects(3)
        integer(IK) :: iobj, itry, oshape(3)
        objects = ["vector", "matrix", "cuboid"]
        assertion = .true._LK

        do iobj = 1, size(objects)
            do itry = 1, 50
                call setUnifRand(oshape, 0_IK, 3_IK)
                if (objects(iobj) == "vector") then
                    vector = getFilled(fill, oshape(1))
                    assertion = assertion .and. all(shape(vector, IK) == oshape(1:iobj))
                    assertion = assertion .and. all(vector IS_EQUAL fill)
                elseif (objects(iobj) == "matrix") then
                    matrix = getFilled(fill, oshape(1), oshape(2))
                    assertion = assertion .and. all(shape(matrix, IK) == oshape(1:iobj))
                    assertion = assertion .and. all(matrix IS_EQUAL fill)
                elseif (objects(iobj) == "cuboid") then
                    cuboid = getFilled(fill, oshape(1), oshape(2), oshape(3))
                    assertion = assertion .and. all(shape(cuboid, IK) == oshape(1:iobj))
                    assertion = assertion .and. all(cuboid IS_EQUAL fill)
                else
                    error stop "Unrecognized object shape." ! LCOV_EXCL_LINE
                end if
                call report()
                call test%assert(assertion, SK_"The shape of the output must be the specified input shape.", int(__LINE__, IK))
                call test%assert(assertion, SK_"The output must be filled with the specified `fill`.", int(__LINE__, IK))
            end do
        end do

    contains

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip
                call disp%show("oshape(1:iobj)")
                call disp%show( oshape(1:iobj) )
                if (objects(iobj) == "vector") then
                    call disp%show("vector")
                    call disp%show( vector )
                elseif (objects(iobj) == "matrix") then
                    call disp%show("matrix")
                    call disp%show( matrix )
                elseif (objects(iobj) == "cuboid") then
                    call disp%show("cuboid")
                    call disp%show( cuboid )
                end if
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#undef IS_EQUAL
#undef GET_SIZE
#undef ALL