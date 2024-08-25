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

!>  \brief This file contains the implementations of the tests of module [pm_distanceEuclid](@ref pm_distanceEuclid).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getDisEuclid_ENABLED || setDisEuclid_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK), parameter :: ntry = 100
        real(TKG), parameter :: EPS = epsilon(0._TKG) * 100
        real(TKG), allocatable :: distance(:,:), distance_ref(:,:)
        real(TKG), allocatable :: point(:,:), ref(:,:)
        real(TKG), allocatable :: diff(:,:)
        type(csp_type), allocatable :: method(:)
        type(css_type), allocatable :: mnames(:)
        integer(IK) :: ndim, npnt, nref
        integer(IK) :: itry, imethod
        logical(LK) :: isorigin

        assertion = .true._LK

        method = [csp_type(euclid), csp_type(euclidu), csp_type(euclidsq)]
        mnames = [css_type("euclid"), css_type("euclidu"), css_type("euclidsq")]

        do itry = 1, ntry

            isorigin = getUnifRand()
            ndim = getUnifRand(1_IK, 10_IK)
            npnt = getUnifRand(1_IK, 20_IK)
            nref = getUnifRand(1_IK, 20_IK)

            diff = getFilled(0._TKG, npnt, nref)
            if (isorigin) then
                ref = getFilled(0._TKG, ndim, nref)
            else
                ref = getUnifRand(-1._TKG, 1._TKG, ndim, nref)
            end if
            point = getUnifRand(-1._TKG, 1._TKG, ndim, npnt)
            call setResized(distance, [npnt, nref])

            do imethod = 1, size(method)
                distance_ref = getDisEuclid_ref(point, ref, method(imethod)%val)
#if             getDisEuclid_ENABLED
                if (isorigin) then
                    ! D1_XX
                    block
                        integer(IK) :: ipnt, iref
                        ipnt = 1_IK; iref = 1_IK
                        distance(ipnt, iref) = getDisEuclid(point(:, ipnt), method(imethod)%val)
                        diff(ipnt, iref) = abs(distance(ipnt, iref) - distance_ref(ipnt, iref))
                        assertion = assertion .and. diff(ipnt, iref) < EPS
                        call report(__LINE__, mnames(imethod)%val, point(:, ipnt), ref(:, iref), distance(ipnt, iref), distance_ref(ipnt, iref), diff(ipnt, iref))
                    end block
                    ! D2_XX
                    block
                        integer(IK) :: iref
                        iref = 1_IK
                        distance(:, iref) = getDisEuclid(point(:,:), method(imethod)%val)
                        diff(:, iref) = abs(distance(:, iref) - distance_ref(:, iref))
                        assertion = assertion .and. all(diff(:, iref) < EPS)
                        call report(__LINE__, mnames(imethod)%val, point(:,:), ref(:, iref), distance(:, iref), distance_ref(:, iref), diff(:, iref))
                    end block
                else
                    ! D1_D1
                    block
                        integer(IK) :: ipnt, iref
                        ipnt = 1_IK; iref = 1_IK
                        distance(ipnt, iref) = getDisEuclid(point(:, ipnt), ref(:, iref), method(imethod)%val)
                        diff(ipnt, iref) = abs(distance(ipnt, iref) - distance_ref(ipnt, iref))
                        assertion = assertion .and. diff(ipnt, iref) < EPS
                        call report(__LINE__, mnames(imethod)%val, point(:, ipnt), ref(:, iref), distance(ipnt, iref), distance_ref(ipnt, iref), diff(ipnt, iref))
                    end block
                    ! D1_D2
                    block
                        integer(IK), allocatable :: ipnt
                        ipnt = 1_IK
                        distance(ipnt, :) = getDisEuclid(point(:, ipnt), ref(:,:), method(imethod)%val)
                        diff(ipnt, :) = abs(distance(ipnt, :) - distance_ref(ipnt, :))
                        assertion = assertion .and. all(diff(ipnt, :) < EPS)
                        call report(__LINE__, mnames(imethod)%val, point(:, ipnt), ref(:,:), distance(ipnt, :), distance_ref(ipnt, :), diff(ipnt, :))
                    end block
                    ! D2_D1
                    block
                        integer(IK) :: iref
                        iref = 1_IK
                        distance(:, iref) = getDisEuclid(point(:,:), ref(:, iref), method(imethod)%val)
                        diff(:, iref) = abs(distance(:, iref) - distance_ref(:, iref))
                        assertion = assertion .and. all(diff(:, iref) < EPS)
                        call report(__LINE__, mnames(imethod)%val, point(:,:), ref(:, iref), distance(:, iref), distance_ref(:, iref), diff(:, iref))
                    end block
                    ! D2_D2
                    block
block
use pm_io
call disp%show([shape(distance), shape(getDisEuclid(point(:,:), ref(:,:), method(imethod)%val)), shape(point), shape(ref)])
end block
                        distance(:,:) = getDisEuclid(point(:,:), ref(:,:), method(imethod)%val)
                        diff(:,:) = abs(distance(:,:) - distance_ref(:,:))
                        assertion = assertion .and. all(diff(:,:) < EPS)
                        call report(__LINE__, mnames(imethod)%val, point(:,:), ref(:,:), distance(:,:), distance_ref(:,:), diff(:,:))
                    end block
                end if
#elif           setDisEuclid_ENABLED
                if (isorigin) then
                    ! D1_XX
                    block
                        integer(IK) :: ipnt, iref
                        ipnt = 1_IK; iref = 1_IK
                        if (same_type_as(method(imethod)%val, euclid)) then
                            call setDisEuclid(distance(ipnt, iref), point(:, ipnt), euclid)
                        elseif (same_type_as(method(imethod)%val, euclidu)) then
                            call setDisEuclid(distance(ipnt, iref), point(:, ipnt), euclidu)
                        elseif (same_type_as(method(imethod)%val, euclidsq)) then
                            call setDisEuclid(distance(ipnt, iref), point(:, ipnt), euclidsq)
                        else
                            error stop "@test_pm_distanceEuclid@test_setDisEuclid(): Internal library error occurred. Unrecognized `method`." ! LCOV_EXCL_LINE
                        end if
                        diff(ipnt, iref) = abs(distance(ipnt, iref) - distance_ref(ipnt, iref))
                        assertion = assertion .and. diff(ipnt, iref) < EPS
                        call report(__LINE__, mnames(imethod)%val, point(:, ipnt), ref(:, iref), distance(ipnt, iref), distance_ref(ipnt, iref), diff(ipnt, iref))
                    end block
                    ! D2_XX
                    block
                        integer(IK) :: iref
                        iref = 1_IK
                        if (same_type_as(method(imethod)%val, euclid)) then
                            call setDisEuclid(distance(:, iref), point(:,:), euclid)
                        elseif (same_type_as(method(imethod)%val, euclidu)) then
                            call setDisEuclid(distance(:, iref), point(:,:), euclidu)
                        elseif (same_type_as(method(imethod)%val, euclidsq)) then
                            call setDisEuclid(distance(:, iref), point(:,:), euclidsq)
                        else
                            error stop "@test_pm_distanceEuclid@test_setDisEuclid(): Internal library error occurred. Unrecognized `method`." ! LCOV_EXCL_LINE
                        end if
                        diff(:, iref) = abs(distance(:, iref) - distance_ref(:, iref))
                        assertion = assertion .and. all(diff(:, iref) < EPS)
                        call report(__LINE__, mnames(imethod)%val, point(:,:), ref(:, iref), distance(:, iref), distance_ref(:, iref), diff(:, iref))
                    end block
                else
                    ! D1_D1
                    block
                        integer(IK) :: ipnt, iref
                        ipnt = 1_IK; iref = 1_IK
                        if (same_type_as(method(imethod)%val, euclid)) then
                            call setDisEuclid(distance(ipnt, iref), point(:, ipnt), ref(:, iref), euclid)
                        elseif (same_type_as(method(imethod)%val, euclidu)) then
                            call setDisEuclid(distance(ipnt, iref), point(:, ipnt), ref(:, iref), euclidu)
                        elseif (same_type_as(method(imethod)%val, euclidsq)) then
                            call setDisEuclid(distance(ipnt, iref), point(:, ipnt), ref(:, iref), euclidsq)
                        else
                            error stop "@test_pm_distanceEuclid@test_setDisEuclid(): Internal library error occurred. Unrecognized `method`." ! LCOV_EXCL_LINE
                        end if
                        diff(ipnt, iref) = abs(distance(ipnt, iref) - distance_ref(ipnt, iref))
                        assertion = assertion .and. diff(ipnt, iref) < EPS
                        call report(__LINE__, mnames(imethod)%val, point(:, ipnt), ref(:, iref), distance(ipnt, iref), distance_ref(ipnt, iref), diff(ipnt, iref))
                    end block
                    ! D1_D2
                    block
                        integer(IK), allocatable :: ipnt
                        ipnt = 1_IK
                        if (same_type_as(method(imethod)%val, euclid)) then
                            call setDisEuclid(distance(ipnt, :), point(:, ipnt), ref(:,:), euclid)
                        elseif (same_type_as(method(imethod)%val, euclidu)) then
                            call setDisEuclid(distance(ipnt, :), point(:, ipnt), ref(:,:), euclidu)
                        elseif (same_type_as(method(imethod)%val, euclidsq)) then
                            call setDisEuclid(distance(ipnt, :), point(:, ipnt), ref(:,:), euclidsq)
                        else
                            error stop "@test_pm_distanceEuclid@test_setDisEuclid(): Internal library error occurred. Unrecognized `method`." ! LCOV_EXCL_LINE
                        end if
                        diff(ipnt, :) = abs(distance(ipnt, :) - distance_ref(ipnt, :))
                        assertion = assertion .and. all(diff(ipnt, :) < EPS)
                        call report(__LINE__, mnames(imethod)%val, point(:, ipnt), ref(:,:), distance(ipnt, :), distance_ref(ipnt, :), diff(ipnt, :))
                    end block
                    ! D2_D1
                    block
                        integer(IK) :: iref
                        iref = 1_IK
                        if (same_type_as(method(imethod)%val, euclid)) then
                            call setDisEuclid(distance(:, iref), point(:,:), ref(:, iref), euclid)
                        elseif (same_type_as(method(imethod)%val, euclidu)) then
                            call setDisEuclid(distance(:, iref), point(:,:), ref(:, iref), euclidu)
                        elseif (same_type_as(method(imethod)%val, euclidsq)) then
                            call setDisEuclid(distance(:, iref), point(:,:), ref(:, iref), euclidsq)
                        else
                            error stop "@test_pm_distanceEuclid@test_setDisEuclid(): Internal library error occurred. Unrecognized `method`." ! LCOV_EXCL_LINE
                        end if
                        diff(:, iref) = abs(distance(:, iref) - distance_ref(:, iref))
                        assertion = assertion .and. all(diff(:, iref) < EPS)
                        call report(__LINE__, mnames(imethod)%val, point(:,:), ref(:, iref), distance(:, iref), distance_ref(:, iref), diff(:, iref))
                    end block
                    ! D2_D2
                    block
                        if (same_type_as(method(imethod)%val, euclid)) then
                            call setDisEuclid(distance(:,:), point(:,:), ref(:,:), euclid)
                        elseif (same_type_as(method(imethod)%val, euclidu)) then
                            call setDisEuclid(distance(:,:), point(:,:), ref(:,:), euclidu)
                        elseif (same_type_as(method(imethod)%val, euclidsq)) then
                            call setDisEuclid(distance(:,:), point(:,:), ref(:,:), euclidsq)
                        else
                            error stop "@test_pm_distanceEuclid@test_setDisEuclid(): Internal library error occurred. Unrecognized `method`." ! LCOV_EXCL_LINE
                        end if
                        diff(:,:) = abs(distance(:,:) - distance_ref(:,:))
                        assertion = assertion .and. all(diff(:,:) < EPS)
                        call report(__LINE__, mnames(imethod)%val, point(:,:), ref(:,:), distance(:,:), distance_ref(:,:), diff(:,:))
                    end block
                end if
#endif
            end do

        end do

    contains

        pure function getDisEuclid_ref(point, ref, method) result(distance_ref)
            real(TKG), intent(in), contiguous :: point(:,:), ref(:,:)
            class(*), intent(in) :: method
            real(TKG) :: distance_ref(size(point, 2, IK), size(ref, 2, IK))
            integer(IK) :: ndim, npnt, nref
            integer(IK) :: ipnt, iref
            ndim = size(point, 1, IK)
            npnt = size(point, 2, IK)
            nref = size(ref, 2, IK)
            do ipnt = 1, npnt
                do iref = 1, nref
                    distance_ref(ipnt, iref) = sum((point(1:ndim, ipnt) - ref(1:ndim, iref))**2)
                end do
            end do
            if (same_type_as(method, euclidsq)) return
            distance_ref = sqrt(distance_ref)
        end function

        subroutine report(line, method, point, ref, distance, distance_ref, diff)
            real(TKG), intent(in) :: point(..), ref(..), distance(..), distance_ref(..), diff(..)
            integer, intent(in) :: line
            character(*, SK) :: method
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[rank(point), rank(ref), rank(diff)]")
                call test%disp%show( [rank(point), rank(ref), rank(diff)] )
                call test%disp%show("[ndim, npnt, nref]")
                call test%disp%show( [ndim, npnt, nref] )
                call test%disp%show("method")
                call test%disp%show( method )
                call display(point, "point")
                call display(ref, "ref")
                call display(distance_ref, "distance_ref")
                call display(distance, "distance")
                call display(diff, "diff")
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The procedure setDisEuclid() must correctly correctly compute the distance.", int(line, IK))
        end subroutine

        ! LCOV_EXCL_START
        subroutine display(object, name)
            real(TKG), intent(in) :: object(..)
            character(*, SK), intent(in) :: name
            select rank(object)
            rank(0)
                call test%disp%show(name)
                call test%disp%show(object)
            rank(1)
                call test%disp%show(name)
                call test%disp%show(object)
            rank(2)
                call test%disp%show(name)
                call test%disp%show(object)
            rank(*)
                error stop "Unrecognized rank for `object`."
            end select
        end subroutine
        ! LCOV_EXCL_STOP

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisMatEuclid_ENABLED || setDisMatEuclid_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK), parameter :: ntry = 100
        real(TKG), parameter :: EPS = epsilon(0._TKG) * 100
        real(TKG), allocatable :: distance(:,:), distance_ref(:,:), point(:,:), diff(:,:)
        type(csp_type), allocatable :: method(:), subset(:), pack(:)
        type(css_type), allocatable :: mnames(:), snames(:), pnames(:)
        integer(IK) :: itry, ipack, isubset, imethod
        integer(IK) :: ndim, npnt

        assertion = .true._LK

        pack = [csp_type(rdpack)]
        pnames = [css_type("rdpack")]

        subset = [csp_type(uppLow), csp_type(uppLowDia)]
        snames = [css_type("uppLow"), css_type("uppLowDia")]

        method = [csp_type(euclid), csp_type(euclidu), csp_type(euclidsq)]
        mnames = [css_type("euclid"), css_type("euclidu"), css_type("euclidsq")]

        do itry = 1, ntry

            ndim = getUnifRand(1_IK, 10_IK)
            npnt = getUnifRand(1_IK, 20_IK)
            point = getUnifRand(-1._TKG, 1._TKG, ndim, npnt)
            do ipack = 1, size(pack)
                do isubset = 1, size(subset)
                    do imethod = 1, size(method)
                        distance_ref = getDisMatEuclid_ref(point, pack(ipack)%val, subset(isubset)%val, method(imethod)%val)
#if                     getDisMatEuclid_ENABLED
                        block
                            if (same_type_as(pack(ipack)%val, rdpack)) then
                                if (same_type_as(subset(isubset)%val, uppLowDia)) then
                                    distance = getDisMatEuclid(rdpack, uppLowDia, point, method(imethod)%val)
                                elseif (same_type_as(subset(isubset)%val, uppLow)) then
                                    distance = getDisMatEuclid(rdpack, uppLow, point, method(imethod)%val)
                                else
                                    error stop "Unrecognized `pack` value." ! LCOV_EXCL_LINE
                                end if
                            else
                                error stop "Unrecognized `pack` value." ! LCOV_EXCL_LINE
                            end if
                            diff = abs(distance - distance_ref)
                            assertion = assertion .and. all(diff < EPS)
                            call report(__LINE__, pnames(ipack)%val, snames(isubset)%val, mnames(imethod)%val, point, distance, distance_ref, diff)
                        end block
#elif                   setDisMatEuclid_ENABLED
                        block
                            if (same_type_as(pack(ipack)%val, rdpack)) then
                                if (same_type_as(subset(isubset)%val, uppLowDia)) then
                                    call setResized(distance, [npnt, npnt])
                                    if (same_type_as(method(imethod)%val, euclid)) then
                                        call setDisMatEuclid(distance, rdpack, uppLowDia, point, euclid)
                                    elseif (same_type_as(method(imethod)%val, euclidu)) then
                                        call setDisMatEuclid(distance, rdpack, uppLowDia, point, euclidu)
                                    elseif (same_type_as(method(imethod)%val, euclidsq)) then
                                        call setDisMatEuclid(distance, rdpack, uppLowDia, point, euclidsq)
                                    else
                                        error stop "Unrecognized `method` value." ! LCOV_EXCL_LINE
                                    end if
                                elseif (same_type_as(subset(isubset)%val, uppLow)) then
                                    call setResized(distance, [npnt - 1, npnt])
                                    if (same_type_as(method(imethod)%val, euclid)) then
                                        call setDisMatEuclid(distance, rdpack, uppLow, point, euclid)
                                    elseif (same_type_as(method(imethod)%val, euclidu)) then
                                        call setDisMatEuclid(distance, rdpack, uppLow, point, euclidu)
                                    elseif (same_type_as(method(imethod)%val, euclidsq)) then
                                        call setDisMatEuclid(distance, rdpack, uppLow, point, euclidsq)
                                    else
                                        error stop "Unrecognized `method` value." ! LCOV_EXCL_LINE
                                    end if
                                else
                                    error stop "Unrecognized `pack` value." ! LCOV_EXCL_LINE
                                end if
                            else
                                error stop "Unrecognized `pack` value." ! LCOV_EXCL_LINE
                            end if
                            diff = abs(distance - distance_ref)
                            assertion = assertion .and. all(diff < EPS)
                            call report(__LINE__, pnames(ipack)%val, snames(isubset)%val, mnames(imethod)%val, point, distance, distance_ref, diff)
                        end block
#endif
                    end do
                end do
            end do

        end do

    contains

        function getDisMatEuclid_ref(point, pack, subset, method) result(distance_ref)
            ! \bug
            ! Intel ifort 2021 passes incorrect size of [0, 0] for distance inside setDisEuclid, call from within getDisEuclid(), called below.
            ! This apparently happens if the contiguous attribute of the `point` argument is missing.
            real(TKG), intent(in), contiguous :: point(:,:)
            class(*), intent(in) :: pack, subset, method
            real(TKG), allocatable :: distance_ref(:,:)
            integer(IK) :: ndim, npnt, ipnt
            ndim = size(point, 1, IK)
            npnt = size(point, 2, IK)
            select type(subset)
            type is (uppLowDia_type)
                if (same_type_as(pack, rdpack)) then
                    call setResized(distance_ref, [npnt, npnt])
                    distance_ref = getDisEuclid(point, point, method)
                else
                    error stop "Unrecognized `pack` value." ! LCOV_EXCL_LINE
                end if
            type is (uppLow_type)
                call setResized(distance_ref, [npnt - 1, npnt])
                if (same_type_as(pack, rdpack)) then
                    do ipnt = 1, npnt - 1
                        ! \bug
                        ! Intel ifort 2024 and older pass a zero length for distance that could not be resolved in any way.
                        ! For now, setDisEuclid() is used as a substitute.
                        !distance_ref(ipnt : npnt - 1, ipnt) = getDisEuclid(point(:, ipnt), point(:, ipnt + 1 : npnt), method)
                        if (same_type_as(method, euclid)) then
                            call setDisEuclid(distance_ref(ipnt : npnt - 1, ipnt), point(:, ipnt), point(:, ipnt + 1 : npnt), euclid)
                        elseif (same_type_as(method, euclidu)) then
                            call setDisEuclid(distance_ref(ipnt : npnt - 1, ipnt), point(:, ipnt), point(:, ipnt + 1 : npnt), euclidu)
                        elseif (same_type_as(method, euclidsq)) then
                            call setDisEuclid(distance_ref(ipnt : npnt - 1, ipnt), point(:, ipnt), point(:, ipnt + 1 : npnt), euclidsq)
                        end if
                        distance_ref(ipnt, ipnt + 1 : npnt) = distance_ref(ipnt : npnt - 1, ipnt)
                    end do
                else
                    error stop "Unrecognized `pack` value." ! LCOV_EXCL_LINE
                end if
            class default
                error stop "Unrecognized `subset` value." ! LCOV_EXCL_LINE
            end select
        end function

        subroutine report(line, pack, subset, method, point, distance, distance_ref, diff)
            real(TKG), intent(in) :: point(..), distance(..), distance_ref(..), diff(..)
            integer, intent(in) :: line
            character(*, SK) :: pack, subset, method
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[ndim, npnt]")
                call test%disp%show( [ndim, npnt] )
                call test%disp%show("pack")
                call test%disp%show( pack )
                call test%disp%show("subset")
                call test%disp%show( subset )
                call test%disp%show("method")
                call test%disp%show( method )
                call display(point, "point")
                call display(distance_ref, "distance_ref")
                call display(distance, "distance")
                call display(diff, "diff")
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The procedure setDisMatEuclid() must correctly correctly compute the distance.", int(line, IK))
        end subroutine

        ! LCOV_EXCL_START
        subroutine display(object, name)
            real(TKG), intent(in) :: object(..)
            character(*, SK), intent(in) :: name
            call test%disp%show(SK_"shape("//name//SK_")")
            call test%disp%show(shape(object))
            select rank(object)
            rank(0)
                call test%disp%show(name)
                call test%disp%show(object)
            rank(1)
                call test%disp%show(name)
                call test%disp%show(object)
            rank(2)
                call test%disp%show(name)
                call test%disp%show(object)
            rank(*)
                error stop "Unrecognized rank for `object`."
            end select
        end subroutine
        ! LCOV_EXCL_STOP

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
