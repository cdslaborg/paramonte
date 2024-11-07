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
!>  This file contains implementations of procedures [pm_clusKmeans](@ref pm_clusKmeans).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     setMember_ENABLED && D0_D1_ENABLED && Def_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: temp
        integer(IK) :: icls
        CHECK_ASSERTION(__LINE__, 0_IK < size(center, 1, IK), SK_"@setMember(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        ! First cluster.
        disq = (sample - center(1))**2
        membership = 1
        ! All the rest of clusters.
        do icls = 2, size(center, 1, IK)
            temp = (sample - center(icls))**2
            if (temp < disq) then
                disq = temp
                membership = icls
            end if
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMember_ENABLED && D1_D1_ENABLED && Def_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: temp
        integer(IK) :: icls, isam, nsam
        nsam = size(sample, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK <  size(center, 1, IK), SK_"@setMember(): The condition `0 < size(center, 2)` must hold. size(center, 1) = "//getStr(size(center, 1, IK)))
        CHECK_ASSERTION(__LINE__, nsam == size(disq, 1, IK), SK_"@setMember(): The condition `size(sample, rank(sample)) == size(disq)` must hold. shape(sample), size(disq) = "//getStr([shape(sample, IK), size(disq, 1, IK)]))
        CHECK_ASSERTION(__LINE__, nsam == size(membership, 1, IK), SK_"@setMember(): The condition `size(sample, rank(sample)) == size(membership)` must hold. shape(sample), size(membership) = "//getStr([shape(sample, IK), size(membership, 1, IK)]))
        ! First cluster.
        do isam = 1, nsam
            disq(isam) = (sample(isam) - center(1))**2
            membership(isam) = 1
        end do
        ! All the rest of clusters.
        do icls = 2, size(center, 1, IK)
            do isam = 1, nsam
                temp = (sample(isam) - center(icls))**2
                if (temp < disq(isam)) then
                    disq(isam) = temp
                    membership(isam) = icls
                end if
            end do
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMember_ENABLED && D1_D2_ENABLED && Def_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: temp
        integer(IK) :: icls, idim, ndim
        ndim = size(sample, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK <  size(center, rank(center), IK), SK_"@setMember(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        CHECK_ASSERTION(__LINE__, ndim == size(center, 1, IK), SK_"@setMember(): The condition `size(sample, 1) == size(center, 1)` must hold. size(sample, 1), size(center, 1) = "//getStr([size(sample, 1, IK), size(center, 1, IK)]))
        ! First cluster.
        !disq = sum((sample - center(1 : ndim, 1))**2)
        disq = 0._RKG
        do idim = 1, ndim
            disq = disq + (sample(idim) - center(idim, 1))**2
        end do
        membership = 1
        ! All the rest of clusters.
        do icls = 2, size(center, 2, IK)
           !temp = sum((sample  - center(1 : ndim, icls))**2)
            temp = 0._RKG
            do idim = 1, ndim
                temp = temp + (sample(idim) - center(idim, icls))**2
            end do
            if (temp < disq) then
                disq = temp
                membership = icls
            end if
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMember_ENABLED && D2_D2_ENABLED && Def_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: temp
        integer(IK) :: icls, idim, ndim, isam, nsam
        ndim = size(sample, 1, IK)
        nsam = size(sample, 2, IK)
        CHECK_ASSERTION(__LINE__, 0_IK <  size(center, 2, IK), SK_"@setMember(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        CHECK_ASSERTION(__LINE__, nsam == size(disq, 1, IK), SK_"@setMember(): The condition `size(sample, rank(sample)) == size(disq)` must hold. shape(sample), size(disq) = "//getStr([shape(sample, IK), size(disq, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ndim == size(center, 1, IK), SK_"@setMember(): The condition `size(sample, 1) == size(center, 1)` must hold. size(sample, 1), size(center, 1) = "//getStr([size(sample, 1, IK), size(center, 1, IK)]))
        CHECK_ASSERTION(__LINE__, nsam == size(membership, 1, IK), SK_"@setMember(): The condition `size(sample, rank(sample)) == size(membership, 1)` must hold. shape(sample), size(membership, 1) = "//getStr([shape(sample, IK), size(membership, 1, IK)]))
        ! First cluster.
        icls = 1
        do isam = 1, nsam
            !disq(isam) = sum((center(1 : ndim, icls) - sample(1 : ndim, isam))**2)
            disq(isam) = 0._RKG
            do idim = 1, ndim
                disq(isam) = disq(isam) + (sample(idim, isam) - center(idim, icls))**2
            end do
            membership(isam) = icls
        end do
        ! All the rest of clusters.
        do icls = 2, size(center, 2, IK)
            do isam = 1, nsam
                !temp = sum((center(1 : ndim, icls) - sample(1 : ndim, isam))**2)
                temp = 0._RKG
                do idim = 1, ndim
                    temp = temp + (sample(idim, isam) - center(idim, icls))**2
                end do
                if (temp < disq(isam)) then
                    disq(isam) = temp
                    membership(isam) = icls
                end if
            end do
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMember_ENABLED && D0_D1_ENABLED && Cng_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: temp
        integer(IK) :: icls, mold ! old membership
        CHECK_ASSERTION(__LINE__, 0_IK < size(center, 1, IK), SK_"@setMember(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        mold = membership
        disq = huge(0._RKG)
        do icls = 1, size(center, 1, IK)
            temp = (sample - center(icls))**2
            if (temp < disq) then
                membership = icls
                disq = temp
            end if
        end do
        changed = mold /= membership

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMember_ENABLED && D1_D1_ENABLED && Cng_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: temp
        integer(IK) :: icls, ncls, isam, nsam, mold ! old membership
        nsam = size(sample, 1, IK)
        ncls = size(center, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK <  size(center, rank(center), IK), SK_"@setMember(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        CHECK_ASSERTION(__LINE__, nsam == size(disq, 1, IK), SK_"@setMember(): The condition `size(sample, rank(sample)) == size(disq)` must hold. shape(sample), size(disq) = "//getStr([shape(sample, IK), size(disq, 1, IK)]))
        CHECK_ASSERTION(__LINE__, nsam == size(membership, 1, IK), SK_"@setMember(): The condition `size(sample, rank(sample)) == size(membership)` must hold. shape(sample), size(membership) = "//getStr([shape(sample, IK), size(membership, 1, IK)]))
        changed = .false._LK
        do isam = 1, nsam
            mold = membership(isam)
            disq(isam) = huge(0._RKG)
            do icls = 1, ncls
                temp = (sample(isam) - center(icls))**2
                if (temp < disq(isam)) then
                    membership(isam) = icls
                    disq(isam) = temp
                end if
            end do
            changed = changed .or. membership(isam) /= mold
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMember_ENABLED && D1_D2_ENABLED && Cng_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: temp
        integer(IK) :: icls, idim, ndim, mold ! old membership
        ndim = size(sample, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK <  size(center, rank(center), IK), SK_"@setMember(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        CHECK_ASSERTION(__LINE__, ndim == size(center, 1, IK), SK_"@setMember(): The condition `size(sample, 1) == size(center, 1)` must hold. size(sample, 1), size(center, 1) = "//getStr([size(sample, 1, IK), size(center, 1, IK)]))
        mold = membership
        disq = huge(0._RKG)
        do icls = 1, size(center, 2, IK)
            temp = 0._RKG
            do idim = 1, ndim
                temp = temp + (sample(idim) - center(idim, icls))**2
            end do
            if (temp < disq) then
                membership = icls
                disq = temp
            end if
        end do
        changed = membership /= mold

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMember_ENABLED && D2_D2_ENABLED && Cng_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: temp
        integer(IK) :: icls, ncls, idim, ndim, isam, nsam, mold ! old membership
        ndim = size(sample, 1, IK)
        nsam = size(sample, 2, IK)
        ncls = size(center, 2, IK)
        CHECK_ASSERTION(__LINE__, 0_IK <  size(center, rank(center), IK), SK_"@setMember(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        CHECK_ASSERTION(__LINE__, nsam == size(disq, 1, IK), SK_"@setMember(): The condition `size(sample, rank(sample)) == size(disq)` must hold. shape(sample), size(disq) = "//getStr([shape(sample, IK), size(disq, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ndim == size(center, 1, IK), SK_"@setMember(): The condition `size(sample, 1) == size(center, 1)` must hold. size(sample, 1), size(center, 1) = "//getStr([size(sample, 1, IK), size(center, 1, IK)]))
        CHECK_ASSERTION(__LINE__, nsam == size(membership, 1, IK), SK_"@setMember(): The condition `size(sample, rank(sample)) == size(membership, 1)` must hold. shape(sample), size(membership, 1) = "//getStr([shape(sample, IK), size(membership, 1, IK)]))
        changed = .false._LK
        do isam = 1, nsam
            mold = membership(isam)
            disq(isam) = huge(0._RKG)
            do icls = 1, ncls
                temp = 0._RKG
                do idim = 1, ndim
                    temp = temp + (sample(idim, isam) - center(idim, icls))**2
                end do
                if (temp < disq(isam)) then
                    membership(isam) = icls
                    disq(isam) = temp
                end if
            end do
            changed = changed .or. membership(isam) /= mold
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCenter_ENABLED && D1_D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: icls, ncls, isam, nsam

        nsam = ubound(sample, rank(sample), IK)
        ncls = ubound(center, rank(center), IK)

        CHECK_ASSERTION(__LINE__, all(0._RKG <= disq), SK_"@setCenter(): The condition `all(0 <= disq)` must hold. disq = "//getStr(disq))
        CHECK_ASSERTION(__LINE__, nsam == ubound(disq, 1, IK), SK_"@setCenter(): The condition `size(sample, rank(sample)) == size(disq)` must hold. shape(sample), size(disq) = "//getStr([shape(sample, IK), ubound(disq, 1, IK)]))
       !check_assertion(__LINE__, nsam >= ubound(center, rank(center), IK), SK_"@setCenter(): The condition `size(sample, rank(sample)) >= size(center, rank(center))` must hold. shape(sample), shape(center) = "//getStr([shape(sample, IK), shape(center, IK)]))
        CHECK_ASSERTION(__LINE__, nsam == ubound(membership, 1, IK), SK_"@setCenter(): The condition `size(sample, rank(sample)) == size(membership)` must hold. shape(sample), size(membership) = "//getStr([shape(sample, IK), ubound(membership, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ncls == ubound(potential, 1, IK), SK_"@setCenter(): The condition `size(center, rank(center)) == size(potential)` must hold. shape(center), size(potential) = "//getStr([shape(center, IK), ubound(potential, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ncls == ubound(size, 1, IK), SK_"@setCenter(): The condition `size(center, rank(center)) == size(size)` must hold. shape(center), size(size) = "//getStr([shape(center, IK), ubound(size, 1, IK)]))
        CHECK_ASSERTION(__LINE__, 0_IK <  ncls, SK_"@setCenter(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        CHECK_ASSERTION(__LINE__, all(0_IK < membership .and. membership <= ncls), SK_"@setCenter(): The condition `all(0 < membership .and. membership <= size(center, rank(center)))` must hold. size(center, rank(center)), membership = "//getStr([ncls, membership]))

        ! Initialize.

        do concurrent(icls = 1 : ncls)
            potential(icls) = 0._RKG
            center(icls) = 0._RKG
            size(icls) = 0_IK
        end do

        ! compute new centers.

        do isam = 1, nsam
            center(membership(isam)) = center(membership(isam)) + sample(isam)
            potential(membership(isam)) = potential(membership(isam)) + disq(isam)
            size(membership(isam)) = size(membership(isam)) + 1
        end do

        ! normalize new centers.

        do concurrent(icls = 1 : ncls)
            if (0_IK < size(icls)) center(icls) = center(icls) / real(size(icls), RKG)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCenter_ENABLED && D2_D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: sizeinv
        integer(IK) :: icls, ncls, ndim, isam, nsam

        ndim = ubound(sample, 1, IK)
        nsam = ubound(sample, 2, IK)
        ncls = ubound(center, 2, IK)

        CHECK_ASSERTION(__LINE__, all(0._RKG <= disq), SK_"@setCenter(): The condition `all(0 <= disq)` must hold. disq = "//getStr(disq))
        CHECK_ASSERTION(__LINE__, nsam == ubound(disq, 1, IK), SK_"@setCenter(): The condition `ubound(sample, rank(sample)) == ubound(disq, 1)` must hold. shape(sample), ubound(disq, 1) = "//getStr([shape(sample, IK), ubound(disq, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ndim == ubound(center, 1, IK), SK_"@setCenter(): The condition `ubound(sample, 1) == ubound(center, 1)` must hold. ubound(sample, 1), ubound(center, 1) = "//getStr([ubound(sample, 1, IK), ubound(center, 1, IK)]))
       !check_assertion(__LINE__, nsam >= ubound(center, rank(center), IK), SK_"@setCenter(): The condition `ubound(sample, rank(sample)) >= ubound(center, rank(center))` must hold. shape(sample), shape(center) = "//getStr([shape(sample, IK), shape(center, IK)]))
        CHECK_ASSERTION(__LINE__, nsam == ubound(membership, 1, IK), SK_"@setCenter(): The condition `ubound(sample, rank(sample)) == ubound(membership, 1)` must hold. shape(sample), ubound(membership, 1) = "//getStr([shape(sample, IK), ubound(membership, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ncls == ubound(potential, 1, IK), SK_"@setCenter(): The condition `ubound(center, rank(center)) == ubound(potential, 1)` must hold. shape(center), ubound(potential, 1) = "//getStr([shape(center, IK), ubound(potential, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ncls == ubound(size, 1, IK), SK_"@setCenter(): The condition `ubound(center, rank(center)) == ubound(size, 1)` must hold. shape(center), ubound(size, 1) = "//getStr([shape(center, IK), ubound(size, 1, IK)]))
        CHECK_ASSERTION(__LINE__, 0_IK <  ncls, SK_"@setCenter(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        CHECK_ASSERTION(__LINE__, all(0 < membership .and. membership <= ncls), SK_"@setCenter(): The condition `all(0 < membership .and. membership <= size(center, rank(center)))` must hold. size(center, rank(center)), membership = "//getStr([ubound(center, 2, IK), membership]))

        ! Initialize.

        do concurrent(icls = 1 : ncls)
            center(1 : ndim, icls) = 0._RKG
            potential(icls) = 0._RKG
            size(icls) = 0_IK
        end do

        ! compute new centers.

        do isam = 1, nsam
            center(1 : ndim, membership(isam)) = center(1 : ndim, membership(isam)) + sample(1 : ndim, isam)
            potential(membership(isam)) = potential(membership(isam)) + disq(isam)
            size(membership(isam)) = size(membership(isam)) + 1
        end do

        ! normalize new centers.

        do concurrent(icls = 1 : ncls)
            if (0_IK < size(icls)) then
                sizeinv = 1._RKG / real(size(icls), RKG)
                center(1 : ndim, icls) = center(1 : ndim, icls) * sizeinv
            end if
        end do

        !%%%%%%%%%%%%%%%%%%
#elif   setKmeansPP_ENABLED
        !%%%%%%%%%%%%%%%%%%

        !!!!
        !!!!    WARNING:
        !!!!    You cannot use the intrinsic `size` here, becuase it is an procedure argument.
        !!!!
#define INTRINSIC_ENABLED 0
        real(RKG)   :: temp
        integer(IK) :: idim, icls, isam, ndim, nsam, cid
       !real(RKG)   :: csdisq(0 : ubound(sample, 2, IK)) ! Sum of distance squared of points up to the given index of the vector.
#if     Optional_ENABLED
        integer(IK) :: ncls
        ncls = ubound(center, 2, IK)
        CHECK_ASSERTION(__LINE__, ubound(sample, 1, IK) == ubound(center, 1, IK), SK_"@setKmeansPP(): The condition `ubound(sample, 1) == ubound(center, 1)` must hold. ubound(sample, 1), ubound(center, 1) = "//getStr([ubound(sample, 1, IK), ubound(center, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ncls == ubound(potential, 1, IK), SK_"@setKmeansPP(): The condition `ncls == ubound(potential, 1)` must hold. ncls, ubound(potential, 1) = "//getStr([ncls, ubound(potential, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ncls == ubound(size, 1, IK), SK_"@setKmeansPP(): The condition `ncls == ubound(size, 1)` must hold. ncls, ubound(size, 1) = "//getStr([ncls, ubound(size, 1, IK)]))
#elif   !Default_ENABLED
#error  "Unrecognized interface."
#endif
        ndim = ubound(sample, 1, IK)
        nsam = ubound(sample, 2, IK)
        CHECK_ASSERTION(__LINE__, 0_IK  < ncls, SK_"@setKmeansPP(): The condition `0 < ncls` must hold. ncls = "//getStr(ncls))
        CHECK_ASSERTION(__LINE__, 0_IK  < ubound(sample, 2, IK), SK_"@setKmeansPP(): The condition `0 < ubound(sample, rank(sample))` must hold. shape(sample) = "//getStr(shape(sample, IK)))
        CHECK_ASSERTION(__LINE__, nsam == ubound(disq, 1, IK), SK_"@setKmeansPP(): The condition `ubound(sample, rank(sample)) == ubound(disq, 1)` must hold. shape(sample), size(disq) = "//getStr([shape(sample, IK), ubound(disq, 1, IK)]))
        CHECK_ASSERTION(__LINE__, nsam == ubound(csdisq, 1, IK), SK_"@setKmeansPP(): The condition `ubound(sample, rank(sample)) == ubound(csdisq, 1)` must hold. shape(sample), size(csdisq) = "//getStr([shape(sample, IK), ubound(csdisq, 1, IK)]))
        CHECK_ASSERTION(__LINE__, nsam == ubound(membership, 1, IK), SK_"@setKmeansPP(): The condition `ubound(sample, rank(sample)) == ubound(membership, 1)` must hold. ubound(sample, rank(sample)), ubound(membership, 1) = "//getStr([nsam, ubound(membership, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ncls <= nsam, SK_"@setKmeansPP(): The condition `ncls < ubound(sample, rank(sample))` must hold. ncls, ubound(sample, rank(sample)) = "//getStr([ncls, nsam]))

        csdisq(0) = 0._RKG ! This element must always remain zero.

        ! Define the first cluster center.

        icls = 1
        call setUnifRand(rng, cid, 1_IK, nsam) ! This random assignment is unnecessary if the `sample` is in random order.
#if     Optional_ENABLED
        center(1 : ndim, icls) = sample(1 : ndim, cid)
#endif
        do isam = 1, nsam
            ! Find the distance-squared of each sample to the #icls cluster cid.
#if         INTRINSIC_ENABLED
            disq(isam) = sum((sample(1 : ndim, isam) - sample(1 : ndim, cid))**2)
#else
            disq(isam) = 0._RKG
            do idim = 1, ndim
                disq(isam) = disq(isam) + (sample(idim, isam) - sample(idim, cid))**2
            end do
#endif
            csdisq(isam) = csdisq(isam - 1) + disq(isam)
            membership(isam) = icls
        end do

        ! Find the second cluster center.

        call setUnifRand(rng, temp)
        temp = temp * csdisq(nsam)
#if     INTRINSIC_ENABLED
        cid = minloc(csdisq, 1, temp < csdisq) - 1 ! `-1` compensates for the 0 starting index of `csdisq`.
        !if (cid < 1) cid = nsam + 1
#else
        do cid = 1, nsam; if (temp < csdisq(cid)) exit; end do
#endif
        !print *, "ndim, nsam, ncls, cid, temp", ndim, nsam, ncls, cid, temp, csdisq(1:nsam)

        ! Take care of the unusual case of only one cluster with a single sample.

        if (1_IK == ncls) then
#if         Optional_ENABLED
            potential(1) = sum(disq)
            size(1) = nsam
#endif
            return
        end if

        ! Find the rest of cluster centers.

        do icls = 2, ncls - 1
#if         Optional_ENABLED
            center(1 : ndim, icls) = sample(1 : ndim, cid)
#endif
            !if (nsam < cid) exit ! only if nsam < ncls.
            do isam = 1, nsam
#if             INTRINSIC_ENABLED
                temp = sum((sample(1 : ndim, isam) - sample(1 : ndim, cid))**2)
#else
                temp = 0._RKG
                do idim = 1, ndim
                    temp = temp + (sample(idim, isam) - sample(idim, cid))**2
                end do
#endif
                if (temp < disq(isam)) then
                    disq(isam) = temp
                    membership(isam) = icls
                end if
                csdisq(isam) = csdisq(isam - 1) + disq(isam)
            end do
            call setUnifRand(rng, temp)
            temp = temp * csdisq(nsam)
#if         INTRINSIC_ENABLED
            cid = minloc(csdisq, 1, temp < csdisq) - 1
            !if (cid < 1) cid = nsam + 1
#else
            do cid = 1, nsam; if (temp < csdisq(cid)) exit; end do
#endif
        end do

        ! Initialize the output.

#if     Optional_ENABLED
        do concurrent(icls = 1 : ncls)
            !center(1 : ndim, icls) = 0._RKG
            potential(icls) = 0._RKG
            size(icls) = 0_IK
        end do
        center(1 : ndim, ncls) = sample(1 : ndim, cid)
#endif

        ! Final re-computation of the cluster centers and sizes and sample memberships.

        do isam = 1, nsam
#if         INTRINSIC_ENABLED
            temp = sum((sample(1 : ndim, isam) - sample(1 : ndim, cid))**2)
#else
            temp = 0._RKG; do idim = 1, ndim; temp = temp + (sample(idim, isam) - sample(idim, cid))**2; end do
#endif
            if (temp < disq(isam)) then
                disq(isam) = temp
                membership(isam) = ncls
            end if
#if         Optional_ENABLED
           !center(1 : ndim, membership(isam)) = center(1 : ndim, membership(isam)) + sample(1 : ndim, isam)
            potential(membership(isam)) = potential(membership(isam)) + disq(isam)
            size(membership(isam)) = size(membership(isam)) + 1
#endif
        end do
#undef  INTRINSIC_ENABLED
        ! Normalize sample. This requires the sample size to be larger than the number of clusters to have `all(size > 0)` hold.

        !do concurrent(icls = 1 : ncls)
        !    temp = 1._RKG / real(size(icls), RKG)
        !    center(1 : ndim, icls) = center(1 : ndim, icls) * temp
        !end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setKmeans_ENABLED && Init_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DEFMET_ENABLED 0
        integer(IK) :: retin, def_minsize, def_nitermax
#if     DEFMET_ENABLED
        integer(IK) :: memberswap(ubound(membership, 1, IK))
#endif
        CHECK_ASSERTION(__LINE__, all(0._RKG <= disq), SK_"@setKmeans(): The condition `all(0 <= disq)` must hold. disq = "//getStr(disq))
        CHECK_ASSERTION(__LINE__, ubound(sample, 2, IK) == ubound(disq, 1, IK), SK_"@setKmeans(): The condition `size(sample, rank(sample)) == size(disq)` must hold. shape(sample), size(disq) = "//getStr([shape(sample, IK), ubound(disq, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ubound(sample, 1, IK) == ubound(center, 1, IK), SK_"@setKmeans(): The condition `size(sample, 1) == size(center, 1)` must hold. size(sample, 1), size(center, 1) = "//getStr([ubound(sample, 1, IK), ubound(center, 1, IK)]))
       !check_assertion(__LINE__, ubound(sample, 2, IK) >= ubound(center, rank(center), IK), SK_"@setKmeans(): The condition `size(sample, rank(sample)) >= size(center, rank(center))` must hold. shape(sample), shape(center) = "//getStr([shape(sample, IK), shape(center, IK)]))
        CHECK_ASSERTION(__LINE__, ubound(sample, 2, IK) == ubound(membership, 1, IK), SK_"@setKmeans(): The condition `size(sample, rank(sample)) == size(membership)` must hold. shape(sample), size(membership) = "//getStr([shape(sample, IK), ubound(membership, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ubound(center, 2, IK) == ubound(potential, 1, IK), SK_"@setKmeans(): The condition `size(center, rank(center)) == size(potential)` must hold. shape(center), size(potential) = "//getStr([shape(center, IK), ubound(potential, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ubound(center, 2, IK) == ubound(size, 1, IK), SK_"@setKmeans(): The condition `size(center, rank(center)) == size(size)` must hold. shape(center), size(size) = "//getStr([shape(center, IK), ubound(size, 1, IK)]))
        CHECK_ASSERTION(__LINE__, ubound(center, 2, IK) >  0_IK, SK_"@setKmeans(): The condition `0 < size(center, rank(center))` must hold. shape(center) = "//getStr(shape(center, IK)))
        CHECK_ASSERTION(__LINE__, all(0 < membership .and. membership <= ubound(center, 2, IK)), SK_"@setKmeans(): The condition `all(0 < membership .and. membership <= size(center, rank(center)))` must hold. size(center, rank(center)), membership = "//getStr([ubound(center, 2, IK), membership]))

        if (present(maxniter)) then
            CHECK_ASSERTION(__LINE__, 0_IK <= maxniter, SK_"@setKmeans(): The condition `0 <= maxniter` must hold. maxniter = "//getStr(maxniter))
            def_nitermax = maxniter
        else
            def_nitermax = 300_IK
        end if
        if (present(minsize)) then
            CHECK_ASSERTION(__LINE__, 0_IK <= minsize, SK_"@setKmeans(): The condition `0 <= minsize` must hold. minsize = "//getStr(minsize))
            def_minsize = minsize
        else
            def_minsize = 1_IK
        end if

        retin = 1
        do

            ! compute the new centers, based on the updated memberships.

            call setCenter(membership, disq, sample, center, size, potential)
            failed = any(size < def_minsize)
            if (failed) exit

            ! compute the new memberships, based on the input initial centers.

#if         DEFMET_ENABLED
            call setMember(memberswap, disq, sample, center)
            if (all(memberswap == membership)) exit
#else
            call setMember(membership, disq, sample, center, changed = failed)
            if (.not. failed) exit
#endif
            failed = def_nitermax < retin
            if (failed) exit

            retin = retin + 1

            ! compute the new centers, based on the updated memberships.

#if         DEFMET_ENABLED
            call setCenter(memberswap, disq, sample, center, size, potential)
#else
            call setCenter(membership, disq, sample, center, size, potential)
#endif
            failed = any(size < def_minsize)
            if (failed) exit

            ! compute the new memberships, based on the input initial centers.

#if         DEFMET_ENABLED
            call setMember(membership, disq, sample, center)
            if (all(memberswap == membership)) exit
#else
            call setMember(membership, disq, sample, center, changed = failed)
            if (.not. failed) exit
#endif
            failed = def_nitermax < retin
            if (failed) exit

            retin = retin + 1

        end do
        if (present(niter)) niter = retin
#undef  DEFMET_ENABLED

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setKmeans_ENABLED && (RNGF_ENABLED || RNGX_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: nfail_def, itry, ncls

        if (present(nfail)) then
            CHECK_ASSERTION(__LINE__, 0_IK <= nfail, SK_"@setKmeans(): The condition `0 <= nfail` must hold. nfail = "//getStr(nfail))
            nfail_def = nfail
        else
            nfail_def = 1
        end if

        ncls = ubound(center, 2, IK)
        do itry = 1, nfail_def
            call setKmeansPP(rng, membership, disq, csdisq, sample, ncls)
            call setKmeans(membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize)
            if (.not. failed) exit
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif