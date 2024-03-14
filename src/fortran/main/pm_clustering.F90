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
!>  This module contains procedures and routines for the computing the Kmeans clustering of a given set of data.
!>
!>  \details
!>
!>  \see
!>  
!>
!>  \test
!>  [test_pm_clustering](@ref test_pm_clustering)<br>
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, April 03, 2017, 2:16 PM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_clustering

    use pm_kind, only: SK, IK, LK
    use pm_distUnif, only: rngf, rngf_type, xoshiro256ssw_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_clustering"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Compute and return the memberships and minimum distances of a set of input points with respect to the an input set of cluster centers.<br>
    !>
    !>  \details
    !>  The `membership` ID of the `i`th sample `sample(:,i)` is determined by computing the distance of the sample from each cluster
    !>  and choosing the cluster ID `j` whose `center(:,j)` has the minimum distance from the sample among all clusters.<br>
    !>  This minimum **squared** distance is output in `disq(:,i)`.<br>
    !>  The metric used within this generic interface is the [Euclidean distance](@ref pm_distanceEuclid).<br>
    !>
    !>  \param[out] membership  :   The output scalar or vector of shape `(1:nsam)` of type `integer` of default kind \IK,
    !>                              containing the membership of each input sample in `sample` from its nearest cluster `center`,
    !>                              such that `cluster(membership(i))` is the nearest cluster center to the `i`th sample `sample(:, i)` at a squared-distance of `disq(i)`.<br>
    !>  \param[out] disq        :   The output scalar or vector of shape `(1:nsam)` of the same type and kind as the input argument `sample`,
    !>                              containing the Euclidean **squared** distance of each input sample in `sample` from its nearest cluster `center`.<br>
    !>  \param[in]  sample      :   The input scalar, vector, or matrix of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sample of `nsam` points in an `ndim`-dimensional space whose
    !>                              memberships and minimum distances with respect to the input `center`s must be computed.<br>
    !>                              <ol>
    !>                                  <li>    If `sample` is a **scalar** and `center` is a **vector** of shape `(1 : ncls)`,
    !>                                          then the input `sample` must be the coordinate of a **single** sample in (univariate space) whose distance from `ncls` cluster `center`s must be computed.<br>
    !>                                  <li>    If `sample` is a **vector** of shape `(1 : ndim)` and `center` is a **matrix** of shape `(1 : ndim, 1 : ncls)`,
    !>                                          then the input `sample` must be a **single** sample (in `ndim`-dimensional space) whose distance from `ncls` cluster `center`s must be computed.<br>
    !>                                  <li>    If `sample` is a **vector** of shape `(1 : nsam)` and `center` is a **vector** of shape `(1 : ncls)`,
    !>                                          then the input `sample` must be a **collection** of `nsam` points (in univariate space) whose distances from `ncls` cluster `center`s must be computed.<br>
    !>                                  <li>    If `sample` is a **matrix** of shape `(1 : ndim, 1 : nsam)` and `center` is a **matrix** of shape `(1 : ndim, 1 : ncls)`,
    !>                                          then the input `sample` must be a **collection** of `nsam` points (in `ndim`-dimensional space) whose distances from `ncls` cluster `center`s must be computed.<br>
    !>                              </ol>
    !>  \param[in]  center      :   The input vector of shape `(1:ncls)` or matrix of shape `(1 : ndim, 1 : ncls)` of the same type and kind as the input argument `sample`,
    !>                              containing the set of `ncls` cluster centers (**centroids**) with respect to which the sample memberships and minimum distances must be computed.<br>
    !>
    !>  \interface{setMember}
    !>  \code{.F90}
    !>
    !>      use pm_clustering, only: setMember
    !>
    !>      call setMember(membership           , disq          , sample                    , center(1 : ncls))
    !>      call setMember(membership(1 : nsam) , disq(1 : nsam), sample(1 : nsam)          , center(1 : ncls))
    !>      call setMember(membership           , disq          , sample(1 : ndim)          , center(1 : ndim, 1 : ncls))
    !>      call setMember(membership(1 : nsam) , disq(1 : nsam), sample(1 : ndim, 1 : nsam), center(1 : ndim, 1 : ncls))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `ubound(center, rank(center)) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, rank(sample)) == ubound(disq, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, rank(sample)) == ubound(membership, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, 1) == ubound(center, 1) .or. rank(sample) == 0 .or. rank(center) == 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  This generic interface implements one of the two major steps involved in Kmeans clustering.<br>
    !>
    !>  \see
    !>  [setKmeans](@ref pm_clustering::setKmeans)<br>
    !>  [setCenter](@ref pm_clustering::setCenter)<br>
    !>  [setMember](@ref pm_clustering::setMember)<br>
    !>  [setKmeansPP](@ref pm_clustering::setKmeansPP)<br>
    !>
    !>  \example{setMember}
    !>  \include{lineno} example/pm_clustering/setMember/main.F90
    !>  \compilef{setMember}
    !>  \output{setMember}
    !>  \include{lineno} example/pm_clustering/setMember/main.out.F90
    !>  \postproc{setMember}
    !>  \include{lineno} example/pm_clustering/setMember/main.py
    !>  \vis{setMember}
    !>  \image html pm_clustering/setMember/setMember.png width=700
    !>
    !>  \test
    !>  [test_pm_clustering](@ref test_pm_clustering)
    !>
    !>  \finmain{setMember}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin
    interface setMember

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMemberEuc_D0_D1_RK5(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D0_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)                    :: sample
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMemberEuc_D0_D1_RK4(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D0_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)                    :: sample
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMemberEuc_D0_D1_RK3(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D0_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)                    :: sample
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMemberEuc_D0_D1_RK2(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D0_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)                    :: sample
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMemberEuc_D0_D1_RK1(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D0_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)                    :: sample
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMemberEuc_D1_D1_RK5(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMemberEuc_D1_D1_RK4(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMemberEuc_D1_D1_RK3(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMemberEuc_D1_D1_RK2(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMemberEuc_D1_D1_RK1(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMemberEuc_D1_D2_RK5(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMemberEuc_D1_D2_RK4(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMemberEuc_D1_D2_RK3(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMemberEuc_D1_D2_RK2(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMemberEuc_D1_D2_RK1(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D1_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)                   :: disq
        integer(IK)             , intent(out)                   :: membership
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMemberEuc_D2_D2_RK5(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D2_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMemberEuc_D2_D2_RK4(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D2_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMemberEuc_D2_D2_RK3(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D2_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMemberEuc_D2_D2_RK2(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D2_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMemberEuc_D2_D2_RK1(membership, disq, sample, center)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMemberEuc_D2_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Compute and return the centers of the clusters corresponding to the input `sample`, cluster `membership` IDs,
    !>  and `sample` distances-squared from their corresponding cluster centers.<br>
    !>
    !>  \details
    !>  This generic interface is the second step in the iterative process of refining a set of
    !>  initial cluster centers toward a final set of clusters and their members.<br>
    !>  As such, this generic interface is of little use with invoking [setMember](@ref pm_clustering::setMember) beforehand.<br>
    !>  The metric used within this generic interface is the [Euclidean distance](@ref pm_distanceEuclid).<br>
    !>
    !>  \param[in]  membership  :   The input vector of shape `(1:nsam)` of type `integer` of default kind \IK,
    !>                              containing the membership of each input sample in `sample` from its nearest cluster `center`,
    !>                              such that `cluster(membership(i))` is the nearest cluster center to the `i`th sample `sample(:, i)` at a squared-distance of `disq(i)`.<br>
    !>  \param[in]  disq        :   The input vector of shape `(1:nsam)` of the same type and kind as the input argument `sample`,
    !>                              containing the Euclidean **squared** distance of each input sample in `sample` from its nearest cluster `center`.<br>
    !>  \param[in]  sample      :   The input vector, or matrix of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sample of `nsam` points in an `ndim`-dimensional space whose corresponding cluster centers must be computed.<br>
    !>                              <ol>
    !>                                  <li>    If `sample` is a **vector** of shape `(1 : nsam)` and `center` is a **vector** of shape `(1 : ncls)`,
    !>                                          then the input `sample` must be a **collection** of `nsam` points (in univariate space).<br>
    !>                                  <li>    If `sample` is a **matrix** of shape `(1 : ndim, 1 : nsam)` and `center` is a **matrix** of shape `(1 : ndim, 1 : ncls)`,
    !>                                          then the input `sample` must be a **collection** of `nsam` points (in `ndim`-dimensional space).<br>
    !>                              </ol>
    !>  \param[out] center      :   The output vector of shape `(1:ncls)` or matrix of shape `(1 : ndim, 1 : ncls)` of the same type and kind as the input argument `sample`,
    !>                              containing the set of `ncls` cluster centers (**centroids**) computed based on the input sample memberships and minimum distances.<br>
    !>  \param[out] size        :   The output vector of shape `(1:ncls)` type `integer` of default kind \IK,
    !>                              containing the sizes (number of members) of the clusters with the corresponding centers output in the argument `center`.<br>
    !>  \param[out] potential   :   The output vector of shape `(1:ncls)` of the same type and kind as the input argument `sample`,
    !>                              the `i`th element of which contains the sum of squared distances of all members of the `i`th cluster from the cluster center as output in the `i`th element of `center`.<br>
    !>
    !>  \interface{setCenter}
    !>  \code{.F90}
    !>
    !>      use pm_clustering, only: setCenter
    !>
    !>      call setCenter(sample(1 : nsam)             , membership(1 : nsam)  , disq(1 : nsam), center(1 : ncls)          , size(1 : ncls), potential(1 : ncls))
    !>      call setCenter(sample(1 : ndim, 1 : nsam)   , membership(1 : nsam)  , disq(1 : nsam), center(1 : ndim, 1 : ncls), size(1 : ncls), potential(1 : ncls))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `ubound(center, rank(center)) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, rank(sample)) == ubound(disq, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(center, rank(center)) == ubound(size, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(center, rank(center)) == ubound(potential, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, rank(sample)) == ubound(membership, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, 1) == ubound(center, 1) .or. (rank(sample) == 1 .and. rank(center) == 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 < membership .and. membership <= ubound(center, rank(center)))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= disq)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  This generic interface implements one of the two major steps involved in Kmeans clustering.<br>
    !>
    !>  \see
    !>  [setKmeans](@ref pm_clustering::setKmeans)<br>
    !>  [setCenter](@ref pm_clustering::setCenter)<br>
    !>  [setMember](@ref pm_clustering::setMember)<br>
    !>  [setKmeansPP](@ref pm_clustering::setKmeansPP)<br>
    !>
    !>  \example{setCenter}
    !>  \include{lineno} example/pm_clustering/setCenter/main.F90
    !>  \compilef{setCenter}
    !>  \output{setCenter}
    !>  \include{lineno} example/pm_clustering/setCenter/main.out.F90
    !>  \postproc{setCenter}
    !>  \include{lineno} example/pm_clustering/setCenter/main.py
    !>  \vis{setCenter}
    !>  \image html pm_clustering/setCenter/setMember.png width=700
    !>  \image html pm_clustering/setCenter/setCenter.png width=700
    !>
    !>  \test
    !>  [test_pm_clustering](@ref test_pm_clustering)
    !>
    !>  \finmain{setCenter}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin
    interface setCenter

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCenterEuc_D1_D1_RK5(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCenterEuc_D1_D1_RK4(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCenterEuc_D1_D1_RK3(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCenterEuc_D1_D1_RK2(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCenterEuc_D1_D1_RK1(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , contiguous    :: sample(:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCenterEuc_D2_D2_RK5(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D2_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCenterEuc_D2_D2_RK4(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D2_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCenterEuc_D2_D2_RK3(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D2_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCenterEuc_D2_D2_RK2(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D2_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCenterEuc_D2_D2_RK1(membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCenterEuc_D2_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(in)    , contiguous    :: disq(:)
        integer(IK)             , intent(in)    , contiguous    :: membership(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Compute and return an asymptotically optimal set of cluster centers for the input `sample`, cluster `membership` IDs,
    !>  and `sample` distances-squared from their corresponding cluster centers.<br>
    !>
    !>  \details
    !>  In data mining, the **k-means++** is an algorithm for choosing the initial values (or **seeds**) for the [k-means clustering algorithm](@ref pm_clustering::setKmeans).<br>
    !>  It was proposed in 2007 by David Arthur and Sergei Vassilvitskii, as an approximation algorithm for the NP-hard k-means problem.<br>
    !>  It offers a way of avoiding the sometimes poor clustering found by the standard k-means algorithm.<br>
    !>
    !>  **Intuition**<br>
    !>
    !>  The intuition behind k-means++ is that spreading out the \f$k\f$ initial cluster centers is a good thing:<br>
    !>  The first cluster center is chosen uniformly at random from the data points that are being clustered,
    !>  after which each subsequent cluster center is chosen from the remaining data points with probability proportional
    !>  to its squared distance from the closest existing cluster center to the point.<br>
    !>
    !>  **Algorithm**<br>
    !>
    !>  The exact algorithm is as follows:<br>
    !>  <ol>
    !>      <li>    Choose one center uniformly at random among the data points.<br>
    !>      <li>    For each data point \f$x\f$ not chosen yet, compute \f$D(x)\f$, the distance between \f$x\f$ and the nearest center that has already been chosen.<br>
    !>      <li>    Choose one new data point at random as a new center, using a weighted probability distribution where a point \f$x\f$ is chosen with probability proportional to \f$D(x)^2\f$.<br>
    !>      <li>    Repeat Steps 2 and 3 until \f$k\f$ centers have been chosen.<br>
    !>      <li>    Now that the initial centers have been chosen, proceed using standard k-means clustering.<br>
    !>  </ol>
    !>
    !>  **Performance improvements**<br>
    !>
    !>  The k-means++ seeding method yields considerable improvement in the final error of k-means algorithm.<br>
    !>  Although the initial selection in the algorithm takes extra time, the k-means part itself converges very quickly after this seeding
    !>  and thus the algorithm actually lowers the computation time.<br>
    !>  Based on the original paper, the method yields typically 2-fold improvements in speed, and for certain datasets, close to 1000-fold improvements in error.<br>
    !>  In these simulations the new method almost always performed at least as well as vanilla k-means in both speed and error.<br>
    !>
    !>  \note
    !>  The metric used within this generic interface is the [Euclidean distance](@ref pm_distanceEuclid).<br>
    !>
    !>  \param[inout]   rng         :   The input/output scalar that can be an object of,
    !>                                  <ol>
    !>                                      <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                              implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                      <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                              implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                                  </ol>
    !>  \param[out]     membership  :   The output vector of shape `(1:nsam)` of type `integer` of default kind \IK,
    !>                                  containing the membership of each input sample in `sample` from its nearest cluster `center`,
    !>                                  such that `cluster(membership(i))` is the nearest cluster center to the `i`th sample `sample(:, i)` at a squared-distance of `disq(i)`.<br>
    !>  \param[out]     disq        :   The output vector of shape `(1:nsam)` of the same type and kind as the input argument `sample`,
    !>                                  containing the Euclidean **squared** distance of each input sample in `sample` from its nearest cluster `center`.<br>
    !>  \param[in]      sample      :   The input vector, or matrix of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the sample of `nsam` points in an `ndim`-dimensional space whose corresponding cluster centers must be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If `sample` is a **vector** of shape `(1 : nsam)` and `center` is a **vector** of shape `(1 : ncls)`,
    !>                                              then the input `sample` must be a **collection** of `nsam` points (in univariate space).<br>
    !>                                      <li>    If `sample` is a **matrix** of shape `(1 : ndim, 1 : nsam)` and `center` is a **matrix** of shape `(1 : ndim, 1 : ncls)`,
    !>                                              then the input `sample` must be a **collection** of `nsam` points (in `ndim`-dimensional space).<br>
    !>                                  </ol>
    !>  \param[in]      ncls        :   The input scalar of type `integer` of default kind \IK,
    !>                                  containing the number of the desired clusters to be identified in the sample.<br>
    !>                                  (**optional**, default = `ubound(center, 2)`. It **must** be present **if and only if** the output arguments `center`, `size`, and `potential` are all missing.)
    !>  \param[out]     center      :   The output vector of shape `(1:ncls)` or matrix of shape `(1 : ndim, 1 : ncls)` of the same type and kind as the input argument `sample`,
    !>                                  containing the set of `ncls` **unique** random cluster centers (**centroids**) selected from the input sample
    !>                                  based on the computed memberships and minimum sample-cluster distances `disq`.<br>
    !>                                  (**optional**. If missing, no cluster `center` information will be output.)
    !>  \param[out]     size        :   The output vector of shape `(1:ncls)` type `integer` of default kind \IK,
    !>                                  containing the sizes (number of members) of the clusters with the corresponding centers output in the argument `center`.<br>
    !>                                  (**optional**. If missing, no cluster `size` information will be output.)
    !>  \param[out]     potential   :   The output vector of shape `(1:ncls)` of the same type and kind as the input argument `sample`,
    !>                                  the `i`th element of which contains the sum of squared distances of all members of the `i`th cluster from the cluster center as output in the `i`th element of `center`.<br>
    !>                                  (**optional**. If missing, no cluster `potential` information will be output.)
    !>
    !>  \interface{setKmeansPP}
    !>  \code{.F90}
    !>
    !>      use pm_clustering, only: setKmeansPP
    !>
    !>      call setKmeansPP(rng, membership(1 : nsam)  , disq(1 : nsam), sample(1 : ndim, 1 : nsam), ncls)
    !>      call setKmeansPP(rng, membership(1 : nsam)  , disq(1 : nsam), sample(1 : ndim, 1 : nsam), center(1 : ndim, 1 : ncls), size(1 : ncls), potential(1 : ncls))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `ubound(center, rank(center)) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, rank(sample)) == ubound(disq, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(center, rank(center)) == ubound(size, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(center, rank(center)) == ubound(potential, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, rank(sample)) == ubound(membership, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(center, rank(center)) <= ubound(sample, rank(sample))` must hold for the corresponding input arguments (the number of clusters must be less than or equal to the sample size).<br>
    !>  The condition `ubound(sample, 1) == ubound(center, 1)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  By definition, the number of points in the input `sample` must be larger than the specified number of clusters.<br>
    !>
    !>  \warnpure
    !>  The procedures of this generic interface are always `impure` when the input `rng` argument is set to an object of type [rngf_type](@ref pm_distUnif::rngf_type).<br>
    !>
    !>  \remark
    !>  The functionality of this generic interface is similar to [setMember](@ref pm_clustering::setMember), with the major difference being that
    !>  [setKmeansPP](@ref pm_clustering::setKmeansPP) simultaneously computes the new cluster centers and sample memberships, whereas
    !>  [setMember](@ref pm_clustering::setMember) computes the new sample memberships based on a given set of cluster centers.<br>
    !>
    !>  \remark
    !>  The output of [setKmeansPP](@ref pm_clustering::setKmeansPP) can be directly passed to
    !>  [setCenter](@ref pm_clustering::setCenter) to the learn the new updated cluster centers and their sizes.<br>
    !>
    !>  \note
    !>  Dropping the optional arguments can aid runtime performance.<br>
    !>  This is particularly relevant when the output of this generic interface is directly passed to the [k-means algorithm](@ref pm_clustering::setKmeans).<br>
    !>
    !>  \see
    !>  [setKmeans](@ref pm_clustering::setKmeans)<br>
    !>  [setCenter](@ref pm_clustering::setCenter)<br>
    !>  [setMember](@ref pm_clustering::setMember)<br>
    !>  [setKmeansPP](@ref pm_clustering::setKmeansPP)<br>
    !>  [Arthur, D.; Vassilvitskii, S. (2007). k-means++: the advantages of careful seeding](http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf)<br>
    !>
    !>  \example{setKmeansPP}
    !>  \include{lineno} example/pm_clustering/setKmeansPP/main.F90
    !>  \compilef{setKmeansPP}
    !>  \output{setKmeansPP}
    !>  \include{lineno} example/pm_clustering/setKmeansPP/main.out.F90
    !>  \postproc{setKmeansPP}
    !>  \include{lineno} example/pm_clustering/setKmeansPP/main.py
    !>  \vis{setKmeansPP}
    !>  \image html pm_clustering/setKmeansPP/setKmeansPP.png width=700
    !>
    !>  \test
    !>  [test_pm_clustering](@ref test_pm_clustering)
    !>
    !>  \finmain{setKmeansPP}
    !>  If you use or redistribute ideas based on this generic interface implementation, you should also cite the original k-means++ article:<br>
    !>  [Arthur, D.; Vassilvitskii, S. (2007). k-means++: the advantages of careful seeding](http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf)<br>
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin
    interface setKmeansPP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setKmeansPPDRNGF_RK5(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGF_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setKmeansPPDRNGF_RK4(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGF_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setKmeansPPDRNGF_RK3(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGF_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setKmeansPPDRNGF_RK2(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGF_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setKmeansPPDRNGF_RK1(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGF_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setKmeansPPDRNGX_RK5(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGX_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setKmeansPPDRNGX_RK4(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGX_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setKmeansPPDRNGX_RK3(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGX_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setKmeansPPDRNGX_RK2(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGX_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setKmeansPPDRNGX_RK1(rng, membership, disq, sample, ncls)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPDRNGX_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)                    :: ncls
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setKmeansPPORNGF_RK5(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGF_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setKmeansPPORNGF_RK4(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGF_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setKmeansPPORNGF_RK3(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGF_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setKmeansPPORNGF_RK2(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGF_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setKmeansPPORNGF_RK1(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGF_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setKmeansPPORNGX_RK5(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGX_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setKmeansPPORNGX_RK4(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGX_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setKmeansPPORNGX_RK3(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGX_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setKmeansPPORNGX_RK2(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGX_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setKmeansPPORNGX_RK1(rng, membership, disq, sample, center, size, potential)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansPPORNGX_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Compute and return an iteratively-refined set of cluster centers given the input `sample` using the k-means approach.<br>
    !>
    !>  \details
    !>  The **k-means** clustering is a method of vector quantization, originally from signal processing,
    !>  that aims to partition n observations into \f$k\f$ clusters in which each observation belongs to
    !>  the cluster with the nearest mean (cluster centers or cluster centroid), serving as a prototype of the cluster.<br>
    !>  This results in a partitioning of the data space into Voronoi cells.<br>
    !>  The k-means clustering minimizes within-cluster variances (squared Euclidean distances),
    !>  but not regular Euclidean distances, which would be the more difficult Weber problem:<br>
    !>  the mean optimizes squared errors, whereas only the geometric median minimizes Euclidean distances.<br>
    !>  For instance, better Euclidean solutions can be found using **k-medians** and **k-medoids**.<br>
    !>  The problem is computationally difficult (NP-hard); however, efficient heuristic algorithms converge quickly to a local optimum.<br>
    !>  These are usually similar to the expectation-maximization algorithm for mixtures of Gaussian distributions
    !>  via an iterative refinement approach employed by both k-means and Gaussian mixture modeling.<br>
    !>  They both use cluster centers to model the data; however, k-means clustering tends to find clusters of comparable spatial extent,
    !>  while the Gaussian mixture model allows clusters to have different shapes.<br>
    !>
    !>  **Algorithm**<br>
    !>
    !>  Given a set of observations \f$(x_1, x_2, \ldots, x_n)\f$, where each observation is a \f$d\f$-dimensional real vector,
    !>  the k-means clustering aims to partition the \f$n\f$ observations into \f$k\f$ (\f$\leq n\f$) sets
    !>  \f$S = \{S_1, S_2, \ldots, S_k\}\f$ so as to minimize the within-cluster sum of squares (WCSS) (i.e. variance).<br>
    !>  Formally, the objective is to find:<br>
    !>  \f{equation}{
    !>      \underset{\mathbf{S}}{\up{arg\,min}}
    !>      \sum_{i=1}^{k} \sum_{\mathbf{x} \in S_{i}} \left\|\mathbf{x} -{\boldsymbol{\mu}}_{i}\right\|^{2}
    !>      = {\underset{\mathbf{S}}{\up{arg\,min}}}\sum_{i=1}^{k}|S_{i}|\up{Var}S_{i} ~,
    !>  \f}
    !>  where \f$\mu_i\f$ is the mean (also called **centroid**) of points in \f$S_{i}\f$, i.e.
    !>  \f{equation}{
    !>  {\boldsymbol {\mu_{i}}} = {\frac{1}{|S_{i}|}} \sum_{\mathbf{x} \in S_{i}} \mathbf{x} ~,
    !>  \f}
    !>  where \f$|S_{i}|\f$ is the size of \f$S_{i}\f$, and \f$\|\cdot\|\f$ is the \f$L^2\f$-norm.<br>
    !>  This is equivalent to minimizing the pairwise squared deviations of points in the same cluster:<br>
    !>  \f{equation}{
    !>      \underset{\mathbf{S}}{\up{arg\,min}} \sum_{i=1}^{k}\,{\frac {1}{|S_{i}|}}\,\sum_{\mathbf{x}, \mathbf{y} \in S_{i}}\left\|\mathbf{x} - \mathbf{y} \right\|^{2} ~,
    !>  \f}
    !>  The equivalence can be deduced from identity
    !>  \f{equation}{
    !>      |S_{i}|\sum_{\mathbf{x} \in S_{i}}\left\|\mathbf{x} -{\boldsymbol{\mu}}_{i}\right\|^{2} =
    !>      {\frac{1}{2}}\sum _{\mathbf {x} ,\mathbf {y} \in S_{i}}\left\|\mathbf {x} -\mathbf {y} \right\|^{2} ~.
    !>  \f}
    !>  Since the total variance is constant, this is equivalent to maximizing the sum of squared deviations between points in different clusters.<br>
    !>  This deterministic relationship is also related to the law of total variance in probability theory.<br>
    !>
    !>  **Performance improvements**<br>
    !>
    !>  The k-means++ seeding method yields considerable improvement in the final error of k-means algorithm.<br>
    !>  For more information, see the documentation of [setKmeansPP](@ref pm_clustering::setKmeansPP).<br>
    !>
    !>  \note
    !>  The metric used within this generic interface is the [Euclidean distance](@ref pm_distanceEuclid).<br>
    !>
    !>  \param[inout]   rng         :   The input/output scalar that can be an object of,
    !>                                  <ol>
    !>                                      <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                              implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                      <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                              implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                                  </ol>
    !>                                  (**optional**,
    !>                                  If this argument is present, then all `intent(inout)` arguments below have `intent(out)` argument
    !>                                  and will be initialized using the [k-means++](@ref pm_clustering::setKmeansPP) algorithm.<br>
    !>                                  If this argument is missing, then the user must initialize all of the following arguments with `intent(inout)`
    !>                                  via [k-means++](@ref pm_clustering::setKmeansPP) or any other method before passing them to this generic interface.)
    !>  \param[inout]   membership  :   The input/output vector of shape `(1:nsam)` of type `integer` of default kind \IK,
    !>                                  containing the membership of each input sample in `sample` from its nearest cluster `center`,
    !>                                  such that `cluster(membership(i))` is the nearest cluster center to the `i`th sample `sample(:, i)` at a squared-distance of `disq(i)`.<br>
    !>                                  If the argument `rng` is missing, then `membership` has `intent(inout)` and must be properly initialized before calling this routine.<br>
    !>                                  If the argument `rng` is present, then `membership` has `intent(out)` and will contain the final cluster memberships on output.<br>
    !>  \param[inout]   disq        :   The input/output vector of shape `(1:nsam)` of the same type and kind as the input argument `sample`,
    !>                                  containing the Euclidean **squared** distance of each input sample in `sample` from its nearest cluster `center`.<br>
    !>                                  If the argument `rng` is missing, then `disq` has `intent(inout)` and must be properly initialized before calling this routine.<br>
    !>                                  If the argument `rng` is present, then `disq` has `intent(out)` and will contain the final squared-distances from cluster centers on output.<br>
    !>  \param[in]      sample      :   The input vector, or matrix of,
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the sample of `nsam` points in an `ndim`-dimensional space whose corresponding cluster centers must be computed.<br>
    !>                                  <ol>
    !>                                      <li>    If `sample` is a **vector** of shape `(1 : nsam)` and `center` is a **vector** of shape `(1 : ncls)`,
    !>                                              then the input `sample` must be a **collection** of `nsam` points (in univariate space).<br>
    !>                                      <li>    If `sample` is a **matrix** of shape `(1 : ndim, 1 : nsam)` and `center` is a **matrix** of shape `(1 : ndim, 1 : ncls)`,
    !>                                              then the input `sample` must be a **collection** of `nsam` points (in `ndim`-dimensional space).<br>
    !>                                  </ol>
    !>  \param[inout]   center      :   The input/output vector of shape `(1:ncls)` or matrix of shape `(1 : ndim, 1 : ncls)` of the same type and kind as the input argument `sample`,
    !>                                  containing the set of `ncls` **unique** cluster centers (**centroids**) computed based on the input sample memberships and minimum distances.<br>
    !>                                  If the argument `rng` is missing, then `center` has `intent(inout)` and must be properly initialized before calling this routine.<br>
    !>                                  If the argument `rng` is present, then `center` has `intent(out)` and will contain the final cluster centers on output.<br>
    !>  \param[inout]   size        :   The input/output vector of shape `(1:ncls)` type `integer` of default kind \IK,
    !>                                  containing the sizes (number of members) of the clusters with the corresponding centers output in the argument `center`.<br>
    !>                                  If the argument `rng` is missing, then `size` has `intent(inout)` and must be properly initialized before calling this routine.<br>
    !>                                  If the argument `rng` is present, then `size` has `intent(out)` and will contain the final cluster sizes (member counts) on output.<br>
    !>  \param[inout]   potential   :   The input/output vector of shape `(1:ncls)` of the same type and kind as the input argument `sample`,
    !>                                  the `i`th element of which contains the sum of squared distances of all members of the `i`th cluster from the cluster center as output in the `i`th element of `center`.<br>
    !>                                  If the argument `rng` is missing, then `potential` has `intent(inout)` and must be properly initialized before calling this routine (although its values are not explicitly referenced).<br>
    !>                                  If the argument `rng` is present, then `potential` has `intent(out)` and will contain the final cluster potentials (sums of squared distances from cluster centers) on output.<br>
    !>  \param[out]     failed      :   The output scalar of type `logical` of default kind \LK that is `.true.` **if and only if** the algorithm fails to converge within the user-specified or default criteria for convergence.<br>
    !>                                  Failure occurs only if `any(size < minsize) .or. maxniter < niter`.<br>
    !>  \param[out]     niter       :   The output scalar of type `integer` of default kind \IK, containing the number of refinement iterations performed within the algorithm to achieve convergence.<br>
    !>                                  An output `niter` value larger than the input `maxniter` implies lack of convergence before return.<br>
    !>                                  (**optional**. If missing, the number of refinement iterations will not be output.)
    !>  \param[in]      maxniter    :   The input non-negative scalar of type `integer` of default kind \IK, containing the maximum number of refinement iterations allowed within the algorithm to achieve convergence.<br>
    !>                                  If convergence does not occur within the maximum specified value, the output arguments can be passed again as is to the generic interface
    !>                                  (**without** the optional `rng` argument) to continue the refinement iterations until convergence.<br>
    !>                                  A reasonable choice can be `300` or comparable values.<br>
    !>                                  (**optional**, default = `300`.)
    !>  \param[in]      minsize     :   The input non-negative scalar of type `integer` of default kind \IK, containing the minimum allowed size of each cluster.<br>
    !>                                  If any cluster has any number of members below the specified `minsize`, the algorithm will return without achieving convergence.<br>
    !>                                  The situation can be detected of any element of the output `size` is smaller than the specified `minsize`.<br>
    !>                                  A reasonable choice can be `ndim = ubound(sample, 1)` or comparable values although any non-negative value including zero is possible.<br>
    !>                                  (**optional**, default = `1`.)
    !>  \param[in]      nfail       :   The input non-negative scalar of type `integer` of default kind \IK, containing the number of times the k-means algorithm is allowed to fail before returning without convergence.<br>
    !>                                  (**optional**. It **can** be present **only if** the argument `rng` is also present, allowing random initializations in case of failures.)
    !>
    !>  \interface{setKmeans}
    !>  \code{.F90}
    !>
    !>      use pm_clustering, only: setKmeans
    !>
    !>      call setKmeansPP(membership(1 : nsam)  , disq(1 : nsam), sample(1 : nsam)          , center(1 : ncls)          , size(1 : ncls), potential(1 : ncls), failed, niter = niter, maxniter = maxniter, minsize = minsize)
    !>      call setKmeansPP(membership(1 : nsam)  , disq(1 : nsam), sample(1 : ndim, 1 : nsam), center(1 : ndim, 1 : ncls), size(1 : ncls), potential(1 : ncls), failed, niter = niter, maxniter = maxniter, minsize = minsize)
    !>
    !>      call setKmeansPP(rng, membership(1 : nsam)  , disq(1 : nsam), sample(1 : nsam)          , center(1 : ncls)          , size(1 : ncls), potential(1 : ncls), failed, niter = niter, maxniter = maxniter, minsize = minsize, nfail = nfail)
    !>      call setKmeansPP(rng, membership(1 : nsam)  , disq(1 : nsam), sample(1 : ndim, 1 : nsam), center(1 : ndim, 1 : ncls), size(1 : ncls), potential(1 : ncls), failed, niter = niter, maxniter = maxniter, minsize = minsize, nfail = nfail)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  If the argument `rng` is present, then all arguments associated with [setKmeansPP](@ref pm_clustering::setKmeansPP) equally apply to this generic interface.<br>
    !>  The condition `ubound(center, rank(center)) > 0` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, rank(sample)) == ubound(disq, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(center, rank(center)) == ubound(size, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(center, rank(center)) == ubound(potential, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(sample, rank(sample)) == ubound(membership, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `ubound(center, rank(center)) <= ubound(sample, rank(sample))` must hold for the corresponding input arguments (the number of clusters must be less than or equal to the sample size).<br>
    !>  The condition `ubound(sample, 1) == ubound(center, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= maxniter` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= minsize` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= nfail` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures of this generic interface are always `impure` when the input `rng` argument is set to an object of type [rngf_type](@ref pm_distUnif::rngf_type).<br>
    !>
    !>  \remark
    !>  The functionality of this generic interface is highly similar to [setCenter](@ref pm_clustering::setCenter), with the major difference being that
    !>  [setKmeans](@ref pm_clustering::setKmeans) simultaneously computes the new cluster centers and sample memberships, whereas
    !>  [setCenter](@ref pm_clustering::setCenter) computes the new cluster centers based on given sample membership.<br>
    !>
    !>  \see
    !>  [setKmeans](@ref pm_clustering::setKmeans)<br>
    !>  [setCenter](@ref pm_clustering::setCenter)<br>
    !>  [setMember](@ref pm_clustering::setMember)<br>
    !>  [setKmeansPP](@ref pm_clustering::setKmeansPP)<br>
    !>  [k-means clustering](https://en.wikipedia.org/wiki/K-means_clustering)<br>
    !>
    !>  \example{setKmeans}
    !>  \include{lineno} example/pm_clustering/setKmeans/main.F90
    !>  \compilef{setKmeans}
    !>  \output{setKmeans}
    !>  \include{lineno} example/pm_clustering/setKmeans/main.out.F90
    !>  \postproc{setKmeans}
    !>  \include{lineno} example/pm_clustering/setKmeans/main.py
    !>  \vis{setKmeans}
    !>  \image html pm_clustering/setKmeans/setKmeans.png width=700
    !>
    !>  \test
    !>  [test_pm_clustering](@ref test_pm_clustering)
    !>
    !>  \finmain{setKmeans}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin
    interface setKmeans

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setKmeansInit_RK5(membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansInit_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(inout) , contiguous    :: disq(:)
        real(RKC)               , intent(inout) , contiguous    :: center(:,:)
        real(RKC)               , intent(inout) , contiguous    :: potential(:)
        integer(IK)             , intent(inout) , contiguous    :: membership(:)
        integer(IK)             , intent(inout) , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setKmeansInit_RK4(membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansInit_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(inout) , contiguous    :: disq(:)
        real(RKC)               , intent(inout) , contiguous    :: center(:,:)
        real(RKC)               , intent(inout) , contiguous    :: potential(:)
        integer(IK)             , intent(inout) , contiguous    :: membership(:)
        integer(IK)             , intent(inout) , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setKmeansInit_RK3(membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansInit_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(inout) , contiguous    :: disq(:)
        real(RKC)               , intent(inout) , contiguous    :: center(:,:)
        real(RKC)               , intent(inout) , contiguous    :: potential(:)
        integer(IK)             , intent(inout) , contiguous    :: membership(:)
        integer(IK)             , intent(inout) , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setKmeansInit_RK2(membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansInit_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(inout) , contiguous    :: disq(:)
        real(RKC)               , intent(inout) , contiguous    :: center(:,:)
        real(RKC)               , intent(inout) , contiguous    :: potential(:)
        integer(IK)             , intent(inout) , contiguous    :: membership(:)
        integer(IK)             , intent(inout) , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setKmeansInit_RK1(membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansInit_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(inout) , contiguous    :: disq(:)
        real(RKC)               , intent(inout) , contiguous    :: center(:,:)
        real(RKC)               , intent(inout) , contiguous    :: potential(:)
        integer(IK)             , intent(inout) , contiguous    :: membership(:)
        integer(IK)             , intent(inout) , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module subroutine setKmeansRNGF_RK5(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGF_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK4_ENABLED
    impure module subroutine setKmeansRNGF_RK4(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGF_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK3_ENABLED
    impure module subroutine setKmeansRNGF_RK3(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGF_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK2_ENABLED
    impure module subroutine setKmeansRNGF_RK2(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGF_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK1_ENABLED
    impure module subroutine setKmeansRNGF_RK1(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGF_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(rngf_type)         , intent(in)                    :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setKmeansRNGX_RK5(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGX_RK5
#endif
        use pm_kind, only: RKC => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setKmeansRNGX_RK4(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGX_RK4
#endif
        use pm_kind, only: RKC => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setKmeansRNGX_RK3(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGX_RK3
#endif
        use pm_kind, only: RKC => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setKmeansRNGX_RK2(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGX_RK2
#endif
        use pm_kind, only: RKC => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setKmeansRNGX_RK1(rng, membership, disq, sample, center, size, potential, failed, niter, maxniter, minsize, nfail)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKmeansRNGX_RK1
#endif
        use pm_kind, only: RKC => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        integer(IK)             , intent(in)    , optional      :: minsize
        integer(IK)             , intent(in)    , optional      :: maxniter
        integer(IK)             , intent(in)    , optional      :: nfail
        real(RKC)               , intent(in)    , contiguous    :: sample(:,:)
        real(RKC)               , intent(out)   , contiguous    :: disq(:)
        real(RKC)               , intent(out)   , contiguous    :: center(:,:)
        real(RKC)               , intent(out)   , contiguous    :: potential(:)
        integer(IK)             , intent(out)   , contiguous    :: membership(:)
        integer(IK)             , intent(out)   , contiguous    :: size(:)
        integer(IK)             , intent(out)   , optional      :: niter
        logical(LK)             , intent(out)                   :: failed
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_clustering