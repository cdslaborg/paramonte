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
!>  This module contains procedures and generic interfaces for computing the nearest neighbor statistics of random samples.<br>
!>
!>  \details
!>  The k-nearest neighbors algorithm (k-NN) is a non-parametric supervised learning
!>  method first developed by Evelyn Fix and Joseph Hodges in 1951 and later expanded by Thomas Cover.<br>
!>  It is used for classification and regression and in both cases, the input consists of the `k` closest training examples in a data set.<br>
!>  The output depends on whether k-NN is used for classification or regression:<br>
!>  <ol>
!>      <li>    In k-NN classification, the output is a class membership.<br>
!>              An object is classified by a plurality vote of its neighbors, with the object being assigned to the
!>              class most common among its `k` nearest neighbors (k is a positive integer, typically small).<br>
!>              If `k = 1`, then the object is simply assigned to the class of that single nearest neighbor.<br>
!>      <li>    In k-NN regression, the output is the property value for the object.<br>
!>              This value is the average of the values of `k` nearest neighbors.<br>
!>              If `k = 1`, then the output is simply assigned to the value of that single nearest neighbor.<br>
!>  </ol>
!>  k-NN is a type of classification where the function is only approximated locally and all computation is deferred until function evaluation.<br>
!>  Since this algorithm relies on distance for classification, if the features represent different physical units or come in
!>  vastly different scales then normalizing the training data can improve its accuracy dramatically.<br>
!>
!>  Algorithm
!>  ---------
!>
!>  The training examples are vectors in a multidimensional feature space, each with a class label.<br>
!>  The training phase of the algorithm consists only of storing the feature vectors and class labels of the training samples.<br>
!>  In the classification phase, `k` is a user-defined constant, and an unlabeled vector (a query or test point)
!>  is classified by assigning the label which is most frequent among the k training samples nearest to that query point.<br>
!>  A commonly used distance metric for continuous variables is Euclidean distance.<br>
!>  For discrete variables, such as for text classification, another metric can be used, such as the overlap metric (or Hamming distance).<br>
!>  In the context of gene expression microarray data, for example, k-NN has been employed with correlation coefficients, such as Pearson and Spearman, as a metric.<br>
!>  Often, the classification accuracy of k-NN can be improved significantly if the distance metric is learned with specialized
!>  algorithms such as Large Margin Nearest Neighbor or Neighborhood components analysis.<br>
!>
!>  Drawbacks
!>  ---------
!>
!>  A major drawback of the basic *majority voting* classification occurs when the class distribution is skewed.<br>
!>  That is, examples of a more frequent class tend to dominate the prediction of the new example,
!>  because they tend to be common among the `k` nearest neighbors due to their large number.<br>
!>  One way to overcome this problem is to weight the classification, taking into account the distance from the test point to each of its k nearest neighbors.<br>
!>  The class (or value, in regression problems) of each of the `k` nearest points is multiplied by a weight proportional to
!>  the inverse of the distance from that point to the test point.<br>
!>  Another way to overcome skew is by abstraction in data representation.<br>
!>  For example, in a self-organizing map (SOM), each node is a representative (a center) of a cluster of similar points,
!>  regardless of their density in the original training data. K-NN can then be applied to the SOM.<br>
!>
!>  \see
!>  [pm_distanceMahal](@ref pm_distanceMahal)<br>
!>  [pm_distanceEuclid](@ref pm_distanceEuclid)<br>
!>  [pm_distanceHellinger](@ref pm_distanceHellinger)<br>
!>  [pm_distanceManhattan](@ref pm_distanceManhattan)<br>
!>  [pm_distanceMinkowski](@ref pm_distanceMinkowski)<br>
!>
!>  \test
!>  [test_pm_distanceEuclid](@ref test_pm_distanceEuclid)<br>
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Thursday 8:40 PM, July 20, 2023, Dallas, TX
!>  \AmirShahmoradi, Saturday 1:00 AM, September, 1, 2018, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_knn

    use pm_kind, only: SK, IK
    use pm_container

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_knn"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the input distance matrix whose columns are sorted in ascending order on output, optionally only up to the `k`th row of each column,
    !>  such that the `k`th row in the `i`th column is the `k`th nearest neighbor to the \f$i^{th}\f$ reference point.<br>
    !>
    !>  \details
    !>  The input matrix of shape `(1:npnt, 1:nref)` represents the distances of `nref` reference points from `npnt` neighbor points.<br>
    !>  The matrix entries can be any measure of distance as long they are non-negative.<br>
    !>
    !>  \param[inout]   distance    :   The input or input/output `contiguous` array of rank `2` of shape `(1:npnt, 1:nref)` of
    !>                                  <ol>
    !>                                      <li>    type `real` of kind \RKALL.<br>
    !>                                  </ol>
    !>                                  each column `i` of which contains the distances of `npnt` neighbor
    !>                                  points stored by the column rows from the `i`th reference point.<br>
    !>                                  For example, such distance matrix can be obtained via:<br>
    !>                                  <ol>
    !>                                      <li>    [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid) or [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>                                  </ol>
    !>                                  or any other distance metric module of the library.<br>
    !>                                  <ol>
    !>                                      <li>    If the input argument `rank` is missing, then `distance` has `intent(inout)`.<br>
    !>                                              On output, the rows of each column of `distance` are sorted separately in ascending order
    !>                                              such that the `k`th row in column `i` contains the distance of the `k`th nearest neighbor
    !>                                              of `reference(1:ndim, i)` represented by column `distance(1:npnt, i)`.<br>
    !>                                      <li>    If the input argument `rank` is present, then `distance` has `intent(in)`.<br>
    !>                                              On output, the rows of each column of `rank` are filled with the rank of
    !>                                              the corresponding element of `distance` ranked by its nearness to the
    !>                                              reference represented by column `distance(1:npnt, i)`.<br>
    !>                                  </ol>
    !>  \param[out]     rank        :   The output `contiguous` array of the same rank and shape as `distance` of type `integer` of default kind \IK
    !>                                  each element of which contains the nearest neighbor rank of the corresponding element of `distance`.<br>
    !>                                  The neighborhood ranking is with respect to the point represented by the column of the element.<br>
    !>                                  (**optional**. It can be present **if and only if** the input argument `k` is missing.)
    !>  \param[in]      k           :   The output scalar `integer` of default kind \IK, representing the row in each column of `distance`
    !>                                  up to which the distances in each column will be sorted in ascending order.<br>
    !>                                  The rest of the rows beyond the `k`th row remain unsorted.<br>
    !>                                  This argument is useful when sorting the entire distance matrix or computing the entire `rank` matrix is unnecessary.<br>
    !>                                  (**optional**, default = `size(distance, 1)`. It can be present **if and only if** the input argument `rank` is missing.)
    !>
    !>  \interface{setKnnSorted}
    !>  \code{.F90}
    !>
    !>      use pm_knn, only: setKnnSorted
    !>
    !>      call setKnnSorted(distance(1:npnt, 1:nref))
    !>      call setKnnSorted(distance(1:npnt, 1:nref), k)
    !>      call setKnnSorted(distance(1:npnt, 1:nref), rank(1:npnt, 1:nref))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(shape(distance) == shape(rank))` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < k .and. k <= size(distance, 1)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  To sort the neighbors of an individual point (corresponding to a single column vector),
    !>  use the generic interfaces of [pm_arrayRank](@ref pm_arrayRank) or [pm_arraySort](@ref pm_arraySort).<br>
    !>
    !>  \see
    !>  [getDisEuclid](@ref pm_distanceEuclid::getDisEuclid)<br>
    !>  [setDisEuclid](@ref pm_distanceEuclid::setDisEuclid)<br>
    !>  [getDisMatEuclid](@ref pm_distanceEuclid::getDisMatEuclid)<br>
    !>  [setDisMatEuclid](@ref pm_distanceEuclid::setDisMatEuclid)<br>
    !>
    !>  \example{setKnnSorted}
    !>  \include{lineno} example/pm_knn/setKnnSorted/main.F90
    !>  \compilef{setKnnSorted}
    !>  \output{setKnnSorted}
    !>  \include{lineno} example/pm_knn/setKnnSorted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_knn](@ref test_pm_knn)
    !>
    !>  \todo
    !>  \phigh
    !>  This generic interface should be extended to support discrete distance matrices.<br>
    !>
    !>  \final{setKnnSorted}
    !>
    !>  \author
    !>  \FatemehBagheri, Thursday 8:40 PM, July 20, 2023, Dallas, TX
    !>  \AmirShahmoradi, Saturday 1:00 AM, September, 1, 2018, Dallas, TX

    ! rdpack, euclid

    interface setKnnSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setKnnSortedVal_RK5(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedVal_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setKnnSortedVal_RK4(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedVal_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setKnnSortedVal_RK3(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedVal_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setKnnSortedVal_RK2(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedVal_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setKnnSortedVal_RK1(distance)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedVal_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setKnnSortedKth_RK5(distance, k)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedKth_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
        integer(IK)         , intent(in)                    :: k
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setKnnSortedKth_RK4(distance, k)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedKth_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
        integer(IK)         , intent(in)                    :: k
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setKnnSortedKth_RK3(distance, k)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedKth_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
        integer(IK)         , intent(in)                    :: k
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setKnnSortedKth_RK2(distance, k)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedKth_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
        integer(IK)         , intent(in)                    :: k
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setKnnSortedKth_RK1(distance, k)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedKth_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: distance(:,:)
        integer(IK)         , intent(in)                    :: k
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setKnnSortedInd_RK5(distance, rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedInd_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: distance(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rank(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setKnnSortedInd_RK4(distance, rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedInd_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: distance(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rank(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setKnnSortedInd_RK3(distance, rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedInd_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: distance(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rank(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setKnnSortedInd_RK2(distance, rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedInd_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: distance(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rank(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setKnnSortedInd_RK1(distance, rank)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setKnnSortedInd_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: distance(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rank(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \cond excluded
#if 0
    type :: minDisEdge_type!(IK, RK)
       !integer , kind  :: IK = IK_DEF
       !integer , kind  :: RK = RK_DEF
        type(IV)        :: index
        type(RV)        :: lenSq
        integer(IK)     :: npnt
    end type minDisEdge_type

    interface minDisEdge_type
        module procedure :: getMinDisEdge
    end interface minDisEdge_type

    !>  \todo
    !>  As of gfortran-10, default values for kind-type-parameters in PDT result in strange compile-time errors
    !>  when the PDT is `use`d in another module. Example error messgae:
    !>  \code
    !>      Error: Cannot convert TYPE(Pdthub_type_4_8) to TYPE(hub_type) at (1)
    !>  \endcode
    type :: hub_type!(IK, RK)
       !integer(IK) , kind          :: IK!= IK_DEF  !< \public integer kind type parameter.
       !integer(IK) , kind          :: RK!= RK_DEF  !< \public real kind type parameter.
        integer(IK)                 :: nh           !< \public number of hubs
        integer(IK) , allocatable   :: nodeIndex(:) !< \public length of `nh`
        integer(IK) , allocatable   :: edgeCount(:) !< \public length of `nh`
        type(IV)    , allocatable   :: edgeIndex(:) !< \public length of `nh`
        type(RV)    , allocatable   :: edgeLenSq(:) !< \public length of `nh`
        type(minDisEdge_type)       :: minDisEdge
    end type

    interface hub_type
        module procedure :: getHub
    end interface hub_type
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if 0
    !> \warning
    !> The input value `npnt` must be larger than 1 at all times. A value of `npnt = 1` is meaningless.
    !> In addition, a value of `npnt = 1` is also meaningless, although this procedure can handle `npnt = 2` gracefully.
    pure function getMinDisEdge(pairDisSq) result(minDisEdge)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinDisEdge
#endif
        real(RK), intent(in) :: pairDisSq(:,:)
        type(minDisEdge_type)           :: minDisEdge

        integer(IK) :: pindex1
        integer(IK) :: pindex2
        integer(IK) :: npnt, jp

        npnt = size(pairDisSq, 1, IK)
        CHECK_ASSERTION(__LINE__, size(pairDisSq, 1, IK) == size(pairDisSq, 2, IK), SK_"@setKnnSorted(): The condition `size(pairDisSq, 1) == size(pairDisSq, 2)` must hold. shape(pairDisSq) = "//getStr([shape(pairDisSq, IK)]))

        allocate(minDisEdge%index%val(npnt), minDisEdge%lenSq%val(npnt))
        minDisEdge%npnt = npnt

        jp = 1_IK
        pindex2 = minloc(pairDisSq(2:npnt, jp), dim = 1, kind = IK) + 1_IK
        minDisEdge%index%val(jp) = pindex2
        minDisEdge%lenSq%val(jp) = pairDisSq(pindex2,jp)

        jp = npnt
        pindex1 = minloc(pairDisSq(1:jp-1, jp), dim = 1, kind = IK)
        minDisEdge%index%val(jp) = pindex1
        minDisEdge%lenSq%val(jp) = pairDisSq(pindex1,jp)

        do concurrent(jp = 2 : npnt - 1)
            pindex1 = minloc(pairDisSq(1 : jp - 1, jp), dim = 1, kind = IK)
            pindex2 = minloc(pairDisSq(jp + 1 : npnt, jp), dim = 1, kind = IK) + jp
            if (pairDisSq(pindex1,jp) < pairDisSq(pindex2,jp)) then
                minDisEdge%index%val(jp) = pindex1
                minDisEdge%lenSq%val(jp) = pairDisSq(pindex1,jp)
            else
                minDisEdge%index%val(jp) = pindex2
                minDisEdge%lenSq%val(jp) = pairDisSq(pindex2,jp)
            end if
        end do

    end function getMinDisEdge

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE function getHub(npnt, pairDisSq) result(hub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHub
#endif
        use pm_kind, only: IK, RK
        use pm_arrayUnique, only: setUnique

        implicit none

        integer(IK) , intent(in)    :: npnt
        real(RK)    , intent(in)    :: pairDisSq(npnt,npnt)
        type(hub_type)              :: hub

        character(*, SK), parameter     :: PROCEDURE_NAME = MODULE_NAME//SK_"@getHub()"
        integer(IK)                 :: ih

        hub%minDisEdge = minDisEdge_type(npnt,pairDisSq)
        call setUnique  ( Array = hub%minDisEdge%index%val & ! LCOV_EXCL_LINE
                        , Unique = hub%nodeIndex & ! LCOV_EXCL_LINE
                        , count = hub%edgeCount & ! LCOV_EXCL_LINE
                        , index = hub%edgeIndex & ! LCOV_EXCL_LINE
                        , order = -1_IK & ! LCOV_EXCL_LINE
                        )
        hub%nh = size(hub%nodeIndex, kind = IK)
        allocate(hub%edgeLenSq(hub%nh))
        do concurrent(ih = 1:hub%nh)
            hub%edgeLenSq(ih)%val = hub%minDisEdge%lenSq%val(hub%edgeIndex(ih)%val)
        end do

    end function getHub
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !> Return the estimated mean and standard deviation of the volume occupied by a set of `npnt`
!    !> points uniformly distributed in an `nd` dimensional domain of arbitrary shape and size.
!    !>
!    !> \param[in]       nd      :   The number of dimensions.
!    !> \param[in]       npnt      :   The number of input points.
!    !> \param[in]       sample   :   The input array of points of shape `(nd,npnt)`.
!    !>
!    !> \return
!    !> `innestIndex` : The index of the point with the smallest sum of distances-squared from all other points.
!    !>
!    !> \remark
!    !>  This method relies on minimizing the KS-test statistic. Specifically, the point with
!    !>
!    !> \todo: The code performance and memory usage could be improved by return a vector instead of a Symmetric matrix.
!    pure function getMinSumPairDistSqPointIndex(nd, npnt, sample) result(innestIndex)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getMinSumPairDistSqPointIndex
!#endif
!        use pm_knn, only: getPairDistSq ! LCOV_EXCL_LINE
!        use pm_kind, only: IK, RK ! LCOV_EXCL_LINE
!        implicit none
!        integer(IK) , intent(in)    :: nd, npnt
!        real(RK)    , intent(in)    :: sample(nd,npnt)
!        integer(IK)                 :: innestIndex
!        innestIndex = minloc(sum(getPairDistSq(nd,npnt,sample), dim = 1_IK), dim = 1_IK)
!    end function getMinSumPairDistSqPointIndex

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!>  \endcond excluded
end module pm_knn ! LCOV_EXCL_LINE