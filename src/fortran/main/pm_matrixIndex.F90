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
!>  This module contains procedures and generic interfaces for
!>  converting the indices of matrix elements between different packing and storage formats.
!>
!>  \details
!>  see the documentation of [pm_matrixPack](@ref pm_matrixPack) and [pm_matrixSubset](@ref pm_matrixSubset) for further details.
!>
!>  \see
!>  [pm_matrixCopy](@ref pm_matrixCopy)<br>
!>  [pm_matrixPack](@ref pm_matrixPack)<br>
!>  [pm_matrixSubset](@ref pm_matrixSubset)<br>
!>
!>  \test[test_pm_matrixIndex](@ref test_pm_matrixIndex)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixIndex

    use pm_kind, only: SK, IK
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: lowDia, lowDia_type
    use pm_matrixPack, only: rdpack, rdpack_type
    use pm_matrixPack, only: lfpack, lfpack_type
    use pm_matrixPack, only: rfpack, rfpack_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_matrixPack"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the index of a specific element of a matrix of a specific storage in
    !>  specific packing format from the corresponding index in an alternative packing format.
    !>
    !>  \details
    !>  Specifically, this generic interface allows the conversion of indices of upper and lower triangular matrices in
    !>  <ol>
    !>      <li>    [Rectangular Default Package](@ref pm_matrixPack::rdpack_type) to [Linear Packed Package](@ref pm_matrixPack::lfpack_type) and vice versa.
    !>      <li>    [Rectangular Default Package](@ref pm_matrixPack::rdpack_type) to [Rectangular Full Package](@ref pm_matrixPack::rfpack_type) and vice versa.
    !>  </ol>
    !>
    !>  \param[in]  dpack   :   The input scalar that can be,
    !>                          <ol>
    !>                              <li>    the constant [lfpack](@ref pm_matrixPack::lfpack)
    !>                                      **only if** the input argument `spack` is [rdpack](@ref pm_matrixPack::rdpack),
    !>                              <li>    the constant [rfpack](@ref pm_matrixPack::rfpack)
    !>                                      **only if** the input argument `spack` is [rdpack](@ref pm_matrixPack::rdpack),
    !>                              <li>    the constant [rdpack](@ref pm_matrixPack::rdpack)
    !>                                      **only if** the input argument `spack` is [rfpack](@ref pm_matrixPack::rfpack) or [lfpack](@ref pm_matrixPack::lfpack),
    !>                          </ol>
    !>                          representing the target packing format to which the input `sindex` must be converted.<br>
    !>  \param[in]  spack   :   The input scalar that can be,
    !>                          <ol>
    !>                              <li>    the constant [lfpack](@ref pm_matrixPack::lfpack)
    !>                                      **only if** the input argument `dpack` is [rdpack](@ref pm_matrixPack::rdpack),
    !>                              <li>    the constant [rfpack](@ref pm_matrixPack::rfpack)
    !>                                      **only if** the input argument `dpack` is [rdpack](@ref pm_matrixPack::rdpack),
    !>                          </ol>
    !>                          representing the original packing format from which the input `sindex` must be converted.<br>
    !>  \param[in]  sindex  :   The input scalar or vector of size `(2)` of type `integer` of default kind \IK,
    !>                          representing the index of an element of the source matrix in `spack` packing format
    !>                          to be converted to the target output index in `dpack` packing format.<br>
    !>                          It must be scalar **if and only if** the input argument `spack` is set to
    !>                          [lfpack](@ref pm_matrixPack::lfpack), otherwise it must be a vector.<br>
    !>  \param[in]  subset  :   The input scalar that can be,
    !>                          <ol>
    !>                              <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia)
    !>                                      indicating the upper-triangular storage of the original matrix in `spack` packing.
    !>                              <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia)
    !>                                      indicating the lower-triangular storage of the original matrix in `spack` packing.
    !>                          </ol>
    !>  \param[in]  shape   :   The input vector of size `1` or `2` representing the shape of the original
    !>                          upper/lower triangular matrix in the **default Rectangular Standard Packing format**.<br>
    !>                          This value is readily returned by the Fortran intrinsic `shape()` if the original/target
    !>                          matrix already exists in the default Rectangular Standard Packing.<br>
    !>  \param[in]  doff    :   The input scalar of type `integer` of default kind \IK, representing the offset of the
    !>                          diagonal of the upper/lower triangular matrix with respect to its top-left corner element.<br>
    !>                          By definition,
    !>                          <ol>
    !>                              <li>    the diagonal offset is non-positive for upper-triangular matrix storage.
    !>                              <li>    the diagonal offset is non-negative for lower-triangular matrix storage.
    !>                              <li>    the diagonal offset is always zero for matrices in Rectangular Full Packing (RFP) format.
    !>                              <li>    the diagonal offset is always zero for upper/lower square matrices.
    !>                          </ol>
    !>                          (**optional**, default = `0`. It can be present **if and only if**
    !>                          neither the original or the target packing is Rectangular Full Package.)
    !>
    !>  \return
    !>  `dindex`            :   The output scalar or vector of size `(2)` of type `integer` of default kind \IK,
    !>                          representing the converted index of the corresponding element of the source
    !>                          matrix in `spack` packing to the target `dpack` packing format.<br>
    !>                          It is a scalar **if and only if** the input argument `dpack` is set to
    !>                          [lfpack](@ref pm_matrixPack::lfpack), otherwise it is a vector of size `2`.<br>
    !>
    !>  \interface{getMatIndex}
    !>  \code{.F90}
    !>
    !>      use pm_matrixIndex, only: getMatIndex
    !>      use pm_matrixSubset, only: uppDia, lowDia
    !>      use pm_matrixPack, only: rdpack, rfpack, lfpack
    !>
    !>      integer(IK) :: dindex, sindex
    !>
    !>      dindex = getMatIndex(dpack, spack, sindex, subset, shape, doff = doff)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0 <= shape)` must hold for the corresponding input arguments.<br>
    !>  The condition `all([0 < sindex])` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= shape(1) + doff .and. doff <= 0` must hold for the corresponding input arguments when `subset = uppDia`.<br>
    !>  The condition `0 <= shape(2) - doff .and. 0 <= doff` must hold for the corresponding input arguments when `subset = lowDia`.<br>
    !>  The condition `all(sindex <= shape)` must hold for the corresponding input arguments when `spack = rdpack`.<br>
    !>  The condition `sindex(1) <= sindex(2) - doff` must hold for the corresponding input arguments when `spack = rdpack, subset = uppDia`.<br>
    !>  The condition `sindex(2) <= sindex(1) + doff` must hold for the corresponding input arguments when `spack = rdpack, subset = lowDia`.<br>
    !>  The specified input `sindex` must correspond to an element within the specified storage `subset` and packing `spack` of the original matrix.<br>
    !>  The condition `shape(1) == shape(2)` must hold for the corresponding input arguments when `spack = rfpack` or `dpack = rfpack`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [pm_matrixCopy](@ref pm_matrixCopy)<br>
    !>  [pm_matrixPack](@ref pm_matrixPack)<br>
    !>  [pm_matrixSubset](@ref pm_matrixSubset)<br>
    !>
    !>  \example{getMatIndex}
    !>  \include{lineno} example/pm_matrixIndex/getMatIndex/main.F90
    !>  \compilef{getMatIndex}
    !>  \output{getMatIndex}
    !>  \include{lineno} example/pm_matrixIndex/getMatIndex/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixIndex](@ref test_pm_matrixIndex)
    !>
    !>  \todo
    !>  \phigh
    !>  This generic interface should be extended to convert indices of matrices from
    !>  Rectangular Standard Packing to Rectangular Band Packing and vice versa.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  The implementation of the procedures for converting indices from LFP to RDP could be likely improved for better performance.<br>
    !>
    !>  \final{getMatIndex}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getMatIndex

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getMatIndex_D0_LFP_RDP_UXD_AIO(dpack, spack, sindex, subset, shape, doff) result(dindex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatIndex_D0_LFP_RDP_UXD_AIO
#endif
        use pm_kind, only: IKG => IK
        type(lfpack_type)   , intent(in)                    :: dpack
        type(rdpack_type)   , intent(in)                    :: spack
        type(uppDia_type)   , intent(in)                    :: subset
        integer(IKG)        , intent(in)                    :: sindex(2)
        integer(IKG)        , intent(in)                    :: shape(2)
        integer(IKG)        , intent(in)    , optional      :: doff
        integer(IKG)                                        :: dindex
    end function

    PURE module function getMatIndex_D0_LFP_RDP_XLD_AIO(dpack, spack, sindex, subset, shape, doff) result(dindex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatIndex_D0_LFP_RDP_XLD_AIO
#endif
        use pm_kind, only: IKG => IK
        type(lfpack_type)   , intent(in)                    :: dpack
        type(rdpack_type)   , intent(in)                    :: spack
        type(lowDia_type)   , intent(in)                    :: subset
        integer(IKG)        , intent(in)                    :: sindex(2)
        integer(IKG)        , intent(in)                    :: shape(2)
        integer(IKG)        , intent(in)    , optional      :: doff
        integer(IKG)                                        :: dindex
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getMatIndex_D0_RDP_LFP_UXD_AIO(dpack, spack, sindex, subset, shape, doff) result(dindex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatIndex_D0_RDP_LFP_UXD_AIO
#endif
        use pm_kind, only: IKG => IK
        type(rdpack_type)   , intent(in)                    :: dpack
        type(lfpack_type)   , intent(in)                    :: spack
        type(uppDia_type)   , intent(in)                    :: subset
        integer(IKG)        , intent(in)                    :: sindex
        integer(IKG)        , intent(in)                    :: shape(2)
        integer(IKG)        , intent(in)    , optional      :: doff
        integer(IKG)                                        :: dindex(2)
    end function

    PURE module function getMatIndex_D0_RDP_LFP_XLD_AIO(dpack, spack, sindex, subset, shape, doff) result(dindex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatIndex_D0_RDP_LFP_XLD_AIO
#endif
        use pm_kind, only: IKG => IK
        type(rdpack_type)   , intent(in)                    :: dpack
        type(lfpack_type)   , intent(in)                    :: spack
        type(lowDia_type)   , intent(in)                    :: subset
        integer(IKG)        , intent(in)                    :: sindex
        integer(IKG)        , intent(in)                    :: shape(2)
        integer(IKG)        , intent(in)    , optional      :: doff
        integer(IKG)                                        :: dindex(2)
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getMatIndex_D0_RFP_RDP_UXD_AIO(dpack, spack, sindex, subset, shape) result(dindex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatIndex_D0_RFP_RDP_UXD_AIO
#endif
        use pm_kind, only: IKG => IK
        type(rfpack_type)   , intent(in)                    :: dpack
        type(rdpack_type)   , intent(in)                    :: spack
        type(uppDia_type)   , intent(in)                    :: subset
        integer(IKG)        , intent(in)                    :: sindex(2)
        integer(IKG)        , intent(in)                    :: shape(2)
        integer(IKG)                                        :: dindex(2)
    end function

    PURE module function getMatIndex_D0_RFP_RDP_XLD_AIO(dpack, spack, sindex, subset, shape) result(dindex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatIndex_D0_RFP_RDP_XLD_AIO
#endif
        use pm_kind, only: IKG => IK
        type(rfpack_type)   , intent(in)                    :: dpack
        type(rdpack_type)   , intent(in)                    :: spack
        type(lowDia_type)   , intent(in)                    :: subset
        integer(IKG)        , intent(in)                    :: sindex(2)
        integer(IKG)        , intent(in)                    :: shape(2)
        integer(IKG)                                        :: dindex(2)
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getMatIndex_D0_RDP_RFP_UXD_AIO(dpack, spack, sindex, subset, shape) result(dindex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatIndex_D0_RDP_RFP_UXD_AIO
#endif
        use pm_kind, only: IKG => IK
        type(rdpack_type)   , intent(in)                    :: dpack
        type(rfpack_type)   , intent(in)                    :: spack
        type(uppDia_type)   , intent(in)                    :: subset
        integer(IKG)        , intent(in)                    :: sindex(2)
        integer(IKG)        , intent(in)                    :: shape(2)
        integer(IKG)                                        :: dindex(2)
    end function

    PURE module function getMatIndex_D0_RDP_RFP_XLD_AIO(dpack, spack, sindex, subset, shape) result(dindex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatIndex_D0_RDP_RFP_XLD_AIO
#endif
        use pm_kind, only: IKG => IK
        type(rdpack_type)   , intent(in)                    :: dpack
        type(rfpack_type)   , intent(in)                    :: spack
        type(lowDia_type)   , intent(in)                    :: subset
        integer(IKG)        , intent(in)                    :: sindex(2)
        integer(IKG)        , intent(in)                    :: shape(2)
        integer(IKG)                                        :: dindex(2)
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixIndex