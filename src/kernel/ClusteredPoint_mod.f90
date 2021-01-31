!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief This module contains the classes and procedures for generating collecitons of clustered points with various distributions.
!> \author Amir Shahmoradi

module ClusteredPoint_mod

    use Constants_mod, only: RK, IK
    use Err_mod, only: Err_type

    implicit none

    character(len=*), parameter :: MODULE_NAME = "@ClusteredPoint_mod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: ClusteredPoint_type
        integer(IK)                 :: nd
        integer(IK)                 :: np
        integer(IK)                 :: nc
        integer(IK)                 :: ndmin
        integer(IK)                 :: ndmax
        integer(IK)                 :: ncmin
        integer(IK)                 :: ncmax
        integer(IK)                 :: sizeMin
        integer(IK)                 :: sizeMax
        real(RK)                    :: etaMin
        real(RK)                    :: etaMax
        real(RK)                    :: stdMin
        real(RK)                    :: stdMax
        real(RK)                    :: centerMin
        real(RK)                    :: centerMax
        real(RK)    , allocatable   :: Eta(:)
        real(RK)    , allocatable   :: Std(:,:)
        real(RK)    , allocatable   :: Point(:,:)
        real(RK)    , allocatable   :: Center(:,:)
        real(RK)    , allocatable   :: ChoDia(:,:)
        real(RK)    , allocatable   :: LogVolume(:)
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)
        integer(IK) , allocatable   :: Membership(:)
        integer(IK) , allocatable   :: Size(:)
        character(:), allocatable   :: dist
        type(Err_type)              :: Err
    contains
        procedure, pass :: get => getClusteredPoint
        procedure, pass :: write => writeClusteredPoint
    end type ClusteredPoint_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getClusteredPoint( self & ! LCOV_EXCL_LINE
                                , nd, ndmin, ndmax & ! LCOV_EXCL_LINE
                                , nc, ncmin, ncmax & ! LCOV_EXCL_LINE
                                , Size, sizeMin, sizeMax & ! LCOV_EXCL_LINE
                                , Center, centerMin, centerMax & ! LCOV_EXCL_LINE
                                , Std, stdmin, stdmax & ! LCOV_EXCL_LINE
                                , Eta, etaMin, etaMax & ! LCOV_EXCL_LINE
                                , dist & ! LCOV_EXCL_LINE
                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getClusteredPoint
#endif
        use Matrix_mod, only: getCholeskyFactor
        use Matrix_mod, only: getInvMatFromCholFac
        use Statistics_mod, only: isInsideEllipsoid
        use Statistics_mod, only: getRandCorMat
        use Statistics_mod, only: getRandMVU
        use Statistics_mod, only: getRandInt
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str !, getLowerCase
        use Math_mod, only: getCumSum
        implicit none
        class(ClusteredPoint_type), intent(inout)       :: self
        integer(IK) , intent(in), optional              :: nd, ndmin, ndmax
        integer(IK) , intent(in), optional              :: nc, ncmin, ncmax
        integer(IK) , intent(in), optional              :: sizeMin, sizeMax
        real(RK)    , intent(in), optional              :: stdMin, stdMax
        real(RK)    , intent(in), optional              :: etaMin, etaMax
        real(RK)    , intent(in), optional              :: centerMin, centerMax
        real(RK)    , intent(in), optional              :: Center(:,:), Std(:,:), Eta(:)
        integer(IK) , intent(in), optional              :: Size(:)
        character(*), intent(in), optional              :: dist
        integer(IK) , allocatable                       :: CumSumSize(:)
        real(RK)    , allocatable                       :: NormedPoint(:), InvCovMat(:,:,:)
        real(RK)    , allocatable                       :: CumSumVolNormed(:)
        real(RK)    , allocatable                       :: VolNormed(:)
        logical                                         :: isUniformSuperposed
        logical                                         :: isUniform, isMember
        logical                                         :: isNormal
        real(RK)                                        :: dummy
        integer(IK)                                     :: i, j, ic, ip, membershipCount

        self%Err%occurred = .false.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! ndmin, ndmax, nd
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(ndmin)) then
            self%ndmin = ndmin
        else
            self%ndmin = 2_IK
        endif

        if (present(ndmax)) then
            self%ndmax = ndmax
        else
            self%ndmax = self%ndmin * 10_IK
        endif

        if (present(nd)) then
            self%nd = nd
        else
            self%nd = getRandInt(self%ndmin, self%ndmax)
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! ncmin, ncmax, nc
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(ncmin)) then
            self%ncmin = ncmin
        else
            self%ncmin = 1_IK
        endif

        if (present(ncmax)) then
            self%ncmax = ncmax
        else
            self%ncmax = self%ncmin * 10_IK
        endif

        if (present(nc)) then
            self%nc = nc
        else
            self%nc = getRandInt(self%ncmin, self%ncmax)
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! sizeMin, sizeMax, Size
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(sizeMin)) then
            self%sizeMin = sizeMin
        else
            self%sizeMin = self%nd + 1
        endif

        if (present(sizeMax)) then
            self%sizeMax = sizeMax
        else
            self%sizeMax = self%sizeMin * 1000_IK
        endif

        if (present(Size)) then
            self%Size = Size
        else
            if (allocated(self%Size)) deallocate(self%Size) ! LCOV_EXCL_LINE ! GFortran crashes without this
            if (allocated(CumSumSize)) deallocate(CumSumSize) ! LCOV_EXCL_LINE ! GFortran crashes without this
            allocate(self%Size(self%nc), CumSumSize(0:self%nc))
            CumSumSize(0) = 0._RK
            do ic = 1, self%nc
                self%Size(ic) = getRandInt(self%sizeMin, self%sizeMax)
                CumSumSize(ic) = CumSumSize(ic-1) + self%Size(ic)
            end do
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! centerMin, centerMax, Size
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(centerMin)) then
            self%centerMin = centerMin
        else
            self%centerMin = 0._RK
        endif

        if (present(centerMax)) then
            self%centerMax = centerMax
        else
            self%centerMax = 1._RK
        endif

        if (present(Center)) then
            self%Center = Center
        else
            if (allocated(self%Center)) deallocate(self%Center) ! LCOV_EXCL_LINE ! GFortran crashes without this
            allocate(self%Center(self%nd,self%nc))
            call random_number(self%Center)
            do concurrent(ic = 1:self%nc)
                self%Center(:,ic) = self%centerMin + self%Center(:,ic) * (self%centerMax - self%centerMin)
            end do
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! stdMin, stdMax, Size
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(stdMin)) then
            self%stdMin = stdMin
        else
            self%stdMin = 0._RK
        endif

        if (present(stdMax)) then
            self%stdMax = stdMax
        else
            self%stdMax = 1._RK
        endif

        if (present(Std)) then
            self%Std = Std
        else
            if (allocated(self%Std)) deallocate(self%Std) ! LCOV_EXCL_LINE ! GFortran crashes without this
            allocate(self%Std(self%nd,self%nc))
            call random_number(self%Std)
            do concurrent(ic = 1:self%nc)
                self%Std(:,ic) = self%stdMin + self%Std(:,ic) * (self%stdMax - self%stdMin)
            end do
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! etaMin, etaMax, eta
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(etaMin)) then
            self%etaMin = etaMin
        else
            self%etaMin = 0.1_RK
        endif

        if (present(etaMax)) then
            self%etaMax = etaMax
        else
            self%etaMax = self%etaMin * (1._RK/self%etaMin)**2
        endif

        if (present(Eta)) then
            self%Eta = Eta
        else
            if (allocated(self%Eta)) deallocate(self%Eta) ! LCOV_EXCL_LINE ! GFortran crashes without this
            allocate(self%Eta(self%nc))
            call random_number(self%Eta)
            do concurrent(ic = 1:self%nc)
                self%Eta(ic) = self%etaMin + self%Eta(ic) * (self%etaMax - self%etaMin)
            end do
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate the representative matrices of the clusters
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (allocated(self%ChoDia)) deallocate(self%ChoDia) ! LCOV_EXCL_LINE ! GFortran crashes without this
        if (allocated(self%LogVolume)) deallocate(self%LogVolume) ! LCOV_EXCL_LINE ! GFortran crashes without this
        if (allocated(self%ChoLowCovUpp)) deallocate(self%ChoLowCovUpp) ! LCOV_EXCL_LINE ! GFortran crashes without this
        if (allocated(CumSumVolNormed)) deallocate(CumSumVolNormed) ! LCOV_EXCL_LINE ! GFortran crashes without this

        allocate(self%LogVolume(self%nc))
        allocate(self%ChoDia(self%nd,self%nc))
        allocate(self%ChoLowCovUpp(self%nd,self%nd,self%nc))
        allocate(CumSumVolNormed(self%nc))

        do ic = 1, self%nc

            ! Generate random correlation matrix.

            self%ChoLowCovUpp(:,:,ic) = getRandCorMat(self%nd, self%Eta(ic))
            do j = 1, self%nd
                do i = 1, self%nd
                    self%ChoLowCovUpp(i,j,ic) = self%ChoLowCovUpp(i,j,ic) * self%Std(i,ic) * self%Std(j,ic)
                end do
            end do

            ! Compute the Cholesky factorization.

            call getCholeskyFactor  ( nd = self%nd & ! LCOV_EXCL_LINE
                                    , PosDefMat = self%ChoLowCovUpp(1:self%nd,1:self%nd,ic) & ! LCOV_EXCL_LINE
                                    , Diagonal = self%ChoDia(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                    )
            if (self%ChoDia(1,ic) < 0._RK) then
                self%Err%occurred = .true.
                self%Err%msg = "Singular Covariance matrix detected."
                write(*,"(A)") "ChoLowCovUpp:"
                write(*,"("//num2str(self%nd)//"(F15.8,:,' '))") self%ChoLowCovUpp(1:self%nd,1:self%nd,ic)
                error stop
            end if

            self%LogVolume(ic) = sum(log(self%ChoDia(1:self%nd,ic)))

        end do

        VolNormed = exp( self%LogVolume - maxval(self%LogVolume) )
        CumSumVolNormed(1:self%nc) = getCumSum(self%nc, VolNormed)
        CumSumVolNormed(1:self%nc) = CumSumVolNormed / CumSumVolNormed(self%nc)


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Set the distribution
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(dist)) then
            self%dist = dist
        else
            self%dist = "uniform-mixture"
        end if

        isUniform = self%dist == "uniform"
        isNormal = self%dist == "normal-mixture"
        isUniformSuperposed = self%dist == "uniform-mixture"

        if (.not. (isUniformSuperposed .or. isUniform)) then
            self%Err%occurred = .true.
            self%Err%msg = "No point distribution other than uniform is currently supported."
            return
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate the cluster members
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        self%np = sum(self%Size)

        if (allocated(self%Point)) deallocate(self%Point) ! LCOV_EXCL_LINE ! GFortran crashes without this
        if (allocated(self%Membership)) deallocate(self%Membership) ! LCOV_EXCL_LINE ! GFortran crashes without this

        allocate(self%Membership(self%np))
        allocate(self%Point(self%nd,self%np))

        if (isUniformSuperposed) then

            do ic = 1, self%nc
                do ip = CumSumSize(ic-1) + 1, CumSumSize(ic)
                    self%Membership(ip) = ic
                    self%Point(1:self%nd,ip) = getRandMVU   ( nd = self%nd & ! LCOV_EXCL_LINE
                                                            , MeanVec = self%Center(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            , CholeskyLower = self%ChoLowCovUpp(1:self%nd,1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            , Diagonal = self%ChoDia(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            )
                end do
            end do

        elseif (isUniform) then

            if (allocated(InvCovMat)) deallocate(InvCovMat) ! LCOV_EXCL_LINE ! GFortran crashes without this
            if (allocated(NormedPoint)) deallocate(NormedPoint) ! LCOV_EXCL_LINE ! GFortran crashes without this
            allocate(NormedPoint(self%nd), InvCovMat(self%nd,self%nd,self%nc))

            do ic = 1, self%nc
                InvCovMat(1:self%nd,1:self%nd,ic) = getInvMatFromCholFac( nd = self%nd & ! LCOV_EXCL_LINE
                                                                        , CholeskyLower = self%ChoLowCovUpp(1:self%nd,1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                                        , CholeskyDiago = self%ChoDia(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                                        )
            end do

            self%Size = 0_IK ! all cluster sizes will be determined after all point generations.

            do ip = 1, self%np

                ! Check for the distribution type.

                loopAddClusterMember: do

                    call random_number(dummy)
                    ic = minloc(CumSumVolNormed, dim = 1, mask = CumSumVolNormed > dummy)
                    !ic = getRandInt(lowerBound = 1_IK, upperBound = self%nc)

                    self%Point(1:self%nd,ip) = getRandMVU   ( nd = self%nd & ! LCOV_EXCL_LINE
                                                            , MeanVec = self%Center(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            , CholeskyLower = self%ChoLowCovUpp(1:self%nd,1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            , Diagonal = self%ChoDia(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            )

                    membershipCount = 0_IK
                    do i = 1, self%nc
                        NormedPoint(1:self%nd) = self%Point(1:self%nd,ip) - self%Center(1:self%nd,i)
                        isMember = isInsideEllipsoid(nd = self%nd, NormedPoint = NormedPoint, InvRepMat = InvCovMat(1:self%nd,1:self%nd,i))
                        if (isMember) membershipCount = membershipCount + 1_IK
                    end do

                    call random_number(dummy)
                    if (dummy < 1._RK / membershipCount) then
                        self%Size(ic) = self%Size(ic) + 1_IK
                        self%Membership(ip) = ic
                        exit loopAddClusterMember
                    end if

                    cycle loopAddClusterMember

                end do loopAddClusterMember

            end do

        end if

    end subroutine getClusteredPoint

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writeClusteredPoint(self, fileUnit)
        use Constants_mod, only: IK, RK
        implicit none
        class(ClusteredPoint_type)  , intent(in)    :: self
        integer(IK)                 , intent(in)    :: fileUnit
        character(*), parameter                     :: fileFormat = "(*(g0.8,:,','))"
        integer(IK)                                 :: i, j, ic

        write(fileUnit,"(A)") "dist"
        write(fileUnit,fileFormat) self%dist

        write(fileUnit,"(A)") "nd, np, nc"
        write(fileUnit,fileFormat) self%nd, self%np, self%nc

        write(fileUnit,"(A)") "Size"
        write(fileUnit,fileFormat) self%Size

        write(fileUnit,"(A)") "Center"
        write(fileUnit,fileFormat) self%Center

        write(fileUnit,"(A)") "LogVolume"
        write(fileUnit,fileFormat) self%LogVolume

        write(fileUnit,"(A)") "CholeskyLower"
        write(fileUnit,fileFormat) ((self%ChoDia(j,ic), (self%ChoLowCovUpp(i,j,ic), i=j+1,self%nd), j=1,self%nd), ic=1,self%nc)

        write(fileUnit,"(A)") "Point"
        write(fileUnit,fileFormat) self%Point

        write(fileUnit,"(A)") "Membership"
        write(fileUnit,fileFormat) self%Membership

    end subroutine writeClusteredPoint

end module ClusteredPoint_mod ! LCOV_EXCL_LINE