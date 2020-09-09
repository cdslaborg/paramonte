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
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SpecMCMC_ProposalStartCovMat_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_ProposalStartCovMat_mod"

    real(RK), allocatable           :: proposalStartCovMat(:,:) ! namelist input

    type                            :: ProposalStartCovMat_type
        logical                     :: isPresent
        real(RK), allocatable       :: Def(:,:)
        real(RK), allocatable       :: Val(:,:)
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setProposalStartCorMat, checkForSanity, nullifyNameListVar
    end type ProposalStartCovMat_type

    interface ProposalStartCovMat_type
        module procedure            :: constructProposalStartCovMat
    end interface ProposalStartCovMat_type

    private :: constructProposalStartCovMat, setProposalStartCorMat, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructProposalStartCovMat(nd,methodName) result(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructProposalStartCovMat
#endif
        use Constants_mod, only: IK, NULL_RK
        use String_mod, only: num2str
        implicit none
        integer(IK), intent(in)         :: nd
        character(*), intent(in)        :: methodName
        type(ProposalStartCovMat_type)  :: self
        integer(IK)                     :: i
        self%isPresent = .false.
        allocate( self%Def(nd,nd) )
        self%Def = 0._RK
        do i = 1,nd
            self%Def(i,i) = 1._RK
        end do
        self%null   = NULL_RK
        self%desc   = &
        "proposalStartCovMat is a real-valued positive-definite matrix of size (ndim,ndim), where ndim is the dimension of the &
        &sampling space. It serves as the best-guess starting covariance matrix of the proposal distribution. &
        &To bring the sampling efficiency of " // methodName // " to within the desired requested range, the covariance matrix will &
        &be adaptively updated throughout the simulation, according to the user's requested schedule. If proposalStartCovMat &
        &is not provided by the user or it is completely missing from the input file, its value will be automatically computed &
        &via the input variables proposalStartCorMat and proposalStartStdVec (or via their default values, if not provided). &
        &The default value of proposalStartCovMat is an ndim-by-ndim Identity matrix."
    end function constructProposalStartCovMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(self,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(ProposalStartCovMat_type), intent(in) :: self
        integer(IK), intent(in)                     :: nd
        if (allocated(proposalStartCovMat)) deallocate(proposalStartCovMat)
        allocate(proposalStartCovMat(nd,nd))
        proposalStartCovMat = self%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setProposalStartCorMat(self, proposalStartStdVec, proposalStartCorMat, proposalStartCovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setProposalStartCorMat
#endif
        use Statistics_mod, only: getCovMatFromCorMat
        use Constants_mod, only: RK, IK
        implicit none
        class(ProposalStartCovMat_type), intent(inout)  :: self
        real(RK), intent(in)                            :: proposalStartCorMat(:,:), proposalStartStdVec(:)
        real(RK), intent(in), optional                  :: proposalStartCovMat(:,:)
        integer(IK)                                     :: i, j, nd

        if (present(proposalStartCovMat)) then
            self%val = proposalStartCovMat
        else
            self%val = self%null
        end if

        self%isPresent = .false.
        nd = size(proposalStartCorMat(:,1))
        do i = 1, nd
            do j = 1, nd
                if (self%val(j,i)==self%null) then
                    self%val(j,i) = self%Def(j,i)
                else
                    self%isPresent = .true.
                end if
            end do
        end do

        if (self%isPresent) return

        self%val = getCovMatFromCorMat  ( nd = nd &
                                        , StdVec = proposalStartStdVec &
                                        , CorMat = proposalStartCorMat &
                                        )

    end subroutine setProposalStartCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(self,Err,methodName,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use Matrix_mod, only: isPosDef
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        class(ProposalStartCovMat_type), intent(in) :: self
        integer(IK), intent(in)                     :: nd
        character(*), intent(in)                    :: methodName
        type(Err_type), intent(inout)               :: Err
        character(*), parameter                     :: PROCEDURE_NAME = "@checkForSanity()"
        if (.not.isPosDef(nd,self%val)) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested proposalStartCovMat for the proposal of " // methodName // &
                        " is not a positive-definite matrix.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecMCMC_ProposalStartCovMat_mod