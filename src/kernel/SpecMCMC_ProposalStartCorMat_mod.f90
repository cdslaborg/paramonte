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

!> \brief
!> This module contains the classes and procedures for setting up the `proposalStartCorMat` attribute of samplers of class [ParaMCMC_type](@ref paramcmc_mod::paramcmc_type).
!> For more information, see the description of this attribute in the body of the module.
!> \author Amir Shahmoradi

module SpecMCMC_ProposalStartCorMat_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_ProposalStartCorMat_mod"

    real(RK), allocatable           :: proposalStartCorMat(:,:) ! namelist input

    type                            :: ProposalStartCorMat_type
        real(RK), allocatable       :: val(:,:)
        real(RK), allocatable       :: def(:,:)
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setProposalStartCorMat, checkForSanity, nullifyNameListVar
    end type ProposalStartCorMat_type

    interface ProposalStartCorMat_type
        module procedure            :: constructProposalStartCorMat
    end interface ProposalStartCorMat_type

    private :: constructProposalStartCorMat, setProposalStartCorMat, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructProposalStartCorMat(nd,methodName) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructProposalStartCorMat
#endif
        use Constants_mod, only: IK, NULL_RK
        use String_mod, only: num2str
        use Matrix_mod, only: getEye
        implicit none
        integer(IK), intent(in)         :: nd
        character(*), intent(in)        :: methodName
        type(ProposalStartCorMat_type)  :: self
        allocate( self%Def(nd,nd) )
        self%Def    = getEye(nd,nd)
        self%null   = NULL_RK
        self%desc   = &
        "proposalStartCorMat is a real-valued positive-definite matrix of size (ndim,ndim), where ndim is the dimension of the &
        &sampling space. It serves as the best-guess starting correlation matrix of the proposal distribution used by " &
        // methodName // ". &
        &It is used (along with the input vector ProposalStartStdVec) to construct the covariance matrix of the proposal &
        &distribution when the input covariance matrix is missing in the input list of variables. &
        &If the covariance matrix is given as input to "//methodName//", any input values for proposalStartCorMat, &
        &as well as ProposalStartStdVec, will be automatically ignored by "//methodName//". As input to " // methodName // &
        ", the variable proposalStartCorMat along with ProposalStartStdVec is especially useful in situations where &
        &obtaining the best-guess covariance matrix is not trivial. The default value of proposalStartCorMat is an ndim-by-ndim Identity matrix."
    end function constructProposalStartCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(self,nd)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(ProposalStartCorMat_type), intent(in) :: self
        integer(IK), intent(in)                     :: nd
        if (allocated(proposalStartCorMat)) deallocate(proposalStartCorMat)
        allocate(proposalStartCorMat(nd,nd))
        proposalStartCorMat = self%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setProposalStartCorMat(self,proposalStartCorMat)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: setProposalStartCorMat
#endif
        use Constants_mod, only: RK
        implicit none
        class(ProposalStartCorMat_type), intent(inout)  :: self
        real(RK), intent(in)                            :: proposalStartCorMat(:,:)
        self%Val = proposalStartCorMat
        where (self%Val==self%null)
            self%Val = self%Def
        end where
    end subroutine setProposalStartCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! There is no need to check for eyeness of the input correlation matrix. Only positive definiteness is enough.
    ! If the input correlation matrix is problematic, it will eventually lead to a non-positive-definite covariance matrix.
    subroutine checkForSanity(self,Err,methodName,nd)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK, SPR
        use Matrix_mod, only: isPosDef
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        class(ProposalStartCorMat_type), intent(in) :: self
        integer(IK), intent(in)                     :: nd
        character(*), intent(in)                    :: methodName
        type(Err_type), intent(inout)               :: Err
        character(*), parameter                     :: PROCEDURE_NAME = "@checkForSanity()"
       !integer(IK)                                 :: i, j
        if (.not.isPosDef(nd,self%Val)) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. The input requested proposalStartCorMat &
                        &for the proposal of " // methodName // " is not a positive-definite matrix.\n\n"
        end if
        !do j = 1, nd
        !    if (abs(proposalStartCorMat(j,j) - 1._RK) > 1.e-10_RK) then
        !        Err%occurred = .true.
        !        Err%msg =   Err%msg // &
        !                    MODULE_NAME // PROCEDURE_NAME // ": Error occurred. The input requested element &
        !                    &proposalStartCorMat("//num2str(j)//","//num2str(j)//") = "//num2str(proposalStartCorMat(j,j))// &
        !                    " must be, by definition, equal to one.\n\n"
        !    end if
        !    do i = 1, j-1
        !        if ( abs(proposalStartCorMat(i,j)) >= 1._RK ) then
        !            Err%occurred = .true.
        !            Err%msg =   Err%msg // &
        !                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. The input requested element &
        !                        &proposalStartCorMat("//num2str(i)//","//num2str(j)//") = "//num2str(proposalStartCorMat(i,j))// &
        !                        " must be, by definition, bounded within the open range (-1,1).\n\n"
        !        end if
        !    end do
        !end do
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecMCMC_ProposalStartCorMat_mod ! LCOV_EXCL_LINE