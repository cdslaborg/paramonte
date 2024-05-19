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
!>  This include file contains procedure implementations of the tests of [pm_distUnif::setUnifCDF](@ref pm_distUnif::setUnifCDF).
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getUnifCDF_ENABLED || setUnifCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define GET_REAL(object) object%re
#define GET_IMAG(object) object%im
#else
#define GET_REAL(object) object
#define GET_IMAG(object) object
#endif
        integer(IK) :: i
#if     SK_ENABLED
        character(1, SK)                :: lower, upper
        character(5, SK), allocatable   :: Sample(:), cdf(:), cdf_ref(:)
#elif   IK_ENABLED
        real(RKG)       , parameter     :: TOL = epsilon(0._RKG) * 100._RKG
        real(RKG)       , allocatable   :: cdf(:), cdf_ref(:), diff(:)
        integer(IKG)    , allocatable   :: Sample(:)
        integer(IKG)                    :: lower, upper
#elif   CK_ENABLED
        complex(CKG)     , parameter    :: TOL = epsilon(0._CKG) * 100._CKG
        complex(CKG)     , allocatable  :: Sample(:), cdf(:), cdf_ref(:), diff(:)
        complex(CKG)                    :: lower, upper
#elif   RK_ENABLED
        real(RKG)       , parameter     :: TOL = epsilon(0._RKG) * 100._RKG
        real(RKG)       , allocatable   :: Sample(:), cdf(:), cdf_ref(:), diff(:)
        real(RKG)                       :: lower, upper
#else
#error  "Unrecognized interface."
#endif
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        lower = achar(1_IK)
        upper = achar(127_IK)
#elif   IK_ENABLED
        lower = 0_IKG
        upper = 1_IKG
        Sample = [-1_IKG, 0_IKG, 1_IKG, 2_IKG]
        cdf_ref = [0._RKG, 0.5_RK, 1._RKG, 1._RKG]
#elif   RK_ENABLED
        lower = 0._RKG
        upper = 1._RKG
        Sample = [-1._RKG, 0._RKG, 1._RKG, 2._RKG]
        cdf_ref = [0._RKG, 0._RKG, 1._RKG, 1._RKG]
#elif   CK_ENABLED
        lower = (0._CKG,0._CKG)
        upper = (1._CKG,1._CKG)
        Sample = [(-1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,1._CKG), (2._CKG,2._CKG)]
        cdf_ref = [(0._CKG,0._CKG), (0._CKG,0._CKG), (1._CKG,1._CKG), (1._CKG,1._CKG)]
#else
#error  "Unrecognized interface."
#endif
        allocate(cdf(size(Sample)), diff(size(Sample)))

        do i = 1_IK, size(Sample, kind = IK)
#if         getUnifCDF_ENABLED
            cdf(i) = getUnifCDF(Sample(i))
#elif       setUnifCDF_ENABLED
            call setUnifCDF(cdf(i), Sample(i))
#else
#error      "Unrecognized interface."
#endif
            diff(i) = cdf(i) - cdf_ref(i)
            assertion = assertion .and. abs(GET_REAL(diff(i))) <= GET_REAL(TOL) .and. abs(GET_IMAG(diff(i))) <= GET_IMAG(TOL)
            call report()
            call test%assert(assertion, desc = "The CDF must be computed correctly for the default range.")
        end do

#if     getUnifCDF_ENABLED
        cdf = getUnifCDF(Sample)
#elif   setUnifCDF_ENABLED
        call setUnifCDF(cdf, Sample)
#endif
        diff = cdf - cdf_ref
        assertion = assertion .and. all(abs(GET_REAL(diff)) <= GET_REAL(TOL)) .and. all(abs(GET_IMAG(diff)) <= GET_IMAG(TOL))
        call report()
        call test%assert(assertion, desc = "The CDF must be computed correctly for the default range in 1D.")

        do i = 1_IK, size(Sample, kind = IK)
#if         getUnifCDF_ENABLED
            cdf(i) = getUnifCDF(Sample(i), lower, upper)
#elif       setUnifCDF_ENABLED
            call setUnifCDF(cdf(i), Sample(i), lower, upper)
#endif
            diff(i) = cdf(i) - cdf_ref(i)
            assertion = assertion .and. abs(GET_REAL(diff(i))) <= GET_REAL(TOL) .and. abs(GET_IMAG(diff(i))) <= GET_IMAG(TOL)
            call report()
            call test%assert(assertion, desc = "The CDF must be computed correctly for when the range is specified to be the default.")
        end do

#if         getUnifCDF_ENABLED
            cdf = getUnifCDF(Sample, lower, upper)
#elif       setUnifCDF_ENABLED
            call setUnifCDF(cdf, Sample, lower, upper)
#endif
        diff = cdf - cdf_ref
        assertion = assertion .and. all(abs(GET_REAL(diff)) <= GET_REAL(TOL)) .and. all(abs(GET_IMAG(diff)) <= GET_IMAG(TOL))
        call report()
        call test%assert(assertion, desc = "The CDF must be computed correctly for when the range is specified to be the default in 1D.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     SK_ENABLED
        lower = achar(1_IK)
        upper = achar(127_IK)
#elif   IK_ENABLED
        lower = -2_IKG
        upper = +3_IKG
        Sample = [-4_IKG, -3_IKG, 2_IKG, 1_IKG, 3_IKG, 4_IKG]
#elif   RK_ENABLED
        lower = -2._RKG
        upper = +3._RKG
        Sample = [-1._RKG, 0._RKG, 1._RKG, 2._RKG]
#elif   CK_ENABLED
        lower = (-2._CKG,+7._CKG)
        upper = (+3._CKG,+5._CKG)
        Sample = [ (-4._CKG,-3._CKG), (-4._CKG,+3._CKG), (-4._CKG,+2._CKG), (-4._CKG,+5._CKG), (-4._CKG,+6._CKG) & ! LCOV_EXCL_LINE
                , (-2._CKG,-3._CKG), (-2._CKG,+3._CKG), (-2._CKG,+2._CKG), (-2._CKG,+5._CKG), (-2._CKG,+6._CKG) & ! LCOV_EXCL_LINE
                , (+0._CKG,-3._CKG), (+0._CKG,+3._CKG), (+0._CKG,+2._CKG), (+0._CKG,+5._CKG), (+0._CKG,+6._CKG) & ! LCOV_EXCL_LINE
                , (+1._CKG,-3._CKG), (+1._CKG,+3._CKG), (+1._CKG,+2._CKG), (+1._CKG,+5._CKG), (+1._CKG,+6._CKG) & ! LCOV_EXCL_LINE
                , (+7._CKG,-3._CKG), (+7._CKG,+3._CKG), (+7._CKG,+2._CKG), (+7._CKG,+5._CKG), (+7._CKG,+6._CKG) & ! LCOV_EXCL_LINE
                , (+9._CKG,-3._CKG), (+9._CKG,+3._CKG), (+9._CKG,+2._CKG), (+9._CKG,+5._CKG), (+9._CKG,+6._CKG) & ! LCOV_EXCL_LINE
                ]
#else
#error  "Unrecognized interface."
#endif
        cdf_ref = getCDF(Sample,lower,upper)
        allocate(cdf, diff, mold = cdf_ref)

        do i = 1_IK, size(Sample, kind = IK)
#if         getUnifCDF_ENABLED
            cdf(i) = getUnifCDF(Sample(i), lower, upper)
#elif       setUnifCDF_ENABLED
            call setUnifCDF(cdf(i), Sample(i), lower, upper)
#endif
            diff(i) = cdf(i) - cdf_ref(i)
            assertion = assertion .and. abs(GET_REAL(diff(i))) <= GET_REAL(TOL) .and. abs(GET_IMAG(diff(i))) <= GET_IMAG(TOL)
            call report()
            call test%assert(assertion, desc = "The CDF must be computed correctly for the non-default range.")
        end do

#if         getUnifCDF_ENABLED
            cdf = getUnifCDF(Sample, lower, upper)
#elif       setUnifCDF_ENABLED
            call setUnifCDF(cdf, Sample, lower, upper)
#endif
        diff = cdf - cdf_ref
        assertion = assertion .and. all(abs(GET_REAL(diff)) <= GET_REAL(TOL)) .and. all(abs(GET_IMAG(diff)) <= GET_IMAG(TOL))
        call report()
        call test%assert(assertion, desc = "The CDF must be computed correctly for the non-default range in 1D.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure elemental function getCDF(x,lower,upper) result(cdf)
#if     IK_ENABLED
            integer(IKG), intent(in) :: x, lower, upper
            real(RKG) :: cdf
            if (x < lower) then
                cdf = 0._RKG
            elseif (x < upper) then
                cdf = real(x-lower+1_IKG,RKG)/ real(upper-lower+1_IKG,RKG)
            else
                cdf = 1._RKG
            end if
#elif       RK_ENABLED
            real(RKG), intent(in) :: x, lower, upper
            real(RKG) :: cdf
            if (x < lower) then
                cdf = 0._RKG
            elseif (x < upper) then
                cdf = (x-lower) / (upper-lower)
            else
                cdf = 1._RKG
            end if
#elif       CK_ENABLED
            complex(CKG), intent(in) :: x, lower, upper
            complex(CKG) :: cdf
            if (x%re < lower%re) then
                cdf%re = 0._CKG
            elseif (x%re < upper%re) then
                cdf%re = (x%re-lower%re) / (upper%re-lower%re)
            else
                cdf%re = 1._CKG
            end if
            if (x%im < lower%im) then
                cdf%im = 0._CKG
            elseif (x%im < upper%im) then
                cdf%im = (x%im-lower%im) / (upper%im-lower%im)
            else
                cdf%im = 1._CKG
            end if
#endif
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(cdf)) deallocate(cdf)
            if (allocated(diff)) deallocate(diff)
            if (allocated(Sample)) deallocate(Sample)
            if (allocated(cdf_ref)) deallocate(cdf_ref)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Sample ", Sample
                write(test%disp%unit,"(*(g0,:,', '))") "CDF_ref", CDF_ref
                write(test%disp%unit,"(*(g0,:,', '))") "CDF    ", CDF
                write(test%disp%unit,"(*(g0,:,', '))") "diff   ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL    ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#elif   getUnifRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define CHECK_LOWER(object) object%re >= lower%re .and. object%im >= lower%im
#define CHECK_UPPER(object) object%re <= upper%re - TOL%re .and. object%im <= upper%im - TOL%im
#elif   SK_ENABLED
        use pm_str, only: getCharVec
!#define CHECK_LOWER(object) getCharVec(reshape(object,shape=[size(object)]))>=lower
!#define CHECK_UPPER(object) getCharVec(reshape(object,shape=[size(object)]))<=upper
#else
#define CHECK_LOWER(object) object>=lower
#define CHECK_UPPER(object) object<=upper-TOL
#endif
        integer(IK) , parameter :: NP = 80000_IK
#if     SK_ENABLED
        character(5,SK)         :: lower, upper
        character(5,SK)         :: Sample(NP)
#elif   LK_ENABLED
        logical(LK)             :: Sample(NP)
#elif   IK_ENABLED
        integer(IKG), parameter :: TOL = 0_IKG
        integer(IKG)            :: lower, upper
        integer(IKG)            :: Sample(NP)
#elif   RK_ENABLED
        real(RK)    , parameter :: TOL = epsilon(1._RK)
        real(RK)                :: lower, upper
        real(RK)                :: Sample(NP)
#elif   CK_ENABLED
        complex(CK) , parameter :: TOL = epsilon(1._CK)
        complex(CK)             :: lower, upper
        complex(CK)             :: Sample(NP)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: i
        assertion = .true._LK

#if     LK_ENABLED

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Sample(1) = getUnifRand()
        call test%assert(assertion, SK_"The procedure must be able to output a scalar logical harvest of default kind.")

        Sample = getUnifRand(size = size(Sample, kind = IK))
        call test%assert(assertion, SK_"The procedure must be able to output 1D array of logical harvest of default kind.")

        !>  \todo
        !>  The following test should be improved by ensuring uniformity of the resulting harvest values.
        Sample = getUnifRand(size = size(Sample, kind = IK))
        assertion = assertion .and. any(Sample)
        call test%assert(assertion, SK_"The procedure must be able to output 1D array of logical harvest of default kind with all values equally sampled.")

        !>  \todo
        !>  The following test should be improved by ensuring uniformity of the resulting harvest values.
        Sample = getUnifRand(size = size(Sample, kind = IK))
        assertion = assertion .and. .not. all(Sample)
        call test%assert(assertion, SK_"The procedure must be able to output 1D array of logical harvest of default kind with all values equally sampled.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! SK_ENABLED
#else

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED
        lower = repeat(char(50_IK),5)
        upper = repeat(char(93_IK),5)
#elif   IK_ENABLED
        lower = -5_IKG
        upper = +5_IKG
#elif   RK_ENABLED
        lower = -10._RK
        upper = +10._RK
#elif   CK_ENABLED
        lower = (-10._CK,-10._CK)
        upper = (+10._CK,+10._CK)
#else
#error  "Unrecognized interface."
#endif

        do i = 1, size(Sample, kind = IK)
            Sample(i) = getUnifRand(lower, upper)
        end do
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output scalar harvest within the specified limits.")

        Sample = getUnifRand(lower, upper, size = size(Sample, kind = IK))
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output 1D array harvest within the specified limits.")

        Sample = getUnifRand(lower, [( upper, i = 1, size(Sample, kind = IK))])
        call report()
        call test%assert(assertion, SK_"The procedure must be able to elementally output 1D array harvest within the specified limits.")

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        subroutine report()
            do i = 1, size(Sample)
#if             SK_ENABLED
                assertion = assertion .and. Sample(i) >= lower .and. Sample(i) <= upper
#elif           IK_ENABLED
                assertion = assertion .and. Sample(i) >= lower .and. Sample(i) <= upper
#elif           RK_ENABLED
                assertion = assertion .and. Sample(i) >= lower .and. Sample(i) < upper
#elif           CK_ENABLED
                assertion = assertion .and. Sample(i)%re >= lower%re .and. Sample(i)%re < upper%re
                assertion = assertion .and. Sample(i)%im >= lower%im .and. Sample(i)%im < upper%im
#else
#error          "Unrecognized interface."
#endif
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    !Sample(1:NPALL) => Sample
                    !do i = 1, NPALL
                        !if (.not. (CHECK_LOWER(Sample(i)) .or. CHECK_UPPER(Sample(i)))) then
                            write(test%disp%unit,"(*(g0,:,', '))")
                            write(test%disp%unit,"(*(g0,:,', '))") "Sample(i)   ", Sample(i)
                            write(test%disp%unit,"(*(g0,:,', '))")
#if                         SK_ENABLED
                            block
                            integer :: j, index(len(Sample))
                            do j = 1, len(Sample)
                                index(j) = ichar(Sample(i)(j:j))
                            end do
                            write(test%disp%unit,"(*(g0,:,', '))")
                            write(test%disp%unit,"(*(g0,:,', '))") "ichar(Sample(i))    ", index
                            write(test%disp%unit,"(*(g0,:,', '))") "ichar(lower), ichar(upper) ", ichar(lower(1:1)), ichar(upper(1:1))
                            write(test%disp%unit,"(*(g0,:,', '))")
                            end block
#endif
                            write(test%disp%unit,"(*(g0,:,', '))")
                            write(test%disp%unit,"(*(g0,:,', '))") "lower, upper", lower, upper
                            write(test%disp%unit,"(*(g0,:,', '))")
                            exit
                        !end if
                    !end do
                    !nullify(Sample)
                    ! LCOV_EXCL_STOP
                end if
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif
! SK_ENABLED

        !%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define CHECK_LOWER(object) object%re>=lower%re.and.object%im>=lower%im
#define CHECK_UPPER(object) object%re<=upper%re-TOL%re.and.object%im<=upper%im-TOL%im
#elif   SK_ENABLED
        use pm_str, only: getCharVec
!#define CHECK_LOWER(object) getCharVec(reshape(object,shape=[size(object)]))>=lower
!#define CHECK_UPPER(object) getCharVec(reshape(object,shape=[size(object)]))<=upper
#else
#define CHECK_LOWER(object) object>=lower
#define CHECK_UPPER(object) object<=upper-TOL
#endif
        integer(IK) , parameter :: NP1 = 20_IK
        integer(IK) , parameter :: NP2 = 20_IK
        integer(IK) , parameter :: NP3 = 20_IK
        integer(IK) , parameter :: NP4 = 20_IK
        integer(IK) , parameter :: NPALL = NP1*NP2*NP3*NP4

#if     SK_ENABLED
        character(5,SK)         :: lower, upper
        character(5,SK), target :: Sample(NP1,NP2,NP3,NP4)
        character(5,SK), pointer:: SamplePointer(:) => null()
#elif   LK_ENABLED
        logical(LK)             :: lower, upper
        logical(LK) , target    :: Sample(NP1,NP2,NP3,NP4)
        logical(LK) , pointer   :: SamplePointer(:) => null()
#elif   IK_ENABLED
        integer(IKG), parameter :: TOL = 0_IKG
        integer(IKG)            :: lower, upper
        integer(IKG), target    :: Sample(NP1,NP2,NP3,NP4)
        integer(IKG), pointer   :: SamplePointer(:) => null()
#elif   RK_ENABLED
        real(RK)    , parameter :: TOL = epsilon(1._RK)
        real(RK)                :: lower, upper
        real(RK)    , target    :: Sample(NP1,NP2,NP3,NP4)
        real(RK)    , pointer   :: SamplePointer(:) => null()
#elif   CK_ENABLED
        complex(CK) , parameter :: TOL = epsilon(1._CK)
        complex(CK)             :: lower, upper
        complex(CK) , target    :: Sample(NP1,NP2,NP3,NP4)
        complex(CK) , pointer   :: SamplePointer(:) => null()
#else
#error  "Unrecognized interface."
#endif

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED
        lower = repeat(char(1_IK, kind = SK), 5)
        upper = repeat(char(127_IK, kind = SK), 5)
#elif   IK_ENABLED
        lower = 0_IKG
        upper = 1_IKG
#elif   LK_ENABLED
        lower = .false._LK
        upper = .true._LK
#elif   RK_ENABLED
        lower = 0._RK
        upper = 1._RK
#elif   CK_ENABLED
        lower = (0._CK,0._CK)
        upper = (1._CK,1._CK)
#else
#error  "Unrecognized interface."
#endif

        call setUnifRand(Sample(1,1,1,1))
        SamplePointer(1:1) => Sample(1:1,1,1,1)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output scalar harvest.")

        call setUnifRand(Sample(:,1,1,1))
        SamplePointer(1:NP1) => Sample(1:NP1,1,1,1)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output 1D array harvest.")

        call setUnifRand(Sample(:,:,1,1))
        SamplePointer(1:NP1*NP2) => Sample(:,:,1,1)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output 2D array harvest.")

        call setUnifRand(Sample(:,:,:,1))
        SamplePointer(1:NP1*NP2*NP3) => Sample(:,:,:,1)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output 3D array harvest.")

        call setUnifRand(Sample(:,:,:,:))
        SamplePointer(1:NP1*NP2*NP3*NP4) => Sample(:,:,:,:)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output 4D array harvest.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !LK_ENABLED

#if     SK_ENABLED
        lower = repeat(char(50_IK),5)
        upper = repeat(char(93_IK),5)
#elif   IK_ENABLED
        lower = -5_IKG
        upper = +5_IKG
#elif   RK_ENABLED
        lower = -10._RK
        upper = +10._RK
#elif   CK_ENABLED
        lower = (-10._CK,-10._CK)
        upper = (+10._CK,+10._CK)
#else
#error  "Unrecognized interface."
#endif

#if     !LK_ENABLED
        call setUnifRand(Sample(1,1,1,1), lower, upper)
        SamplePointer(1:1) => Sample(1:1,1,1,1)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output scalar harvest within the specified limits.")

        call setUnifRand(Sample(:,1,1,1), lower, upper)
        SamplePointer(1:NP1) => Sample(1:NP1,1,1,1)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output 1D array harvest within the specified limits.")

        call setUnifRand(Sample(:,:,1,1), lower, upper)
        SamplePointer(1:NP1*NP2) => Sample(:,:,1,1)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output 2D array harvest within the specified limits.")

        call setUnifRand(Sample(:,:,:,1), lower, upper)
        SamplePointer(1:NP1*NP2*NP3) => Sample(:,:,:,1)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output 3D array harvest within the specified limits.")

        call setUnifRand(Sample(:,:,:,:), lower, upper)
        SamplePointer(1:NP1*NP2*NP3*NP4) => Sample(:,:,:,:)
        call report()
        call test%assert(assertion, SK_"The procedure must be able to output 4D array harvest within the specified limits.")
#endif

#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#if     SK_ENABLED
!        lower = char(50_IK)
!        upper = char(93_IK)
!        call setUnifRand(Sample(1,1,1,1), ichar(lower,kind=IK), ichar(upper,kind=IK))
!        assertion = assertion .and. all(ichar(getCharVec(Sample(1,1,1,1)),kind=IK) >= ichar(lower,kind=IK)) .and. all(ichar(getCharVec(Sample(1,1,1,1)),kind=IK) <= ichar(upper,kind=IK))
!        call test%assert(assertion, SK_"The procedure must be able to output scalar harvest within the specified limits.")
!#endif

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        subroutine report()
            integer(IK) :: i
            do i = 1, size(SamplePointer)
#if             SK_ENABLED
                assertion = assertion .and. SamplePointer(i) >= lower .and. SamplePointer(i) <= upper
#elif           IK_ENABLED
                assertion = assertion .and. SamplePointer(i) >= lower .and. SamplePointer(i) <= upper
#elif           LK_ENABLED
                assertion = assertion .and. (SamplePointer(i) .eqv. lower) .or. (SamplePointer(i) .eqv. upper)
#elif           RK_ENABLED
                assertion = assertion .and. SamplePointer(i) >= lower .and. SamplePointer(i) < upper
#elif           CK_ENABLED
                assertion = assertion .and. SamplePointer(i)%re >= lower%re .and. SamplePointer(i)%re < upper%re
                assertion = assertion .and. SamplePointer(i)%im >= lower%im .and. SamplePointer(i)%im < upper%im
#else
#error          "Unrecognized interface."
#endif
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    !SamplePointer(1:NPALL) => Sample
                    !do i = 1, NPALL
                        !if (.not. (CHECK_LOWER(SamplePointer(i)) .or. CHECK_UPPER(SamplePointer(i)))) then
                            write(test%disp%unit,"(*(g0,:,', '))")
                            write(test%disp%unit,"(*(g0,:,', '))") "SamplePointer(i)   ", SamplePointer(i)
                            write(test%disp%unit,"(*(g0,:,', '))")
#if                         SK_ENABLED
                            block
                            integer :: j, index(len(Sample))
                            do j = 1, len(Sample)
                                index(j) = ichar(SamplePointer(i)(j:j))
                            end do
                            write(test%disp%unit,"(*(g0,:,', '))")
                            write(test%disp%unit,"(*(g0,:,', '))") "ichar(SamplePointer(i))    ", index
                            write(test%disp%unit,"(*(g0,:,', '))") "ichar(lower), ichar(upper) ", ichar(lower(1:1)), ichar(upper(1:1))
                            write(test%disp%unit,"(*(g0,:,', '))")
                            end block
#elif                       !LK_ENABLED
                            write(test%disp%unit,"(*(g0,:,', '))")
                            write(test%disp%unit,"(*(g0,:,', '))") "lower, upper", lower, upper
                            write(test%disp%unit,"(*(g0,:,', '))")
                            exit
#endif
                        !end if
                    !end do
                    !nullify(SamplePointer)
                    ! LCOV_EXCL_STOP
                end if
            end do
        end subroutine
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_LOWER
#undef  CHECK_UPPER
#undef  GET_REAL
#undef  GET_IMAG