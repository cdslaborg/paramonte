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
!>  This include file contains the implementations of the tests of procedures in [pm_matrixInitDia](@ref pm_matrixInitDia).
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_arraySpace, only: setLinSpace

        integer(IK)                 :: nrow, ncol
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        integer(IK)                 :: i
        integer(IKC), parameter     :: ZERO = 0_IKC
        integer(IKC), allocatable   :: Eye_ref(:,:)
        integer(IKC), allocatable   :: diff(:,:)
        integer(IKC), allocatable   :: Eye(:,:)
        integer(IKC), allocatable   :: Diag(:)
        integer(IKC)                :: offdiag
        integer(IKC)                :: upper
        integer(IKC)                :: lower
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        complex(CK) , parameter     :: ZERO = (0._CK, 0._CK)
        complex(CK) , allocatable   :: Eye_ref(:,:)
        complex(CK) , allocatable   :: diff(:,:)
        complex(CK) , allocatable   :: Eye(:,:)
        complex(CK) , allocatable   :: Diag(:)
        complex(CK)                 :: offdiag
        complex(CK)                 :: upper
        complex(CK)                 :: lower
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        real(RK)    , parameter     :: ZERO = 0._RK
        real(RK)    , allocatable   :: Eye_ref(:,:)
        real(RK)    , allocatable   :: diff(:,:)
        real(RK)    , allocatable   :: Eye(:,:)
        real(RK)    , allocatable   :: Diag(:)
        real(RK)                    :: offdiag
        real(RK)                    :: upper
        real(RK)                    :: lower
#else
#error  "Unrecognized interface."
#endif
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 1_IK
        ncol = 1_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = 1_IKC
        lower = 0_IKC
        upper = 0_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        Diag = (1._CK, 0._CK)
        lower = (0._CK, 0._CK)
        upper = (0._CK, 0._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        Diag = 1._RK
        lower = 0._RK
        upper = 0._RK
#endif

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag(1))
#elif   setMatEye_ENABLED
        call setMatEye(Eye)
        call report()
        call test%assert(assertion, SK_"setMatEye() must generate an identity matrix [nrow, ncol] = [1,1].")
        call setMatEye(Eye, Diag(1))
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate an identity matrix with scalar Diag and [nrow, ncol] = [1,1].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 3_IK
        ncol = 3_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = 2_IKC
        lower = 1_IKC
        upper = 1_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        Diag = (2._CK, 2._CK)
        lower = (1._CK, 1._CK)
        upper = (1._CK, 1._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        Diag = 2._RK
        lower = 1._RK
        upper = 1._RK
#endif
        offdiag = lower

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag(1), offdiag)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag(1), offdiag)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate an identity matrix with scalar Diag and offdiag and [nrow, ncol] = [3,3].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 4_IK
        ncol = 6_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = 2_IKC
        lower = 1_IKC
        upper = 1_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        Diag = (2._CK, 2._CK)
        lower = (1._CK, 1._CK)
        upper = (1._CK, 1._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        Diag = 2._RK
        lower = 1._RK
        upper = 1._RK
#endif
        offdiag = lower

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag(1), offdiag)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag(1), offdiag)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate an identity matrix with scalar Diag and offdiag and [nrow, ncol] = [4,6].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 6_IK
        ncol = 4_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = 2_IKC
        lower = 1_IKC
        upper = 1_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        Diag = (2._CK, 2._CK)
        lower = (1._CK, 1._CK)
        upper = (1._CK, 1._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        Diag = 2._RK
        lower = 1._RK
        upper = 1._RK
#endif
        offdiag = lower

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag(1), offdiag)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag(1), offdiag)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate an identity matrix with scalar Diag and offdiag and [nrow, ncol] = [6,4].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 4_IK
        ncol = 6_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = 2_IKC
        lower = 1_IKC
        upper = 1_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        Diag = (2._CK, 2._CK)
        lower = (1._CK, 1._CK)
        upper = (1._CK, 1._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        Diag = 2._RK
        lower = 1._RK
        upper = 1._RK
#endif

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag(1), lower, upper)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag(1), lower, upper)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate an identity matrix with scalar Diag and lower and upper and [nrow, ncol] = [4,6].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 6_IK
        ncol = 4_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = 2_IKC
        lower = 1_IKC
        upper = 1_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        Diag = (2._CK, 2._CK)
        lower = (1._CK, 1._CK)
        upper = (1._CK, 1._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        Diag = 2._RK
        lower = 1._RK
        upper = 1._RK
#endif

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag(1), lower, upper)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag(1), lower, upper)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate an identity matrix with scalar Diag and lower and upper and [nrow, ncol] = [6,4].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 3_IK
        ncol = 3_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = [( int(i,IKC), i = 1, min(nrow, ncol) )]
        lower = 0_IKC
        upper = 0_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        call setLinSpace(Diag, (1._CK, 1._CK), (3._CK, 3._CK))
        lower = (0._CK, 0._CK)
        upper = (0._CK, 0._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        call setLinSpace(Diag, 1._RK, 3._RK)
        lower = 0._RK
        upper = 0._RK
#endif

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate a diagonal matrix with vector Diag and [nrow, ncol] = [3,3].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 2_IK
        ncol = 3_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = [( int(i,IKC), i = 1, min(nrow, ncol) )]
        lower = 0_IKC
        upper = 0_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        call setLinSpace(Diag, (1._CK, 1._CK), (3._CK, 3._CK))
        lower = (0._CK, 0._CK)
        upper = (0._CK, 0._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        call setLinSpace(Diag, 1._RK, 3._RK)
        lower = 0._RK
        upper = 0._RK
#endif

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate a diagonal matrix with vector Diag and [nrow, ncol] = [2,3].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 5_IK
        ncol = 2_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = [( int(i,IKC), i = 1, min(nrow, ncol) )]
        lower = 0_IKC
        upper = 0_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        call setLinSpace(Diag, (1._CK, 1._CK), (3._CK, 3._CK))
        lower = (0._CK, 0._CK)
        upper = (0._CK, 0._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        call setLinSpace(Diag, 1._RK, 3._RK)
        lower = 0._RK
        upper = 0._RK
#endif

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate a diagonal matrix with vector Diag and [nrow, ncol] = [1,1].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 5_IK
        ncol = 2_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = [( int(i,IKC), i = 1, min(nrow, ncol) )]
        lower = 0_IKC
        upper = 0_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        call setLinSpace(Diag, (1._CK, 1._CK), (3._CK, 3._CK))
        lower = (0._CK, 0._CK)
        upper = (0._CK, 0._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        call setLinSpace(Diag, 1._RK, 3._RK)
        lower = 0._RK
        upper = 0._RK
#endif
        offdiag = lower

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag, offdiag)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag, offdiag)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate a diagonal matrix with vector Diag and offdiag and [nrow, ncol] = [5,2].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 3_IK
        ncol = 6_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = [( int(i,IKC), i = 1, min(nrow, ncol) )]
        lower = 0_IKC
        upper = 0_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        call setLinSpace(Diag, (1._CK, 1._CK), (3._CK, 3._CK))
        lower = (0._CK, 0._CK)
        upper = (0._CK, 0._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        call setLinSpace(Diag, 1._RK, 3._RK)
        lower = 0._RK
        upper = 0._RK
#endif
        offdiag = lower

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag, offdiag)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag, offdiag)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate a diagonal matrix with vector Diag and offdiag and [nrow, ncol] = [3,6].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 5_IK
        ncol = 2_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = [( int(i,IKC), i = 1, min(nrow, ncol) )]
        lower = 10_IKC
        upper = 20_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        call setLinSpace(Diag, (1._CK, 1._CK), (3._CK, 3._CK))
        lower = (10._CK, 10._CK)
        upper = (20._CK, 20._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        call setLinSpace(Diag, 1._RK, 3._RK)
        lower = 10._RK
        upper = 20._RK
#endif

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag, lower, upper)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag, lower, upper)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate a diagonal matrix with vector Diag and lower and upper and [nrow, ncol] = [5,2].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 2_IK
        ncol = 5_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = [( int(i,IKC), i = 1, min(nrow, ncol) )]
        lower = 10_IKC
        upper = 20_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        call setLinSpace(Diag, (1._CK, 1._CK), (3._CK, 3._CK))
        lower = (10._CK, 10._CK)
        upper = (20._CK, 20._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        call setLinSpace(Diag, 1._RK, 3._RK)
        lower = 10._RK
        upper = 20._RK
#endif

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag, lower, upper)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag, lower, upper)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate a diagonal matrix with vector Diag and lower and upper and [nrow, ncol] = [2,5].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nrow = 3_IK
        ncol = 3_IK
        call reset()
#if     getMatEye_IK_ENABLED || setMatEye_IK_ENABLED
        Diag = [( int(i,IKC), i = 1, min(nrow, ncol) )]
        lower = 10_IKC
        upper = 20_IKC
#elif   getMatEye_CK_ENABLED || setMatEye_CK_ENABLED
        call setLinSpace(Diag, (1._CK, 1._CK), (3._CK, 3._CK))
        lower = (10._CK, 10._CK)
        upper = (20._CK, 20._CK)
#elif   getMatEye_RK_ENABLED || setMatEye_RK_ENABLED
        call setLinSpace(Diag, 1._RK, 3._RK)
        lower = 10._RK
        upper = 20._RK
#endif

        call setMatEye_ref()

#if     getMatEye_ENABLED
        Eye = getMatEye(nrow, ncol, Diag, lower, upper)
#elif   setMatEye_ENABLED
        call setMatEye(Eye, Diag, lower, upper)
#else
#error  "Unrecognized interface."
#endif

        call report()
        call test%assert(assertion, SK_"The procedure must generate a diagonal matrix with vector Diag and lower and upper and [nrow, ncol] = [3,3].")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine setMatEye_ref()
            integer(IK) :: i
            do i = 1, min(nrow, ncol)
                Eye_ref(1:i-1,i) = upper
                Eye_ref(i,i) = Diag(i)
                Eye_ref(i+1:nrow,i) = lower
            end do
            do i = min(nrow, ncol) + 1, ncol
                Eye_ref(1:nrow,i) = upper
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(Eye_ref)) deallocate(Eye_ref)
            if (allocated(Diag)) deallocate(Diag)
            if (allocated(Eye)) deallocate(Eye)
            allocate(Diag(min(nrow, ncol)))
            allocate(Eye_ref(nrow, ncol))
            allocate(Eye(nrow, ncol))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            integer(IK) :: i,j
            assertion = assertion .and. all(Eye - Eye_ref == ZERO)
            if (test%traceable .and. .not. assertion) then
                diff = Eye - Eye_ref
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Eye, Eye_ref"
                write(test%disp%unit,"(4(g0,:,', '))")  ((Eye(i,j), Eye_ref(i,j), i = 1, size(Eye,1)), j = 1, size(Eye,2))
                write(test%disp%unit,"(*(g0,:,', '))") "diff"
                write(test%disp%unit,"(2(g0,:,', '))")  ((diff(i,j), i = 1, size(diff,1)), j = 1, size(diff,2))
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
