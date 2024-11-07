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

!>  \brief This file contains the implementations of the tests of module [pm_distanceMahal](@ref pm_distanceMahal).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_matrixInit, only: setMatDia
        use pm_distUnif, only: setUnifRand

#if     CK_ENABLED
        real(CK)    , parameter     :: EPS = epsilon(0._CK) * 100
        complex(CK) , allocatable   :: MahalSq_ref(:)
        complex(CK) , allocatable   :: invCov(:,:)
        complex(CK) , allocatable   :: mahalSq(:)
        complex(CK) , allocatable   :: diff(:)
        complex(CK) , allocatable   :: mean(:)
        complex(CK) , allocatable   :: X(:,:)
#elif   RK_ENABLED
        real(RK)    , parameter     :: EPS = epsilon(0._RK) * 100
        real(RK)    , allocatable   :: MahalSq_ref(:)
        real(RK)    , allocatable   :: invCov(:,:)
        real(RK)    , allocatable   :: mahalSq(:)
        real(RK)    , allocatable   :: diff(:)
        real(RK)    , allocatable   :: mean(:)
        real(RK)    , allocatable   :: X(:,:)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: ndim, npnt, i

        interface getDisMahalSq_ref
            procedure :: getDisMahalSq_ref_D0, getDisMahalSq_ref_D1, getDisMahalSq_ref_D2
        end interface

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ndim = 1_IK
        npnt = 5_IK
        call reset()
#if     CK_ENABLED
        call setMatDia(invCov, 1._RK, ndim, 0_IK, 0_IK)
        mean = (0._CK, 0._CK)
#elif   RK_ENABLED
        call setMatDia(invCov, 1._RK, ndim, 0_IK, 0_IK)
        mean = 0._RK
#endif

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(1,i), invCov(1,1), mean(1))
            mahalSq(i) = getDisMahalSq(X(1,i), invCov(1,1))
        end do
        call report()
        call test%assert(assertion, SK_"ndim = 1: The procedure must compute mahalSq correctly for scalar values.")

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(1,i), invCov(1,1), mean(1))
            mahalSq(i) = getDisMahalSq(X(1,i), invCov(1,1), mean(1))
        end do
        call report()
        call test%assert(assertion, SK_"ndim = 1: The procedure must compute mahalSq correctly for scalar values with zero mean.")

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(:,i), invCov, mean)
            mahalSq(i) = getDisMahalSq(X(:,i), invCov)
        end do
        call report()
        call test%assert(assertion, SK_"ndim = 1: The procedure must compute mahalSq correctly for vector values.")

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(:,i), invCov, mean)
            mahalSq(i) = getDisMahalSq(X(:,i), invCov, mean)
        end do
        call report()
        call test%assert(assertion, SK_"ndim = 1: The procedure must compute mahalSq correctly for vector values with zero mean.")

        call setUnifRand(X)
        MahalSq_ref = getDisMahalSq_ref(X, invCov, mean)
        mahalSq = getDisMahalSq(X, invCov)
        call report()
        call test%assert(assertion, SK_"ndim = 1: The procedure must compute mahalSq correctly for matrix values.")

        call setUnifRand(X)
        MahalSq_ref = getDisMahalSq_ref(X, invCov, mean)
        mahalSq = getDisMahalSq(X, invCov, mean)
        call report()
        call test%assert(assertion, SK_"ndim = 1: The procedure must compute mahalSq correctly for matrix values with zero mean.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ndim = 3_IK
        npnt = 5_IK
        call reset()
#if     CK_ENABLED
        call setMatDia(invCov, 1._RK, ndim, 0_IK, 0_IK)
        mean = (0._CK, 0._CK)
#elif   RK_ENABLED
        call setMatDia(invCov, 1._RK, ndim, 0_IK, 0_IK)
        mean = 0._RK
#endif

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(1,i), invCov(1,1), mean(1))
            mahalSq(i) = getDisMahalSq(X(1,i), invCov(1,1))
        end do
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for scalar values.")

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(1,i), invCov(1,1), mean(1))
            mahalSq(i) = getDisMahalSq(X(1,i), invCov(1,1), mean(1))
        end do
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for scalar values with zero mean.")

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(:,i), invCov, mean)
            mahalSq(i) = getDisMahalSq(X(:,i), invCov)
        end do
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for vector values.")

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(:,i), invCov, mean)
            mahalSq(i) = getDisMahalSq(X(:,i), invCov, mean)
        end do
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for vector values with zero mean.")

        call setUnifRand(X)
        MahalSq_ref = getDisMahalSq_ref(X, invCov, mean)
        mahalSq = getDisMahalSq(X, invCov)
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for matrix values.")

        call setUnifRand(X)
        MahalSq_ref = getDisMahalSq_ref(X, invCov, mean)
        mahalSq = getDisMahalSq(X, invCov, mean)
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for matrix values with zero mean.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ndim = 3_IK
        npnt = 5_IK
        call reset()
#if     CK_ENABLED
        call setMatDia(invCov, (+0.5_CK, -0.5_CK), ndim, 0_IK, 0_IK)
        mean = (0._CK, 0._CK)
#elif   RK_ENABLED
        call setMatDia(invCov, +0.5_RK, ndim, 0_IK, 0_IK)
        mean = 0._RK
#endif

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(1,i), invCov(1,1), mean(1))
            mahalSq(i) = getDisMahalSq(X(1,i), invCov(1,1))
        end do
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for scalar values.")

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(1,i), invCov(1,1), mean(1))
            mahalSq(i) = getDisMahalSq(X(1,i), invCov(1,1), mean(1))
            call report()
            call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for scalar values with zero mean.")
        end do

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(:,i), invCov, mean)
            mahalSq(i) = getDisMahalSq(X(:,i), invCov)
            call report()
            call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for vector values.")
        end do

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(:,i), invCov, mean)
            mahalSq(i) = getDisMahalSq(X(:,i), invCov, mean)
            call report()
            call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for vector values with zero mean.")
        end do

        call setUnifRand(X)
        MahalSq_ref = getDisMahalSq_ref(X, invCov, mean)
        mahalSq = getDisMahalSq(X, invCov)
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for matrix values.")

        call setUnifRand(X)
        MahalSq_ref = getDisMahalSq_ref(X, invCov, mean)
        mahalSq = getDisMahalSq(X, invCov, mean)
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for matrix values with zero mean.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ndim = 3_IK
        npnt = 5_IK
        call reset()
#if     CK_ENABLED
        call setMatDia(invCov, (+0.8_CK, -0.5_CK), ndim, 0_IK, 0_IK)
        mean = (2._CK, -5._CK)
#elif   RK_ENABLED
        call setMatDia(invCov, +0.8_RK, ndim, 0_IK, 0_IK)
        mean = 2._RK
#endif

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(1,i), invCov(1,1), mean(1))
                mahalSq(i) =     getDisMahalSq(X(1,i), invCov(1,1), mean(1))
        end do
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for scalar values with zero mean.")

        call setUnifRand(X)
        do i = 1, npnt
            MahalSq_ref(i) = getDisMahalSq_ref(X(:,i), invCov, mean)
            mahalSq(i) = getDisMahalSq(X(:,i), invCov, mean)
        end do
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for vector values with zero mean.")

        call setUnifRand(X)
        MahalSq_ref = getDisMahalSq_ref(X, invCov, mean)
        mahalSq = getDisMahalSq(X, invCov, mean)
        call report()
        call test%assert(assertion, SK_"The procedure must compute mahalSq correctly for matrix values with zero mean.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function getDisMahalSq_ref_D0(X, invCov, mean) result(mahalSq)
#if         RK_ENABLED
            real(RK)    , intent(in)            :: X, invCov, mean
            real(RK)                            :: mahalSq
#elif       CK_ENABLED
            complex(CK) , intent(in)            :: X, invCov, mean
            complex(CK)                         :: mahalSq
#else
#error      "Unrecognized interface."
#endif
            mahalSq = (x - mean)**2 * invCov
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function getDisMahalSq_ref_D1(X, invCov, mean) result(mahalSq)
#if         RK_ENABLED
            real(RK)    , intent(in)            :: X(:), invCov(:,:), mean(:)
            real(RK)                            :: mahalSq
#elif       CK_ENABLED
            complex(CK) , intent(in)            :: X(:), invCov(:,:), mean(:)
            complex(CK)                         :: mahalSq
#else
#error      "Unrecognized interface."
#endif
            mahalSq = dot_product(X - mean, matmul(invCov, X - mean))
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function getDisMahalSq_ref_D2(X, invCov, mean) result(mahalSq)
#if         RK_ENABLED
            real(RK)    , intent(in)            :: X(:,:), invCov(:,:), mean(:)
            real(RK)                            :: mahalSq(size(X,2))
#elif       CK_ENABLED
            complex(CK) , intent(in)            :: X(:,:), invCov(:,:), mean(:)
            complex(CK)                         :: mahalSq(size(X,2))
#else
#error      "Unrecognized interface."
#endif
            integer :: i
            do i = 1, size(X,2)
                mahalSq(i) = dot_product(X(:,i) - mean, matmul(invCov, X(:,i) - mean))
            end do
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#if 0
!        subroutine getDisMahalSq_ref(X, invCov, mean)
!#if         RK_ENABLED
!            real(RK)    , intent(in)            :: X(..), invCov(..), mean(..)
!#elif       CK_ENABLED
!            complex(CK) , intent(in)            :: X(..), invCov(..), mean(..)
!#else
!#error      "Unrecognized interface."
!#endif
!            integer(IK) :: i
!            select rank(X)
!                rank(0)
!                    select rank(invCov)
!                        rank(0)
!                            select rank(mean)
!                                rank(0)
!                                    MahalSq_ref(1) = (X - mean)**2 * invCov
!                                rank default
!                                    error stop ! LCOV_EXCL_LINE
!                            end select
!                        rank default
!                            error stop ! LCOV_EXCL_LINE
!                    end select
!                rank(1)
!                    select rank(invCov)
!                        rank(2)
!                            select rank(mean)
!                                rank(1)
!                                    MahalSq_ref(i) = dot_product(X - mean, matmul(invCov, (X - mean)))
!                                rank default
!                                    error stop ! LCOV_EXCL_LINE
!                            end select
!                        rank default
!                            error stop ! LCOV_EXCL_LINE
!                    end select
!                rank(2)
!                    select rank(invCov)
!                        rank(2)
!                            select rank(mean)
!                                rank(1)
!                                    MahalSq_ref(i) = dot_product(X(:,i) - mean, matmul(invCov, (X(:,i) - mean)))
!                                rank default
!                                    error stop ! LCOV_EXCL_LINE
!                            end select
!                        rank default
!                            error stop ! LCOV_EXCL_LINE
!                    end select
!                rank default
!                    error stop ! LCOV_EXCL_LINE
!            end select
!        end subroutine
!#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(MahalSq_ref)) deallocate(MahalSq_ref)
            if (allocated(invCov)) deallocate(invCov)
            if (allocated(mahalSq)) deallocate(mahalSq)
            if (allocated(mean)) deallocate(mean)
            if (allocated(X)) deallocate(X)
            allocate(invCov(ndim,ndim))
            allocate(MahalSq_ref(npnt))
            allocate(mahalSq(npnt))
            allocate(X(ndim,npnt))
            allocate(mean(ndim))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            integer(IK) :: i
            diff = mahalSq - MahalSq_ref
            assertion = assertion .and. all(abs(diff) <= EPS)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "mahalSq, MahalSq_ref"
#if             getDisMahalSq_CK_ENABLED
                write(test%disp%unit,"(4(g0,:,', '))") (mahalSq(i), MahalSq_ref(i), i = 1, size(mahalSq))
                write(test%disp%unit,"(*(g0,:,', '))") "diff"
                write(test%disp%unit,"(2(g0,:,', '))") (diff(i), i = 1, size(diff))
#elif           getDisMahalSq_RK_ENABLED
                write(test%disp%unit,"(2(g0,:,', '))") (mahalSq(i), MahalSq_ref(i), i = 1, size(mahalSq))
                write(test%disp%unit,"(*(g0,:,', '))") "diff"
                write(test%disp%unit,"(1(g0,:,', '))") (diff(i), i = 1, size(diff))
#endif
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
