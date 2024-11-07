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
!>  This module contains implementations of the tests of the procedures under the generic interfaces
!>  [getPadded](@ref pm_arrayPad::getPadded),
!>  [setPadded](@ref pm_arrayPad::setPadded).
!>
!>  \todo
!>  \phigh The tests in this file still benefit from expansion and improvement.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getPadded_ENABLED || setPadded_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif
#if     SK_ENABLED && D0_ENABLED
#define GET_LBOUND(Array) 1_IK
#define GEN_UBOUND(Array) len(Array, kind = IK)
#define GET_SIZE(Array) len(Array, kind = IK)
#define GEN_LBOLD(lb) 1_IK
#define GEN_LBNEW(lb) 1_IK
#elif   getPadded_ENABLED
#define GET_LBOUND(Array) 1_IK
#define GEN_UBOUND(Array) size(Array, kind = IK)
#define GET_SIZE(Array) size(Array, kind = IK)
#define GEN_LBOLD(lb) 1_IK
#define GEN_LBNEW(lb) 1_IK
#else
#define GET_LBOUND(Array) lbound(Array, dim = 1, kind = IK)
#define GEN_UBOUND(Array) ubound(Array, dim = 1, kind = IK)
#define GET_SIZE(Array) size(Array, kind = IK)
#define GEN_LBOLD(lb) lb
#define GEN_LBNEW(lb) lb
#endif

#if     SK_ENABLED && D0_ENABLED
#define ALL
        character(:,SKG), allocatable   :: Array, arrayPadded
        character(1,SKG), parameter     :: lpfill = SKG_"/"
        character(1,SKG), parameter     :: rpfill = SKG_"*"
        character(1,SKG), parameter     :: lmfill = SKG_"-"
        character(1,SKG), parameter     :: rmfill = SKG_"+"
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: Array, arrayPadded
        character(2,SKG), parameter                 :: lpfill = SKG_"//"
        character(2,SKG), parameter                 :: rpfill = SKG_"**"
        character(2,SKG), parameter                 :: lmfill = SKG_"--"
        character(2,SKG), parameter                 :: rmfill = SKG_"++"
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: Array, arrayPadded
        integer(IKG)    , parameter                 :: lpfill = huge(1_IKG)
        integer(IKG)    , parameter                 :: rpfill = huge(1_IKG)
        integer(IKG)    , parameter                 :: lmfill = huge(1_IKG)
        integer(IKG)    , parameter                 :: rmfill = huge(1_IKG)
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: Array, arrayPadded
        logical(LKG)    , parameter                 :: lpfill = .false._LKG
        logical(LKG)    , parameter                 :: rpfill = .false._LKG
        logical(LKG)    , parameter                 :: lmfill = .false._LKG
        logical(LKG)    , parameter                 :: rmfill = .false._LKG
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: Array, arrayPadded
        complex(CKG)    , parameter                 :: lpfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
        complex(CKG)    , parameter                 :: rpfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
        complex(CKG)    , parameter                 :: lmfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
        complex(CKG)    , parameter                 :: rmfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: Array, arrayPadded
        real(RKG)       , parameter                 :: lpfill = huge(0._RKG)
        real(RKG)       , parameter                 :: rpfill = huge(0._RKG)
        real(RKG)       , parameter                 :: lmfill = huge(0._RKG)
        real(RKG)       , parameter                 :: rmfill = huge(0._RKG)
#else
#error  "Unrecognized interface."
#endif
        !integer(IK) :: sizepadded
        !integer(IK) :: sizeold, sizenew
        !integer(IK) :: lpsize, rpsize
        !integer(IK) :: lmsize, rmsize
        !integer(IK) :: lbcold, ubcold
        !integer(IK) :: lbcnew, ubcnew
        !integer(IK) :: lbold, ubold
        !integer(IK) :: lbnew, ubnew
        !logical(LK) :: menabled
        integer(IK) :: i, j, k

        !>  \bug
        !>  Avoid zero margin and pad sizes in the following because of the GNU gfortran bug as of 10.3.
        integer(IK) , parameter :: SizePad(2,3) = reshape ( [ 1_IK, 3_IK &
                                                            , 2_IK, 2_IK &
                                                            , 3_IK, 1_IK &
                                                            ], shape = shape(SizePad) )
        integer(IK) , parameter :: SizeMarg(2,3) = reshape( [ 1_IK, 3_IK &
                                                            , 2_IK, 2_IK &
                                                            , 2_IK, 1_IK &
                                                            ], shape = shape(SizePad) )
        integer(IK) , parameter :: SizeArray(3) =   [ 1_IK &
                                                    , 2_IK &
                                                    , 3_IK &
                                                    ] ! Avoid zero-sized arrays in the following because it messes up with the array lower bounds and resets it to 1 which causes the tests to wrongly fail.
#if     setPadded_ENABLED
        logical(LK) :: failed
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        do i = 1, size(SizeArray,1,IK)
            do j = 1, size(SizePad,2,IK)
                do k = 1, size(SizeMarg,2,IK)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill, lmsize = SizeMarg(1,k), rmsize = SizeMarg(2,k))
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill, lmsize = SizeMarg(1,k), rmsize = SizeMarg(2,k), lmfill = lmfill)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill, lmsize = SizeMarg(1,k), rmsize = SizeMarg(2,k), rmfill = rmfill)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill, lmsize = SizeMarg(1,k), rmsize = SizeMarg(2,k), lmfill = lmfill, rmfill = rmfill)
#if                 setPadded_ENABLED
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill, failed = failed)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill, lmsize = SizeMarg(1,k), rmsize = SizeMarg(2,k), failed = failed)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill, lmsize = SizeMarg(1,k), rmsize = SizeMarg(2,k), lmfill = lmfill, failed = failed)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill, lmsize = SizeMarg(1,k), rmsize = SizeMarg(2,k), rmfill = rmfill, failed = failed)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(1,j), rpsize = SizePad(2,j), lpfill = lpfill, rpfill = rpfill, lmsize = SizeMarg(1,k), rmsize = SizeMarg(2,k), lmfill = lmfill, rmfill = rmfill, failed = failed)
#endif
                end do
            end do
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(sizeOld, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)

            integer(IK)     , intent(in)            :: sizeOld
            integer(IK)     , intent(in)            :: lpsize, rpsize
            integer(IK)     , intent(in), optional  :: lmsize, rmsize
            logical(LK)                 , optional  :: failed
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG), intent(in)            :: lpfill, rpfill
            character(1,SKG), intent(in), optional  :: lmfill, rmfill
#elif       SK_ENABLED && D1_ENABLED
            character(2,SKG), intent(in)            :: lpfill, rpfill
            character(2,SKG), intent(in), optional  :: lmfill, rmfill
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in)            :: lpfill, rpfill
            integer(IKG)    , intent(in), optional  :: lmfill, rmfill
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in)            :: lpfill, rpfill
            logical(LKG)    , intent(in), optional  :: lmfill, rmfill
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in)            :: lpfill, rpfill
            complex(CKG)    , intent(in), optional  :: lmfill, rmfill
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in)            :: lpfill, rpfill
            real(RKG)       , intent(in), optional  :: lmfill, rmfill
#else
#error      "Unrecognized interface."
#endif
            logical(LK) :: menabled
            integer(IK) :: sizeNew, lmsize_def, rmsize_def
            integer(IK) :: lbp, ubp
            type :: OldNew_type
                integer(IK) :: old, new
            end type OldNew_type
            type(OldNew_type) :: lb, ub, lbc, ubc

            if (present(lmsize) .neqv. present(rmsize)) error stop "Internal ParaMonte Testing error occurred: `lmsize` and `rmsize` must be both present or both missing."
            menabled = present(lmsize) .and. present(rmsize)

            !>  \bug
            !>  GNU Fortran 10.3 cannot concatenate empty character array of length 2 with a non-empty character array of the same length.
            !>  Fortran runtime error: Different CHARACTER lengths (0/2) in array constructor
            if (present(lmsize) .and. present(rmsize)) then
                if (lmsize == 0_IK .and. rmsize == 0_IK) error stop "Internal ParaMonte Testing error occurred: GNU bug exception."
            end if

            lmsize_def = getOption(0_IK, lmsize)
            rmsize_def = getOption(0_IK, rmsize)

            assertion = .true._LK

            ! Enlarge and pad and empty array.

            call reset()

            call setUnifRand(lb%old, -10_IK, 10_IK)
            lb%old = GEN_LBOLD(lb%old)
            ub%old = lb%old + sizeOld - 1_IK
            lbc%old = lb%old
            ubc%old = ub%old

            sizeNew = sizeOld + lmsize_def + lpsize + rpsize + rmsize_def
            lb%new = lb%old
            ub%new = lb%new + sizeNew - 1_IK
            lbp = lb%new + lmsize_def
            ubp = ub%new - rmsize_def
            lbc%new = lbp + lpsize
            ubc%new = ubp - rpsize

#if         SK_ENABLED && D0_ENABLED
            allocate(character(sizeOld,SKG) :: Array)
            call setUnifRand(Array, repeat(SKG_"A",len(Array)), repeat(SKG_"Z",len(Array)))
            arrayPadded = genRepeat(lmsize_def,lmfill)//genRepeat(lpsize,lpfill)//Array//genRepeat(rpsize,rpfill)//genRepeat(rmsize_def,rmfill)
#else
            allocate(Array(lb%old : ub%old))
#if         SK_ENABLED && D1_ENABLED
            call setUnifRand(Array, SKG_"AA", SKG_"ZZ")
#elif       IK_ENABLED && D1_ENABLED
            call setUnifRand(Array, -100_IKG, +100_IKG)
#elif       LK_ENABLED && D1_ENABLED
            call setUnifRand(Array)
#elif       CK_ENABLED && D1_ENABLED
            call setUnifRand(Array, (-100._CKG,-500._CKG), (+100._CKG,+500._CKG))
#elif       RK_ENABLED && D1_ENABLED
            call setUnifRand(Array, -100._RKG, +100._RKG)
#endif
            allocate(arrayPadded(lb%new : ub%new))
            !>  \bug
            !>  Bypass the GNU 10.3 bug for concatenation of zero-sized character arrays.
            if (lmsize_def > 0_IK .and. rmsize_def > 0_IK) then
                arrayPadded(:) = [genRepeat(lmsize_def,lmfill), genRepeat(lpsize,lpfill), Array, genRepeat(rpsize,rpfill), genRepeat(rmsize_def,rmfill)]
            elseif (lmsize_def > 0_IK) then
                arrayPadded(:) = [genRepeat(lmsize_def,lmfill), genRepeat(lpsize,lpfill), Array, genRepeat(rpsize,rpfill)]
            elseif (rmsize_def > 0_IK) then
                arrayPadded(:) = [genRepeat(lpsize,lpfill), Array, genRepeat(rpsize,rpfill), genRepeat(rmsize_def,rmfill)]
            else
                arrayPadded(:) = [genRepeat(lpsize,lpfill), Array, genRepeat(rpsize,rpfill)]
            end if
#endif

#if         setPadded_ENABLED
            if (menabled) then
                call setPadded(Array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill, failed)
            else
                call setPadded(Array, lpsize, rpsize, lpfill, rpfill, failed)
            end if
#elif       getPadded_ENABLED
            if (menabled) then
                Array = getPadded(Array, lpsize, rpsize, lpfill, rpfill, lmsize, rmsize, lmfill, rmfill)
            else
                Array = getPadded(Array, lpsize, rpsize, lpfill, rpfill)
            end if
#else
#error      "Unrecognized interface."
#endif

            if (present(failed)) then
                assertion = assertion .and. .not. failed
                call report()
                call test%assert(assertion, SK_"Call to setPadded() must happen without failure.")
            end if

            assertion = assertion .and. GET_SIZE(Array) == GET_SIZE(arrayPadded)
            call report()
            call test%assert(assertion, SK_"Call to setPadded()/getPadded() must yield an array of proper size, with present(lmfill), present(rmfill) = "//getStr([present(lmfill), present(rmfill)]), int(__LINE__, IK))

            assertion = assertion .and. GET_LBOUND(Array) == GET_LBOUND(arrayPadded)
            call report()
            call test%assert(assertion, SK_"Call to setPadded()/getPadded() must properly set the lower bound of the output array, with present(lmfill), present(rmfill), present(failed) = "//getStr([present(lmfill), present(rmfill), present(failed)]), int(__LINE__, IK))

            assertion = assertion .and. GEN_UBOUND(Array) == GEN_UBOUND(arrayPadded)
            call report()
            call test%assert(assertion, SK_"Call to setPadded()/getPadded() must properly set the upper bound of the output array, with present(lmfill), present(rmfill), present(failed) = "//getStr([present(lmfill), present(rmfill), present(failed)]), int(__LINE__, IK))

            assertion = assertion .and. ALL(Array(lbp : ubp) IS_EQUAL arrayPadded(lbp : ubp))
            call report()
            call test%assert(assertion, SK_"Call to setPadded()/getPadded() must properly set the contents of the output array, with present(lmfill), present(rmfill) = "//getStr([present(lmfill), present(rmfill)]), int(__LINE__, IK))

            if (menabled .and. present(lmfill)) then
                assertion = assertion .and. ALL(Array(lb%new : lbp - 1_IK) IS_EQUAL arrayPadded(lb%new : lbp - 1_IK))
                call report()
                call test%assert(assertion, SK_"Call to setPadded()/getPadded() must properly fill the new left margin elements with `lmfill`", int(__LINE__, IK))
            end if

            if (menabled .and. present(rmfill)) then
                assertion = assertion .and. ALL(Array(ubp + 1_IK : ub%new) IS_EQUAL arrayPadded(ubp + 1_IK : ub%new))
                call report()
                call test%assert(assertion, SK_"Call to setPadded()/getPadded() must properly fill the new right margin elements with `lmfill`", int(__LINE__, IK))
            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(Array)) deallocate(Array)
            if (allocated(arrayPadded)) deallocate(arrayPadded)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure function genRepeat(count,fill) result(Array)
            integer(IK)     , intent(in)            :: count
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG), intent(in), optional  :: fill
            character(count,SKG)                    :: Array
            if (present(fill)) Array(:) = repeat(fill, count)
#else
#if         SK_ENABLED && D1_ENABLED
            character(2,SKG), intent(in), optional  :: fill
            character(2,SKG) :: Array(count)
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in), optional  :: fill
            integer(IKG)                            :: Array(count)
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in), optional  :: fill
            logical(LKG)                            :: Array(count)
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in), optional  :: fill
            complex(CKG)                            :: Array(count)
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in), optional  :: fill
            real(RKG)                               :: Array(count)
#else
#error      "Unrecognized interface."
#endif
            if (present(fill)) Array(:) = fill
#endif
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Array                      ", Array
                write(test%disp%unit,"(*(g0,:,', '))") "arrayPadded                ", arrayPadded
                write(test%disp%unit,"(*(g0,:,', '))") "GET_LBOUND(Array      )    ", GET_LBOUND(Array      )
                write(test%disp%unit,"(*(g0,:,', '))") "GET_LBOUND(arrayPadded)    ", GET_LBOUND(arrayPadded)
                write(test%disp%unit,"(*(g0,:,', '))") "GEN_UBOUND(Array      )    ", GEN_UBOUND(Array      )
                write(test%disp%unit,"(*(g0,:,', '))") "GEN_UBOUND(arrayPadded)    ", GEN_UBOUND(arrayPadded)
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GEN_UBOUND
#undef GET_LBOUND
#undef GEN_LBOLD
#undef GEN_LBNEW
#undef IS_EQUAL
#undef GET_SIZE
#undef ALL

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPaddedl_ENABLED || setPaddedl_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     SK_ENABLED && D0_ENABLED
#define GET_LBOUND(Array) 1_IK
#define GEN_UBOUND(Array) len(Array, kind = IK)
#define GET_SIZE(Array) len(Array, kind = IK)
#define GEN_LBOLD(lb) 1_IK
#define GEN_LBNEW(lb) 1_IK
#elif   getPaddedl_ENABLED
#define GET_LBOUND(Array) 1_IK
#define GEN_UBOUND(Array) size(Array, kind = IK)
#define GET_SIZE(Array) size(Array, kind = IK)
#define GEN_LBOLD(lb) 1_IK
#define GEN_LBNEW(lb) 1_IK
#else
#define GET_LBOUND(Array) lbound(Array, dim = 1, kind = IK)
#define GEN_UBOUND(Array) ubound(Array, dim = 1, kind = IK)
#define GET_SIZE(Array) size(Array, kind = IK)
#define GEN_LBOLD(lb) lb
#define GEN_LBNEW(lb) lb
#endif

#if     SK_ENABLED && D0_ENABLED
#define ALL
        character(:,SKG), allocatable   :: Array, arrayPadded
        character(1,SKG), parameter     :: lpfill = SKG_"/"
        character(1,SKG), parameter     :: lmfill = SKG_"-"
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: Array, arrayPadded
        character(2,SKG), parameter                 :: lpfill = SKG_"//"
        character(2,SKG), parameter                 :: lmfill = SKG_"--"
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: Array, arrayPadded
        integer(IKG)    , parameter                 :: lpfill = huge(1_IKG)
        integer(IKG)    , parameter                 :: lmfill = huge(1_IKG)
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: Array, arrayPadded
        logical(LKG)    , parameter                 :: lpfill = .false._LKG
        logical(LKG)    , parameter                 :: lmfill = .false._LKG
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: Array, arrayPadded
        complex(CKG)    , parameter                 :: lpfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
        complex(CKG)    , parameter                 :: lmfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: Array, arrayPadded
        real(RKG)       , parameter                 :: lpfill = huge(0._RKG)
        real(RKG)       , parameter                 :: lmfill = huge(0._RKG)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: i, j, k

        !>  \bug
        !>  Avoid zero margin and setPaddedl sizes in the following because of the GNU gfortran bug as of 10.3.
        integer(IK) , parameter :: SizePad(3) = [ 1_IK &
                                                , 2_IK &
                                                , 3_IK &
                                                ]
        integer(IK) , parameter :: SizeMarg(3)= [ 1_IK &
                                                , 2_IK &
                                                , 2_IK &
                                                ]
        integer(IK) , parameter :: SizeArray(3) =   [ 1_IK &
                                                    , 2_IK &
                                                    , 3_IK &
                                                    ] ! Avoid zero-sized arrays in the following because it messes up with the array lower bounds and resets it to 1 which causes the tests to wrongly fail.
#if     setPaddedl_ENABLED
        logical(LK) :: failed
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        do i = 1, size(SizeArray,1,IK)
            do j = 1, size(SizePad,1,IK)
                do k = 1, size(SizeMarg,1,IK)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(j), lpfill = lpfill)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(j), lpfill = lpfill, lmsize = SizeMarg(k))
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(j), lpfill = lpfill, lmsize = SizeMarg(k), lmfill = lmfill)
#if                 setPaddedl_ENABLED
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(j), lpfill = lpfill, failed = failed)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(j), lpfill = lpfill, lmsize = SizeMarg(k), failed = failed)
                    call runTestsWith(sizeOld = SizeArray(i), lpsize = SizePad(j), lpfill = lpfill, lmsize = SizeMarg(k), lmfill = lmfill, failed = failed)
#endif
                end do
            end do
        end do

    contains

        subroutine runTestsWith(sizeOld, lpsize, lpfill, lmsize, lmfill, failed)

            integer(IK)     , intent(in)            :: sizeOld
            integer(IK)     , intent(in)            :: lpsize
            integer(IK)     , intent(in), optional  :: lmsize
            logical(LK)                 , optional  :: failed
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG), intent(in)            :: lpfill
            character(1,SKG), intent(in), optional  :: lmfill
#elif       SK_ENABLED && D1_ENABLED
            character(2,SKG), intent(in)            :: lpfill
            character(2,SKG), intent(in), optional  :: lmfill
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in)            :: lpfill
            integer(IKG)    , intent(in), optional  :: lmfill
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in)            :: lpfill
            logical(LKG)    , intent(in), optional  :: lmfill
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in)            :: lpfill
            complex(CKG)    , intent(in), optional  :: lmfill
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)        , intent(in)            :: lpfill
            real(RKG)        , intent(in), optional  :: lmfill
#else
#error      "Unrecognized interface."
#endif
            integer(IK) :: sizeNew, lmsize_def
            integer(IK) :: lbp, ubp
            type :: OldNew_type
                integer(IK) :: old, new
            end type OldNew_type
            type(OldNew_type) :: lb, ub, lbc, ubc

            !>  \bug
            !>  GNU Fortran 10.3 cannot concatenate empty character array of length 2 with a non-empty character array of the same length.
            !>  Fortran runtime error: Different CHARACTER lengths (0/2) in array constructor
            if (present(lmsize)) then
                if (lmsize == 0_IK) error stop "Internal ParaMonte Testing error occurred: GNU bug exception."
            end if

            lmsize_def = getOption(0_IK, lmsize)

            assertion = .true._LK

            ! Enlarge and setPaddedl and empty array

            call reset()

            call setUnifRand(lb%old, -10_IK, 10_IK)
            lb%old = GEN_LBOLD(lb%old)
            ub%old = lb%old + sizeOld - 1_IK
            lbc%old = lb%old
            ubc%old = ub%old

            sizeNew = sizeOld + lmsize_def + lpsize
            lb%new = lb%old
            ub%new = lb%new + sizeNew - 1_IK
            lbp = lb%new + lmsize_def
            ubp = ub%new
            lbc%new = lbp + lpsize
            ubc%new = ubp

#if         SK_ENABLED && D0_ENABLED
            allocate(character(sizeOld,SKG) :: Array)
            call setUnifRand(Array, repeat(SKG_"A",len(Array)), repeat(SKG_"Z",len(Array)))
            arrayPadded = genRepeat(lmsize_def,lmfill)//genRepeat(lpsize,lpfill)//Array
#else
            allocate(Array(lb%old : ub%old))
#if         SK_ENABLED && D1_ENABLED
            call setUnifRand(Array, SKG_"AA", SKG_"ZZ")
#elif       LK_ENABLED && D1_ENABLED
            call setUnifRand(Array)
#elif       IK_ENABLED && D1_ENABLED
            call setUnifRand(Array, -100_IKG, +100_IKG)
#elif       CK_ENABLED && D1_ENABLED
            call setUnifRand(Array, (-100._CKG,-500._CKG), (+100._CKG,+500._CKG))
#elif       RK_ENABLED && D1_ENABLED
            call setUnifRand(Array, -100._RKG, +100._RKG)
#endif
            allocate(arrayPadded(lb%new : ub%new))
            !>  \bug
            !>  Bypass the GNU 10.3 bug for concatenation of zero-sized character arrays.
            if (lmsize_def > 0_IK) then
                arrayPadded(:) = [genRepeat(lmsize_def,lmfill), genRepeat(lpsize,lpfill), Array]
            else
                arrayPadded(:) = [genRepeat(lpsize,lpfill), Array]
            end if
#endif

#if         setPaddedl_ENABLED
            if (present(lmsize)) then
                call setPaddedl(Array, lpsize, lpfill, lmsize, lmfill, failed)
            else
                call setPaddedl(Array, lpsize, lpfill, failed)
            end if
#elif       getPaddedl_ENABLED
            if (present(lmsize)) then
                Array = getPaddedl(Array, lpsize, lpfill, lmsize, lmfill)
            else
                Array = getPaddedl(Array, lpsize, lpfill)
            end if
#else
#error      "Unrecognized interface."
#endif

            if (present(failed)) then
                assertion = assertion .and. .not. failed
                call report()
                call test%assert(assertion, desc = "Call to setPaddedl() must happen without failure.")
            end if

            assertion = assertion .and. GET_SIZE(Array) == GET_SIZE(arrayPadded)
            call report()
            call test%assert(assertion, desc = "Call to setPaddedl()/getPaddedl() must yield an array of proper size, with present(lmfill) = "//getStr([present(lmfill)]))

            assertion = assertion .and. GET_LBOUND(Array) == GET_LBOUND(arrayPadded)
            call report()
            call test%assert(assertion, desc = "Call to setPaddedl()/getPaddedl() must properly set the lower bound of the output array, with present(lmfill), present(failed) = "//getStr([present(lmfill), present(failed)]))

            assertion = assertion .and. GEN_UBOUND(Array) == GEN_UBOUND(arrayPadded)
            call report()
            call test%assert(assertion, desc = "Call to setPaddedl()/getPaddedl() must properly set the upper bound of the output array, with present(lmfill), present(failed) = "//getStr([present(lmfill), present(failed)]))

            assertion = assertion .and. ALL(Array(lbp : ubp) IS_EQUAL arrayPadded(lbp : ubp))
            call report()
            call test%assert(assertion, desc = "Call to setPaddedl()/getPaddedl() must properly set the contents of the output array, with present(lmfill) = "//getStr([present(lmfill)]))

            if (present(lmsize) .and. present(lmfill)) then
                assertion = assertion .and. ALL(Array(lb%new : lbp - 1_IK) IS_EQUAL arrayPadded(lb%new : lbp - 1_IK))
                call report()
                call test%assert(assertion, desc = "Call to setPaddedl()/getPaddedl() must properly fill the new left margin elements with `lmfill`")
            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(Array)) deallocate(Array)
            if (allocated(arrayPadded)) deallocate(arrayPadded)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure function genRepeat(count,fill) result(Array)
            integer(IK)     , intent(in)            :: count
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG), intent(in), optional  :: fill
            character(count,SKG) :: Array
            if (present(fill)) Array(:) = repeat(fill, count)
#else
#if         SK_ENABLED && D1_ENABLED
            character(2,SKG), intent(in), optional  :: fill
            character(2,SKG) :: Array(count)
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in), optional  :: fill
            logical(LKG)                            :: Array(count)
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in), optional  :: fill
            integer(IKG)                            :: Array(count)
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in), optional  :: fill
            complex(CKG)                            :: Array(count)
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in), optional  :: fill
            real(RKG)                               :: Array(count)
#else
#error      "Unrecognized interface."
#endif
            if (present(fill)) Array(:) = fill
#endif
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Array                      ", Array
                write(test%disp%unit,"(*(g0,:,', '))") "arrayPadded                ", arrayPadded
                write(test%disp%unit,"(*(g0,:,', '))") "GET_LBOUND(Array      )    ", GET_LBOUND(Array      )
                write(test%disp%unit,"(*(g0,:,', '))") "GET_LBOUND(arrayPadded)    ", GET_LBOUND(arrayPadded)
                write(test%disp%unit,"(*(g0,:,', '))") "GEN_UBOUND(Array      )    ", GEN_UBOUND(Array      )
                write(test%disp%unit,"(*(g0,:,', '))") "GEN_UBOUND(arrayPadded)    ", GEN_UBOUND(arrayPadded)
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GEN_UBOUND
#undef GET_LBOUND
#undef GEN_LBOLD
#undef GEN_LBNEW
#undef IS_EQUAL
#undef GET_SIZE
#undef ALL

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPaddedr_ENABLED || setPaddedr_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     SK_ENABLED && D0_ENABLED
#define GET_LBOUND(Array) 1_IK
#define GEN_UBOUND(Array) len(Array, kind = IK)
#define GET_SIZE(Array) len(Array, kind = IK)
#define GEN_LBOLD(lb) 1_IK
#define GEN_LBNEW(lb) 1_IK
#elif   getPaddedr_ENABLED
#define GET_LBOUND(Array) 1_IK
#define GEN_UBOUND(Array) size(Array, kind = IK)
#define GET_SIZE(Array) size(Array, kind = IK)
#define GEN_LBOLD(lb) 1_IK
#define GEN_LBNEW(lb) 1_IK
#else
#define GET_LBOUND(Array) lbound(Array, dim = 1, kind = IK)
#define GEN_UBOUND(Array) ubound(Array, dim = 1, kind = IK)
#define GET_SIZE(Array) size(Array, kind = IK)
#define GEN_LBOLD(lb) lb
#define GEN_LBNEW(lb) lb
#endif

#if     SK_ENABLED && D0_ENABLED
#define ALL
        character(:,SKG), allocatable   :: Array, arrayPadded
        character(1,SKG), parameter     :: rpfill = SKG_"/"
        character(1,SKG), parameter     :: rmfill = SKG_"-"
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: Array, arrayPadded
        character(2,SKG), parameter                 :: rpfill = SKG_"//"
        character(2,SKG), parameter                 :: rmfill = SKG_"--"
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: Array, arrayPadded
        integer(IKG)    , parameter                 :: rpfill = huge(1_IKG)
        integer(IKG)    , parameter                 :: rmfill = huge(1_IKG)
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: Array, arrayPadded
        logical(LKG)    , parameter                 :: rpfill = .false._LKG
        logical(LKG)    , parameter                 :: rmfill = .false._LKG
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: Array, arrayPadded
        complex(CKG)    , parameter                 :: rpfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
        complex(CKG)    , parameter                 :: rmfill = cmplx(huge(0._CKG), huge(0._CKG), kind = CKG)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: Array, arrayPadded
        real(RKG)       , parameter                 :: rpfill = huge(0._RKG)
        real(RKG)       , parameter                 :: rmfill = huge(0._RKG)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: i, j, k

        !>  \bug
        !>  Avoid zero margin and setPaddedr sizes in the following because of the GNU gfortran bug as of 10.3.
        integer(IK) , parameter :: SizePad(3) = [ 1_IK &
                                                , 2_IK &
                                                , 3_IK &
                                                ]
        integer(IK) , parameter :: SizeMarg(3)= [ 1_IK &
                                                , 2_IK &
                                                , 2_IK &
                                                ]
        integer(IK) , parameter :: SizeArray(3) =   [ 1_IK &
                                                    , 2_IK &
                                                    , 3_IK &
                                                    ] ! Avoid zero-sized arrays in the following because it messes up with the array lower bounds and resets it to 1 which causes the tests to wrongly fail.
#if     setPaddedr_ENABLED
        logical(LK) :: failed
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        do i = 1, size(SizeArray,1,IK)
            do j = 1, size(SizePad,1,IK)
                do k = 1, size(SizeMarg,1,IK)
                    call runTestsWith(sizeOld = SizeArray(i), rpsize = SizePad(j), rpfill = rpfill)
                    call runTestsWith(sizeOld = SizeArray(i), rpsize = SizePad(j), rpfill = rpfill, rmsize = SizeMarg(k))
                    call runTestsWith(sizeOld = SizeArray(i), rpsize = SizePad(j), rpfill = rpfill, rmsize = SizeMarg(k), rmfill = rmfill)
#if                 setPaddedr_ENABLED
                    call runTestsWith(sizeOld = SizeArray(i), rpsize = SizePad(j), rpfill = rpfill, failed = failed)
                    call runTestsWith(sizeOld = SizeArray(i), rpsize = SizePad(j), rpfill = rpfill, rmsize = SizeMarg(k), failed = failed)
                    call runTestsWith(sizeOld = SizeArray(i), rpsize = SizePad(j), rpfill = rpfill, rmsize = SizeMarg(k), rmfill = rmfill, failed = failed)
#endif
                end do
            end do
        end do

    contains

        subroutine runTestsWith(sizeOld, rpsize, rpfill, rmsize, rmfill, failed)

            integer(IK)     , intent(in)            :: sizeOld
            integer(IK)     , intent(in)            :: rpsize
            integer(IK)     , intent(in), optional  :: rmsize
            logical(LK)                 , optional  :: failed
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG), intent(in)            :: rpfill
            character(1,SKG), intent(in), optional  :: rmfill
#elif       SK_ENABLED && D1_ENABLED
            character(2,SKG), intent(in)            :: rpfill
            character(2,SKG), intent(in), optional  :: rmfill
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in)            :: rpfill
            integer(IKG)    , intent(in), optional  :: rmfill
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in)            :: rpfill
            logical(LKG)    , intent(in), optional  :: rmfill
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in)            :: rpfill
            complex(CKG)    , intent(in), optional  :: rmfill
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in)            :: rpfill
            real(RKG)       , intent(in), optional  :: rmfill
#else
#error      "Unrecognized interface."
#endif
            integer(IK) :: sizeNew, rmsize_def
            integer(IK) :: lbp, ubp
            type :: OldNew_type
                integer(IK) :: old, new
            end type OldNew_type
            type(OldNew_type) :: lb, ub, lbc, ubc

            !>  \bug
            !>  GNU Fortran 10.3 cannot concatenate empty character array of length 2 with a non-empty character array of the same length.
            !>  Fortran runtime error: Different CHARACTER lengths (0/2) in array constructor
            if (present(rmsize)) then
                if (rmsize == 0_IK) error stop "Internal ParaMonte Testing error occurred: GNU bug exception."
            end if

            rmsize_def = getOption(0_IK, rmsize)

            assertion = .true._LK

            ! Enlarge and setPaddedr and empty array.

            call reset()

            call setUnifRand(lb%old, -10_IK, 10_IK)
            lb%old = GEN_LBOLD(lb%old)
            ub%old = lb%old + sizeOld - 1_IK
            lbc%old = lb%old
            ubc%old = ub%old

            sizeNew = sizeOld + rmsize_def + rpsize
            lb%new = lb%old
            ub%new = lb%new + sizeNew - 1_IK
            lbp = lb%new
            ubp = ub%new - rmsize_def
            lbc%new = lbp
            ubc%new = ubp - rpsize

#if         SK_ENABLED && D0_ENABLED
            allocate(character(sizeOld,SKG) :: Array)
            call setUnifRand(Array, repeat(SKG_"A",len(Array,IK)), repeat(SKG_"Z",len(Array,IK)))
            arrayPadded = Array//genRepeat(rpsize,rpfill)//genRepeat(rmsize_def,rmfill)
#else
            allocate(Array(lb%old : ub%old))
#if         SK_ENABLED && D1_ENABLED
            call setUnifRand(Array, SKG_"AA", SKG_"ZZ")
#elif       LK_ENABLED && D1_ENABLED
            call setUnifRand(Array)
#elif       IK_ENABLED && D1_ENABLED
            call setUnifRand(Array, -100_IKG, +100_IKG)
#elif       CK_ENABLED && D1_ENABLED
            call setUnifRand(Array, (-100._CKG,-500._CKG), (+100._CKG,+500._CKG))
#elif       RK_ENABLED && D1_ENABLED
            call setUnifRand(Array, -100._RKG, +100._RKG)
#endif
            allocate(arrayPadded(lb%new : ub%new))
            !>  \bug
            !>  Bypass the GNU 10.3 bug for concatenation of zero-sized character arrays.
            if (rmsize_def > 0_IK) then
                arrayPadded(:) = [Array, genRepeat(rpsize,rpfill), genRepeat(rmsize_def,rmfill)]
            else
                arrayPadded(:) = [Array, genRepeat(rpsize,rpfill)]
            end if
#endif

#if         setPaddedr_ENABLED
            if (present(rmsize)) then
                call setPaddedr(Array, rpsize, rpfill, rmsize, rmfill, failed)
            else
                call setPaddedr(Array, rpsize, rpfill, failed)
            end if
#elif       getPaddedr_ENABLED
            if (present(rmsize)) then
                Array = getPaddedr(Array, rpsize, rpfill, rmsize, rmfill)
            else
                Array = getPaddedr(Array, rpsize, rpfill)
            end if
#else
#error      "Unrecognized interface."
#endif

            if (present(failed)) then
                assertion = assertion .and. .not. failed
                call report()
                call test%assert(assertion, desc = "Call to setPaddedr() must happen without failure.")
            end if

            assertion = assertion .and. GET_SIZE(Array) == GET_SIZE(arrayPadded)
            call report()
            call test%assert(assertion, desc = "Call to setPaddedr()/getPaddedr() must yield an array of proper size, with present(rmfill) = "//getStr([present(rmfill)]))

            assertion = assertion .and. GET_LBOUND(Array) == GET_LBOUND(arrayPadded)
            call report()
            call test%assert(assertion, desc = "Call to setPaddedr()/getPaddedr() must properly set the lower bound of the output array, with present(rmfill), present(failed) = "//getStr([present(rmfill), present(failed)]))

            assertion = assertion .and. GEN_UBOUND(Array) == GEN_UBOUND(arrayPadded)
            call report()
            call test%assert(assertion, desc = "Call to setPaddedr()/getPaddedr() must properly set the upper bound of the output array, with present(rmfill), present(failed) = "//getStr([present(rmfill), present(failed)]))

            assertion = assertion .and. ALL(Array(lbp : ubp) IS_EQUAL arrayPadded(lbp : ubp))
            call report()
            call test%assert(assertion, desc = "Call to setPaddedr()/getPaddedr() must properly set the contents of the output array, with present(rmfill) = "//getStr([present(rmfill)]))

            if (present(rmsize) .and. present(rmfill)) then
                assertion = assertion .and. ALL(Array(lb%new : lbp - 1_IK) IS_EQUAL arrayPadded(lb%new : lbp - 1_IK))
                call report()
                call test%assert(assertion, desc = "Call to setPaddedr()/getPaddedr() must properly fill the new right margin elements with `rmfill`")
            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(Array)) deallocate(Array)
            if (allocated(arrayPadded)) deallocate(arrayPadded)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure function genRepeat(count,fill) result(Array)
            integer(IK)     , intent(in)            :: count
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG), intent(in), optional  :: fill
            character(count,SKG) :: Array
            if (present(fill)) Array(:) = repeat(fill, count)
#else
#if         SK_ENABLED && D1_ENABLED
            character(2,SKG), intent(in), optional  :: fill
            character(2,SKG) :: Array(count)
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in), optional  :: fill
            logical(LKG)                            :: Array(count)
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in), optional  :: fill
            integer(IKG)                            :: Array(count)
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in), optional  :: fill
            complex(CKG)                            :: Array(count)
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in), optional  :: fill
            real(RKG)                               :: Array(count)
#else
#error      "Unrecognized interface."
#endif
            if (present(fill)) Array(:) = fill
#endif
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Array                      ", Array
                write(test%disp%unit,"(*(g0,:,', '))") "arrayPadded                ", arrayPadded
                write(test%disp%unit,"(*(g0,:,', '))") "GET_LBOUND(Array      )    ", GET_LBOUND(Array      )
                write(test%disp%unit,"(*(g0,:,', '))") "GET_LBOUND(arrayPadded)    ", GET_LBOUND(arrayPadded)
                write(test%disp%unit,"(*(g0,:,', '))") "GEN_UBOUND(Array      )    ", GEN_UBOUND(Array      )
                write(test%disp%unit,"(*(g0,:,', '))") "GEN_UBOUND(arrayPadded)    ", GEN_UBOUND(arrayPadded)
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GEN_UBOUND
#undef GET_LBOUND
#undef GEN_LBOLD
#undef GEN_LBNEW
#undef IS_EQUAL
#undef GET_SIZE
#undef ALL

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif