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
!>  This include file contains the implementation of procedures in [pm_distUnif](@ref pm_distUnif).
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define array elements.
#if     setUnifCDF_ENABLED && D0_ENABLED
#define GET_ELEMENT(ARRAY) ARRAY
#elif   setUnifCDF_ENABLED && D1_ENABLED
#define GET_ELEMENT(ARRAY) ARRAY(i)
#elif   setUnifCDF_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define the optional RNG argument for the RNG procedures.
#if     getUnifRand_ENABLED || setUnifRand_ENABLED
#if     RNGD_ENABLED || RNGF_ENABLED
#define RNG
#elif   SM64_ENABLED || X256SSG_ENABLED || X256SSW_ENABLED
#define RNG rng,
#else
#error  "Unrecognized interface."
#endif
#endif
        !%%%%%%%%%%%%%%%%%
#if     getUnifCDF_ENABLED
        !%%%%%%%%%%%%%%%%%

#if     DD_ENABLED
        call setUnifCDF(cdf, x)
#elif   LU_ENABLED
        call setUnifCDF(cdf, x, lower, upper)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%
#elif   setUnifCDF_ENABLED
        !%%%%%%%%%%%%%%%%%

        ! Set default bounds.
#if     IK_ENABLED && DD_ENABLED
        integer(IKG), parameter :: lower = 0_IKG, upper = 1_IKG
#elif   CK_ENABLED && DD_ENABLED
        complex(CKG), parameter :: lower = 0._CKG, upper = 1._CKG
#elif   RK_ENABLED && DD_ENABLED
        real(RKG), parameter :: lower = 0._RKG, upper = 1._RKG
#elif   !LU_ENABLED
#error  "Unrecognized interface."
#endif
        ! Define the normalization constant for vector computations.
#if     D1_ENABLED
        integer(IK) :: i
#if     IK_ENABLED && LU_ENABLED
        real(RKG) :: inverseUpperMinusLowerPlusOne
        inverseUpperMinusLowerPlusOne = 1._RKG / real(upper - lower + 1_IKG, RKG)
#elif   CK_ENABLED && LU_ENABLED
        complex(CKG) :: inverseUpperMinusLowerPlusOne
        inverseUpperMinusLowerPlusOne%re = 1._CKG / (upper%re - lower%re)
        inverseUpperMinusLowerPlusOne%im = 1._CKG / (upper%im - lower%im)
#elif   RK_ENABLED && LU_ENABLED
        real(RKG) :: inverseUpperMinusLowerPlusOne
        inverseUpperMinusLowerPlusOne = 1._RKG / (upper - lower)
#elif   !DD_ENABLED
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, size(cdf, kind = IK) == size(x, kind = IK), SK_"setUnifCDF(): The condition `size(cdf) == size(x)` must hold. size(cdf), size(x) = "//getStr([size(cdf, kind = IK), size(x, kind = IK)])) ! fpp
#elif   !D0_ENABLED
#error  "Unrecognized interface."
#endif
#if     LU_ENABLED && CK_ENABLED
        CHECK_ASSERTION(__LINE__, lower%re <= upper%re .and. lower%im <= upper%im, SK_"setUnifCDF(): The conditions `lower%re <= upper%re .and. lower%im <= upper%im` must hold. lower, upper = "//getStr([lower, upper])) ! fpp
#elif   LU_ENABLED
        CHECK_ASSERTION(__LINE__, lower <= upper, SK_"setUnifCDF(): The condition `lower <= upper` must hold. lower, upper = "//getStr([lower, upper])) ! fpp
#elif   !DD_ENABLED
#error  "Unrecognized interface."
#endif
        ! Begin the computation.
#if     D1_ENABLED
        do concurrent(i = 1 : size(x, kind = IK))
#endif
            ! integer.
#if         IK_ENABLED
            if (GET_ELEMENT(x) < lower) then
                GET_ELEMENT(cdf) = 0._RKG
            elseif (GET_ELEMENT(x) < upper) then
#if             DD_ENABLED
                GET_ELEMENT(cdf) = real(GET_ELEMENT(x) + 1_IKG, RKG) * 0.5_RKG
#elif           LU_ENABLED && D0_ENABLED
                GET_ELEMENT(cdf) = real(GET_ELEMENT(x) + 1_IKG - lower, RKG) / real(upper - lower + 1_IKG, RKG)
#elif           LU_ENABLED && D1_ENABLED
                GET_ELEMENT(cdf) = real(GET_ELEMENT(x) + 1_IKG - lower, RKG) * inverseUpperMinusLowerPlusOne
#else
#error          "Unrecognized interface."
#endif
            else
                GET_ELEMENT(cdf) = 1._RKG
            end if
            ! real.
#elif       RK_ENABLED
            if (GET_ELEMENT(x) < lower) then
                GET_ELEMENT(cdf) = 0._RKG
            elseif (GET_ELEMENT(x) < upper) then
#if             DD_ENABLED
                GET_ELEMENT(cdf) = GET_ELEMENT(x)
#elif           LU_ENABLED && D0_ENABLED
                GET_ELEMENT(cdf) = (GET_ELEMENT(x) - lower) / (upper - lower)
#elif           LU_ENABLED && D1_ENABLED
                GET_ELEMENT(cdf) = (GET_ELEMENT(x) - lower) * inverseUpperMinusLowerPlusOne
#else
#error          "Unrecognized interface."
#endif
            else
                GET_ELEMENT(cdf) = 1._RKG
            end if
            ! complex.
#elif       CK_ENABLED
            ! real part.
            if (GET_ELEMENT(x)%re < real(lower,CKG)) then
                GET_ELEMENT(cdf)%re = 0._CKG ! fpp
            elseif (GET_ELEMENT(x)%re < real(upper,CKG)) then
#if             DD_ENABLED
                GET_ELEMENT(cdf)%re = GET_ELEMENT(x)%re
#elif           LU_ENABLED && D0_ENABLED
                GET_ELEMENT(cdf)%re = (GET_ELEMENT(x)%re - real(lower,CKG)) / (real(upper,CKG) - real(lower,CKG))
#elif           LU_ENABLED && D1_ENABLED
                GET_ELEMENT(cdf)%re = (GET_ELEMENT(x)%re - real(lower,CKG)) * inverseUpperMinusLowerPlusOne%re
#else
#error          "Unrecognized interface."
#endif
            else
                GET_ELEMENT(cdf)%re = 1._CKG ! fpp
            end if
            ! imaginary part.
            if (GET_ELEMENT(x)%im < aimag(lower)) then
                GET_ELEMENT(cdf)%im = 0._CKG ! fpp
            elseif (GET_ELEMENT(x)%im < aimag(upper)) then
#if             DD_ENABLED
                GET_ELEMENT(cdf)%im = GET_ELEMENT(x)%im
#elif           LU_ENABLED && D0_ENABLED
                GET_ELEMENT(cdf)%im = (GET_ELEMENT(x)%im - aimag(lower)) / (aimag(upper) - aimag(lower))
#elif           LU_ENABLED && D1_ENABLED
                GET_ELEMENT(cdf)%im = (GET_ELEMENT(x)%im - aimag(lower)) * inverseUpperMinusLowerPlusOne%im
#else
#error          "Unrecognized interface."
#endif
            else
                GET_ELEMENT(cdf)%im = 1._CKG ! fpp
            end if
#else
#error  "Unrecognized interface."
#endif
#if     D1_ENABLED
        end do
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   splitmix64_typer_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK64) :: count
        if (present(seed)) then
            rng%state = seed
        else ! By default, the seed is randomly initialized for every new instance of the RNG.
            call system_clock(count)
            rng%state = 324108011427370141_IK64 ! This must be present, otherwise GNU 10.3 uninitliazation warning bug.
            rng%state = ieor(rng%state, count)
        end if
        if (present(imageID)) then
            CHECK_ASSERTION(__LINE__, 0_IK < imageID, \
            SK_"@splitmix64_typer(): The condition `0 < imageID` must hold. imageID = "//getStr(imageID))
            rng%state = ieor(rng%state, int(imageID, IK64))
        else
            rng%state = ieor(rng%state, 1_IK64)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setStateNext_ENABLED && SM64_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Equivalent to unsigned hexadecimal integers
        ! Z"9e3779b97f4a7c15", Z"bf58476d1ce4e5b9", Z"94d049bb133111eb"
        integer(IK64), parameter :: TRIPLE(3) = [ -7046029254386353131_IK64 &
                                                , -4658895280553007687_IK64 &
                                                , -7723592293110705685_IK64 ]
        if (rng%state < 0_IK64) then
            if (rng%state < -huge(0_IK64) - 1_IK64 - TRIPLE(1)) then
                rng%state = rng%state + huge(0_IK64) + 1_IK64
                rng%state = rng%state + TRIPLE(1)
                rng%state = rng%state + huge(0_IK64) + 1_IK64
            else
                rng%state = rng%state + TRIPLE(1)
            end if
        else
            rng%state = rng%state + TRIPLE(1)
        end if
        !rng%state = rng%state + TRIPLE(1)
        rng%stream = rng%state
        rng%stream = ieor(rng%stream, shiftr(rng%stream, 30)) * TRIPLE(2)
        rng%stream = ieor(rng%stream, shiftr(rng%stream, 27)) * TRIPLE(3)
        rng%stream = ieor(rng%stream, shiftr(rng%stream, 31))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   xoshiro256ssg_typer_ENABLED || xoshiro256ssw_typer_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_kind, only: IKG => IK64, RKG => RK64
        integer(IKG) :: ijump
        type(splitmix64_type) :: rngsm
        rngsm = splitmix64_type(seed = seed)
        call setUnifRand(rngsm, rng%state)
        call setStateJump(rng)
        if (present(imageID)) then
            CHECK_ASSERTION(__LINE__, 0_IK < imageID, SK_"@xoshiro256ssw_typer(): The condition `0 < imageID` must hold. imageID = "//getStr(imageID))
            do ijump = 2_IKG, imageID
                call setStateJump(rng)
            end do
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setStateNext_ENABLED && (X256SSG_ENABLED || X256SSW_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK64) :: dum
        rng%stream = ishftc(rng%state(2) * 5_IK64, 7) * 9_IK64
        dum = shiftl(rng%state(2), 17)
        rng%state(3) =   ieor(rng%state(3), rng%state(1))
        rng%state(4) =   ieor(rng%state(4), rng%state(2))
        rng%state(2) =   ieor(rng%state(2), rng%state(3))
        rng%state(1) =   ieor(rng%state(1), rng%state(4))
        rng%state(3) =   ieor(rng%state(3), dum)
        rng%state(4) = ishftc(rng%state(4), 45)
#if     X256SSG_ENABLED
        rng%pos = 0_IK
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setStateJump_ENABLED && (X256SSG_ENABLED || X256SSW_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK64)   :: state(4)
        integer(IK)     :: istate, ibit
#if     DJ_ENABLED
#define JUMP xoshiro256ssJump128
#elif   AJ_ENABLED
        CHECK_ASSERTION(__LINE__, size(jump, 1, IK) == size(rng%state, 1, IK), \
        SK_": The condition `size(jump, 1) == size(rng%state, 1)` must hold. size(jump), size(rng%state) = "\
        //getStr([size(jump, 1, IK), size(rng%state, 1, IK)]))
#else
#error  "Unreocgnized interface."
#endif
        state = 0_IK64
        ! Jump 2^128 (or 2^64) steps ahead from the specified `jump`.
        do istate = 1_IK, size(rng%state, 1, IK)
            do ibit = 0_IK, int(bit_size(rng%stream), IK) - 1_IK ! 63_IK
                if (btest(JUMP(istate), ibit)) state = ieor(state, rng%state) ! fpp
                call setStateNext(rng)
            end do
        end do
        rng%state = state
#undef  JUMP

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getUnifRand_ENABLED && D0_ENABLED && LK_ENABLED && (RNGD_ENABLED || RNGF_ENABLED) && DD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real :: dummy
        call random_number(dummy)
        rand = logical(dummy < .5, LKG)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getUnifRand_ENABLED && D0_ENABLED && LK_ENABLED && (SM64_ENABLED || X256SSW_ENABLED) && DD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rand = logical(int(0, kind(rng%stream)) < rng%stream, LKG)
        call setStateNext(rng)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getUnifRand_ENABLED && D0_ENABLED && LK_ENABLED && X256SSG_ENABLED && DD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rand = logical(btest(rng%stream, rng%pos), LKG)
        rng%pos = rng%pos + 1_IK
        if (rng%pos == xoshiro256ssStreamBitSize) call setStateNext(rng)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getUnifRand_ENABLED && (RNGD_ENABLED || RNGF_ENABLED || SM64_ENABLED || X256SSG_ENABLED || X256SSW_ENABLED) && LU_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setUnifRand(RNG rand, lb, ub)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && D0_ENABLED && SK_ENABLED && DD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, randint
        do i = 1_IK, len(rand, IK)
            call setUnifRand(RNG randint, 1_IK, 127_IK)
            rand(i:i) = char(randint)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && D0_ENABLED && SK_ENABLED && LU_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lbLen, ubLen
        integer(IK) :: i, randint
        integer(IK) :: lbi, ubi
        lbLen = len(lb, IK)
        ubLen = len(ub, IK)
        CHECK_ASSERTION(__LINE__, lb <= ub, SK_"@setUnifRand(): The condition `lb <= ub` must hold. lb, ub = "//getStr([lb, ub]))
        CHECK_ASSERTION(__LINE__, lbLen == ubLen .or. lbLen == 1_IK .or. ubLen == 1_IK, SK_"@setUnifRand(): The condition `len(lb) == len(ub) .or. len(lb) == 1 .or. len(ub) == 1` must hold. len(lb), len(ub) = "//getStr([lbLen, ubLen]))
        CHECK_ASSERTION(__LINE__, lbLen == len(rand, IK) .or. lbLen == 1_IK, SK_"@setUnifRand(): The condition `len(lb) == len(rand) .or. len(lb) == 1`  must hold. len(lb), len(rand)  = "//getStr([lbLen, len(rand , IK)]))
        CHECK_ASSERTION(__LINE__, ubLen == len(rand, IK) .or. ubLen == 1_IK, SK_"@setUnifRand(): The condition `len(ub) == len(rand) .or. len(ub) == 1`  must hold. len(ub), len(rand)  = "//getStr([ubLen, len(rand , IK)]))
        if (1_IK < lbLen .and. 1_IK < ubLen) then
            do i = 1_IK, len(rand, IK)
                call setUnifRand(RNG randint, ichar(lb(i:i), IK), ichar(ub(i:i), IK))
                rand(i:i) = char(randint)
            end do
        elseif (1_IK == lbLen .and. 1_IK == ubLen) then
            lbi = ichar(lb, IK)
            ubi = ichar(ub, IK)
            do i = 1_IK, len(rand, IK)
                call setUnifRand(RNG randint, lbi, ubi)
                rand(i:i) = char(randint)
            end do
        elseif (1_IK == lbLen) then
            lbi = ichar(lb, IK)
            do i = 1_IK, len(rand, IK)
                call setUnifRand(RNG randint, lbi, ichar(ub(i:i), IK))
                rand(i:i) = char(randint)
            end do
        elseif (1_IK == ubLen) then
            ubi = ichar(ub, IK)
            do i = 1_IK, len(rand, IK)
                call setUnifRand(RNG randint, ichar(lb(i:i), IK), ubi)
                rand(i:i) = char(randint)
            end do
        elseif (len(rand, IK) /= 0_IK) then
            error stop "@setUnifRand(): Invalid user-specified input arguments. "& ! LCOV_EXCL_LINE
            //"Recompile & rerun with macro CHECK_ENABLED=1 for more information." ! LCOV_EXCL_LINE
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && D0_ENABLED && IK_ENABLED && (RNGD_ENABLED || RNGF_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Convert significant bit count to decimal precision.
        ! RKT: target real kind with at least `digits(rand)` bits.
        ! This is to ensure full coverage of the range of the specific-kind integer.
        integer, parameter :: RKT = selected_real_kind(floor(digits(rand) * log10(2._RKB)))
        integer, parameter :: RKG = merge(RKB, RKT, -6 < RKT .and. RKT < 0)
        real(RKG) :: temp
#if     DD_ENABLED
        !real(RKD) :: temp
        !call random_number(temp)
        ! \bug GNU Fortran compiler 10.3
        ! The following comment is a bug with `IKG = integer_kinds(5)`.
        !rand = floor(.5_RKD + temp, kind = IKG) ! rand = nint(temp, kind = IKG)
        integer(IKG), parameter :: lb = -huge(0_IKG), ub = +huge(0_IKG)
#elif   LU_ENABLED
        CHECK_ASSERTION(__LINE__, lb <= ub, SK_"@setUnifRand(): The condition `lb <= ub` must hold. lb, ub = "//getStr([lb, ub]))
#else
#error  "Unrecognized interface."
#endif
        if (lb + 1_IKG < ub) then
            call random_number(temp)
            ! early conversion to `real` avoids possible overflow with `huge` limits.
            ! rand = lb + int(temp * real(ub - lb + 1_IKG, kind(temp)), kind = IKG)
            rand = floor((1._RKG - temp) * real(lb, RKG) + temp * real(ub, RKG) + temp, kind = IKG)
        elseif (lb == ub) then
            rand = lb
        else
            call random_number(temp)
            if (temp < 0.5_RKG) then
                rand = lb
            else
                rand = ub
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && D0_ENABLED && IK_ENABLED && (SM64_ENABLED || X256SSG_ENABLED || X256SSW_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer     , parameter :: IKS = kind(rng%stream)
        integer(IK) , parameter :: streamBitSize = bit_size(rng%stream) ! this has to remain generic, supportive of all RNGs.
#if     DD_ENABLED
        integer(IK) , parameter :: randBitSize = bit_size(rand)
        integer(IK) , parameter :: streamBitExcess = streamBitSize - randBitSize
        integer(IK) , parameter :: randStreamBitSizeRatio = int(real(randBitSize) / real(streamBitSize), IKG) + 1_IK
        integer(IKS)            :: buffer(randStreamBitSizeRatio)
        integer(IK)             :: ibuf
#if     X256SSG_ENABLED
        ibuf = rng%pos + randBitSize
        if (ibuf < streamBitSize) then
            rand = int(ibits(rng%stream, rng%pos, randBitSize), IKG)
            rng%pos = ibuf
            return
        elseif (ibuf == streamBitSize) then
            rand = int(ibits(rng%stream, rng%pos, randBitSize), IKG)
            call setStateNext(rng)
            return
        end if
        ! Here, the `rand` kind is higher than kind(rng%stream).
        ! There is currently no neat solution in Fortran for greedy storage of the remaining bits.
        ! For now, accept the wasteful approach below for the greedy algorithm too.
        if (0 < rng%pos) call setStateNext(rng)
#else
        if (0_IK < streamBitExcess) then
            ! Both shifting and transfer below are possible solutions.
            ! The `abs` below, though redundant, bypasses the standard constraint.
            !rand = int(shiftr(rng%stream, abs(streamBitExcess)), IKG)
            rand = transfer(source = rng%stream, mold = rand)
            call setStateNext(rng)
            return
        elseif (0_IK == streamBitExcess) then
            rand = int(rng%stream, IKG)
            call setStateNext(rng)
            return
        end if
#endif
        do ibuf = 1_IK, randStreamBitSizeRatio
            buffer(ibuf) = rng%stream
            call setStateNext(rng)
        end do
        rand = transfer(source = buffer, mold = rand)
#elif   LU_ENABLED
        integer(IKG), parameter :: HUGE_IKG = huge(0_IKG)
        integer(IKG) :: scale, lbb, ubb, nzeros, nsignif, imask, temp, diff
        CHECK_ASSERTION(__LINE__, lb <= ub, SK_"@setUnifRand(): The condition `lb <= ub` must hold. lb, ub = "//getStr([lb, ub]))
        if (lb /= ub) then
            if (lb + 1_IKG == ub) then ! impossible range overflow.
#if             X256SSG_ENABLED
                if (btest(rng%stream, rng%pos)) then
                    rand = lb
                else
                    rand = ub
                end if
                rng%pos = rng%pos + 1_IK
                if (rng%pos == streamBitSize) call setStateNext(rng)
#else
                if (0_IKS < rng%stream) then
                    rand = lb
                else
                    rand = ub
                end if
                call setStateNext(rng)
#endif
                return
            end if
            if (lb < 0_IKG .and. 0_IKG < ub) then ! possible range overflow.
                if (HUGE_IKG + lb < ub) then ! overflowed.
                    lbb = -lb - HUGE_IKG - 1_IKG
                    ubb = +ub - HUGE_IKG - 1_IKG
                    scale = ubb + lbb
                    nzeros = 0_IKG
                else
                    scale = ub - lb
                    nzeros = int(leadz(scale), IKG)
                end if
                !scale = ub .uadd. -lb
            else ! impossible range overflow (assuming `lb < ub`).
                scale = ub - lb
                nzeros = int(leadz(scale), IKG)
            end if
            nsignif = bit_size(scale) - nzeros
            imask = shiftr(not(0_IKG), nzeros)
            loopTry: do
                call setUnifRand(rng, temp)
                rand = iand(temp, imask)
                if(rand <= scale) exit loopTry
                diff = nzeros
                loopReject: do
                    if(diff < nsignif) exit loopReject
                    temp = shiftr(temp, nsignif)
                    rand = iand(temp, imask)
                    if(rand <= scale) exit loopTry
                    diff = diff - nsignif
                end do loopReject
            end do loopTry
            rand = rand + lb
        else
            rand = lb
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && D0_ENABLED && LK_ENABLED && (RNGD_ENABLED || RNGF_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     DD_ENABLED
        real :: dummy
        call random_number(dummy)
        rand = logical(dummy < .5, kind = LKG)
#elif   LU_ENABLED
        real :: dummy
        if (lb .neqv. ub) then
            call random_number(dummy)
            rand = logical(dummy < .5, kind = LKG)
        else
            rand = lb
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && D0_ENABLED && LK_ENABLED && (SM64_ENABLED || X256SSG_ENABLED || X256SSW_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     DD_ENABLED && X256SSG_ENABLED
        integer(IK) , parameter :: streamBitSize = bit_size(rng%stream)
        rand = logical(btest(rng%stream, rng%pos), LKG)
        rng%pos = rng%pos + 1_IK
        if (rng%pos == streamBitSize) call setStateNext(rng)
#elif   DD_ENABLED && (SM64_ENABLED || X256SSW_ENABLED)
        rand = logical(int(0, kind(rng%stream)) < rng%stream, kind = LKG)
        call setStateNext(rng)
#elif   LU_ENABLED
        if (lb .neqv. ub) then
            call setUnifRand(rng, rand)
        else
            rand = lb
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && D0_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Curse you Intel bug.
#if     __INTEL_COMPILER
        real(CKG) :: temp(2)
#if     DD_ENABLED
        call setUnifRand(RNG temp)
#elif   LU_ENABLED
        if (lb /= ub) then
            call setUnifRand(RNG temp(1), lb%re, ub%re)
            call setUnifRand(RNG temp(2), lb%im, ub%im)
            !call setUnifRand(RNG temp(1), real(lb, CKG), real(ub, CKG))
            !call setUnifRand(RNG temp(2), aimag(lb), aimag(ub))
        else
            rand = lb
        end if
#else
#error  "Unrecognized interface."
#endif
        rand = cmplx(temp(1), temp(2), CKG)
#else
#if     DD_ENABLED
        call setUnifRand(RNG rand%re)
        call setUnifRand(RNG rand%im)
#elif   LU_ENABLED
        if (lb /= ub) then
            call setUnifRand(RNG rand%re, lb%re, ub%re)
            call setUnifRand(RNG rand%im, lb%im, ub%im)
            !call setUnifRand(RNG rand%re, real(lb, CKG), real(ub, CKG))
            !call setUnifRand(RNG rand%im, aimag(lb), aimag(ub))
        else
            rand = lb
        end if
#else
#error  "Unrecognized interface."
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && D0_ENABLED && RK_ENABLED && (RNGD_ENABLED || RNGF_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     DD_ENABLED
        call random_number(rand)
#elif   LU_ENABLED
        if (lb /= ub) then
            CHECK_ASSERTION(__LINE__, lb < ub, SK_"@setUnifRand(): The condition `lb <= ub` must hold. lb, ub = "//getStr([lb, ub]))
            ! The following loop ensures the random number falls in the half-open range `[lb, ub)`,
            ! even in the extreme case when `ub = lb + spacing(lb)`.
            do
                call random_number(rand)
                !rand = lb + rand * (ub - lb)
                ! The product expansion, although more expensive, avoids possible overflow with `huge` limits.
                rand = (1._RKG - rand) * lb + rand * ub
                ! The equality (instead of <) ensures NAN cases are handled gracefully without infinite loop.
                ! Hard lesson learned.
                !if (rand < ub) return
                if (ub <= rand) cycle
                exit
            end do
        else
            rand = lb
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && D0_ENABLED && RK_ENABLED && (SM64_ENABLED || X256SSG_ENABLED || X256SSW_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     DD_ENABLED
        integer     , parameter :: IKS = kind(rng%stream)
        integer     , parameter :: DIGITS_RKG = digits(rand)
        integer     , parameter :: DIGITS_RNG = digits(rng%stream)
!#if     X256SSG_ENABLED
        ! The following original approach is too costly.
        ! We will instead use the wasteful approach for the greedy.
        !real(RKG)   , parameter :: INVPDIG_RKG = 1._RKG / 2._RKG**DIGITS_RKG
        !integer                 :: remaining, remainders
        !remainders = DIGITS_RNG - DIGITS_RKG - rng%pos
        !if (0 < remainders) then
        !    rand = real(ibits(rng%stream, rng%pos, DIGITS_RKG), RKG) * INVPDIG_RKG
        !    rng%pos = rng%pos + DIGITS_RKG
        !else
        !    remainders = DIGITS_RNG - rng%pos
        !    rand = real(ibits(rng%stream, rng%pos, remainders), RKG) * 0.5_RKG**remainders
        !    call setStateNext(rng)
        !    remaining = DIGITS_RKG - remainders
        !    do
        !        remainders = remaining - DIGITS_RNG
        !        if (0 < remainders) then ! use full stream and cycle
        !            rand = (rand + real(abs(rng%stream), RKG)) * 0.5_RKG**DIGITS_RNG
        !            remaining = remainders
        !            call setStateNext(rng)
        !        else ! use part of stream and exit.
        !            rng%pos = -remainders
        !            rand = (rand + real(ibits(rng%stream, 0_IK, rng%pos), RKG)) * 0.5_RKG**rng%pos
        !            return
        !        end if
        !    end do
        !end if
!#else
        integer                 :: i
        integer     , parameter :: NUMBIT_RNG = bit_size(rng%stream) ! DIGITS_RNG + 1 ! The `1` makes up for the missing sign bit in the counting.
        integer     , parameter :: REPEAT_RNG = int(real(DIGITS_RKG) / real(DIGITS_RNG))
        integer     , parameter :: REMAINDERS = DIGITS_RKG - DIGITS_RNG * REPEAT_RNG
        ! bit-shifting or integer exponentiation overflows for exponent 63. Use real exponentation.
        real(RKG)   , parameter :: INV_POWREM = 1._RKG / 2._RKG**REMAINDERS ! real(shiftl(1_IKS, REMAINDERS), RKG) ! 1. / 2_IKS**REMAINDERS
        real(RKG)   , parameter :: INV_POWDIG = 1._RKG / 2._RKG**DIGITS_RNG ! real(shiftl(1_IKS, DIGITS_RNG), RKG) ! 1. / 2_IKS**DIGITS_RNG
#if     X256SSG_ENABLED
        if (DIGITS_RNG < REMAINDERS + rng%pos) call setStateNext(rng)
        rand = real(shiftr(rng%stream, NUMBIT_RNG - REMAINDERS + rng%pos), RKG) * INV_POWREM
#else
        rand = real(shiftr(rng%stream, NUMBIT_RNG - REMAINDERS), RKG) * INV_POWREM
#endif
        call setStateNext(rng)
        do i = 1, REPEAT_RNG ! exists only for higher-than-double precision real kinds.
            rand = (rand + real(shiftr(rng%stream, NUMBIT_RNG - DIGITS_RNG), RKG)) * INV_POWDIG
            call setStateNext(rng)
        end do
!#endif
#elif   LU_ENABLED
        if (lb /= ub) then
            CHECK_ASSERTION(__LINE__, lb < ub, SK_"@setUnifRand(): The condition `lb <= ub` must hold. lb, ub = "//getStr([lb, ub]))
            ! The following loop ensures the random number falls in the half-open range `[lb, ub)`,
            ! even in the extreme case when `ub = lb + spacing(lb)`.
            do
                call setUnifRand(rng, rand)
                ! The product expansion, although more expensive, avoids possible overflow with `huge` limits.
                rand = (1._RKG - rand) * lb + rand * ub
                !rand = lb + rand * (ub - lb)
                if (ub <= rand) cycle
                exit
            end do
        else
            rand = lb
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setUnifRand_ENABLED && !D0_ENABLED && (RNGD_ENABLED || RNGF_ENABLED || SM64_ENABLED || X256SSG_ENABLED || X256SSW_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the bounds.
#if     DD_ENABLED
#define BOUNDS
#elif   LU_ENABLED
#define BOUNDS, lb, ub
#else
#error  "Unrecognized interface."
#endif
#if     D1_ENABLED
        integer(IK) :: irow, nrow
        nrow = size(rand, 1, IK)
        do irow = 1, nrow
            call setUnifRand(RNG rand(irow) BOUNDS)
        end do
#elif   D2_ENABLED
        integer(IK) :: irow, nrow, icol, ncol
        nrow = size(rand, 1, IK)
        ncol = size(rand, 2, IK)
        do icol = 1, size(rand, 2, IK)
            do irow = 1, nrow
                call setUnifRand(RNG rand(irow, icol) BOUNDS)
            end do
        end do
#elif   D3_ENABLED
        integer(IK) :: irow, nrow, icol, ncol, idim, ndim
        nrow = size(rand, 1, IK)
        ncol = size(rand, 2, IK)
        ndim = size(rand, 3, IK)
        do idim = 1, ndim
            do icol = 1, ncol
                do irow = 1, nrow
                    call setUnifRand(RNG rand(irow, icol, idim) BOUNDS)
                end do
            end do
        end do
#else
#error  "Unrecognized interface."
#endif
#undef  BOUNDS

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_ELEMENT
#undef  RNG