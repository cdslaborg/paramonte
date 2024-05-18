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
!>  This include file contains procedure implementations of [pm_polynomial](@ref pm_polynomial).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 3:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define zero.
#if     CK_ENABLED || CK_CK_ENABLED
#define TYPE_KIND complex(TKC)
        complex(TKC), parameter :: ZERO = (0._TKC, 0._TKC), ONE = (1._TKC, 0._TKC)
#elif   RK_ENABLED || RK_RK_ENABLED || RK_CK_ENABLED
#define TYPE_KIND real(TKC)
        real(TKC), parameter :: ZERO = 0._TKC, ONE = 1._TKC
#else
#error  "Unrecognized interface."
#endif
#if     setPolyRoot_ENABLED && Jen_ENABLED
!       Bypass the non-existence of `SIND` and `COSD` for compilers (other than Intel/GNU): !(INTEL_ENABLED || GNU_ENABLED)
#if     1
#define COSD(x) cos(x * acos(real(-1, kind(x))) / real(180, kind(x)))
#define SIND(x) sin(x * acos(real(-1, kind(x))) / real(180, kind(x)))
#endif
!       Enable robust calculations.
!       \todo `hypot` argument still needs work.
#if     0
        use pm_complexDiv, only: getDiv
#define GET_RATIO(x, y) getDiv(x, y)
#define GET_ABS(x) hypot(x%re, x%im)
#else
#define GET_RATIO(x, y) x / y
#define GET_ABS(x) abs(x)
#endif
#endif

        !%%%%%%%%%%%%%%%%%
#if     getPolyVal_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: i, degree
        degree = size(coef, 1, IK) - 1_IK
        CHECK_ASSERTION(__LINE__, 0_IK <= degree, SK_"@setPolyRoot(): The condition `0_IK < size(coef)` must hold. size(coef) = "//getStr(degree + 1_IK))
        poly = coef(degree)
        do i = degree - 1_IK, 0_IK, -1_IK
            poly = poly * x + coef(i)
        end do

        !%%%%%%%%%%%%%%%%%
#elif   getPolyAdd_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: lenLhs, lenRhs, lenMax
        lenLhs = size(lhs, 1, IK)
        lenRhs = size(rhs, 1, IK)
        do
            if (lenLhs == 0_IK) exit
            if (lhs(lenLhs) /= ZERO) exit
            lenLhs = lenLhs - 1_IK
        end do
        do
            if (lenRhs == 0_IK) exit
            if (rhs(lenRhs) /= ZERO) exit
            lenRhs = lenRhs - 1_IK
        end do
        if (0_IK < lenLhs .and. 0_IK < lenRhs) then
            lenMax = max(lenLhs, lenRhs)
            call setPolyAdd(add(1 : lenMax), lhs(1 : lenLhs), rhs(1 : lenRhs))
            add(lenMax + 1 :) = ZERO
        elseif (0_IK < lenLhs) then
            add(1 : lenLhs) = lhs(1 : lenLhs)
            add(lenLhs + 1:) = ZERO
        else
            add(1 : lenRhs) = rhs(1 : lenRhs)
            add(lenRhs + 1:) = ZERO
        end if

        !%%%%%%%%%%%%%%%%%
#elif   setPolyAdd_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: lenRhs, lenLhs
        lenLhs = size(lhs, 1, IK)
        lenRhs = size(rhs, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK < lenLhs, SK_"@setPolyAdd(): The condition `0 < size(lhs)` must hold. size(lhs) = "//getStr(lenLhs))
        CHECK_ASSERTION(__LINE__, 0_IK < lenRhs, SK_"@setPolyAdd(): The condition `0 < size(rhs)` must hold. size(rhs) = "//getStr(lenRhs))
        CHECK_ASSERTION(__LINE__, size(add, 1, IK) == max(lenLhs, lenRhs), \
        SK_"@setPolyAdd(): The condition `size(add) == max(lenLhs, lenRhs)` must hold. size(add), size(lhs), size(rhs) = "//\
        getStr([size(add, 1, IK), lenLhs, lenRhs]))
        CHECK_ASSERTION(__LINE__, lhs(lenLhs) /= ZERO, SK_"@setPolyAdd(): The condition `lhs(size(lhs)) /= 0.` must hold. lhs = "//getStr(lhs))
        CHECK_ASSERTION(__LINE__, rhs(lenRhs) /= ZERO, SK_"@setPolyAdd(): The condition `rhs(size(rhs)) /= 0.` must hold. rhs = "//getStr(rhs))
        if (lenLhs < lenRhs) then
            add(1 : lenLhs) = lhs(1 : lenLhs) + rhs(1 : lenLhs)
            add(lenLhs + 1 : lenRhs) = rhs(lenLhs + 1 : lenRhs)
        else
            add(1 : lenRhs) = lhs(1 : lenRhs) + rhs(1 : lenRhs)
            add(lenRhs + 1 : lenLhs) = lhs(lenRhs + 1 : lenLhs)
        end if

        !%%%%%%%%%%%%%%%%%
#elif   getPolySub_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: lenLhs, lenRhs, lenMax
        lenLhs = size(lhs, 1, IK)
        lenRhs = size(rhs, 1, IK)
        do
            if (lenLhs == 0_IK) exit
            if (lhs(lenLhs) /= ZERO) exit
            lenLhs = lenLhs - 1_IK
        end do
        do
            if (lenRhs == 0_IK) exit
            if (rhs(lenRhs) /= ZERO) exit
            lenRhs = lenRhs - 1_IK
        end do
        if (0_IK < lenLhs .and. 0_IK < lenRhs) then
            lenMax = max(lenLhs, lenRhs)
            call setPolySub(sub(1 : lenMax), lhs(1 : lenLhs), rhs(1 : lenRhs))
            sub(lenMax + 1 :) = ZERO
        elseif (0_IK < lenLhs) then
            sub(1 : lenLhs) = lhs(1 : lenLhs)
            sub(lenLhs + 1:) = ZERO
        else
            sub(1 : lenRhs) = rhs(1 : lenRhs)
            sub(lenRhs + 1:) = ZERO
        end if

        !%%%%%%%%%%%%%%%%%
#elif   setPolySub_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: lenRhs, lenLhs
        lenLhs = size(lhs, 1, IK)
        lenRhs = size(rhs, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK < lenLhs, SK_"@setPolySub(): The condition `0 < size(lhs)` must hold. size(lhs) = "//getStr(lenLhs))
        CHECK_ASSERTION(__LINE__, 0_IK < lenRhs, SK_"@setPolySub(): The condition `0 < size(rhs)` must hold. size(rhs) = "//getStr(lenRhs))
        CHECK_ASSERTION(__LINE__, size(sub, 1, IK) == max(lenLhs, lenRhs), \
        SK_"@setPolySub(): The condition `size(sub) == max(lenLhs, lenRhs)` must hold. size(sub), size(lhs), size(rhs) = "//\
        getStr([size(sub, 1, IK), lenLhs, lenRhs]))
        CHECK_ASSERTION(__LINE__, lhs(lenLhs) /= ZERO, SK_"@setPolySub(): The condition `lhs(size(lhs)) /= 0.` must hold. lhs = "//getStr(lhs))
        CHECK_ASSERTION(__LINE__, rhs(lenRhs) /= ZERO, SK_"@setPolySub(): The condition `rhs(size(rhs)) /= 0.` must hold. rhs = "//getStr(rhs))
        if (lenLhs < lenRhs) then
            sub(1 : lenLhs) =  lhs(1 : lenLhs) - rhs(1 : lenLhs)
            sub(lenLhs + 1 : lenRhs) = -rhs(lenLhs + 1 : lenRhs)
        else
            sub(1 : lenRhs) =  lhs(1 : lenRhs) + rhs(1 : lenRhs)
            sub(lenRhs + 1 : lenLhs) = lhs(lenRhs + 1 : lenLhs)
        end if

        !%%%%%%%%%%%%%%%%%%
#elif   getPolyMul_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenLhs, lenRhs
        lenLhs = size(lhs, 1, IK)
        lenRhs = size(rhs, 1, IK)
        do
            if (lenLhs == 0_IK) exit
            if (lhs(lenLhs) /= ZERO) exit
            lenLhs = lenLhs - 1_IK
        end do
        do
            if (lenRhs == 0_IK) exit
            if (rhs(lenRhs) /= ZERO) exit
            lenRhs = lenRhs - 1_IK
        end do
        if (0_IK < lenLhs .and. 0_IK < lenRhs) then
            call setPolyMul(mul(1 : lenLhs + lenRhs - 1), lhs(1 : lenLhs), rhs(1 : lenRhs))
        else
            mul = ZERO
        end if

        !%%%%%%%%%%%%%%%%%
#elif   setPolyMul_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: i, j
        integer(IK) :: lhsdeg, rhsdeg
        lhsdeg = size(lhs, 1, IK) - 1_IK
        rhsdeg = size(rhs, 1, IK) - 1_IK
        CHECK_ASSERTION(__LINE__, 0_IK <= lhsdeg, SK_"@setPolyMul(): The condition `0 < size(lhs)` must hold. size(lhs) = "//getStr(lhsdeg + 1_IK))
        CHECK_ASSERTION(__LINE__, 0_IK <= rhsdeg, SK_"@setPolyMul(): The condition `0 < size(rhs)` must hold. size(rhs) = "//getStr(rhsdeg + 1_IK))
        CHECK_ASSERTION(__LINE__, size(mul, 1, IK) == lhsdeg + rhsdeg + 1_IK, \
        SK_"@setPolyMul(): The condition `size(mul) == size(lhs) + size(rhs) - 1` must hold. size(mul), size(lhs), size(rhs) = "//\
        getStr([size(mul, 1, IK), lhsdeg + 1_IK, rhsdeg + 1_IK]))
        CHECK_ASSERTION(__LINE__, lhs(lhsdeg) /= ZERO, SK_"@setPolyMul(): The condition `lhs(size(lhs)) /= 0.` must hold. lhs = "//getStr(lhs))
        CHECK_ASSERTION(__LINE__, rhs(rhsdeg) /= ZERO, SK_"@setPolyMul(): The condition `rhs(size(rhs)) /= 0.` must hold. rhs = "//getStr(rhs))
        mul = ZERO
        do i = 0_IK, lhsdeg
            do j = 0_IK, rhsdeg
                mul(i + j) = mul(i + j) + lhs(i) * rhs(j)
            end do
        end do

        !%%%%%%%%%%%%%%%%%
#elif   setPolyDiv_ENABLED
        !%%%%%%%%%%%%%%%%%

        TYPE_KIND :: remainder, normfac
        TYPE_KIND, allocatable :: DividendTerm(:), DivisorTerm(:), QuotientTerm(:)
        integer(IK) :: i, degDividend, degDivisor, lenDivisor, lenDividend
        lenDividend = size(dividend, 1, IK)
        lenDivisor = size(divisor, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK < lenDivisor, SK_"@setPolyDiv(): The condition `0 < size(divisor)` must hold. size(divisor) = "//getStr(lenDivisor))
        CHECK_ASSERTION(__LINE__, 0_IK < lenDividend, SK_"@setPolyDiv(): The condition `0 < size(dividend)` must hold. size(dividend) = "//getStr(lenDividend))
        CHECK_ASSERTION(__LINE__, size(quorem, 1, IK) == lenDividend, SK_"@setPolyDiv(): The condition `size(quorem) == size(dividend)` must hold. size(quorem), size(dividend) = "//getStr([size(quorem, 1, IK), lenDividend]))
        CHECK_ASSERTION(__LINE__, dividend(lenDividend) /= ZERO, SK_"@setPolyDiv(): The condition `dividend(size(dividend)) /= 0.` must hold. dividend = "//getStr(dividend))
        CHECK_ASSERTION(__LINE__, divisor(lenDivisor) /= ZERO, SK_"@setPolyDiv(): The condition `divisor(size(divisor)) /= 0.` must hold. divisor = "//getStr(divisor))
        if (lenDivisor <= lenDividend) then
            if (lenDivisor > 2_IK) then
                allocate(DivisorTerm(lenDividend), QuotientTerm(lenDividend), source = ZERO)
                DivisorTerm(1:lenDivisor) = divisor
                degDividend = lenDividend - 1_IK
                degDivisor = lenDivisor - 1_IK
                DividendTerm = dividend
                lenQuo = 0_IK
                do
                    DivisorTerm = eoshift(DivisorTerm, degDivisor - degDividend)
                    QuotientTerm(degDividend - degDivisor + 1) = DividendTerm(degDividend + 1) / DivisorTerm(degDividend + 1)
                    DividendTerm = DividendTerm - DivisorTerm * QuotientTerm(degDividend - degDivisor + 1)
                    lenQuo = max(lenQuo, degDividend - degDivisor)
                    do
                        degDividend = degDividend - 1_IK
                        if (DividendTerm(degDividend + 1) /= ZERO) exit
                    end do
                    DivisorTerm(1 : lenDivisor) = divisor
                    DivisorTerm(lenDivisor + 1 : lenDividend) = ZERO
                    if (degDividend < degDivisor) exit
                end do
                lenQuo = lenQuo + 1_IK
                quorem(1 : lenQuo) = QuotientTerm(1 : lenQuo)
                quorem(lenQuo + 1 : lenDividend) = DividendTerm(1 : degDividend + 1)
                ! \bug gfortran as of version 12 fails to automatically deallocate heap allocations.
                deallocate(DividendTerm, DivisorTerm, QuotientTerm)
            elseif (lenDivisor == 2_IK) then
                if (divisor(2) /= ONE) then
                    normfac = ONE / divisor(2)
                    remainder = dividend(lenDividend) * normfac
                    do i = lenDividend - 1_IK, 1_IK, -1_IK
                        quorem(i) = remainder
                        remainder = (dividend(i) - divisor(1) * remainder) * normfac
                    end do
                else
                    remainder = dividend(lenDividend)
                    do i = lenDividend - 1_IK, 1_IK, -1_IK
                        quorem(i) = remainder
                        remainder = dividend(i) - divisor(1) * remainder
                    end do
                end if
                quorem(lenDividend) = remainder
                lenQuo = lenDividend - 1_IK
            else ! lenDivisor == 1_IK
                normfac = ONE / divisor(1)
                quorem = dividend * normfac
            end if
        else
            lenQuo = 0_IK
            quorem = dividend
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getPolyDiff_ENABLED || setPolyDiff_ENABLED) && Def_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        integer(IK), parameter :: order = 1
#if     setPolyDiff_ENABLED
        CHECK_ASSERTION(__LINE__, size(diff, 1, IK) == max(0_IK, size(coef, 1, IK) - order),\
        SK_"@setPolyDiff(): The condition `size(diff) == max(0, size(coef) - order)` must hold. size(diff), size(coef), order = "//\
        getStr([size(diff, 1, IK), size(coef, 1, IK), order]))
#elif   !getPolyDiff_ENABLED
#error  "Unrecognized interface."
#endif
        do concurrent(i = order : size(coef, 1, IK) - 1_IK)
            diff(i) = coef(i) * i
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getPolyDiff_ENABLED || setPolyDiff_ENABLED) && Ord_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenCoef
        lenCoef = size(coef, 1, IK)
#if     setPolyDiff_ENABLED
        CHECK_ASSERTION(__LINE__, size(diff, 1, IK) == max(0_IK, lenCoef - order), \
        SK_"@setPolyDiff(): The condition `size(diff) == max(0, size(coef) - order)` must hold. size(diff), size(coef), order = "//\
        getStr([size(diff, 1, IK), lenCoef, order]))
#elif   !getPolyDiff_ENABLED
#error  "Unrecognized interface."
#endif
        if (0_IK < order) then
            do concurrent(i = order : lenCoef - 1_IK)
                diff(i) = coef(i) * real(i, TKC)
            end do
        else
            CHECK_ASSERTION(__LINE__, 0_IK <= order, SK_"@setPolyDiff(): The condition `0 <= order` must hold. order = "//getStr(order))
            diff = coef
        end if

        !%%%%%%%%%%%%%%%%%
#elif   getPolyStr_ENABLED
        !%%%%%%%%%%%%%%%%%

        !> The maximum possible length of the string output from the functions under the generic interface `getStr`.
        !> Note that `NUMERIC_MAXLEN = 127` is very generous given that,
        !>  +   The maximum length of a complex of kind `real128` including signs, exponentiation, comma, and parentheses is `93` characters.
        !>      `(-1.189731495357231765085759326628007016E+4932,-1.189731495357231765085759326628007016E+4932)`
        !>  +   The maximum length of a real of kind `real128` including signs and exponentiation is `44` characters.
        !>      `-1.18973149535723176508575932662800702E+4932`
        !>  +   The maximum length of an integer of kind `int64` including signs `20` characters.
        !>      `-9223372036854775807`
        character(127, SK)          :: iomsg
        integer(IK)                 :: i, iostat, lenCoef, strlen
        integer(IK)     , parameter :: NUMERIC_MAXLEN = 127_IK
        character(*,SKC), parameter :: PRODSYM = SKC_""
        character(*,SKC), parameter :: VARNAME = SKC_"x"
        character(*,SKC), parameter :: POWSYM = SKC_"^"
        character(*,SKC), parameter :: FIXED = PRODSYM//VARNAME//POWSYM
#if     CK_ENABLED
#define GET_POLY_TERM(i) SKC_"(", coef(i)%re, SKC_",", coef(i)%im, SKC_")", FIXED, i
#define POLY_TERM_CONST GET_POLY_TERM(0)
#define GET_SIGNSYM(i) SKC_" + "
#elif   RK_ENABLED
        character(*,SKC), parameter :: SIGNSYM(-1:1) = [SKC_" - ", SKC_"   ", SKC_" + "]
#define POLY_TERM_CONST coef(0), FIXED, 0
#define GET_POLY_TERM(i) abs(coef(i)), FIXED, i
#define GET_SIGNSYM(i) SIGNSYM(int(sign(1.1_TKC, coef(i))))
#else
#error  "Unrecognized interface."
#endif
        lenCoef = size(coef, 1, IK)
        if (0 < lenCoef) then
            strlen = NUMERIC_MAXLEN * lenCoef
            allocate(character(strlen,SKC) :: str)
            loopExpandStr: do
                write(str, "(ss,*(g0))", iostat = iostat, iomsg = iomsg) POLY_TERM_CONST, (GET_SIGNSYM(i), GET_POLY_TERM(i), i = 1_IK, lenCoef - 1_IK)
                !write(*,*); write(*,*) is_iostat_eor(iostat), iostat, iomsg; write(*,*)
                if (iostat == 0_IK) then
                    str = getTrimmedTZ(str)
                    return
                elseif (is_iostat_eor(iostat)) then
                    call setResized(str) ! LCOV_EXCL_LINE
                    cycle loopExpandStr ! LCOV_EXCL_LINE
                else
                    error stop MODULE_NAME//& ! LCOV_EXCL_LINE
                    SK_"@getPolyStr(): A runtime error occurred while converting the polynomial coefficients to a polynomial expression: "//trim(iomsg) ! LCOV_EXCL_LINE
                end if
            end do loopExpandStr
        else
            str = SKC_""
        end if
#undef  POLY_TERM_CONST
#undef  GET_POLY_TERM
#undef  GET_SIGNSYM

        !%%%%%%%%%%%%%%%%%%
#elif   getPolyRoot_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: count
#if     Def_ENABLED
        type(eigen_type), parameter :: method = eigen_type()
#elif   !(Eig_ENABLED || Jen_ENABLED || Lag_ENABLED)
#error  "Unrecognized interface."
#endif
        allocate(root(size(coef, 1, IK) - 1_IK))
        call setPolyRoot(root, count, coef, method)
        root = root(1 : count)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPolyRoot_ENABLED && Eig_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC), parameter :: EPS = epsilon(0._TKC), BASE = radix(1._TKC), BASESQ = BASE * BASE
#if     RK_CK_ENABLED
#define GET_ABS(x) abs(x)
        real(TKC) :: p, q, t, w, x, y, zz
        integer(IK) :: k, l, m, ll, low, mm, mp2, na, niter, enm2
        logical(LK) :: notlas
#elif   CK_CK_ENABLED
#define GET_ABS(x) abs(x%re) + abs(x%im)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: i, j, degree
        logical(LK) :: normReductionEnabled
        TYPE_KIND :: workspace(size(root, 1, IK), size(root, 1, IK))
        real(TKC) :: r, s, c, scaleFac, g
        degree = size(coef, 1, IK) - 1_IK
        CHECK_ASSERTION(__LINE__, 0_IK < degree, SK_"@setPolyRoot(): The condition `1 < size(coef)` must hold. size(coef) = "//getStr(degree + 1_IK))
        CHECK_ASSERTION(__LINE__, size(root, 1, IK) == degree, SK_"@setPolyRoot(): The condition `size(root) == size(coef) - 1` must hold. size(root), size(coef) - 1 = "//getStr([size(root, 1, IK), degree]))
        CHECK_ASSERTION(__LINE__, all(shape(workspace, kind = IK) == degree), SK_"@setPolyRoot(): The condition `all(shape(workspace) == size(coef) - 1)` must hold. shape(workspace), size(coef) = "//getStr([shape(workspace, kind = IK), degree + 1_IK]))
        CHECK_ASSERTION(__LINE__, coef(degree) /= ZERO, SK_"@setPolyRoot(): The condition `coef(size(coef)) /= 0.` must hold. degree = "//getStr(coef(degree)))

        ! Gracefully capture the error and return if the highest polynomial coefficient is zero.

        if (coef(degree) == ZERO) then
            count = 0_IK ! degree + 1_IK
            return
        end if

        ! Build the count row of the upper Hessenberg Companion matrix.

        do concurrent(i = 1_IK : degree)
            workspace(1, i) = -coef(degree - i) / coef(degree)
        end do

        ! Extract any exact ZERO roots and set degree = degree of remaining polynomial.

        do j = degree, 1_IK, -1_IK
            if (workspace(1, j) /= ZERO) exit
            root(j) = cmplx(ZERO, kind = TKC) ! xxx
            degree = degree - 1_IK
        end do

        ! The special case of `degree == 0` or `degree == 1`.

        if (degree == 0_IK) then ! Gracefully capture the error and return if the degree of polynomial is zero.
            count = 0_IK
            return
        elseif (degree == 1_IK) then
            root(1) = workspace(1,1)
            count = 1_IK
            return
        end if

        ! Build rows 2 through `degree` of the Companion matrix.

        do i = 2_IK, degree
            do j = 1_IK, degree
                workspace(i, j) = ZERO
            end do
            workspace(i, i - 1) = ONE
        end do

        !   Balance the matrix.
        !   This is an adaption of the EISPACK subroutine `balanc` to the special case of a Companion matrix.
        !   The EISPACK balance is a translation of the Algol procedure balance, num. math. 13, 293-304(1969) by parlett and reinsch. handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

        !   Iterative loop for norm reduction

        loopNormReduction: do

            normReductionEnabled = .false.

            do i = 1_IK, degree

                !   r = sum of magnitudes in row i skipping diagonal.
                !   c = sum of magnitudes in col i skipping diagonal.

                if (i == 1_IK) then
                    r = GET_ABS(workspace(1, 2))
                    do j = 3, degree
                        r = r + GET_ABS(workspace(1, j))
                    end do
                    c = GET_ABS(workspace(2, 1))
                else
                    r = GET_ABS(workspace(i, i - 1))
                    c = GET_ABS(workspace(1, i))
                    if (i /= degree) c = c + GET_ABS(workspace(i + 1, i))
                end if

                ! Determine column scale factor, `scaleFac`.

                g = r / BASE
                scaleFac = real(ONE, TKC)
                s = c + r
                do
                    if (c < g) then
                        scaleFac = scaleFac * BASE
                        c = c * BASESQ
                        cycle
                    end if
                    exit
                end do
                g = r * BASE
                do
                    if (c >= g) then
                        scaleFac = scaleFac / BASE
                        c = c / BASESQ
                        cycle
                    end if
                    exit
                end do

                ! Will the factor `scaleFac` have a significant effect?

                if ((c + r) / scaleFac < .95_TKC * s) then ! yes, so do the scaling.

                    g = real(ONE, TKC) / scaleFac
                    normReductionEnabled = .true.

                    ! scale row `i`

                    if (i == 1_IK) then
                        do j = 1_IK, degree
                            workspace(1, j) = workspace(1, j) * g
                        end do
                    else
                        workspace(i, i - 1) = workspace(i, i - 1) * g
                    end if

                    ! Scale column `i`.

                    workspace(1, i) = workspace(1, i) * scaleFac
                    if (i /= degree) workspace(i + 1, i) = workspace(i + 1, i) * scaleFac

                end if

            end do

            if (normReductionEnabled) cycle loopNormReduction
            exit loopNormReduction

        end do loopNormReduction

#if     RK_CK_ENABLED
        !   ***************** QR eigenvalue algorithm ***********************
        !   This is the EISPACK subroutine hqr that uses the QR
        !   algorithm to compute all eigenvalues of an upper
        !   hessenberg matrix. original algol code was due to Martin,
        !   Peters, and Wilkinson, numer. math., 14, 219-231(1970).
        count = degree
        low = 1_IK
        t = ZERO
        !   ********** search for next eigenvalues **********
        loopNextEigen: do
            if (count < low) exit loopNextEigen
            niter = 0_IK
            na = count - 1_IK
            enm2 = na - 1_IK
            !   ********** look for single small sub-diagonal element for `l = count step -1` until `low` do -- **********
            do
                do ll = low, count
                    l = count + low - ll
                    if (l == low) exit
                    if (abs(workspace(l, l - 1)) <= EPS * (abs(workspace(l - 1, l - 1)) + abs(workspace(l, l)))) exit
                end do
                !   ********** form shift **********
                x = workspace(count, count)
                if (l == count) exit
                y = workspace(na,na)
                w = workspace(count,na) * workspace(na,count)
                if (l == na) then ! two roots found
                    p = (y - x) * .5_TKC
                    q = p * p + w
                    zz = sqrt(abs(q))
                    x = x + t
                    if (q < ZERO) then ! complex pair
                        root(na) = cmplx(x + p, zz, TKC)
                        root(count) = cmplx(x + p, -zz, TKC)
                    else ! pair of reals
                        zz = p + sign(zz, p)
                        root(na) = cmplx(x + zz, ZERO, TKC)
                        root(count) = root(na)
                        if (zz /= ZERO) root(count) = cmplx(x - w / zz, ZERO, TKC)
                    end if
                    count = enm2
                    cycle loopNextEigen
                end if
                ! increased iterations for quad-precision.
               !if (niter == 30_IK) exit loopNextEigen ! failed. No convergence to an eigenvalue after 30 iterations.
                if (niter == 90_IK) exit loopNextEigen ! failed. No convergence to an eigenvalue after 30 iterations.
                if (niter == 10_IK .or. niter == 20_IK) then ! form exceptional shift.
                    t = t + x
                    do i = low, count
                        workspace(i,i) = workspace(i,i) - x
                    end do
                    s = abs(workspace(count, na)) + abs(workspace(na, enm2))
                    x = .75_TKC * s
                    y = x
                    w = -.4375_TKC * s * s
                end if
                niter = niter + 1_IK
                ! Look for two consecutive small sub-diagonal elements. for `m = count - 2 step -1` until l do.
                do mm = l, enm2
                    m = enm2 + l - mm
                    zz = workspace(m, m)
                    r = x - zz
                    s = y - zz
                    p = (r * s - w) / workspace(m + 1,m) + workspace(m, m + 1)
                    q = workspace(m + 1, m + 1) - zz - r - s
                    r = workspace(m + 2, m + 1)
                    s = abs(p) + abs(q) + abs(r)
                    p = p / s
                    q = q / s
                    r = r / s
                    if (m == l) exit
                    if (abs(workspace(m, m - 1)) * (abs(q) + abs(r)) <= EPS * abs(p) * (abs(workspace(m - 1, m - 1)) + abs(zz) + abs(workspace(m + 1, m + 1)))) exit
                end do
                mp2 = m + 2_IK
                do i = mp2, count
                    workspace(i, i - 2) = ZERO
                    if (i == mp2) cycle
                    workspace(i, i - 3) = ZERO
                end do
                ! Double QR step involving rows l to count and columns m to count.
                do k = m, na
                    notlas = k /= na
                    if (k /= m) then
                        p = workspace(k, k - 1)
                        q = workspace(k + 1, k - 1)
                        r = ZERO
                        if (notlas) r = workspace(k + 2, k - 1)
                        x = abs(p) + abs(q) + abs(r)
                        if (x == ZERO) cycle
                        p = p / x
                        q = q / x
                        r = r / x
                    end if
                    s = sign(sqrt(p * p + q * q + r * r), p)
                    if (k /= m) then
                        workspace(k, k - 1) = -s * x
                    elseif (l /= m) then
                        workspace(k, k - 1) = -workspace(k, k - 1)
                    end if
                    p = p + s
                    x = p / s
                    y = q / s
                    zz = r / s
                    q = q / p
                    r = r / p
                    ! row modification.
                    do j = k, count
                        p = workspace(k, j) + q * workspace(k + 1, j)
                        if (notlas) then
                            p = p + r * workspace(k + 2, j)
                            workspace(k + 2, j) = workspace(k + 2, j) - p * zz
                        end if
                        workspace(k + 1, j) = workspace(k + 1, j) - p * y
                        workspace(k, j) = workspace(k, j) - p * x
                    end do
                    j = min(count, k + 3)
                    ! column modification.
                    do i = l, j
                        p = x * workspace(i,k) + y * workspace(i, k + 1)
                        if (notlas) then
                            p = p + zz * workspace(i, k + 2)
                            workspace(i, k + 2) = workspace(i, k + 2) - p * r
                        end if
                        workspace(i, k + 1) = workspace(i, k + 1) - p * q
                        workspace(i,k) = workspace(i,k) - p
                    end do
                end do
            end do
            !   ********** ONE root found **********
            root(count) = cmplx(x + t, ZERO, TKC)
            count = na
        end do loopNextEigen
        count = count + 1_IK
#elif   CK_CK_ENABLED
        call dcomqr(1_IK, degree, workspace(:, 1 : degree), root(1 : degree), count)
#else
#error  "Unrecognized interface."
#endif
        !   This block reverses the ordering of the identified roots, such that `root(1:count)` contains it all.
        !   \todo
        !   This ordering reversal could be merged with the rest of the algorithm to avoid the extra final reversal copy below.
        !   However, given that most polynomials are of low degree, the final reversal copy cost should generally be quite negligible.
        !   Because of the complexity of the rest of the algorithm, a final copy was the preferred approach over merging the reversal with the algorithm, for now.
        block
            complex(TKC) :: temp
            integer(IK) :: lenRoot
            lenRoot = size(root, 1, IK)
            do i = lenRoot, count, -1_IK
                temp = root(i)
                root(i) = root(lenRoot - i + 1_IK)
                root(lenRoot - i + 1_IK) = temp
            end do
            count = lenRoot - count + 1_IK
        end block
#if CK_CK_ENABLED
    contains
        !   This subroutine is a translation of a unitary analogue of the algol
        !   procedure comlr, num. math. 12, 369-376(1968) by Martin and Wilkinson.
        !   Handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
        !   the unitary analogue substitutes the QR algorithm of Francis
        !   (comp. jour. 4, 332-345(1962)) for the LR algorithm.
        !   This subroutine finds the eigenvalues of a complex
        !   upper hessenberg matrix by the QR method.
        !   `low` and `igh` are integers determined by the balancing subroutine `cbal`.
        !   if `cbal` has not been used, set `low = 1, igh = order`.
        !   `hessen` contains the complex upper hessenberg matrix.
        !   The lower triangle of `hessen` below the subdiagonal contains
        !   information about the unitary transformations used in the reduction by `corth`, if performed.
        !   ON OUTPUT:
        !       The upper hessenberg portions of `hessen` are destroyed.
        !       Therefore, they must be saved before calling `comqr`
        !       if subsequent calculation of eigenvectors is to be performed.
        !
        !       `root` contains the eigenvalues.
        !       If an error exit is made, the eigenvalues should be correct for indices `first+1,...,order`,
        !       `first` is set to,
        !           zero       for normal return,
        !           j          if the j-th eigenvalue has not been
        !                      determined after 30 iterations.
        !
        !   arithmetic is real except for the replacement of the algol
        !   procedure cdiv by complex division and use of the subroutines
        !   sqrt and cmplx in computing complex square roots.
        pure subroutine dcomqr(low, igh, hessen, eigen, first)
            complex(TKC)    , intent(out)   :: eigen(:)
            complex(TKC)    , intent(inout) :: hessen(:,:) ! hessen%im(nm, order),hessen%re(nm,order)
            integer(IK)     , intent(in)    :: low, igh
            integer(IK)     , intent(out)   :: first
            complex(TKC)                    :: s, t, x, y, z3, zz
            integer(IK)                     :: enm1, i, niter, j, l, ll, lp1, order
            real(TKC)                       :: norm
            order = size(eigen, 1, IK) ! The order of the Hessenberg matrix.
            if (low /= igh) then ! 180
                !   create real subdiagonal elements
                l = low + 1_IK
                do i = l, igh
                    if (hessen(i, i - 1)%im /= 0._TKC) then
                        norm = abs(hessen(i, i - 1))
                        y = hessen(i, i - 1) / norm
                        hessen(i, i - 1) = cmplx(norm, 0._TKC, TKC)
                        do j = i, igh
                            s%im = y%re * hessen(i, j)%im - y%im * hessen(i, j)%re
                            hessen(i, j)%re = y%re * hessen(i, j)%re + y%im * hessen(i, j)%im
                            hessen(i, j)%im = s%im
                        end do
                        do j = low, min(i + 1_IK, igh)
                            s%im = y%re * hessen(j, i)%im + y%im * hessen(j, i)%re
                            hessen(j, i)%re = y%re * hessen(j, i)%re - y%im * hessen(j, i)%im
                            hessen(j, i)%im = s%im
                        end do
                    end if
                end do
            end if
            !   Store roots isolated by cbal.
            do i = 1_IK, order
                if (low <= i .and. i <= igh) cycle
                eigen(i) = hessen(i, i)
            end do
            first = igh
            t = (0._TKC, 0._TKC)
            !   Search for next eigenvalue.
            loopNextEigen: do ! 220
                if (first < low) then ! convergence achieved.
                    first = low
                    return
                end if
                niter = 0_IK
                enm1 = first - 1_IK
                !   Look for single small sub-diagonal element for `l = first` step `-1` until `low`.
                loop240: do
                    do ll = low, first
                        l = first + low - ll
                        if (l == low) exit
                        if (abs(hessen(l, l - 1)%re) <= EPS * & ! LCOV_EXCL_LINE
                                                            ( abs(hessen(l - 1, l - 1)%re) & ! LCOV_EXCL_LINE
                                                            + abs(hessen(l - 1, l - 1)%im) & ! LCOV_EXCL_LINE
                                                            + abs(hessen(l, l)%re) & ! LCOV_EXCL_LINE
                                                            + abs(hessen(l, l)%im) & ! LCOV_EXCL_LINE
                                                            ) ) exit
                    end do
                    !   Form shift. 300
                    if (l == first) then ! 660 : a root found.
                        ! A root found ! 660
                        eigen(first) = hessen(first, first) + t
                        first = enm1
                        cycle loopNextEigen
                    end if
                    ! amir: increased maximum iteration to allow for quad-precision convergence.
                    if (niter == 90_IK) then ! No convergence to an eigenvalue after 30 iterations.
                        first = first + 1_IK
                        return
                    end if
                    if (niter == 10_IK .or. niter == 20_IK) then ! 320 Form exceptional shift.
                        s = cmplx(abs(hessen(first, enm1)%re) + abs(hessen(enm1, first - 2)%re), 0._TKC, TKC)
                    else
                        s = hessen(first, first)
                        x%re = hessen(enm1, first)%re * hessen(first, enm1)%re
                        x%im = hessen(enm1, first)%im * hessen(first, enm1)%re
                        if (x%re /= 0._TKC .or. x%im /= 0._TKC) then ! 340 ! xx z3 and zz can be replaced with a single variable `z`.
                            y = 0.5_TKC * (hessen(enm1, enm1) - s)
                            z3 = sqrt(cmplx(y%re**2 - y%im**2 + x%re, 2._TKC * y%re * y%im + x%im, TKC))
                            zz = z3
                            if (y%re * zz%re + y%im * zz%im < 0._TKC) zz = -zz ! 310
                            z3 = x / (y + zz)
                            s = s - z3
                        end if
                    end if
                    do i = low, first ! 340
                        hessen(i,i) = hessen(i,i) - s
                    end do
                    t = t + s
                    niter = niter + 1_IK
                    ! Reduce to triangle (rows).
                    lp1 = l + 1_IK
                    do i = lp1, first ! 500
                        s%re = hessen(i, i - 1)%re
                        hessen(i, i - 1)%re = 0._TKC
                        norm = sqrt(s%re**2 + hessen(i - 1, i - 1)%re**2 + hessen(i - 1, i - 1)%im**2)
                        x = hessen(i - 1, i - 1) / norm
                        eigen(i - 1) = x
                        hessen(i - 1, i - 1) = cmplx(norm, 0._TKC, TKC)
                        hessen(i, i - 1)%im = s%re / norm
                        do j = i, first ! 490
                            y = hessen(i - 1, j)
                            zz = hessen(i, j)
                            hessen(i - 1, j)%re = x%re * y%re + x%im * y%im + hessen(i, i - 1)%im * zz%re
                            hessen(i - 1, j)%im = x%re * y%im - x%im * y%re + hessen(i, i - 1)%im * zz%im
                            hessen(i, j)%re = x%re * zz%re - x%im * zz%im - hessen(i, i - 1)%im * y%re
                            hessen(i, j)%im = x%re * zz%im + x%im * zz%re - hessen(i, i - 1)%im * y%im
                        end do
                    end do
                    s%im = hessen(first, first)%im
                    if (s%im /= 0._TKC) then
                        norm = abs(cmplx(hessen(first, first)%re, s%im, TKC))
                        s%re = hessen(first, first)%re / norm
                        s%im = s%im / norm
                        hessen(first, first) = cmplx(norm, 0._TKC, TKC)
                    end if
                    ! Inverse operation (columns).
                    do j = lp1, first ! 540
                        x = eigen(j - 1)
                        do i = l, j
                            y = cmplx(hessen(i, j - 1)%re, 0._TKC, TKC)
                            zz = hessen(i, j)
                            if (i /= j) then ! 560
                                y%im = hessen(i, j - 1)%im
                                hessen(i, j - 1)%im = x%re * y%im + x%im * y%re + hessen(j, j - 1)%im * zz%im
                            end if
                            hessen(i, j - 1)%re = x%re * y%re - x%im * y%im + hessen(j, j - 1)%im * zz%re
                            hessen(i, j)%re = x%re * zz%re + x%im * zz%im - hessen(j, j - 1)%im * y%re
                            hessen(i, j)%im = x%re * zz%im - x%im * zz%re - hessen(j, j - 1)%im * y%im
                        end do
                    end do
                    if (s%im == 0._TKC) cycle loop240
                    do i = l, first
                        y = hessen(i,first)
                        hessen(i,first)%re = s%re * y%re - s%im * y%im
                        hessen(i,first)%im = s%re * y%im + s%im * y%re
                    end do
                end do loop240
            end do loopNextEigen
        end subroutine
#endif
#undef  GET_ABS

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPolyRoot_ENABLED && Jen_ENABLED && CK_CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! the program has been written to reduce the chance of overflow occurring.
        ! if it does occur, there is still a possibility that the zerofinder will work
        ! provided the overflowed quantity is replaced by a large number.

        ! Global variables.
        real(TKC)   , parameter     :: BASE = radix(0._TKC), EPS = epsilon(0._TKC), TIN = tiny(0._TKC), INF = huge(0._TKC)
        real(TKC)   , parameter     :: SINR = SIND(94._TKC), COSR = COSD(94._TKC) ! rotation by 94 degrees.
        real(TKC)   , parameter     :: MRE = 2._TKC * sqrt(2._TKC) * EPS ! error bound on complex multiplication.
        real(TKC)   , parameter     :: ARE = EPS ! error bound on complex addition.
        real(TKC)   , parameter     :: HALF_SQRT2 = 0.5_TKC * sqrt(2._TKC) ! cos(45d)
        complex(TKC), allocatable   :: h(:), qp(:), qh(:), sh(:), workspace(:)
        complex(TKC)                :: s, t, pv
        integer(IK)                 :: degree
        integer(IK)                 :: nn

        ! Local variables.
        complex(TKC)                :: temp
        real(TKC)   , parameter     :: HIGH = sqrt(INF), LOW = TIN / EPS, INV_LOG_BASE = 1._TKC / log(BASE)
        real(TKC)   , parameter     :: NEG_HALF_INV_LOG_BASE = -.5_TKC * INV_LOG_BASE
        real(TKC)                   :: xx, yy, xxx, bnd, max, min
        real(TKC)                   :: x, xm, f, dx, df, xni ! cauchy, noshft routines
        integer(IK)                 :: cnt1, cnt2, i, n, j, jj, nm1
        logical(LK)                 :: converged

        nn = size(coef, 1, IK)
        degree = nn - 1_IK
        xx = +HALF_SQRT2
        yy = -HALF_SQRT2

        CHECK_ASSERTION(__LINE__, 1_IK < nn, SK_"@setPolyRoot(): The condition `1 < size(coef)` must hold. size(coef) = "//getStr(nn))
        CHECK_ASSERTION(__LINE__, coef(degree) /= ZERO, SK_"@setPolyRoot(): The condition `coef(size(coef)) /= 0.` must hold. coef = "//getStr(coef))
        CHECK_ASSERTION(__LINE__, size(root, 1, IK) == degree, SK_"@setPolyRoot(): The condition `size(root) == size(coef) - 1` must hold. size(root), size(coef) - 1 = "//getStr([size(root, 1, IK), degree]))

        ! Algorithm fails if the leading coefficient is zero.
        if (coef(degree) == ZERO) then
            count = 0_IK
            return
        end if

        ! Allocate arrays.
#if     __GFORTRAN__
        if (allocated(workspace)) deallocate(workspace)
        if (allocated(qp)) deallocate(qp)
        if (allocated(qh)) deallocate(qh)
        if (allocated(sh)) deallocate(sh)
        if (allocated(h)) deallocate(h)
#endif
        allocate(h(nn), qp(nn), qh(nn), sh(nn), workspace(nn)) ! xxx this should be moved down once cleaning is complete.

        ! Extract any exact ZERO roots and set count = degree of remaining polynomial.
        do count = 1_IK, degree ! 10
            if (coef(count - 1) /= ZERO) exit
            root(count) = ZERO
            nn = nn - 1_IK
        end do

        ! Make a copy of the coefficients.
        do i = 1_IK, nn
            workspace(i) = coef(nn - i)
            sh(i)%re = GET_ABS(workspace(i)) ! modulus of coefficients of `p`.
        end do

        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute a scale factor to multiply the coefficients of the polynomial.
        ! The scaling is done to avoid overflow and to avoid undetected underflow interfering with the convergence criterion.
        ! The factor is a power of the BASE.

        ! Find largest and smallest moduli of coefficients.
        min = INF
        max = 0._TKC
        do i = 1_IK, nn
            if (sh(i)%re > max) max = sh(i)%re
            if (sh(i)%re /= 0._TKC .and. sh(i)%re < min) min = sh(i)%re
        end do

        ! Scale only if there are very large or very small components.
        if (LOW <= min .and. max <= HIGH) then
            bnd = 1._TKC
        else
            bnd = LOW / min
            if (bnd <= 1._TKC) then
                bnd = BASE ** int((log(max) + log(min)) * NEG_HALF_INV_LOG_BASE + .5_TKC)
            elseif (max < INF / bnd) then
                bnd = 1._TKC
            else
                bnd = BASE ** int(log(bnd) * INV_LOG_BASE + .5_TKC)
            end if
            if (bnd /= 1._TKC) then
                ! Scale the polynomial.
                do i = 1_IK, nn
                    workspace(i) = bnd * workspace(i)
                end do
            end if
        end if

        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Start the algorithm for one zero.
        loopIter: do ! 40

            if (nn <= 2) then
                ! Calculate the final zero and return.
                count = degree
                root(degree) = GET_RATIO(-workspace(2), workspace(1))
                return
            end if

            ! Calculate bnd, a lower bound on the modulus of the zeros.
            do i = 1, nn
                sh(i)%re = GET_ABS(workspace(i)) ! hypot
            end do

            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! call cauchy(nn, sh, bnd)
            ! Cauchy computes a lower bound on the moduli of the zeros of a polynomial -
            ! the real component of `sh` is the modulus of the coefficients.
            ! Compute upper estimate of bound.
            n = nn - 1_IK
            sh(nn)%re = -sh(nn)%re
            x = exp((log(-sh(nn)%re) - log(sh(1)%re)) / n)
            if (sh(n)%re /= 0._TKC) then
                ! if newton step at the origin is better, use it.
                xm = -sh(nn)%re / sh(n)%re
                if (xm < x) x = xm
            end if
            ! Chop the interval (0,x) until `f <= 0`.
            do ! 10
                xm = x * .1_TKC
                f = sh(1)%re
                do i = 2_IK, nn
                    f = f * xm + sh(i)%re
                end do
                if (f > 0._TKC) then
                    x = xm
                    cycle ! 10
                end if
                exit
            end do
            dx = x
            ! Do newton iteration until x converges to two decimal places.
            do ! 30
                if (abs(dx / x) > .005_TKC) then
                    sh(1)%im = sh(1)%re
                    do i = 2_IK, nn
                        sh(i)%im = sh(i - 1)%im * x + sh(i)%re
                    end do
                    f = sh(nn)%im
                    df = sh(1)%im
                    do i = 2_IK, n
                        df = df * x + sh(i)%im
                    end do
                    dx = f / df
                    x = x - dx
                    cycle
                end if
                exit
            end do
            bnd = x
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! The outer loop to control two major passes with different sequences of shifts.
            do cnt1 = 1_IK, 2_IK
                ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! First stage calculation, no shift.
                ! Compute the derivative polynomial as the initial `h`
                ! polynomial and compute `l1` no-shift `h` polynomials.
                ! call noshft(5_IK)
                n = nn - 1_IK
                nm1 = n - 1_IK
                do i = 1_IK, n
                    xni = nn - i
                    h(i) = xni * workspace(i) / real(n, TKC)
                end do
                do jj = 1_IK, 5_IK
                    if (GET_ABS(h(n)) > 10._TKC * EPS * GET_ABS(workspace(n))) then ! hypot
                        t = GET_RATIO(-workspace(nn), h(n))
                        do i = 1, nm1
                            j = nn - i
                            temp = h(j - 1)
                            h(j)%re = t%re * temp%re - t%im * temp%im + workspace(j)%re
                            h(j)%im = t%re * temp%im + t%im * temp%re + workspace(j)%im
                        end do
                        h(1) = workspace(1)
                    else
                        ! If the constant term is essentially zero, shift `h` coefficients.
                        do i = 1, nm1
                            j = nn - i
                            h(j) = h(j - 1)
                        end do
                        h(1) = ZERO
                    end if
                end do
                ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Inner loop to select a shift.
                do cnt2 = 1_IK, 9_IK
                    ! Shift is chosen with modulus `bnd` and amplitude rotated by 94 degrees from the previous shift.
                    xxx = COSR * xx - SINR * yy
                    yy = SINR * xx + COSR * yy
                    xx = xxx
                    s%re = bnd * xx
                    s%im = bnd * yy
                    ! Second stage calculation, fixed shift.
                    count = degree - nn + 2_IK
                    call fxshft(10_IK * cnt2, root(count), converged)
                    if (converged) then
                        ! The second stage jumps directly to the third stage iteration.
                        ! If successful the zero is stored and the polynomial deflated.
                        nn = nn - 1_IK
                        do i = 1_IK, nn
                            workspace(i) = qp(i)
                        end do
                        cycle loopIter ! 40
                    end if
                    ! if the iteration is unsuccessful another shift is chosen.
                end do
                ! if 9 shifts fail, the outer loop is repeated with another sequence of shifts.
            end do
            ! The root finder has failed on two major passes. Return empty handed.
            count = count - 1_IK
            return
        end do loopIter

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine fxshft(l2, z, converged)
            ! Computes l2 fixed-shift h polynomials and tests for convergence.
            ! Initiates a variable-shift iteration and returns with the approximate zero if successful.
            ! l2 - limit of fixed shift steps.
            ! z - the output approximate zero if `converged == .true._LK`.
            ! converged  - logical indicating convergence of stage 3 iteration.
            integer(IK) , intent(in)    :: l2
            complex(TKC), intent(out)   :: z
            logical(LK) , intent(out)   :: converged
            complex(TKC)                :: ot, svs
            logical(LK)                 :: test, pasd, bool
            integer(LK)                 :: i, j, n
            n = nn - 1_IK
            ! Evaluate `p` at `s`.
            call polyev(nn, s, workspace, qp, pv)
            test = .true._LK
            pasd = .false._LK
            ! calculate first t = -p(s)/h(s).
            call calct(bool)
            ! main loop for one second stage step.
            do j = 1_IK, l2
                ot = t
                ! Compute next `h` polynomial and new `t`.
                call nexth(bool)
                call calct(bool)
                z = s + t
                ! Test for convergence unless stage 3 has failed once or this is the last h polynomial.
                if (.not. (bool .or. .not. test .or. j == l2)) then
                    if (GET_ABS(t - ot) < .5_TKC * GET_ABS(z)) then
                        if (pasd) then
                            ! The weak convergence test has been passed twice, start the third stage iteration,
                            ! after saving the current h polynomial and shift.
                            do i = 1, n
                                sh(i) = h(i)
                            end do
                            svs = s
                            call vrshft(10_IK, z, converged)
                            if (converged) return
                            ! The iteration failed to converge. turn off testing and restore `h, s, pv, t`.
                            test = .false._LK
                            do i = 1_IK, n
                                h(i) = sh(i)
                            end do
                            s = svs
                            call polyev(nn, s, workspace, qp, pv)
                            call calct(bool)
                            cycle
                        end if
                        pasd = .true._LK
                    else
                        pasd = .false._LK
                    end if
                end if
            end do
            ! Attempt an iteration with final `h` polynomial from the second stage.
            call vrshft(10_IK, z, converged)
        end subroutine fxshft

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine vrshft(limss3, z, converged)
            ! carries out the third stage iteration.
            ! limss3:       limit of steps in stage 3.
            ! z:            on entry contains the initial iterate, if the iteration converges it contains the final iterate on exit.
            ! converged:    `.true.` if iteration converges.
            integer(IK) , intent(in)        :: limss3
            complex(TKC), intent(inout)     :: z
            logical(LK) , intent(out)       :: converged
            real(TKC)                       :: mp, ms, omp, relstp, r1, r2, tp
            logical(LK)                     :: b, bool
            integer(IK)                     :: i, j
            converged = .false._LK
            b = .false._LK
            s = z
            ! main loop for stage three.
            do i = 1_IK, limss3
                ! Evaluate p at s and test for convergence.
                call polyev(nn, s, workspace, qp, pv)
                mp = GET_ABS(pv)
                ms = GET_ABS(s)
                if (mp <= 20._TKC * getErrHorner(nn, qp, ms, mp)) then
                    ! Polynomial value is smaller in value than a bound on the error in evaluating p, terminate the iteration.
                    converged = .true._LK
                    z = s
                    return
                end if
                if (i /= 1_IK) then
                    if (.not. (b .or. mp < omp .or. relstp >= .05_TKC)) then
                        ! Iteration has stalled. probably a cluster of zeros.
                        ! Do 5 fixed shift steps into the cluster to force one zero to dominate.
                        tp = relstp
                        b = .true._LK
                        if (relstp < EPS) tp = EPS
                        r1 = sqrt(tp)
                        r2 = s%re * (1._TKC + r1) - s%im * r1
                        s%im = s%re * r1 + s%im * (1._TKC + r1)
                        s%re = r2
                        call polyev(nn, s, workspace, qp, pv)
                        do j = 1_IK, 5_IK
                            call calct(bool)
                            call nexth(bool)
                        end do
                        omp = INF
                        goto 20
                    end if
                    ! Exit if polynomial value increases significantly.
                    if (omp < .1_TKC * mp) return
                end if
                omp = mp
                ! Calculate the next iterate.
                20 call calct(bool)
                call nexth(bool)
                call calct(bool)
                if (.not. bool) then
                    relstp = GET_ABS(t) / GET_ABS(s)
                    s = s + t
                end if
            end do
        end subroutine vrshft

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine calct(bool)
            ! Computes  t = -p(s) / h(s).
            ! bool: logical(LK) , set true if h(s) is essentially zero.
            logical(LK) , intent(out)   :: bool
            complex(TKC)                :: hv
            integer(IK)                 :: n
            n = nn - 1_IK
            ! Evaluate h(s).
            call polyev(n, s, h, qh, hv)
            bool = GET_ABS(hv) <= 10._TKC * ARE * GET_ABS(h(n))
            if (.not. bool) then
                t = GET_RATIO(-pv, hv)
                return
            end if
            t = ZERO
        end subroutine calct

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine nexth(bool)
            ! Compute the next shifted `h` polynomial.
            ! bool: logical(LK) , if .true._LK h(s) is essentially zero
            logical(LK) , intent(in)    :: bool
            complex(TKC)                :: temp
            integer(IK)                 :: j,n
            n = nn - 1_IK
            if (.not. bool) then
                do j = 2_IK, n
                    temp = qh(j - 1)
                    h(j)%re = t%re * temp%re - t%im * temp%im + qp(j)%re
                    h(j)%im = t%re * temp%im + t%im * temp%re + qp(j)%im
                end do
                h(1) = qp(1)
                return
            end if
            ! If `h(s) == zero`, replace `h` with `qh`.
            do j = 2_IK, n
                h(j) = qh(j - 1)
            end do
            h(1) = ZERO
        end subroutine nexth

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure subroutine polyev(nn, s, workspace, q, pv)
            ! Evaluate a polynomial  p  at  s  by the horner recurrence placing the partial sums in q and the computed value in pv.
            integer(IK) , intent(in)                    :: nn
            complex(TKC), intent(in)                    :: s
            complex(TKC), intent(in)    , contiguous    :: workspace(:)
            complex(TKC), intent(out)   , contiguous    :: q(:)
            complex(TKC), intent(out)                   :: pv
            real(TKC)   :: tempo
            integer(IK) :: i
            q(1) = workspace(1)
            pv = q(1)
            do i = 2_IK, nn
                tempo = pv%re * s%re - pv%im * s%im + workspace(i)%re
                pv%im = pv%re * s%im + pv%im * s%re + workspace(i)%im
                pv%re = tempo
                q(i) = pv
            end do
        end subroutine polyev

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure function getErrHorner(nn, q, ms, mp) result(errHorner)
            ! Bound the error in evaluating the polynomial by the horner recurrence.
            ! q:    the partial sums.
            ! ms:   modulus of the point.
            ! mp:   modulus of polynomial value.
            integer(IK) , intent(in)    :: nn
            complex(TKC), intent(in)    :: q(:)
            real(TKC)   , intent(in)    :: ms
            real(TKC)   , intent(in)    :: mp
            real(TKC)                   :: errHorner
            real(TKC)                   :: e
            integer                     :: i
            e = GET_ABS(q(1)) * MRE / (ARE + MRE)
            do i = 1_IK, nn
                e = e * ms + GET_ABS(q(i))
            end do
            errHorner = e * (ARE + MRE) - mp * MRE
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPolyRoot_ENABLED && Jen_ENABLED && RK_CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Global variables.

        integer(IK)     , parameter     :: MAX_ITER = 30_IK
        integer         , parameter     :: RKR = TKC ! SELECTED_REAL_KIND(int(precision(0._TKC) / 3.))
        real(RKR)       , parameter     :: BASE = radix(0._RKR)
        real(RKR)       , parameter     :: EPS = epsilon(0._TKC) ! important to ask for the `TKC` precision.
        real(RKR)       , parameter     :: TIN = tiny(0._RKR)
        real(RKR)       , parameter     :: INF = huge(0._RKR)
        real(RKR)       , parameter     :: MRE = EPS ! error bound on complex multiplication.
        real(RKR)       , parameter     :: ARE = EPS ! The error bound on `+` operation.
        real(TKC)                       :: u, v, a, b, c, d, a1, a3, a7, e, f, g, h ! , a2, a6
        real(TKC)       , allocatable   :: qp(:), k(:), qk(:), svk(:), workspace(:)
        complex(TKC)                    :: s, sz, lz
        integer(IK)                     :: n, nn

        ! Local variables.

        real(RKR)                       :: max, min, xx, yy, xxx, x, sc, bnd, xm, ff, df, dx
        real(RKR)       , parameter     :: SINR = SIND(94._RKR), COSR = COSD(94._RKR) ! rotation by 94 degrees.
        real(RKR)       , parameter     :: LOW = TIN / EPS
        real(RKR)       , allocatable   :: pt(:)
        real(TKC)       , allocatable   :: temp(:)
        real(TKC)                       :: t, aa, bb, cc, factor
        integer(IK)                     :: cnt, nz, i, j, jj, nm1
        logical(LK)                     :: zerok, scalingEnabled

        nn = size(coef, 1, IK)
        count = nn - 1_IK ! degree of the polynomial.
        n = count

        CHECK_ASSERTION(__LINE__, 1_IK < nn, SK_"@setPolyRoot(): The condition `1 < size(coef)` must hold. size(coef) = "//getStr(nn))
        CHECK_ASSERTION(__LINE__, coef(count) /= 0._TKC, SK_"@setPolyRoot(): The condition `coef(size(coef)) /= 0.` must hold. coef = "//getStr(coef))
        CHECK_ASSERTION(__LINE__, size(root, 1, IK) == count, SK_"@setPolyRoot(): The condition `size(root) == size(coef) - 1` must hold. size(root), size(coef) - 1 = "//getStr([size(root, 1, IK), count]))

        ! Initialization of constants for shift rotation.
        xx = sqrt(.5_RKR)
        yy = -xx

        ! Algorithm fails if the leading coefficient is zero.
        if (coef(count) == 0._TKC) then
            count = 0_IK
            return
        end if

        ! Extract any exact ZERO roots and set count = degree of remaining polynomial.
        do j = 1_IK, count
            if (coef(j - 1) /= 0._TKC) exit
            root(j) = (0._TKC, 0._TKC)
            nn = nn - 1_IK
            n = n - 1_IK
            !count = count - 1_IK
        end do

        ! bypass gfortran bug for heap allocations.
#if     __GFORTRAN__
        if (allocated(k))   deallocate(k)
        if (allocated(pt))  deallocate(pt)
        if (allocated(qp))  deallocate(qp)
        if (allocated(qk))  deallocate(qk)
        if (allocated(svk)) deallocate(svk)
        if (allocated(temp)) deallocate(temp)
        if (allocated(workspace)) deallocate(workspace)
#endif
        allocate(k(nn), pt(nn), qp(nn), qk(nn), svk(nn), temp(nn), workspace(nn))

        ! Make a copy of the coefficients.
        do j = 1_IK, nn
            workspace(j) = coef(nn - j)
        end do

        !write(*,*) "nn, workspace", nn, workspace

        ! Start the algorithm for one zero.
        loopIter: do ! 30

            if (n <= 2_IK) then
                if (n < 1_IK) return
                ! Compute the final zero or pair of zeros.
                if (n /= 2_IK) then
                    root(count)%re = -workspace(2) / workspace(1)
                    root(count)%im = 0._TKC
                    return
                end if
                call quad(workspace(1), workspace(2), workspace(3), root(count - 1), root(count))
                return
            end if

            ! Find largest and smallest moduli of coefficients.
            max = 0._RKR
            min = INF
            do i = 1_IK, nn
                x = abs(real(workspace(i), RKR))
                if (x > max) max = x
                if (x /= 0._RKR .and. x < min) min = x
            end do

            ! Scale if there are large or very small coefficients.
            ! Computes a scale factor to multiply the coefficients of the polynomial.
            ! The scaling is done to avoid overflow and to avoid undetected underflow interfering with the convergence criterion.
            ! The factor is a power of the `BASE`.
            sc = LOW / min
            if (sc <= 1._RKR) then
                scalingEnabled = logical(max >= 10._RKR, LK)
            else
                scalingEnabled = logical(INF / sc >= max, LK)
            end if
            if (scalingEnabled) then
                if (sc == 0._RKR) sc = TIN
                factor = real(BASE, TKC) ** int(log(sc) / log(BASE) + .5_RKR, IK)
                if (factor /= 1._TKC) then
                    do concurrent(i = 1_IK : nn)
                        workspace(i) = factor * workspace(i)
                    end do
                end if
            end if

            ! Compute lower bound on moduli of zeros.
            do concurrent(i = 1_IK : nn)
                pt(i) = real(abs(workspace(i)), RKR)
            end do
            pt(nn) = -pt(nn)

            ! Compute upper estimate of bound.
            x = exp((log(-pt(nn)) - log(pt(1))) / n)
            if (pt(n) /= 0._RKR) then
                ! If newton step at the origin is better, use it.
                xm = -pt(nn) / pt(n)
                if (xm < x) x = xm
            end if

            ! Chop the interval `(0, x)` until `ff <= 0`.
            do
                xm = x * .1_RKR
                ff = pt(1)
                do i = 2_IK, nn
                    ff = ff * xm + pt(i)
                end do
                if (ff <= 0._RKR) exit
                x = xm
            end do
            dx = x

            ! Do newton iteration until `x` converges to two decimal places.
            do
                if (abs(dx / x) <= .005_RKR) exit
                ff = pt(1)
                df = ff
                do i = 2_IK, n
                    ff = ff * x + pt(i)
                    df = df * x + ff
                end do
                ff = ff * x + pt(nn)
                dx = ff / df
                x = x - dx
            end do
            bnd = x

            ! Compute the derivative as the initial `k` polynomial and do `5` steps with no shift.
            nm1 = n - 1_IK
            do i = 2_IK, n
                k(i) = (nn - i) * workspace(i) / n
            end do
            k(1) = workspace(1)
            aa = workspace(nn)
            bb = workspace(n)
            zerok = logical(k(n) == 0._TKC, LK)
            do jj = 1_IK, 5_IK
                cc = k(n)
                if (.not. zerok) then
                    ! Use scaled form of recurrence if value of k at 0 is nonzero
                    t = -aa / cc
                    do i = 1_IK, nm1
                        j = nn - i
                        k(j) = t * k(j - 1) + workspace(j)
                    end do
                    k(1) = workspace(1)
                    zerok = abs(k(n)) <= abs(bb) * EPS * 10._RKR
                else
                    ! Use unscaled form of recurrence.
                    do i = 1_IK, nm1
                        j = nn - i
                        k(j) = k(j - 1)
                    end do
                    k(1) = 0._TKC
                    zerok = logical(k(n) == 0._TKC, LK)
                end if
            end do

            ! Save k for restarts with new shifts.
            temp(1 : n) = k(1 : n)

            ! Loop to select the quadratic corresponding to each new shift.
            do cnt = 1_IK, MAX_ITER
                ! Quadratic corresponds to a double shift to a non-real point and niter complex conjugate.
                ! The point has modulus bnd and amplitude rotated by 94 degrees from the previous shift.
                xxx = COSR * xx - SINR * yy
                yy = SINR * xx + COSR * yy
                xx = xxx
                s%re = bnd * xx
                s%im = bnd * yy
                u = -2._TKC * s%re
                v = bnd
                ! Second stage calculation, fixed quadratic.
                call fxshfr(MAX_ITER * cnt, nz)
                if (nz /= 0_IK) then
                    ! The second stage jumps directly to one of the third stage iterations and returns here if successful.
                    ! Deflate the polynomial, store the zero or zeros and return to the main algorithm.
                    j = count - n + 1_IK
                    root(j) = sz
                    !write(*,*) "j", j, count, n, nn, nz
                    nn = nn - nz
                    n = nn - 1_IK
                    workspace(1 : nn) = qp(1 : nn)
                    if (nz /= 1_IK) root(j + 1) = lz
                    cycle loopIter
                end if
                ! If the iteration is unsuccessful another quadratic is chosen after restoring `k`.
                k(1 : nn) = temp(1 : nn)
            end do
            ! Return with failure if no convergence with `20` shifts.
            count = count - n
            return
        end do loopIter

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine fxshfr(l2, nz)
            ! Compute up to `l2` fixed shift `k`-polynomials, testing for convergence in the linear or quadratic case.
            ! Initiate one of the variable shift iterations and returns with the number of zeros found.
            ! l2:   limit of fixed shift steps
            ! nz:   number of zeros found
            integer(IK) , intent(in)    :: l2
            integer(IK) , intent(out)   :: nz
            real(TKC)                   :: svu, svv, ui, vi, slocal
            real(RKR)                   :: betas, betav, oss, ovv, ss, vv, ts, tv, ots, otv, tvv, tss
            logical(LK)                 :: vpass, spass, vtry, stry
            integer(IK)                 :: itype, j, iflag
            nz = 0_IK
            betav = .25_RKR
            betas = .25_RKR
            oss = real(s%re, RKR)
            ovv = real(v, RKR)
            ! Evaluate polynomial by synthetic division.
            !write(*,*) "quadsd1, nn, u, v, workspace, qp, a, b", nn, u, v, workspace, qp, a, b
            call quadsd(nn, u, v, workspace, qp, a, b)
            !write(*,*) "quadsd2, nn, u, v, workspace, qp, a, b", nn, u, v, workspace, qp, a, b
            !write(*,*) "calcsc1, itype", itype
            call calcsc(itype)
            !write(*,*) "calcsc2, itype", itype
            do j = 1_IK, l2
                ! Compute the next `k` polynomial and estimate `v`.
                call nextk(itype)
                call calcsc(itype)
                call newest(itype, ui, vi)
                vv = real(vi, RKR)
                ! Estimate `s`.
                ss = 0._RKR
                if (k(n) /= 0._TKC) ss = real(-workspace(nn) / k(n), RKR)
                tv = 1._RKR
                ts = 1._RKR
                !write(*,*) "BLOCK: ", j, itype
                if (j /= 1_IK .and. itype /= 3_IK) then
                    ! Compute relative measures of convergence of `s` and `v` sequences.
                    if (vv /= 0._RKR) tv = abs((vv - ovv) / vv)
                    if (ss /= 0._RKR) ts = abs((ss - oss) / ss)
                    ! If decreasing, multiply two most recent convergence measures.
                    tvv = 1._RKR
                    if (tv < otv) tvv = tv * otv
                    tss = 1._RKR
                    if (ts < ots) tss = ts * ots
                    ! Compare with convergence criteria.
                    vpass = tvv < betav
                    spass = tss < betas
                    if (spass .or. vpass) then
                        ! At least one sequence has passed the convergence test. Store variables before iterating.
                        svu = u
                        svv = v
                        svk(1 : n) = k(1 : n)
                        slocal = ss
                        ! Choose iteration according to the fastest converging sequence.
                        vtry = .false._LK
                        stry = .false._LK
                        if (spass .and. ((.not. vpass) .or. tss < tvv)) goto 40
                        !write(*,*) "quadit1", ui, vi, nz
                        20 call quadit(ui, vi, nz)
                        !write(*,*) "quadit2", ui, vi, nz
                        if (nz > 0_IK) return
                        ! Quadratic iteration has failed. flag that it has been tried and decrease the convergence criterion.
                        vtry = .true._LK
                        betav = betav * .25_RKR
                        ! Try linear iteration if it has not been tried and the `s` sequence is converging.
                        if (stry .or. (.not. spass)) goto 50
                        k(1 : n) = svk(1 : n)
                        40 call realit(slocal, nz, iflag)
                        if (nz > 0_IK) return
                        ! Linear iteration has failed.  flag that it has been tried and decrease the convergence criterion
                        stry = .true.
                        betas = betas * .25_RKR
                        if (iflag /= 0_IK) then
                            ! if linear iteration signals an almost double real zero attempt quadratic interation
                            ui = -(slocal + slocal)
                            vi = slocal * slocal
                            goto 20
                        end if
                        ! restore variables
                        50 u = svu
                        v = svv
                        k(1 : n) = svk(1 : n)
                        ! try quadratic iteration if it has not been tried and the v sequence is converging
                        if (vpass .and. .not. vtry) goto 20
                        ! recompute qp and scalar values to continue the second stage
                        call quadsd(nn, u, v, workspace, qp, a, b)
                        call calcsc(itype)
                    end if
                end if
                ovv = vv
                oss = ss
                otv = tv
                ots = ts
            end do
        end subroutine fxshfr

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine quadit(uu, vv, nz)
            ! Variable-shift `k`-polynomial iteration for a quadratic factor, converges only if the zeros are equimodular or nearly so.
            ! uu, vv:   coefficients of starting quadratic
            ! nz:       number of zero found
            real(TKC)   , intent(in)    :: uu
            real(TKC)   , intent(in)    :: vv
            integer(IK) , intent(out)   :: nz
            real(TKC)                   :: ui, vi
            real(RKR)                   :: mp, omp, ee, relstp, t, zm
            integer(IK)                 :: itype, i, j
            logical(LK)                 :: tried
            nz = 0_IK
            tried = .false._LK
            u = uu
            v = vv
            j = 0_IK
            ! main loop
            10 call quad(1._TKC, u, v, sz, lz)
            ! Return if roots of the quadratic are real and not close to multiple or nearly equal and  of opposite sign.
            if (abs(abs(sz%re) - abs(lz%re)) > .01_TKC * abs(lz%re)) return
            ! Evaluate polynomial by quadratic synthetic division.
            call quadsd(nn, u, v, workspace, qp, a, b)
            mp = real(abs(a - sz%re * b) + abs(sz%im * b), RKR)
            ! Compute a rigorous  bound on the rounding error in evaluating P (`workspace`).
            zm = sqrt(abs(real(v, RKR)))
            ee = 2._RKR * abs(real(qp(1), RKR))
            t = real(-sz%re * b, RKR)
            do i = 2_IK, n
                ee = ee * zm + abs(real(qp(i), RKR))
            end do
            ee = ee * zm + abs(real(a, RKR) + t)
            ee = (5._RKR * MRE + 4._RKR * ARE) * ee - (5._RKR * MRE + 2._RKR * ARE) * (abs(real(a, RKR) + t) + abs(real(b, RKR)) * zm) + 2._RKR * ARE * abs(t)
            ! Iteration has converged sufficiently if the polynomial value is less than 20 times this bound.
            if (mp <= MAX_ITER * ee) then
                nz = 2_IK
                return
            end if
            j = j + 1_IK
            ! Stop iteration after 20 steps.
            if (j > MAX_ITER) return
            if (j >= 2_IK) then
                if (.not. (relstp > .01_RKR .or. mp < omp .or. tried)) then
                    ! A cluster appears to be stalling the convergence. Five fixed shift steps are taken with a `u, v` close to the cluster.
                    if (relstp < EPS) relstp = EPS
                    relstp = sqrt(relstp)
                    u = u - u * relstp
                    v = v + v * relstp
                    call quadsd(nn, u, v, workspace, qp, a, b)
                    do i = 1_IK, 5_IK
                        call calcsc(itype)
                        call nextk(itype)
                    end do
                    tried = .true._LK
                    j = 0_IK
                end if
            end if
            omp = mp
            ! Calculate next k polynomial and new u and v.
            call calcsc(itype)
            call nextk(itype)
            call calcsc(itype)
            call newest(itype, ui, vi)
            ! If vi is zero the iteration is not converging.
            if (vi == 0._TKC) return
            relstp = real(abs((vi - v) / vi), RKR)
            u = ui
            v = vi
            goto 10
        end subroutine quadit

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine realit(sss,nz,iflag)
            ! Variable-shift `h` polynomial iteration for a real zero.
            ! sss  : starting iterate.
            ! nz   : number of zero found.
            ! iflag: flag to indicate a pair of zeros near real axis.
            real(TKC)   , intent(inout) :: sss
            integer(IK) , intent(out)   :: nz, iflag
            real(TKC)                   :: pv, kv, t, s
            real(RKR)                   :: ms, mp, omp, ee
            integer(IK)                 :: i, j
            nz = 0_IK
            s = sss
            iflag = 0_IK
            j = 0_IK
            ! main loop
            do ! 10
                pv = workspace(1)
                ! Evaluate P (`workspace`) at s.
                qp(1) = pv
                do i = 2_IK, nn
                    pv = pv * s + workspace(i)
                    qp(i) = pv
                end do
                mp = real(abs(pv), RKR)
                !write(*,*) "workspace", workspace
                !write(*,*) "mp, pv, abs(pv)", mp, pv, abs(pv)
                ! Compute a rigorous bound on the error in evaluating P (`workspace`).
                ms = abs(real(s, RKR))
                ee = (MRE / (ARE + MRE)) * abs(real(qp(1), RKR))
                do i = 2_IK, nn
                    ee = ee * ms + abs(real(qp(i), RKR))
                end do
                ! Iteration has converged sufficiently if the polynomial value is less than 20 times this bound.
                if (mp <= MAX_ITER * ((ARE + MRE) * ee - MRE * mp)) then
                    sz%re = s
                    sz%im = 0._TKC
                    nz = 1_IK
                    return
                end if
                j = j + 1_IK
                ! Stop iteration after `10` steps.
                if (j > 10_IK) return
                if (j >= 2_IK) then
                    if (abs(t) <= .001_TKC * abs(s - t) .and. mp > omp) then
                        ! A cluster of zeros near the real axis has been encountered, return with iflag set to initiate a quadratic iteration.
                        iflag = 1_IK
                        sss = s
                        return
                    end if
                end if
                ! Return if the polynomial value has increased significantly.
                omp = mp
                ! Compute t, the next polynomial, and the new iterate.
                kv = k(1)
                qk(1) = kv
                do i = 2_IK, n
                    kv = kv * s + k(i)
                    qk(i) = kv
                end do
                if (abs(kv) > abs(k(n)) * EPS * 10._TKC) then
                    ! Use the scaled form of the recurrence if the value of k at s is nonzero.
                    t = -pv / kv
                    k(1) = qp(1)
                    do i = 2_IK, n
                        k(i) = t * qk(i - 1) + qp(i)
                    end do
                else
                    ! Use unscaled form.
                    k(1) = 0._TKC
                    do i = 2_IK, n
                        k(i) = qk(i - 1)
                    end do
                end if
                kv = k(1)
                do i = 2_IK, n
                    kv = kv * s + k(i)
                end do
                t = 0._TKC
                if (abs(kv) > abs(k(n)) * EPS * 10._TKC) t = -pv / kv
                s = s + t
            end do ! goto 10
        end subroutine realit

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine calcsc(itype) ! checked
            ! Compute scalar quantities need for the next `k` polynomial and new estimates of the quadratic coefficients.
            ! itype: integer variable that indicates how the calculations are normalized to avoid overflow.
            integer(IK), intent(out) :: itype
            ! Synthetic division of k by the quadratic 1, u, v.
            call quadsd(n, u, v, k, qk, c, d)
            if (abs(c) <= abs(k(n)) * EPS * 100._TKC) then
                if (abs(d) <= abs(k(n - 1)) * EPS * 100._TKC) then
                    itype = 3_IK ! Indicates the quadratic is almost a factor of `k`.
                    return
                end if
            end if
            if (abs(d) >= abs(c)) then
                itype = 2_IK ! Indicates that all formulas are divided by `d`.
                e = a / d
                f = c / d
                g = u * b
                h = v * b
                a3 = (a + g) * e + h * (b / d)
                a1 = b * f - a
                a7 = (f + u) * a + h
                return
            end if
            itype = 1_IK ! Indicates that all formulas are divided by `c`.
            e = a / c
            f = d / c
            g = u * e
            h = v * b
            a3 = a * e + (h / c + g) * b
            a1 = b - a * (d / c)
            a7 = a + g * d + h * f
            return
        end subroutine calcsc

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine nextk(itype) ! checked
            ! Compute the next `k` polynomials using scalars computed in `calcsc()`.
            integer(IK) , intent(in)    :: itype
            real(TKC)                   :: temp
            integer(IK)                 :: i
            if (itype /= 3_IK) then
                temp = a
                if (itype == 1_IK) temp = b
                if (abs(a1) <= abs(temp) * EPS * 10._TKC) then
                    ! If `a1` is nearly zero then use a special form of the recurrence.
                    k(1) = 0._TKC
                    k(2) = -a7 * qp(1)
                    do i = 3_IK, n
                        k(i) = a3 * qk(i - 2) - a7 * qp(i - 1)
                    end do
                    return
                end if
                ! Use scaled form of the recurrence.
                a7 = a7 / a1
                a3 = a3 / a1
                k(1) = qp(1)
                k(2) = qp(2) - a7 * qp(1)
                do i = 3_IK, n
                    k(i) = a3 * qk(i - 2) - a7 * qp(i - 1) + qp(i)
                end do
                return
            end if
            ! Use unscaled form of the recurrence if itype is 3.
            k(1) = 0._TKC
            k(2) = 0._TKC
            do i = 3_IK, n
                k(i) = qk(i - 2)
            end do
        end subroutine nextk

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine newest(itype, uu, vv)
            ! Compute new estimates of the quadratic coefficients using the scalars computed in calcsc.
            integer(IK) , intent(in)    :: itype
            real(TKC)   , intent(out)   :: uu
            real(TKC)   , intent(out)   :: vv
            real(TKC)                   :: a4, a5, b1, b2, c1, c2, c3, c4, temp
            ! Use formulas appropriate to setting of itype.
            if (itype /= 3_IK) then
                if (itype /= 2_IK) then
                    a4 = a + u * b + h * f
                    a5 = c + (u + v * f) * d
                else
                    a4 = (a + g) * f + h
                    a5 = (f + u) * c + v * d
                end if
                ! Evaluate new quadratic coefficients.
                b1 = -k(n) / workspace(nn)
                b2 = -(k(n - 1) + b1 * workspace(n)) / workspace(nn)
                c1 = v * b2 * a1
                c2 = b1 * a7
                c3 = b1 * b1 * a3
                c4 = c1 - c2 - c3
                temp = a5 + b1 * a4 - c4
                if (temp /= 0._TKC) then
                    uu = u - (u * (c3 + c2) + v * (b1 * a1 + b2 * a7)) / temp
                    vv = v * (1._TKC + c4 / temp)
                    return
                end if
            end if
            ! If itype = 3 the quadratic is zeroed.
            uu = 0._TKC
            vv = 0._TKC
        end subroutine newest

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure subroutine quadsd(nn, u, v, workspace, q, a, b)
            ! Divide P (`workspace`) by the quadratic `1, u, v` placing the quotient in `q` and the remainder in `a, b`.
            integer(IK) , intent(in)    :: nn
            real(TKC)   , intent(in)    :: u, v, workspace(nn)
            real(TKC)   , intent(out)   :: q(nn), a, b
            real(TKC)                   :: c
            integer(IK)                 :: i
            b = workspace(1)
            q(1) = b
            a = workspace(2) - u * b
            q(2) = a
            do i = 3_IK, nn
                c = workspace(i) - u * a - v * b
                q(i) = c
                b = a
                a = c
            end do
        end subroutine quadsd

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure subroutine quad(a, b1, c, sroot, lroot)
            ! Compute the zeros of the quadratic `a * z**2 + b1 * z + c`.
            ! The quadratic formula, modified to avoid overflow, is used to find the larger zero if the zeros are real and both zeros are complex.
            ! The smaller real zero is found directly from the product of the zeros `c / a`.
            real(TKC)   , intent(in)    :: a, b1, c
            complex(TKC), intent(out)   :: sroot, lroot
            complex(TKC), parameter     :: ZERO = (0._TKC, 0._TKC)
            real(TKC)                   :: b, d, e
            if (a == 0._TKC) then
                lroot = ZERO
                sroot = ZERO
                if (b1 /= 0._TKC) sroot%re = -c / b1
                return
            end if
            if (c == 0._TKC) then
                lroot%re = -b1 / a
                lroot%im = 0._TKC
                sroot = ZERO
                return
            end if
            ! Compute discriminant avoiding overflow.
            b = b1 * .5_TKC
            if (abs(b) >= abs(c)) then
                e = 1._TKC - (a / b) * (c / b)
                d = sqrt(abs(e)) * abs(b)
            else
                e = a
                if (c < 0._TKC) e = -a
                e = b * (b / abs(c)) - e
                d = sqrt(abs(e)) * sqrt(abs(c))
            end if
            if (e >= 0._TKC) then
                ! Real zeros.
                if (b >= 0._TKC) d = -d
                lroot%re = (-b + d) / a
                lroot%im = 0._TKC
                sroot = ZERO
                if (lroot%re /= 0._TKC) sroot%re = (c / lroot%re) / a
                return
            end if
            ! Complex conjugate zeros.
            sroot%re = -b / a
            lroot%re = sroot%re
            sroot%im = abs(d / a)
            lroot%im = -sroot%im
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPolyRoot_ENABLED && Lag_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ideg, jdeg, niter, degree
        real(TKC), parameter :: EPS = 2 * epsilon(0._TKC) ! the estimated fractional roundoff error.
        complex(TKC) :: temp, b, c, workspace(size(coef, 1, IK))
        degree = size(coef, 1, IK) - 1_IK
        CHECK_ASSERTION(__LINE__, 0_IK < degree, SK_"@setPolyRoot(): The condition `1 < size(coef)` must hold. size(coef) = "//getStr(size(coef, 1, IK)))
        CHECK_ASSERTION(__LINE__, coef(degree) /= ZERO, SK_"@setPolyRoot(): The condition `coef(size(coef)) /= 0.` must hold. coef = "//getStr(coef))
        CHECK_ASSERTION(__LINE__, size(root, 1, IK) == degree, SK_"@setPolyRoot(): The condition `size(root) == size(coef) - 1` must hold. size(root), size(coef) - 1 = "//getStr([size(root, 1, IK), degree + 1_IK]))
        count = 0_IK
        workspace = coef
        if (degree < 1_IK) return
        do ideg = degree, 1_IK, -1_IK
            temp = ZERO
            call setPolyRootPolished(temp, niter, workspace)
            if (0_IK < niter) then
                if (abs(temp%im) <= 2._TKC * EPS**2 * abs(temp%re)) temp%im = 0._TKC
                count = count + 1_IK
                root(count) = temp
                ! deflate.
                b = workspace(ideg + 1)
                do jdeg = ideg, 1_IK, -1_IK
                    c = workspace(jdeg)
                    workspace(jdeg) = b
                    b = temp * b + c
                end do
            else
                exit
            end if
        end do
        ! Polish the roots using the un-deflated coefficients.
        ideg = 0_IK
        do
            ideg = ideg + 1_IK
            if (count < ideg) exit
            call setPolyRootPolished(root(ideg), niter, coef)
            if (0_IK < niter) cycle
            count = count - 1_IK
            root(ideg : count) = root(ideg + 1 : count + 1)
        end do
        ! Sort roots by their real parts by straight insertion.
        !do jdeg = 2_IK, degree
        !    x = root(jdeg)
        !    do ideg = jdeg - 1, 1, -1
        !        if(root(ideg)%re < x%re) exit
        !        root(ideg + 1) = root(ideg)
        !    end do
        !    root(ideg + 1) = x
        !end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPolyRootPolished_ENABLED && Lag_ENABLED && RK_RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        complex(TKC) :: croot
        croot = root
        call setPolyRootPolished(croot, niter, coef)
        root = croot%re

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPolyRootPolished_ENABLED && Lag_ENABLED && (CK_CK_ENABLED || RK_CK_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Given the complex polynomial coefficients `coef(1 : degree + 1)` of the polynomial in increasing order,
        ! a complex root `x`, this routine improves `x` by the Laguerre method until it converges,
        ! within the achievable roundoff limit, to a root of the given polynomial.<br>
        ! The number of iterations taken is returned as `niter`.<br>
        ! Break (rare) limit cycles with different fractional values once every `nstep` steps, for `nitermax` total allowed iterations.<br>
        real(TKC), parameter :: EPS = 10_IK * epsilon(0._TKC) ! the estimated fractional roundoff error.
        real(TKC), parameter :: fraction(8) = [.5_TKC, .25_TKC, .75_TKC, .13_TKC, .38_TKC, .62_TKC, .88_TKC, 1._TKC] ! Fractions used to break a limit cycle.
        integer(IK) , parameter :: nstep = 10, nitermax = size(fraction, 1, IK) * nstep * 10000 ! amir: Extra 10000 was added to allow convergence for extremely awkward coefficients.
        complex(TKC) :: droot, rootnew, der0, der1, der2, der12, der12sq, hterm, dterm, denom1, denom2
        real(TKC) :: absroot, absdenom1, absdenom2, relerr
        integer(IK) :: degree, id, ifrac, counter
        degree = size(coef, 1, IK) - 1
        CHECK_ASSERTION(__LINE__, 0_IK < degree, SK_"@setPolyRootPolished(): The condition `1 < size(coef)` must hold. size(coef) = "//getStr(size(coef, 1, IK)))
        counter = 0_IK
        ifrac = 0_IK
        niter = 0_IK
        if (degree < 1_IK) return
        do
            niter = niter + 1_IK
            der2 = coef(degree)
            relerr = abs(der2)
            der0 = ZERO
            der1 = ZERO
            absroot = abs(root)
            do id = degree - 1, 0, -1
                der0 = root * der0 + der1
                der1 = root * der1 + der2
                der2 = root * der2 + coef(id)
                relerr = abs(der2) + absroot * relerr
            end do
            relerr = EPS * relerr
            if(abs(der2) <= relerr) return
            ! Use the Laguerre method.
            der12 = der1 / der2
            der12sq = der12 * der12
            hterm = der12sq - 2 * der0 / der2
            dterm = sqrt((degree - 1) * (degree * hterm - der12sq))
            denom1 = der12 + dterm
            denom2 = der12 - dterm
            absdenom1 = abs(denom1)
            absdenom2 = abs(denom2)
            if(absdenom1 < absdenom2) then
                absdenom1 = absdenom2
                denom1 = denom2
            end if
            if (absdenom1 /= 0._TKC) then
                droot = degree / denom1
            else
                droot = exp(cmplx(log(1._TKC + absroot), real(niter, TKC), TKC))
            endif
            rootnew = root - droot
            if(root == rootnew) return
            counter = counter + 1
            if (counter /= nstep) then
                root = rootnew
            else
                counter = 0_IK
                ifrac = ifrac + 1_IK
                root = root - droot * fraction(ifrac)
                if (ifrac == size(fraction)) ifrac = 0
            end if
            if (niter == nitermax) exit
        end do
        niter = -niter

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_KIND
#undef  GET_RATIO
#undef  GET_ABS
#undef  COSD
#undef  SIND
