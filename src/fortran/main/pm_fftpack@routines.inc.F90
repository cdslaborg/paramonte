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
!>  This file contains procedure implementations of [pm_fftpack](@ref pm_fftpack).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     setFactorFFT_ENABLED && DCA_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        allocate(coef(size(data, 1, IK)))
        factor = getFactorFFT(data, coef)

        !%%%%%%%%%%%%%%%%%%%
#elif   setFactorFFT_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenData, i, j, lenFactor, nl, nq, nr, ntry
#if     CK_ENABLED
        integer(IK), parameter :: NUM_TRYH(4) = [integer(IK) :: 3, 4, 2, 5]
#elif   RK_ENABLED
        integer(IK), parameter :: NUM_TRYH(4) = [integer(IK) :: 4, 2, 3, 5]
#else
#error  "Unrecognized interface."
#endif
#if     DCX_ENABLED
        integer(IK) :: ido, ii, is, k1, l1, l2, ld, itmp
        real(TKC), parameter :: TWO_TIMES_PI = 2.0_TKC * acos( - 1.0_TKC)
        real(TKC) :: arg, argh, argld, fi
#endif
        CHECK_ASSERTION(__LINE__, size(data, kind = IK) > 1_IK, SK_"@getFactorFFT(): The condition `size(data) > 1` must hold. size(data) = "//getStr([size(data, kind = IK)]))
        lenData = size(data, kind = IK)
        allocate(factor(15))
        lenFactor = 0_IK
        nl = lenData
        j = 0_IK
        loop1: do
            j = j + 1_IK
            if (j <= 4_IK) then
                ntry = NUM_TRYH(j)
            else
                ntry = ntry + 2_IK
            end if
            loop2: do
                nq = nl / ntry
                nr = nl - ntry * nq
                if (nr /= 0_IK) cycle loop1
                lenFactor = lenFactor + 1_IK
                if (lenFactor > size(factor, kind = IK)) call setResized(factor)
                factor(lenFactor) = ntry
                nl = nq
                if (ntry == 2_IK) then
                    if (lenFactor /= 1_IK) then
                        do i = lenFactor, 2_IK, -1_IK
                            factor(i) = factor(i - 1_IK)
                        end do
                        factor(1_IK) = 2_IK
                    end if
                end if
                if (nl /= 1_IK) cycle loop2
                exit loop1
            end do loop2
        end do loop1
        factor = factor(1:lenFactor)
        !factor(2_IK) = lenFactor
        !factor(1_IK) = lenData
#if     DCX_ENABLED
        argh = TWO_TIMES_PI / real(lenData, TKC)
#if     CK_ENABLED
        i = 1_IK
        l1 = 1_IK
        itmp = lenFactor
        do k1 = 1_IK, itmp
            ld = 0_IK
            l2 = l1 * factor(k1)
            ido = lenData / l2
            itmp = ido + ido + 2_IK
            do j = 1_IK, factor(k1) - 1_IK
                is = i
                coef(i) = (1._TKC, 0._TKC)
                fi = 0._TKC
                ld = ld + l1
                argld = ld * argh
                do ii = 4_IK, itmp, 2_IK
                    i = i + 1_IK
                    fi = fi + 1._TKC
                    arg = fi * argld
                    coef(i) = cmplx(cos(arg), sin(arg), TKC)
                end do
                if (factor(k1) > 5_IK) coef(is) = coef(i)
            end do
            l1 = l2
        end do
#elif   RK_ENABLED
        if (lenFactor == 1_IK) return
        itmp = lenFactor - 1_IK
        is = 0_IK
        l1 = 1_IK
        do k1 = 1_IK, itmp
            ld = 0_IK
            l2 = l1 * factor(k1)
            ido = lenData / l2
            do j = 1_IK, factor(k1) - 1_IK
                i = is
                fi = 0._TKC
                ld = ld + l1
                argld = ld * argh
                do ii = 3_IK, ido, 2_IK
                    i = i + 2_IK
                    fi = fi + 1._TKC
                    arg = fi * argld
                    coef(i) = sin(arg)
                    coef(i - 1_IK) = cos(arg)
                end do
                is = is + ido
            end do
            l1 = l2
        end do
#else
#error  "Unrecognized interface."
#endif
#elif   !DXX_ENABLED
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setFFTF_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenData, idl1, ido, ip, ix1, ix2, ix3, ix4, k1, l1, l2, na, nac, lenFactor
        lenFactor = size(factor, kind = IK)
        lenData = size(data, kind = IK)
        CHECK_ASSERTION(__LINE__, lenData == size(coef, kind = IK), SK_"@setFFTF(): The condition `size(data) == size(coef)` must hold. size(data), size(coef) = "//getStr([lenData, size(coef, kind = IK)]))
        CHECK_ASSERTION(__LINE__, lenData == size(work, kind = IK), SK_"@setFFTF(): The condition `size(data) == size(work)` must hold. size(data), size(work) = "//getStr([lenData, size(work, kind = IK)]))
        inwork = logical(1_IK < lenData, LK)
        if (.not. inwork) return
        na = 0_IK
        l1 = 1_IK
        ix1 = 1_IK
        do k1 = 1_IK, lenFactor
            ip = factor(k1)
            l2 = ip * l1
            ido = lenData / l2
            idl1 = ido * l1
            if (ip == 4_IK) then
                ix2 = ix1 + ido
                ix3 = ix2 + ido
                if (na /= 0_IK) then
                    call passf4(ido, l1, ix1, ix2, ix3, coef, work, data)
                else
                    call passf4(ido, l1, ix1, ix2, ix3, coef, data, work)
                end if
                na = 1_IK - na
            elseif (ip == 2_IK) then
                if (na /= 0_IK) then
                    call passf2(ido, l1, ix1, coef, work, data)
                else
                    call passf2(ido, l1, ix1, coef, data, work)
                end if
                na = 1_IK - na
            elseif (ip == 3_IK) then
                ix2 = ix1 + ido
                if (na /= 0_IK) then
                    call passf3(ido, l1, ix1, ix2, coef, work, data)
                else
                    call passf3(ido, l1, ix1, ix2, coef, data, work)
                end if
                na = 1_IK - na
            elseif (ip /= 5_IK) then
                if (na /= 0_IK) then
                    call passf(nac, ido, ip, l1, idl1, ix1, coef, work, work, work, data, data)
                else
                    call passf(nac, ido, ip, l1, idl1, ix1, coef, data, data, data, work, work)
                end if
                if (nac /= 0_IK) na = 1_IK - na
            else
                ix2 = ix1 + ido
                ix3 = ix2 + ido
                ix4 = ix3 + ido
                if (na /= 0_IK) then
                    call passf5(ido, l1, ix1, ix2, ix3, ix4, coef, work, data)
                else
                    call passf5(ido, l1, ix1, ix2, ix3, ix4, coef, data, work)
                end if
                na = 1_IK - na
            end if
            l1 = l2
            ix1 = ix1 + (ip - 1_IK) * ido
        end do
        inwork = logical(na /= 0_IK, LK)
        !if (inwork) return
        !data(1:lenData) = work(1:lenData) ! \todo these copies can be avoided.

    contains

        pure subroutine passf2(ido, l1, ix1, coef, data, work)
            integer(IK) , intent(in)                    :: ido, l1, ix1
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(in)                    :: data(ido, 2, l1)
            complex(TKC), intent(out)                   :: work(ido, l1, 2)
            integer(IK) :: i, k
            if (ido > 1_IK) then
                do k = 1_IK, l1
                    do i = 1_IK, ido
                        work(i, k, 1) =  data(i, 1, k) + data(i, 2, k)
                        work(i, k, 2) = (data(i, 1, k) - data(i, 2, k)) * conjg(coef(ix1 + i))
                    end do
                end do
            else
                do k = 1_IK, l1
                    work(1, k, 1) = data(1, 1, k) + data(1, 2, k)
                    work(1, k, 2) = data(1, 1, k) - data(1, 2, k)
                end do
            end if
        end subroutine

        pure subroutine passf3(ido, l1, ix1, ix2, coef, data, work)
            integer(IK) , intent(in)                    :: ido, l1, ix1, ix2
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(in)                    :: data(ido, 3, l1)
            complex(TKC), intent(out)                   :: work(ido, l1, 3)
            complex(TKC), parameter                     :: TAU = -cmplx(0.5_TKC, sqrt(3._TKC) / 2._TKC, TKC)
            complex(TKC)                                :: c2, c3, d2, d3, t2
            integer(IK)                                 :: i, k
            if (ido /= 1_IK) then
                do k = 1_IK, l1
                    do i = 1_IK, ido
                        t2 = data(i, 2, k) + data(i, 3, k)
                        work(i, k, 1) = data(i, 1, k) + t2
                        c2 = data(i, 1, k) + TAU%re * t2
                        c3 = (data(i, 2, k) - data(i, 3, k)) * TAU%im
                        d2%re = c2%re - c3%im
                        d2%im = c2%im + c3%re
                        d3%re = c2%re + c3%im
                        d3%im = c2%im - c3%re
                        work(i, k, 2)%im = coef(ix1 + i)%re * d2%im - coef(ix1 + i)%im * d2%re
                        work(i, k, 2)%re = coef(ix1 + i)%re * d2%re + coef(ix1 + i)%im * d2%im
                        work(i, k, 3)%im = coef(ix2 + i)%re * d3%im - coef(ix2 + i)%im * d3%re
                        work(i, k, 3)%re = coef(ix2 + i)%re * d3%re + coef(ix2 + i)%im * d3%im
                    end do
                end do
            else
                do k = 1_IK, l1
                    t2 = data(1, 2, k) + data(1, 3, k)
                    c2 = data(1, 1, k) + TAU%re * t2
                    work(1, k, 1) = data(1, 1, k) + t2
                    c3 = (data(1, 2, k) - data(1, 3, k)) * TAU%im
                    work(1, k, 2)%re = c2%re - c3%im
                    work(1, k, 2)%im = c2%im + c3%re
                    work(1, k, 3)%re = c2%re + c3%im
                    work(1, k, 3)%im = c2%im - c3%re
                end do
            end if
        end subroutine

        pure subroutine passf4(ido, l1, ix1, ix2, ix3, coef, data, work)
            integer(IK) , intent(in)                    :: ido, l1, ix1, ix2, ix3
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(in)                    :: data(ido, 4, l1)
            complex(TKC), intent(out)                   :: work(ido, l1, 4)
            complex(TKC)                                :: c2, c3, c4, t1, t2, t3, t4
            integer(IK)                                 :: i, k
            if (ido /= 1_IK) then
                do k = 1_IK, l1
                    do i = 1_IK, ido
                        t1 = data(i, 1, k) - data(i, 3, k)
                        t2 = data(i, 1, k) + data(i, 3, k)
                        t3 = data(i, 2, k) + data(i, 4, k)
                        t4%re = data(i, 2, k)%im - data(i, 4, k)%im
                        t4%im = data(i, 4, k)%re - data(i, 2, k)%re
                        work(i, k, 1) = t2 + t3
                        c3 = t2 - t3
                        c2 = t1 + t4
                        c4 = t1 - t4
                        work(i, k, 2)%re = coef(ix1 + i)%re * c2%re + coef(ix1 + i)%im * c2%im
                        work(i, k, 2)%im = coef(ix1 + i)%re * c2%im - coef(ix1 + i)%im * c2%re
                        work(i, k, 3)%re = coef(ix2 + i)%re * c3%re + coef(ix2 + i)%im * c3%im
                        work(i, k, 3)%im = coef(ix2 + i)%re * c3%im - coef(ix2 + i)%im * c3%re
                        work(i, k, 4)%re = coef(ix3 + i)%re * c4%re + coef(ix3 + i)%im * c4%im
                        work(i, k, 4)%im = coef(ix3 + i)%re * c4%im - coef(ix3 + i)%im * c4%re
                    end do
                end do
            else
                do k = 1_IK, l1
                    t1 = data(1, 1, k) - data(1, 3, k)
                    t2 = data(1, 1, k) + data(1, 3, k)
                    t3 = data(1, 2, k) + data(1, 4, k)
                    t4%re = data(1, 2, k)%im - data(1, 4, k)%im
                    t4%im = data(1, 4, k)%re - data(1, 2, k)%re
                    work(1, k, 1) = t2 + t3
                    work(1, k, 2) = t1 + t4
                    work(1, k, 3) = t2 - t3
                    work(1, k, 4) = t1 - t4
                end do
            end if
        end subroutine

        PURE subroutine passf5(ido, l1, ix1, ix2, ix3, ix4, coef, data, work)
            integer(IK) , intent(in)                    :: ido, l1, ix1, ix2, ix3, ix4
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(in)                    :: data(ido, 5, l1)
            complex(TKC), intent(out)                   :: work(ido, l1, 5)
            real(TKC)   , parameter                     :: PI = acos( - 1._TKC)
            complex(TKC), parameter                     :: T11 = cmplx(cos(2._TKC * PI / 5._TKC), -sin(2._TKC * PI / 5._TKC), TKC)
            complex(TKC), parameter                     :: T12 = cmplx(cos(4._TKC * PI / 5._TKC), -sin(4._TKC * PI / 5._TKC), TKC)
            complex(TKC)                                :: c2, c3, c4, c5, d2, d3, d4, d5, t2, t3, t4, t5
            integer(IK)                                 :: i, k
            if (ido /= 1_IK) then
                do k = 1_IK, l1
                    do i = 1_IK, ido
                        t2 = data(i, 2, k) + data(i, 5, k)
                        t3 = data(i, 3, k) + data(i, 4, k)
                        t4 = data(i, 3, k) - data(i, 4, k)
                        t5 = data(i, 2, k) - data(i, 5, k)
                        work(i, k, 1) = data(i, 1, k) + t2 + t3
                        c2%re = data(i, 1, k)%re + T11%re * t2%re + T12%re * t3%re
                        c2%im = data(i, 1, k)%im + T11%re * t2%im + T12%re * t3%im
                        c3%re = data(i, 1, k)%re + T12%re * t2%re + T11%re * t3%re
                        c3%im = data(i, 1, k)%im + T12%re * t2%im + T11%re * t3%im
                        c4%re = T12%im * t5%re - T11%im * t4%re
                        c4%im = T12%im * t5%im - T11%im * t4%im
                        c5%re = T11%im * t5%re + T12%im * t4%re
                        c5%im = T11%im * t5%im + T12%im * t4%im
                        d2%re = c2%re - c5%im
                        d2%im = c2%im + c5%re
                        d3%re = c3%re - c4%im
                        d3%im = c3%im + c4%re
                        d4%re = c3%re + c4%im
                        d4%im = c3%im - c4%re
                        d5%re = c2%re + c5%im
                        d5%im = c2%im - c5%re
                        work(i, k, 2)%re = coef(ix1 + i)%re * d2%re + coef(ix1 + i)%im * d2%im
                        work(i, k, 2)%im = coef(ix1 + i)%re * d2%im - coef(ix1 + i)%im * d2%re
                        work(i, k, 3)%re = coef(ix2 + i)%re * d3%re + coef(ix2 + i)%im * d3%im
                        work(i, k, 3)%im = coef(ix2 + i)%re * d3%im - coef(ix2 + i)%im * d3%re
                        work(i, k, 4)%re = coef(ix3 + i)%re * d4%re + coef(ix3 + i)%im * d4%im
                        work(i, k, 4)%im = coef(ix3 + i)%re * d4%im - coef(ix3 + i)%im * d4%re
                        work(i, k, 5)%re = coef(ix4 + i)%re * d5%re + coef(ix4 + i)%im * d5%im
                        work(i, k, 5)%im = coef(ix4 + i)%re * d5%im - coef(ix4 + i)%im * d5%re
                    end do
                end do
            else
                do k = 1_IK, l1
                    t2 = data(1, 2, k) + data(1, 5, k)
                    t3 = data(1, 3, k) + data(1, 4, k)
                    t4 = data(1, 3, k) - data(1, 4, k)
                    t5 = data(1, 2, k) - data(1, 5, k)
                    work(1, k, 1) = data(1, 1, k) + t2 + t3
                    c2%re = data(1, 1, k)%re + T11%re * t2%re + T12%re * t3%re
                    c2%im = data(1, 1, k)%im + T11%re * t2%im + T12%re * t3%im
                    c3%re = data(1, 1, k)%re + T12%re * t2%re + T11%re * t3%re
                    c3%im = data(1, 1, k)%im + T12%re * t2%im + T11%re * t3%im
                    c5%re = T11%im * t5%re + T12%im * t4%re
                    c5%im = T11%im * t5%im + T12%im * t4%im
                    c4%re = T12%im * t5%re - T11%im * t4%re
                    c4%im = T12%im * t5%im - T11%im * t4%im
                    work(1, k, 2)%re = c2%re - c5%im
                    work(1, k, 2)%im = c2%im + c5%re
                    work(1, k, 3)%re = c3%re - c4%im
                    work(1, k, 3)%im = c3%im + c4%re
                    work(1, k, 4)%re = c3%re + c4%im
                    work(1, k, 4)%im = c3%im - c4%re
                    work(1, k, 5)%re = c2%re + c5%im
                    work(1, k, 5)%im = c2%im - c5%re
                    !work(1, k, 2) = c2%re - c5%im
                    !work(2, k, 2) = c2%im + c5%re
                    !work(1, k, 3) = c3%re - c4%im
                    !work(2, k, 3) = c3%im + c4%re
                    !work(1, k, 4) = c3%re + c4%im
                    !work(2, k, 4) = c3%im - c4%re
                    !work(1, k, 5) = c2%re + c5%im
                    !work(2, k, 5) = c2%im - c5%re
                end do
            end if
        end subroutine

        pure subroutine passf(nac, ido, ip, l1, idl1, ix1, coef, data, C1, C2, work, Ch2)
            integer(IK) , intent(out)                   :: nac
            integer(IK) , intent(in)                    :: ido, ip, l1, idl1, ix1
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(inout)                 :: C1(ido,l1,ip), C2(idl1,ip)
            complex(TKC), intent(inout)                 :: data(ido,ip,l1), work(ido,l1,ip), Ch2(idl1,ip)
            integer(IK) :: i, idij, idj, idl, idlj, idp, jk, inc, ipp2, ipph, j, jc, k, l, lc !, nt, idot
            !idot = ido / 2_IK
            !nt = ip * idl1 what in the world is this doing here?
            idp = ip * ido
            ipp2 = ip + 2_IK
            ipph = (ip + 1_IK) / 2_IK
            if (2_IK * ido < l1) then
                do j = 2_IK, ipph
                    jc = ipp2 - j
                    do i = 1_IK, ido
                        do k = 1_IK, l1
                            work(i, k, j) = data(i, j, k) + data(i, jc, k)
                            work(i, k, jc) = data(i, j, k) - data(i, jc, k)
                        end do
                    end do
                end do
                do i = 1_IK, ido
                    do k = 1_IK, l1
                        work(i, k, 1) = data(i, 1, k)
                    end do
                end do
            else
                do j = 2_IK, ipph
                    jc = ipp2 - j
                    do k = 1_IK, l1
                        do i = 1_IK, ido
                            work(i, k, j) = data(i, j, k) + data(i, jc, k)
                            work(i, k, jc) = data(i, j, k) - data(i, jc, k)
                        end do
                    end do
                end do
                do k = 1_IK, l1
                    do i = 1, ido
                        work(i, k, 1) = data(i, 1, k)
                    end do
                end do
            end if
            inc = 0_IK
            idl = 1_IK - ido
            do l = 2_IK, ipph
                lc = ipp2 - l
                idl = idl + ido
                do jk = 1_IK, idl1
                    C2(jk, l) = Ch2(jk, 1) + coef(ix1 + idl)%re * Ch2(jk, 2)
                    C2(jk, lc) = -coef(ix1 + idl)%im * Ch2(jk, ip)
                end do
                idlj = idl
                inc = inc + ido
                do j = 3_IK, ipph
                    jc = ipp2 - j
                    idlj = idlj + inc
                    if (idlj > idp) idlj = idlj - idp
                    do jk = 1_IK, idl1
                        C2(jk, l) = C2(jk, l) + coef(ix1 + idlj)%re * Ch2(jk, j)
                        C2(jk, lc) = C2(jk, lc) - coef(ix1 + idlj)%im * Ch2(jk, jc)
                    end do
                end do
            end do
            do j = 2_IK, ipph
                do jk = 1_IK, idl1
                    Ch2(jk, 1) = Ch2(jk, 1) + Ch2(jk, j)
                end do
            end do
            do j = 2_IK, ipph
                jc = ipp2 - j
                do jk = 1_IK, idl1
                    Ch2(jk, j)%re = C2(jk, j)%re - C2(jk, jc)%im
                    Ch2(jk, j)%im = C2(jk, j)%im + C2(jk, jc)%re
                    Ch2(jk, jc)%re = C2(jk, j)%re + C2(jk, jc)%im
                    Ch2(jk, jc)%im = C2(jk, j)%im - C2(jk, jc)%re
                end do
            end do
            nac = 1_IK
            if (ido == 1_IK) return
            nac = 0_IK
            do jk = 1_IK, idl1
                C2(jk, 1) = Ch2(jk, 1)
            end do
            do j = 2_IK, ip
                do k = 1_IK, l1
                    C1(1, k, j) = work(1, k, j)
                end do
            end do
            if (ido > l1) then
                idj = 1_IK - ido
                do j = 2_IK, ip
                    idj = idj + ido
                    do k = 1_IK, l1
                        idij = idj
                        do i = 2_IK, ido
                            idij = idij + 1_IK
                            C1(i, k, j)%re = coef(ix1 + idij)%re * work(i, k, j)%re + coef(ix1 + idij)%im * work(i, k, j)%im
                            C1(i, k, j)%im = coef(ix1 + idij)%re * work(i, k, j)%im - coef(ix1 + idij)%im * work(i, k, j)%re
                        end do
                    end do
                end do
            else
                idij = 0_IK
                do j = 2_IK, ip
                    idij = idij + 1_IK
                    do i = 2_IK, ido
                        idij = idij + 2_IK
                        do k = 1_IK, l1
                            C1(i, k, j)%re = coef(ix1 + idij)%re * work(i, k, j)%re + coef(ix1 + idij)%im * work(i, k, j)%im
                            C1(i, k, j)%im = coef(ix1 + idij)%re * work(i, k, j)%im - coef(ix1 + idij)%im * work(i, k, j)%re
                        end do
                    end do
                end do
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setFFTF_ENABLED && RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenData, idl1, ido, ip, ix1, ix2, ix3, ix4, k1, l1, l2, na, lenFactor
        lenFactor = size(factor, kind = IK)
        lenData = size(data, kind = IK)
        CHECK_ASSERTION(__LINE__, lenData == size(coef, kind = IK), SK_"@setFFTF(): The condition `size(data) == size(coef)` must hold. size(data), size(coef) = "//getStr([lenData, size(coef, kind = IK)]))
        CHECK_ASSERTION(__LINE__, lenData == size(work, kind = IK), SK_"@setFFTF(): The condition `size(data) == size(work)` must hold. size(data), size(work) = "//getStr([lenData, size(work, kind = IK)]))
        inwork = logical(1_IK < lenData, LK)
        if (.not. inwork) return
        na = 1_IK
        l2 = lenData
        ix1 = lenData
        do k1 = 1_IK, lenFactor
            ip = factor(lenFactor - k1 + 1_IK)
            l1 = l2 / ip
            ido = lenData / l2
            idl1 = ido * l1
            ix1 = ix1 - (ip - 1_IK) * ido
            na = 1_IK - na
            if (ip == 4) then
                ix2 = ix1 + ido
                ix3 = ix2 + ido
                if (na /= 0_IK) then
                    call radf4(ido, l1, ix1, ix2, ix3, coef, work, data)
                else
                    call radf4(ido, l1, ix1, ix2, ix3, coef, data, work)
                end if
            elseif (ip /= 2_IK) then
                if (ip == 3_IK) then
                    ix2 = ix1 + ido
                    if (na /= 0_IK) then
                        call radf3(ido, l1, ix1, ix2, coef, work, data)
                    else
                        call radf3(ido, l1, ix1, ix2, coef, data, work)
                    end if
                elseif (ip /= 5_IK) then
                    if (ido == 1_IK) na = 1_IK - na
                    if (na /= 0_IK) then
                        na = 0_IK
                        call radfg(ido, ip, l1, idl1, ix1, coef, work, work, work, data, data)
                    else
                        call radfg(ido, ip, l1, idl1, ix1, coef, data, data, data, work, work)
                        na = 1_IK
                    end if
                else
                    ix2 = ix1 + ido
                    ix3 = ix2 + ido
                    ix4 = ix3 + ido
                    if (na /= 0_IK) then
                        call radf5(ido, l1, ix1, ix2, ix3, ix4, coef, work, data)
                    else
                        call radf5(ido, l1, ix1, ix2, ix3, ix4, coef, data, work)
                    end if
                end if
            elseif (na /= 0_IK) then
                call radf2(ido, l1, ix1, coef, work, data)
            else
                call radf2(ido, l1, ix1, coef, data, work)
            end if
            l2 = l1
        end do
        inwork = logical(na /= 1_IK, LK)
        !if (na == 1_IK) return
        !data(1:lenData) = work(1:lenData)

    contains

        pure subroutine radf2(ido, l1, ix1, coef, data, work)
            integer(IK) , intent(in)                :: ido, l1, ix1
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(in)                :: data(ido, l1, 2)
            real(TKC)   , intent(out)               :: work(ido, 2, l1)
            integer(IK)                             :: i, ic, idp2, k
            complex(TKC)                            :: t2
            do k = 1_IK, l1
                work(1, 1, k) = data(1, k, 1) + data(1, k, 2)
                work(ido, 2, k) = data(1, k, 1) - data(1, k, 2)
            end do
            if (ido < 2_IK) return
            if (ido /= 2_IK) then
                idp2 = ido + 2_IK
                do k = 1_IK, l1
                    do i = 3_IK, ido, 2_IK
                        ic = idp2 - i
                        t2%re = coef(ix1 + i - 2) * data(i - 1, k, 2) + coef(ix1 + i - 1) * data(i, k, 2)
                        t2%im = coef(ix1 + i - 2) * data(i, k, 2) - coef(ix1 + i - 1) * data(i - 1, k, 2)
                        work(i, 1, k) = data(i, k, 1) + t2%im
                        work(ic, 2, k) = t2%im - data(i, k, 1)
                        work(i - 1, 1, k) = data(i - 1, k, 1) + t2%re
                        work(ic - 1, 2, k) = data(i - 1, k, 1) - t2%re
                    end do
                end do
                if (mod(ido, 2_IK) == 1_IK) return
            end if
            do k = 1_IK, l1
                work(1, 2, k) = -data(ido, k, 2)
                work(ido, 1, k) = data(ido, k, 1)
            end do
        end subroutine

        pure subroutine radf3(ido, l1, ix1, ix2, coef, data, work)
            integer(IK) , intent(in)                :: ido, l1, ix1, ix2
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(in)                :: data(ido, l1, 3)
            real(TKC)   , intent(out)               :: work(ido, 3, l1)
            complex(TKC), parameter                 :: TAU = cmplx(-.5_TKC, sqrt(3._TKC) / 2._TKC, TKC)
            complex(TKC)                            :: c2, d2, d3, t2, t3
            integer(IK)                             :: i, ic, idp2, k
            do k = 1_IK, l1
                c2%re = data(1, k, 2) + data(1, k, 3)
                work(1, 1, k) = data(1, k, 1) + c2%re
                work(1, 3, k) = TAU%im * (data(1, k, 3) - data(1, k, 2))
                work(ido, 2, k) = data(1, k, 1) + TAU%re * c2%re
            end do
            if (ido == 1_IK) return
            idp2 = ido + 2_IK
            do k = 1_IK, l1
                do i = 3_IK, ido, 2_IK
                    ic = idp2 - i
                    d2%re = coef(ix1 + i - 2) * data(i - 1, k, 2) + coef(ix1 + i - 1) * data(i, k, 2)
                    d2%im = coef(ix1 + i - 2) * data(i, k, 2) - coef(ix1 + i - 1) * data(i - 1, k, 2)
                    d3%re = coef(ix2 + i - 2) * data(i - 1, k, 3) + coef(ix2 + i - 1) * data(i, k, 3)
                    d3%im = coef(ix2 + i - 2) * data(i, k, 3) - coef(ix2 + i - 1) * data(i - 1, k, 3)
                    c2%re = d2%re + d3%re
                    c2%im = d2%im + d3%im
                    work(i - 1, 1, k) = data(i - 1, k, 1) + c2%re
                    work(i, 1, k) = data(i, k, 1) + c2%im
                    t2%re = data(i - 1, k, 1) + TAU%re * c2%re
                    t2%im = data(i, k, 1) + TAU%re * c2%im
                    t3%re = TAU%im * (d2%im - d3%im)
                    t3%im = TAU%im * (d3%re - d2%re)
                    work(i - 1, 3, k) = t2%re + t3%re
                    work(ic - 1, 2, k) = t2%re - t3%re
                    work(i, 3, k) = t2%im + t3%im
                    work(ic, 2, k) = t3%im - t2%im
                end do
            end do
        end subroutine

        pure subroutine radf4(ido, l1, ix1, ix2, ix3, coef, data, work)
            integer(IK) , intent(in)                :: ido, l1, ix1, ix2, ix3
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(in)                :: data(ido, l1, 4)
            real(TKC)   , intent(out)               :: work(ido, 4, l1)
            real(TKC)   , parameter                 :: NHSQT2 = -sqrt(2._TKC) / 2._TKC
            complex(TKC)                            :: c2, c3, c4, t1, t2, t3, t4
            integer(IK)                             :: i, ic, idp2, k
            do k = 1_IK, l1
                t1%re = data(1, k, 2) + data(1, k, 4)
                t2%re = data(1, k, 1) + data(1, k, 3)
                work(1, 1, k) = t1%re + t2%re
                work(ido, 4, k) = t2%re - t1%re
                work(ido, 2, k) = data(1, k, 1) - data(1, k, 3)
                work(1, 3, k) = data(1, k, 4) - data(1, k, 2)
            end do
            if (ido < 2_IK) return
            if (ido /= 2_IK) then
               idp2 = ido + 2_IK
               do k = 1_IK, l1
                    do i = 3_IK, ido, 2_IK
                        ic = idp2 - i
                        c2%re = coef(ix1 + i - 2) * data(i - 1, k, 2) + coef(ix1 + i - 1) * data(i, k, 2)
                        c2%im = coef(ix1 + i - 2) * data(i, k, 2) - coef(ix1 + i - 1) * data(i - 1, k, 2)
                        c3%re = coef(ix2 + i - 2) * data(i - 1, k, 3) + coef(ix2 + i - 1) * data(i, k, 3)
                        c3%im = coef(ix2 + i - 2) * data(i, k, 3) - coef(ix2 + i - 1) * data(i - 1, k, 3)
                        c4%re = coef(ix3 + i - 2) * data(i - 1, k, 4) + coef(ix3 + i - 1) * data(i, k, 4)
                        c4%im = coef(ix3 + i - 2) * data(i, k, 4) - coef(ix3 + i - 1) * data(i - 1, k, 4)
                        t1%re = c2%re + c4%re
                        t4%re = c4%re - c2%re
                        t1%im = c2%im + c4%im
                        t4%im = c2%im - c4%im
                        t2%im = data(i, k, 1) + c3%im
                        t3%im = data(i, k, 1) - c3%im
                        t2%re = data(i - 1, k, 1) + c3%re
                        t3%re = data(i - 1, k, 1) - c3%re
                        work(i - 1, 1, k) = t1%re + t2%re
                        work(ic - 1, 4, k) = t2%re - t1%re
                        work(i, 1, k) = t1%im + t2%im
                        work(ic, 4, k) = t1%im - t2%im
                        work(i - 1, 3, k) = t4%im + t3%re
                        work(ic - 1, 2, k) = t3%re - t4%im
                        work(i, 3, k) = t4%re + t3%im
                        work(ic, 2, k) = t4%re - t3%im
                    end do
               end do
               if (mod(ido, 2_IK) == 1_IK) return
            end if
            do k = 1_IK, l1
                t1%im = NHSQT2 * (data(ido, k, 2) + data(ido, k, 4))
                t1%re = NHSQT2 * (data(ido, k, 4) - data(ido, k, 2))
                work(ido, 1, k) = t1%re + data(ido, k, 1)
                work(ido, 3, k) = data(ido, k, 1) - t1%re
                work(1, 2, k) = t1%im - data(ido, k, 3)
                work(1, 4, k) = t1%im + data(ido, k, 3)
            end do
        end subroutine

        pure subroutine radf5(ido, l1, ix1, ix2, ix3, ix4, coef, data, work)
            integer(IK) , intent(in)                :: ido, l1, ix1, ix2, ix3, ix4
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(in)                :: data(ido, l1, 5)
            real(TKC)   , intent(out)               :: work(ido, 5, l1)
            real(TKC)   , parameter                 :: PI = acos(-1._TKC)
            complex(TKC), parameter                 :: T11 = cmplx(cos(2._TKC * PI / 5._TKC), sin(2._TKC * PI / 5._TKC), TKC)
            complex(TKC), parameter                 :: T12 = cmplx(cos(4._TKC * PI / 5._TKC), sin(4._TKC * PI / 5._TKC), TKC)
            complex(TKC)                            :: c2, c3, c4, c5, d2, d3, d4, d5, t2, t3, t4, t5
            integer(IK)                             :: i, ic, idp2, k
            do k = 1_IK, l1
                c2%re = data(1, k, 5) + data(1, k, 2)
                c5%im = data(1, k, 5) - data(1, k, 2)
                c3%re = data(1, k, 4) + data(1, k, 3)
                c4%im = data(1, k, 4) - data(1, k, 3)
                work(1, 1, k) = data(1, k, 1) + c2%re + c3%re
                work(ido, 2, k) = data(1, k, 1) + T11%re * c2%re + T12%re * c3%re
                work(1, 3, k) = T11%im * c5%im + T12%im * c4%im
                work(ido, 4, k) = data(1, k, 1) + T12%re * c2%re + T11%re * c3%re
                work(1, 5, k) = T12%im * c5%im - T11%im * c4%im
            end do
            if (ido == 1_IK) return
            idp2 = ido + 2_IK
            do k = 1_IK, l1
                do i = 3_IK, ido, 2_IK
                    ic = idp2 - i
                    d2%re = coef(ix1 + i - 2) * data(i - 1, k, 2) + coef(ix1 + i - 1) * data(i, k, 2)
                    d2%im = coef(ix1 + i - 2) * data(i, k, 2) - coef(ix1 + i - 1) * data(i - 1, k, 2)
                    d3%re = coef(ix2 + i - 2) * data(i - 1, k, 3) + coef(ix2 + i - 1) * data(i, k, 3)
                    d3%im = coef(ix2 + i - 2) * data(i, k, 3) - coef(ix2 + i - 1) * data(i - 1, k, 3)
                    d4%re = coef(ix3 + i - 2) * data(i - 1, k, 4) + coef(ix3 + i - 1) * data(i, k, 4)
                    d4%im = coef(ix3 + i - 2) * data(i, k, 4) - coef(ix3 + i - 1) * data(i - 1, k, 4)
                    d5%re = coef(ix4 + i - 2) * data(i - 1, k, 5) + coef(ix4 + i - 1) * data(i, k, 5)
                    d5%im = coef(ix4 + i - 2) * data(i, k, 5) - coef(ix4 + i - 1) * data(i - 1, k, 5)
                    c2%re = d2%re + d5%re
                    c5%im = d5%re - d2%re
                    c5%re = d2%im - d5%im
                    c2%im = d2%im + d5%im
                    c3%re = d3%re + d4%re
                    c4%im = d4%re - d3%re
                    c4%re = d3%im - d4%im
                    c3%im = d3%im + d4%im
                    work(i - 1, 1, k) = data(i - 1, k, 1) + c2%re + c3%re
                    work(i, 1, k) = data(i, k, 1) + c2%im + c3%im
                    t2%re = data(i - 1, k, 1) + T11%re * c2%re + T12%re * c3%re
                    t2%im = data(i, k, 1) + T11%re * c2%im + T12%re * c3%im
                    t3%re = data(i - 1, k, 1) + T12%re * c2%re + T11%re * c3%re
                    t3%im = data(i, k, 1) + T12%re * c2%im + T11%re * c3%im
                    t5%re = T11%im * c5%re + T12%im * c4%re
                    t5%im = T11%im * c5%im + T12%im * c4%im
                    t4%re = T12%im * c5%re - T11%im * c4%re
                    t4%im = T12%im * c5%im - T11%im * c4%im
                    work(i - 1, 3, k) = t2%re + t5%re
                    work(ic - 1, 2, k) = t2%re - t5%re
                    work(i, 3, k) = t2%im + t5%im
                    work(ic, 2, k) = t5%im - t2%im
                    work(i - 1, 5, k) = t3%re + t4%re
                    work(ic - 1, 4, k) = t3%re - t4%re
                    work(i, 5, k) = t3%im + t4%im
                    work(ic, 4, k) = t4%im - t3%im
                end do
            end do
        end subroutine radf5

        pure subroutine radfg(ido, ip, l1, idl1, ix1, coef, data, C1, C2, work, Ch2)
            integer(IK) , intent(in)                :: ido, ip, l1, idl1, ix1
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(inout)             :: C1(ido,l1,ip), C2(idl1,ip)
            real(TKC)   , intent(inout)             :: data(ido,ip,l1), work(ido,l1,ip), Ch2(idl1,ip)
            real(TKC)   , parameter                 :: TWO_PI = 2._TKC * acos(-1._TKC)
            real(TKC)                               :: ar1h, ar2h, arg, dc2, dcp, ds2, dsp
            integer(IK)                             :: i, ic, idij, idp2, jk, ipp2, ipph, is, j, j2, jc, k, l, lc, nbd
            complex(TKC)                            :: a1, a2
            arg = TWO_PI / real(ip, TKC)
            dcp = cos(arg)
            dsp = sin(arg)
            ipph = (ip + 1_IK) / 2_IK
            ipp2 = ip + 2_IK
            idp2 = ido + 2_IK
            nbd = (ido - 1_IK) / 2_IK
            if (ido == 1_IK) then
                C2(1:idl1, 1) = Ch2(1:idl1, 1)
            else
                Ch2(1:idl1, 1) = C2(1:idl1, 1)
                work(1, 1:l1, 2:ip) = C1(1, 1:l1, 2:ip)
                if (nbd > l1) then
                    is = -ido
                    do j = 2_IK, ip
                        is = is + ido
                        do k = 1_IK, l1
                            idij = is
                            do i = 3_IK, ido, 2_IK
                                idij = idij + 2_IK
                                work(i - 1, k, j) = coef(ix1 + idij - 1) * C1(i - 1, k, j) + coef(ix1 + idij) * C1(i, k, j)
                                work(i, k, j) = coef(ix1 + idij - 1) * C1(i, k, j) - coef(ix1 + idij) * C1(i - 1, k, j)
                            end do
                        end do
                    end do
                else
                    is = -ido
                    do j = 2_IK, ip
                        is = is + ido
                        idij = is
                        do i = 3_IK, ido, 2_IK
                            idij = idij + 2_IK
                            do k = 1_IK, l1
                                work(i - 1, k, j) = coef(ix1 + idij - 1) * C1(i - 1, k, j) + coef(ix1 + idij) * C1(i, k, j)
                                work(i, k, j) = coef(ix1 + idij - 1) * C1(i, k, j) - coef(ix1 + idij) * C1(i - 1, k, j)
                            end do
                        end do
                    end do
                end if
                if (nbd < l1) then
                    do j = 2_IK, ipph
                        jc = ipp2 - j
                        do i = 3_IK, ido, 2_IK
                            do k = 1_IK, l1
                                C1(i - 1, k, j) = work(i - 1, k, j) + work(i - 1, k, jc)
                                C1(i - 1, k, jc) = work(i, k, j) - work(i, k, jc)
                                C1(i, k, j) = work(i, k, j) + work(i, k, jc)
                                C1(i, k, jc) = work(i - 1, k, jc) - work(i - 1, k, j)
                            end do
                        end do
                    end do
                else
                    do j = 2_IK, ipph
                        jc = ipp2 - j
                        do k = 1_IK, l1
                            do i = 3_IK, ido, 2_IK
                                C1(i - 1, k, j) = work(i - 1, k, j) + work(i - 1, k, jc)
                                C1(i - 1, k, jc) = work(i, k, j) - work(i, k, jc)
                                C1(i, k, j) = work(i, k, j) + work(i, k, jc)
                                C1(i, k, jc) = work(i - 1, k, jc) - work(i - 1, k, j)
                            end do
                        end do
                    end do
                end if
            end if
            do j = 2_IK, ipph
                jc = ipp2 - j
                do k = 1_IK, l1
                    C1(1, k, j) = work(1, k, j) + work(1, k, jc)
                    C1(1, k, jc) = work(1, k, jc) - work(1, k, j)
                end do
            end do
            a1%re = 1._TKC
            a1%im = 0._TKC
            do l = 2_IK, ipph
                lc = ipp2 - l
                ar1h = dcp * a1%re - dsp * a1%im
                a1%im = dcp * a1%im + dsp * a1%re
                a1%re = ar1h
                do jk = 1, idl1
                    Ch2(jk, l) = C2(jk, 1) + a1%re * C2(jk, 2)
                    Ch2(jk, lc) = a1%im * C2(jk, ip)
                end do
                dc2 = a1%re
                ds2 = a1%im
                a2%re = a1%re
                a2%im = a1%im
                do j = 3_IK, ipph
                    jc = ipp2 - j
                    ar2h = dc2 * a2%re - ds2 * a2%im
                    a2%im = dc2 * a2%im + ds2 * a2%re
                    a2%re = ar2h
                    do jk = 1_IK, idl1
                        Ch2(jk, l) = Ch2(jk, l) + a2%re * C2(jk, j)
                        Ch2(jk, lc) = Ch2(jk, lc) + a2%im * C2(jk, jc)
                    end do
                end do
            end do
            do j = 2_IK, ipph
                do jk = 1_IK, idl1
                    Ch2(jk, 1) = Ch2(jk, 1) + C2(jk, j)
                end do
            end do
            if (ido < l1) then
                do i = 1, ido
                    do k = 1, l1
                        data(i, 1, k) = work(i, k, 1)
                    end do
                end do
            else
                do k = 1_IK, l1
                    do i = 1, ido
                        data(i, 1, k) = work(i, k, 1)
                    end do
                end do
            end if
            do j = 2_IK, ipph
                jc = ipp2 - j
                j2 = j + j
                do k = 1, l1
                    data(ido, j2 - 2, k) = work(1, k, j)
                    data(1, j2 - 1, k) = work(1, k, jc)
                end do
            end do
            if (ido == 1_IK) return
            if (nbd < l1) then
                do j = 2_IK, ipph
                    jc = ipp2 - j
                    j2 = j + j
                    do i = 3_IK, ido, 2_IK
                        ic = idp2 - i
                        do k = 1, l1
                            data(i - 1, j2 - 1, k) = work(i - 1, k, j) + work(i - 1, k, jc)
                            data(ic - 1, j2 - 2, k) = work(i - 1, k, j) - work(i - 1, k, jc)
                            data(i, j2 - 1, k) = work(i, k, j) + work(i, k, jc)
                            data(ic, j2 - 2, k) = work(i, k, jc) - work(i, k, j)
                        end do
                    end do
                end do
            else
                do j = 2_IK, ipph
                    jc = ipp2 - j
                    j2 = j + j
                    do k = 1_IK, l1
                        do i = 3_IK, ido, 2_IK
                            ic = idp2 - i
                            data(i - 1, j2 - 1, k) = work(i - 1, k, j) + work(i - 1, k, jc)
                            data(ic - 1, j2 - 2, k) = work(i - 1, k, j) - work(i - 1, k, jc)
                            data(i, j2 - 1, k) = work(i, k, j) + work(i, k, jc)
                            data(ic, j2 - 2, k) = work(i, k, jc) - work(i, k, j)
                        end do
                    end do
                end do
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setFFTR_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenData, idl1, ido, ip, ix1, ix2, ix3, ix4, k1, l1, l2, na, nac, lenFactor
        lenFactor = size(factor, kind = IK)
        lenData = size(data, kind = IK)
        CHECK_ASSERTION(__LINE__, lenData == size(coef, kind = IK), SK_"@setFFTR(): The condition `size(data) == size(coef)` must hold. size(data), size(coef) = "//getStr([lenData, size(coef, kind = IK)]))
        CHECK_ASSERTION(__LINE__, lenData == size(work, kind = IK), SK_"@setFFTR(): The condition `size(data) == size(work)` must hold. size(data), size(work) = "//getStr([lenData, size(work, kind = IK)]))
        inwork = logical(1_IK < lenData, LK)
        if (.not. inwork) return
        na = 0_IK
        l1 = 1_IK
        ix1 = 1_IK
        do k1 = 1_IK, lenFactor
            ip = factor(k1)
            l2 = ip * l1
            ido = lenData / l2
            !idot = ido + ido
            idl1 = ido * l1
            if (ip == 4_IK) then
                ix2 = ix1 + ido
                ix3 = ix2 + ido
                if (na /= 0_IK) then
                    call passb4(ido, l1, ix1, ix2, ix3, coef, work, data)
                else
                    call passb4(ido, l1, ix1, ix2, ix3, coef, data, work)
                end if
                na = 1_IK - na
            elseif (ip == 2_IK) then
                if (na /= 0_IK) then
                    call passb2(ido, l1, ix1, coef, work, data)
                else
                    call passb2(ido, l1, ix1, coef, data, work)
                end if
                na = 1_IK - na
            elseif (ip == 3_IK) then
                ix2 = ix1 + ido
                if (na /= 0_IK) then
                    call passb3(ido, l1, ix1, ix2, coef, work, data)
                else
                    call passb3(ido, l1, ix1, ix2, coef, data, work)
                end if
                na = 1_IK - na
            elseif (ip /= 5_IK) then
                if (na /= 0_IK) then
                    call passb(nac, ido, ip, l1, idl1, ix1, coef, work, work, work, data, data)
                else
                    call passb(nac, ido, ip, l1, idl1, ix1, coef, data, data, data, work, work)
                end if
                if (nac /= 0_IK) na = 1_IK - na
            else
                ix2 = ix1 + ido
                ix3 = ix2 + ido
                ix4 = ix3 + ido
                if (na /= 0_IK) then
                    call passb5(ido, l1, ix1, ix2, ix3, ix4, coef, work, data)
                else
                    call passb5(ido, l1, ix1, ix2, ix3, ix4, coef, data, work)
                end if
                na = 1_IK - na
            end if
            l1 = l2
            ix1 = ix1 + (ip - 1_IK) * ido
        end do
        inwork = logical(na /= 0_IK, LK)
        !if (na == 0_IK) return
        !data(1 : lenData) = work(1 : lenData)

    contains

        pure subroutine passb2(ido, l1, ix1, coef, data, work)
            integer(IK) , intent(in)                    :: ido, l1, ix1
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(in)                    :: data(ido, 2, l1)
            complex(TKC), intent(out)                   :: work(ido, l1, 2)
            complex(TKC)                                :: t2
            integer(IK)                                 :: i, k
            if (ido > 1_IK) then
                do k = 1_IK, l1
                    do i = 1_IK, ido
                        work(i, k, 1) = data(i, 1, k) + data(i, 2, k)
                        t2 = data(i, 1, k) - data(i, 2, k)
                        work(i, k, 2)%im = coef(ix1 + i)%re * t2%im + coef(ix1 + i)%im * t2%re
                        work(i, k, 2)%re = coef(ix1 + i)%re * t2%re - coef(ix1 + i)%im * t2%im
                    end do
                end do
            else
                do k = 1_IK, l1
                    work(1, k, 1) = data(1, 1, k) + data(1, 2, k)
                    work(1, k, 2) = data(1, 1, k) - data(1, 2, k)
                end do
            end if
        end subroutine

        pure subroutine passb3(ido, l1, ix1, ix2, coef, data, work)
            integer(IK) , intent(in)                    :: ido, l1, ix1, ix2
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(in)                    :: data(ido, 3, l1)
            complex(TKC), intent(out)                   :: work(ido, l1, 3)
            complex(TKC), parameter                     :: TAU = cmplx(-.5_TKC, sqrt(3._TKC) / 2._TKC, TKC)
            complex(TKC)                                :: c2, c3, d2, d3, t2
            integer(IK)                                 :: i, k
            if (ido /= 1_IK) then
                do k = 1_IK, l1
                    do i = 1_IK, ido
                        t2 = data(i, 2, k) + data(i, 3, k)
                        c2 = data(i, 1, k) + TAU%re * t2
                        work(i, k, 1) = data(i, 1, k) + t2
                        c3 = TAU%im * (data(i, 2, k) - data(i, 3, k))
                        d2%re = c2%re - c3%im
                        d2%im = c2%im + c3%re
                        d3%re = c2%re + c3%im
                        d3%im = c2%im - c3%re
                        work(i, k, 2)%re = coef(ix1 + i)%re * d2%re - coef(ix1 + i)%im * d2%im
                        work(i, k, 2)%im = coef(ix1 + i)%re * d2%im + coef(ix1 + i)%im * d2%re
                        work(i, k, 3)%re = coef(ix2 + i)%re * d3%re - coef(ix2 + i)%im * d3%im
                        work(i, k, 3)%im = coef(ix2 + i)%re * d3%im + coef(ix2 + i)%im * d3%re
                    end do
                end do
            else
                do k = 1_IK, l1
                    t2 = data(1, 2, k) + data(1, 3, k)
                    c2 = data(1, 1, k) + TAU%re * t2
                    work(1, k, 1) = data(1, 1, k) + t2
                    c3 = TAU%im * (data(1, 2, k) - data(1, 3, k))
                    work(1, k, 2)%re = c2%re - c3%im
                    work(1, k, 2)%im = c2%im + c3%re
                    work(1, k, 3)%re = c2%re + c3%im
                    work(1, k, 3)%im = c2%im - c3%re
                end do
            end if
        end subroutine

        pure subroutine passb4(ido, l1, ix1, ix2, ix3, coef, data, work)
            integer(IK) , intent(in)                    :: ido, l1, ix1, ix2, ix3
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(in)                    :: data(ido, 4, l1)
            complex(TKC), intent(out)                   :: work(ido, l1, 4)
            complex(TKC)                                :: c2, c3, c4, t1, t2, t3, t4
            integer(IK)                                 :: i, k
            if (ido /= 1_IK) then
                do k = 1_IK, l1
                    do i = 1_IK, ido, 1_IK
                        t1 = data(i, 1, k) - data(i, 3, k)
                        t2 = data(i, 1, k) + data(i, 3, k)
                        t3 = data(i, 2, k) + data(i, 4, k)
                        t4%re = data(i, 4, k)%im - data(i, 2, k)%im
                        t4%im = data(i, 2, k)%re - data(i, 4, k)%re
                        work(i, k, 1) = t2 + t3
                        c3 = t2 - t3
                        c2 = t1 + t4
                        c4 = t1 - t4
                        work(i, k, 2)%re = coef(ix1 + i)%re * c2%re - coef(ix1 + i)%im * c2%im
                        work(i, k, 2)%im = coef(ix1 + i)%re * c2%im + coef(ix1 + i)%im * c2%re
                        work(i, k, 3)%re = coef(ix2 + i)%re * c3%re - coef(ix2 + i)%im * c3%im
                        work(i, k, 3)%im = coef(ix2 + i)%re * c3%im + coef(ix2 + i)%im * c3%re
                        work(i, k, 4)%re = coef(ix3 + i)%re * c4%re - coef(ix3 + i)%im * c4%im
                        work(i, k, 4)%im = coef(ix3 + i)%re * c4%im + coef(ix3 + i)%im * c4%re
                    end do
                end do
            else
                do k = 1_IK, l1
                    t1 = data(1, 1, k) - data(1, 3, k)
                    t2 = data(1, 1, k) + data(1, 3, k)
                    t3 = data(1, 2, k) + data(1, 4, k)
                    t4%re = data(1, 4, k)%im - data(1, 2, k)%im
                    t4%im = data(1, 2, k)%re - data(1, 4, k)%re
                    work(1, k, 1) = t2 + t3
                    work(1, k, 3) = t2 - t3
                    work(1, k, 2) = t1 + t4
                    work(1, k, 4) = t1 - t4
                end do
            end if
        end subroutine passb4

        pure subroutine passb5(ido, l1, ix1, ix2, ix3, ix4, coef, data, work)
            integer(IK) , intent(in)                    :: ido, l1, ix1, ix2, ix3, ix4
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(in)                    :: data(ido, 5, l1)
            complex(TKC), intent(out)                   :: work(ido, l1, 5)
            complex(TKC)                                :: c2, c3, c4, c5, d2, d3, d4, d5, t2, t3, t4, t5
            real(TKC)   , parameter                     :: pi = acos( - 1._TKC)
            complex(TKC), parameter                     :: T11 = cmplx(cos(2._TKC * pi / 5._TKC), sin(2._TKC * pi / 5._TKC), TKC)
            complex(TKC), parameter                     :: T12 = cmplx(cos(4._TKC * pi / 5._TKC), sin(4._TKC * pi / 5._TKC), TKC)
            integer(IK)                                 :: i, k
            if (ido /= 1_IK) then
                do k = 1_IK, l1
                    do i = 1_IK, ido
                        t2 = data(i, 2, k) + data(i, 5, k)
                        t3 = data(i, 3, k) + data(i, 4, k)
                        t4 = data(i, 3, k) - data(i, 4, k)
                        t5 = data(i, 2, k) - data(i, 5, k)
                        work(i, k, 1) = data(i, 1, k) + t2 + t3
                        c2%re = data(i, 1, k)%re + T11%re * t2%re + T12%re * t3%re
                        c2%im = data(i, 1, k)%im + T11%re * t2%im + T12%re * t3%im
                        c3%re = data(i, 1, k)%re + T12%re * t2%re + T11%re * t3%re
                        c3%im = data(i, 1, k)%im + T12%re * t2%im + T11%re * t3%im
                        c5%re = T11%im * t5%re + T12%im * t4%re
                        c5%im = T11%im * t5%im + T12%im * t4%im
                        c4%re = T12%im * t5%re - T11%im * t4%re
                        c4%im = T12%im * t5%im - T11%im * t4%im
                        d3%re = c3%re - c4%im
                        d3%im = c3%im + c4%re
                        d4%re = c3%re + c4%im
                        d4%im = c3%im - c4%re
                        d5%re = c2%re + c5%im
                        d5%im = c2%im - c5%re
                        d2%re = c2%re - c5%im
                        d2%im = c2%im + c5%re
                        work(i, k, 2)%re = coef(ix1 + i)%re * d2%re - coef(ix1 + i)%im * d2%im
                        work(i, k, 2)%im = coef(ix1 + i)%re * d2%im + coef(ix1 + i)%im * d2%re
                        work(i, k, 3)%re = coef(ix2 + i)%re * d3%re - coef(ix2 + i)%im * d3%im
                        work(i, k, 3)%im = coef(ix2 + i)%re * d3%im + coef(ix2 + i)%im * d3%re
                        work(i, k, 4)%re = coef(ix3 + i)%re * d4%re - coef(ix3 + i)%im * d4%im
                        work(i, k, 4)%im = coef(ix3 + i)%re * d4%im + coef(ix3 + i)%im * d4%re
                        work(i, k, 5)%re = coef(ix4 + i)%re * d5%re - coef(ix4 + i)%im * d5%im
                        work(i, k, 5)%im = coef(ix4 + i)%re * d5%im + coef(ix4 + i)%im * d5%re
                    end do
                end do
            else
                do k = 1_IK, l1
                    t2 = data(1, 2, k) + data(1, 5, k)
                    t3 = data(1, 3, k) + data(1, 4, k)
                    t4 = data(1, 3, k) - data(1, 4, k)
                    t5 = data(1, 2, k) - data(1, 5, k)
                    work(1, k, 1) = data(1, 1, k) + t2 + t3
                    c2%re = data(1, 1, k)%re + T11%re * t2%re + T12%re * t3%re
                    c2%im = data(1, 1, k)%im + T11%re * t2%im + T12%re * t3%im
                    c3%re = data(1, 1, k)%re + T12%re * t2%re + T11%re * t3%re
                    c3%im = data(1, 1, k)%im + T12%re * t2%im + T11%re * t3%im
                    c5%re = T11%im * t5%re + T12%im * t4%re
                    c5%im = T11%im * t5%im + T12%im * t4%im
                    c4%re = T12%im * t5%re - T11%im * t4%re
                    c4%im = T12%im * t5%im - T11%im * t4%im
                    work(1, k, 2)%re = c2%re - c5%im
                    work(1, k, 2)%im = c2%im + c5%re
                    work(1, k, 3)%re = c3%re - c4%im
                    work(1, k, 3)%im = c3%im + c4%re
                    work(1, k, 4)%re = c3%re + c4%im
                    work(1, k, 4)%im = c3%im - c4%re
                    work(1, k, 5)%re = c2%re + c5%im
                    work(1, k, 5)%im = c2%im - c5%re
                end do
            end if
        end subroutine passb5

        pure subroutine passb(nac, ido, ip, l1, idl1, ix1, coef, data, C1, C2, work, Ch2)
            integer(IK) , intent(out)                   :: nac
            integer(IK) , intent(in)                    :: ido, ip, l1, idl1, ix1
            complex(TKC), intent(in)    , contiguous    :: coef(2:)
            complex(TKC), intent(in)                    :: data(ido, ip, l1)
            complex(TKC), intent(inout)                 :: C1(ido, l1, ip)
            complex(TKC), intent(inout)                 :: work(ido, l1, ip), C2(idl1,ip), Ch2(idl1,ip)
            integer(IK)                                 :: i, idij, idj, idl, idlj, idp, jk, inc, ipp2, ipph, j, jc, k, l, lc!, nt
            !idot = ido / 2_IK
            !nt = ip * idl1 what in the world is this doing here?
            idp = ip * ido
            ipp2 = ip + 2_IK
            ipph = (ip + 1_IK) / 2_IK
            if (2_IK * ido < l1) then
                do j = 2_IK, ipph
                    jc = ipp2 - j
                    do i = 1_IK, ido
                        do k = 1_IK, l1
                            work(i, k, j) = data(i, j, k) + data(i, jc, k)
                            work(i, k, jc) = data(i, j, k) - data(i, jc, k)
                        end do
                    end do
                end do
                do i = 1_IK, ido
                    do k = 1_IK, l1
                        work(i, k, 1) = data(i, 1, k)
                    end do
                end do
            else
                do j = 2_IK, ipph
                    jc = ipp2 - j
                    do k = 1_IK, l1
                        do i = 1_IK, ido
                            work(i, k, j) = data(i, j, k) + data(i, jc, k)
                            work(i, k, jc) = data(i, j, k) - data(i, jc, k)
                        end do
                    end do
                end do
                do k = 1_IK, l1
                    do i = 1_IK, ido
                        work(i, k, 1) = data(i, 1, k)
                    end do
                end do
            end if
            inc = 0_IK
            idl = 1_IK - ido
            do l = 2_IK, ipph
                lc = ipp2 - l
                idl = idl + ido
                do jk = 1_IK, idl1
                    C2(jk, l) = Ch2(jk, 1) + coef(ix1 + idl)%re * Ch2(jk, 2)
                    C2(jk, lc) = coef(ix1 + idl)%im * Ch2(jk, ip)
                end do
                idlj = idl
                inc = inc + ido
                do j = 3_IK, ipph
                    jc = ipp2 - j
                    idlj = idlj + inc
                    if (idlj > idp) idlj = idlj - idp
                    do jk = 1_IK, idl1
                        C2(jk, l) = C2(jk, l) + coef(ix1 + idlj)%re * Ch2(jk, j)
                        C2(jk, lc) = C2(jk, lc) + coef(ix1 + idlj)%im * Ch2(jk, jc)
                    end do
                end do
            end do
            do j = 2_IK, ipph
                do jk = 1_IK, idl1
                    Ch2(jk, 1) = Ch2(jk, 1) + Ch2(jk, j)
                end do
            end do
            do j = 2_IK, ipph
                jc = ipp2 - j
                do jk = 1_IK, idl1
                    Ch2(jk, j)%re = C2(jk, j)%re - C2(jk, jc)%im
                    Ch2(jk, j)%im = C2(jk, j)%im + C2(jk, jc)%re
                    Ch2(jk, jc)%re = C2(jk, j)%re + C2(jk, jc)%im
                    Ch2(jk, jc)%im = C2(jk, j)%im - C2(jk, jc)%re
                end do
            end do
            nac = 1_IK
            if (ido == 1_IK) return
            nac = 0_IK
            do jk = 1_IK, idl1
                C2(jk, 1) = Ch2(jk, 1)
            end do
            do j = 2_IK, ip
                do k = 1_IK, l1
                    C1(1, k, j) = work(1, k, j)
                end do
            end do
            if (ido > l1) then
                idj = 1_IK - ido
                do j = 2_IK, ip
                    idj = idj + ido
                    do k = 1_IK, l1
                        idij = idj
                        do i = 2_IK, ido
                            idij = idij + 1_IK
                            C1(i, k, j)%re = coef(ix1 + idij)%re * work(i, k, j)%re - coef(ix1 + idij)%im * work(i, k, j)%im
                            C1(i, k, j)%im = coef(ix1 + idij)%re * work(i, k, j)%im + coef(ix1 + idij)%im * work(i, k, j)%re
                        end do
                    end do
                end do
            else
                idij = 0_IK
                do j = 2_IK, ip
                    idij = idij + 1_IK
                    do i = 2_IK, ido
                        idij = idij + 2_IK
                        do k = 1_IK, l1
                            C1(i, k, j)%re = coef(ix1 + idij)%re * work(i, k, j)%re - coef(ix1 + idij)%im * work(i, k, j)%im
                            C1(i, k, j)%im = coef(ix1 + idij)%re * work(i, k, j)%im + coef(ix1 + idij)%im * work(i, k, j)%re
                        end do
                    end do
                end do
            end if
        end subroutine passb

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setFFTR_ENABLED && RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenData, idl1, ido, ip, ix1, ix2, ix3, ix4, k1, l1, l2, na, lenFactor
        lenFactor = size(factor, kind = IK)
        lenData = size(data, kind = IK)
        CHECK_ASSERTION(__LINE__, lenData == size(coef, kind = IK), SK_"@setFFTR(): The condition `size(data) == size(coef)` must hold. size(data), size(coef) = "//getStr([lenData, size(coef, kind = IK)]))
        CHECK_ASSERTION(__LINE__, lenData == size(work, kind = IK), SK_"@setFFTR(): The condition `size(data) == size(work)` must hold. size(data), size(work) = "//getStr([lenData, size(work, kind = IK)]))
        inwork = logical(1_IK < lenData, LK)
        if (.not. inwork) return
        na = 0_IK
        l1 = 1_IK
        ix1 = 1_IK
        do k1 = 1_IK, lenFactor
            ip = factor(k1)
            l2 = ip * l1
            ido = lenData / l2
            idl1 = ido * l1
            if (ip == 4_IK) then
                ix2 = ix1 + ido
                ix3 = ix2 + ido
                if (na /= 0_IK) then
                    call radb4(ido, l1, ix1, ix2, ix3, coef, work, data)
                else
                    call radb4(ido, l1, ix1, ix2, ix3, coef, data, work)
                end if
                na = 1_IK - na
            elseif (ip == 2_IK) then
                if (na /= 0_IK) then
                    call radb2(ido, l1, ix1, coef, work, data)
                else
                    call radb2(ido, l1, ix1, coef, data, work)
                end if
                na = 1_IK - na
            elseif (ip == 3_IK) then
                ix2 = ix1 + ido
                if (na /= 0_IK) then
                    call radb3(ido, l1, ix1, ix2, coef, work, data)
                else
                    call radb3(ido, l1, ix1, ix2, coef, data, work)
                end if
                na = 1_IK - na
            elseif (ip /= 5_IK) then
                if (na /= 0_IK) then
                    call radbg(ido, ip, l1, idl1, ix1, coef, work, work, work, data, data)
                else
                    call radbg(ido, ip, l1, idl1, ix1, coef, data, data, data, work, work)
                end if
                if (ido == 1_IK) na = 1_IK - na
            else
                ix2 = ix1 + ido
                ix3 = ix2 + ido
                ix4 = ix3 + ido
                if (na /= 0_IK) then
                    call radb5(ido, l1, ix1, ix2, ix3, ix4, coef, work, data)
                else
                    call radb5(ido, l1, ix1, ix2, ix3, ix4, coef, data, work)
                end if
                na = 1_IK - na
            end if
            l1 = l2
            ix1 = ix1 + (ip - 1_IK) * ido
        end do
        inwork = logical(na /= 0_IK, LK)
        !if (na == 0_IK) return
        !data(1 : lenData) = work(1 : lenData)

    contains

        pure subroutine radb2(ido, l1, ix1, coef, data, work)
            integer(IK) , intent(in)                :: ido, l1, ix1
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(in)                :: data(ido, 2, l1)
            real(TKC)   , intent(out)               :: work(ido, l1, 2)
            integer(IK)                             :: i, ic, idp2, k
            complex(TKC)                            :: t2
            do k = 1_IK, l1
                work(1, k, 1) = data(1, 1, k) + data(ido, 2, k)
                work(1, k, 2) = data(1, 1, k) - data(ido, 2, k)
            end do
            if (ido < 2_IK) return
            if (ido /= 2_IK) then
                idp2 = ido + 2_IK
                do k = 1_IK, l1
                    do i = 3_IK, ido, 2_IK
                        ic = idp2 - i
                        work(i - 1, k, 1) = data(i - 1, 1, k) + data(ic - 1, 2, k)
                        t2%re = data(i - 1, 1, k) - data(ic - 1, 2, k)
                        work(i, k, 1) = data(i, 1, k) - data(ic, 2, k)
                        t2%im = data(i, 1, k) + data(ic, 2, k)
                        work(i - 1, k, 2) = coef(ix1 + i - 2) * t2%re - coef(ix1 + i - 1) * t2%im
                        work(i, k, 2) = coef(ix1 + i - 2) * t2%im + coef(ix1 + i - 1) * t2%re
                    end do
                end do
                if (mod(ido, 2_IK) == 1_IK) return
            end if
            do k = 1_IK, l1
                work(ido, k, 1) = data(ido, 1, k) + data(ido, 1, k)
                work(ido, k, 2) = -(data(1, 2, k) + data(1, 2, k))
            end do
        end subroutine

        pure subroutine radb3(ido, l1, ix1, ix2, coef, data, work)
            integer(IK) , intent(in)                :: ido, l1, ix1, ix2
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(in)                :: data(ido, 3, l1)
            real(TKC)   , intent(out)               :: work(ido, l1, 3)
            complex(TKC), parameter                 :: TAU = cmplx(-.5_TKC, sqrt(3._TKC) / 2._TKC, TKC)
            complex(TKC)                            :: c2, c3, d2, d3, t2
            integer(IK)                             :: i, ic, idp2, k
            do k = 1_IK, l1
                t2%re = data(ido, 2, k) + data(ido, 2, k)
                c2%re = data(1, 1, k) + TAU%re * t2%re
                work(1, k, 1) = data(1, 1, k) + t2%re
                c3%im = TAU%im * (data(1, 3, k) + data(1, 3, k))
                work(1, k, 2) = c2%re - c3%im
                work(1, k, 3) = c2%re + c3%im
            end do
            if (ido == 1_IK) return
            idp2 = ido + 2_IK
            do k = 1_IK, l1
                do i = 3_IK, ido, 2_IK
                    ic = idp2 - i
                    t2%re = data(i - 1, 3, k) + data(ic - 1, 2, k)
                    c2%re = data(i - 1, 1, k) + TAU%re * t2%re
                    work(i - 1, k, 1) = data(i - 1, 1, k) + t2%re
                    t2%im = data(i, 3, k) - data(ic, 2, k)
                    c2%im = data(i, 1, k) + TAU%re * t2%im
                    work(i, k, 1) = data(i, 1, k) + t2%im
                    c3%re = TAU%im * (data(i - 1, 3, k) - data(ic - 1, 2, k))
                    c3%im = TAU%im * (data(i, 3, k) + data(ic, 2, k))
                    d2%re = c2%re - c3%im
                    d3%re = c2%re + c3%im
                    d2%im = c2%im + c3%re
                    d3%im = c2%im - c3%re
                    work(i - 1, k, 2) = coef(ix1 + i - 2) * d2%re - coef(ix1 + i - 1) * d2%im
                    work(i, k, 2) = coef(ix1 + i - 2) * d2%im + coef(ix1 + i - 1) * d2%re
                    work(i - 1, k, 3) = coef(ix2 + i - 2) * d3%re - coef(ix2 + i - 1) * d3%im
                    work(i, k, 3) = coef(ix2 + i - 2) * d3%im + coef(ix2 + i - 1) * d3%re
                end do
            end do
        end subroutine

        pure subroutine radb4(ido, l1, ix1, ix2, ix3, coef, data, work)
            integer(IK) , intent(in)                :: ido, l1, ix1, ix2, ix3
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(in)                :: data(ido, 4, l1)
            real(TKC)   , intent(out)               :: work(ido, l1, 4)
            real(TKC)   , parameter                 :: NSQT2 = -sqrt(2._TKC)
            complex(TKC)                            :: c2, c3, c4, t1, t2, t3, t4
            integer(IK)                             :: i, ic, idp2, k
            do k = 1_IK, l1
                t1%re = data(1, 1, k) - data(ido, 4, k)
                t2%re = data(1, 1, k) + data(ido, 4, k)
                t3%re = data(ido, 2, k) + data(ido, 2, k)
                t4%re = data(1, 3, k) + data(1, 3, k)
                work(1, k, 1) = t2%re + t3%re
                work(1, k, 2) = t1%re - t4%re
                work(1, k, 3) = t2%re - t3%re
                work(1, k, 4) = t1%re + t4%re
            end do
            if (ido < 2_IK) return
            if (ido /= 2_IK) then
                idp2 = ido + 2_IK
                do k = 1_IK, l1
                    do i = 3_IK, ido, 2_IK
                        ic = idp2 - i
                        t1%im = data(i, 1, k) + data(ic, 4, k)
                        t2%im = data(i, 1, k) - data(ic, 4, k)
                        t3%im = data(i, 3, k) - data(ic, 2, k)
                        t4%re = data(i, 3, k) + data(ic, 2, k)
                        t1%re = data(i - 1, 1, k) - data(ic - 1, 4, k)
                        t2%re = data(i - 1, 1, k) + data(ic - 1, 4, k)
                        t4%im = data(i - 1, 3, k) - data(ic - 1, 2, k)
                        t3%re = data(i - 1, 3, k) + data(ic - 1, 2, k)
                        work(i - 1, k, 1) = t2%re + t3%re
                        c3%re = t2%re - t3%re
                        work(i, k, 1) = t2%im + t3%im
                        c3%im = t2%im - t3%im
                        c2%re = t1%re - t4%re
                        c4%re = t1%re + t4%re
                        c2%im = t1%im + t4%im
                        c4%im = t1%im - t4%im
                        work(i - 1, k, 2) = coef(ix1 + i - 2) * c2%re - coef(ix1 + i - 1) * c2%im
                        work(i, k, 2) = coef(ix1 + i - 2) * c2%im + coef(ix1 + i - 1) * c2%re
                        work(i - 1, k, 3) = coef(ix2 + i - 2) * c3%re - coef(ix2 + i - 1) * c3%im
                        work(i, k, 3) = coef(ix2 + i - 2) * c3%im + coef(ix2 + i - 1) * c3%re
                        work(i - 1, k, 4) = coef(ix3 + i - 2) * c4%re - coef(ix3 + i - 1) * c4%im
                        work(i, k, 4) = coef(ix3 + i - 2) * c4%im + coef(ix3 + i - 1) * c4%re
                    end do
                end do
                if (mod(ido, 2_IK) == 1_IK) return
            end if
            do k = 1_IK, l1
                t1%im = data(1, 2, k) + data(1, 4, k)
                t2%im = data(1, 4, k) - data(1, 2, k)
                t1%re = data(ido, 1, k) - data(ido, 3, k)
                t2%re = data(ido, 1, k) + data(ido, 3, k)
                work(ido, k, 1) = t2%re + t2%re
                work(ido, k, 2) = NSQT2 * (t1%im - t1%re)
                work(ido, k, 3) = t2%im + t2%im
                work(ido, k, 4) = NSQT2 * (t1%re + t1%im)
            end do
        end subroutine radb4

        pure subroutine radb5(ido, l1, ix1, ix2, ix3, ix4, coef, data, work)
            integer(IK) , intent(in)                :: ido, l1, ix1, ix2, ix3, ix4
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(in)                :: data(ido, 5, l1)
            real(TKC)   , intent(out)               :: work(ido, l1, 5)
            real(TKC)   , parameter                 :: PI = acos(-1._TKC)
            complex(TKC), parameter                 :: T11 = cmplx(cos(2._TKC * PI / 5._TKC), sin(2._TKC * PI / 5._TKC), TKC)
            complex(TKC), parameter                 :: T12 = cmplx(cos(4._TKC * PI / 5._TKC), sin(4._TKC * PI / 5._TKC), TKC)
            complex(TKC)                            :: c2, c3, c4, c5, d2, d3, d4, d5, t2, t3, t4, t5
            integer(IK)                             :: i, ic, idp2, k
            do k = 1_IK, l1
                t5%im = data(1, 3, k) + data(1, 3, k)
                t4%im = data(1, 5, k) + data(1, 5, k)
                t2%re = data(ido, 2, k) + data(ido, 2, k)
                t3%re = data(ido, 4, k) + data(ido, 4, k)
                work(1, k, 1) = data(1, 1, k) + t2%re + t3%re
                c2%re = data(1, 1, k) + T11%re * t2%re + T12%re * t3%re
                c3%re = data(1, 1, k) + T12%re * t2%re + T11%re * t3%re
                c5%im = T11%im * t5%im + T12%im * t4%im
                c4%im = T12%im * t5%im - T11%im * t4%im
                work(1, k, 2) = c2%re - c5%im
                work(1, k, 3) = c3%re - c4%im
                work(1, k, 4) = c3%re + c4%im
                work(1, k, 5) = c2%re + c5%im
            end do
            if (ido == 1_IK) return
            idp2 = ido + 2_IK
            do k = 1_IK, l1
                do i = 3_IK, ido, 2_IK
                    ic = idp2 - i
                    t5%im = data(i, 3, k) + data(ic, 2, k)
                    t2%im = data(i, 3, k) - data(ic, 2, k)
                    t4%im = data(i, 5, k) + data(ic, 4, k)
                    t3%im = data(i, 5, k) - data(ic, 4, k)
                    t5%re = data(i - 1, 3, k) - data(ic - 1, 2, k)
                    t2%re = data(i - 1, 3, k) + data(ic - 1, 2, k)
                    t4%re = data(i - 1, 5, k) - data(ic - 1, 4, k)
                    t3%re = data(i - 1, 5, k) + data(ic - 1, 4, k)
                    work(i - 1, k, 1) = data(i - 1, 1, k) + t2%re + t3%re
                    work(i, k, 1) = data(i, 1, k) + t2%im + t3%im
                    c2%re = data(i - 1, 1, k) + T11%re * t2%re + T12%re * t3%re
                    c2%im = data(i, 1, k) + T11%re * t2%im + T12%re * t3%im
                    c3%re = data(i - 1, 1, k) + T12%re * t2%re + T11%re * t3%re
                    c3%im = data(i, 1, k) + T12%re * t2%im + T11%re * t3%im
                    c5%re = T11%im * t5%re + T12%im * t4%re
                    c5%im = T11%im * t5%im + T12%im * t4%im
                    c4%re = T12%im * t5%re - T11%im * t4%re
                    c4%im = T12%im * t5%im - T11%im * t4%im
                    d3%re = c3%re - c4%im
                    d4%re = c3%re + c4%im
                    d3%im = c3%im + c4%re
                    d4%im = c3%im - c4%re
                    d5%re = c2%re + c5%im
                    d2%re = c2%re - c5%im
                    d5%im = c2%im - c5%re
                    d2%im = c2%im + c5%re
                    work(i - 1, k, 2) = coef(ix1 + i - 2) * d2%re - coef(ix1 + i - 1) * d2%im
                    work(i, k, 2) = coef(ix1 + i - 2) * d2%im + coef(ix1 + i - 1) * d2%re
                    work(i - 1, k, 3) = coef(ix2 + i - 2) * d3%re - coef(ix2 + i - 1) * d3%im
                    work(i, k, 3) = coef(ix2 + i - 2) * d3%im + coef(ix2 + i - 1) * d3%re
                    work(i - 1, k, 4) = coef(ix3 + i - 2) * d4%re - coef(ix3 + i - 1) * d4%im
                    work(i, k, 4) = coef(ix3 + i - 2) * d4%im + coef(ix3 + i - 1) * d4%re
                    work(i - 1, k, 5) = coef(ix4 + i - 2) * d5%re - coef(ix4 + i - 1) * d5%im
                    work(i, k, 5) = coef(ix4 + i - 2) * d5%im + coef(ix4 + i - 1) * d5%re
                end do
            end do
        end subroutine

        pure subroutine radbg(ido, ip, l1, idl1, ix1, coef, data, C1, C2, work, Ch2)
            integer(IK) , intent(in)                :: ido, ip, l1, idl1, ix1
            real(TKC)   , intent(in), contiguous    :: coef(2:)
            real(TKC)   , intent(in)                :: data(ido,ip,l1)
            real(TKC)   , intent(inout)             :: C1(ido,l1,ip)
            real(TKC)   , intent(inout)             :: work(ido,l1,ip), C2(idl1,ip), Ch2(idl1,ip)
            real(TKC)   , parameter                 :: TWO_PI = 2._TKC * acos(-1._TKC)
            real(TKC)                               :: ar1h, ar2h, arg, dc2, dcp, ds2, dsp
            integer(IK)                             :: i, ic, idij, idp2, jk, ipp2, ipph, is, j, j2, jc, k, l, lc, nbd
            complex(TKC)                            :: a1, a2
            arg = TWO_PI / ip
            dcp = cos(arg)
            dsp = sin(arg)
            idp2 = ido + 2_IK
            nbd = (ido - 1_IK) / 2_IK
            ipp2 = ip + 2
            ipph = (ip + 1_IK) / 2_IK
            if (ido < l1) then
                do i = 1_IK, ido
                    do k = 1_IK, l1
                        work(i, k, 1) = data(i, 1, k)
                    end do
                end do
            else
                do k = 1_IK, l1
                    do i = 1_IK, ido
                    work(i, k, 1) = data(i, 1, k)
                    end do
                end do
            end if
            do j = 2_IK, ipph
                jc = ipp2 - j
                j2 = j + j
                do k = 1_IK, l1
                    work(1, k, j) = data(ido, j2 - 2, k) + data(ido, j2 - 2, k)
                    work(1, k, jc) = data(1, j2 - 1, k) + data(1, j2 - 1, k)
                end do
            end do
            if (ido /= 1_IK) then
                if (nbd < l1) then
                    do j = 2_IK, ipph
                        jc = ipp2 - j
                        do i = 3_IK, ido, 2_IK
                            ic = idp2 - i
                            do k = 1_IK, l1
                                work(i - 1, k, j) = data(i - 1, 2 * j - 1, k) + data(ic - 1, 2 * j - 2, k)
                                work(i - 1, k, jc) = data(i - 1, 2 * j - 1, k) - data(ic - 1, 2 * j - 2, k)
                                work(i, k, j) = data(i, 2 * j - 1, k) - data(ic, 2 * j - 2, k)
                                work(i, k, jc) = data(i, 2 * j - 1, k) + data(ic, 2 * j - 2, k)
                            end do
                        end do
                    end do
                else
                    do j = 2_IK, ipph
                        jc = ipp2 - j
                        do k = 1_IK, l1
                            do i = 3_IK, ido, 2_IK
                                ic = idp2 - i
                                work(i - 1, k, j) = data(i - 1, 2 * j - 1, k) + data(ic - 1, 2 * j - 2, k)
                                work(i - 1, k, jc) = data(i - 1, 2 * j - 1, k) - data(ic - 1, 2 * j - 2, k)
                                work(i, k, j) = data(i, 2 * j - 1, k) - data(ic, 2 * j - 2, k)
                                work(i, k, jc) = data(i, 2 * j - 1, k) + data(ic, 2 * j - 2, k)
                            end do
                        end do
                    end do
                end if
            end if
            a1%re = 1._TKC
            a1%im = 0._TKC
            do l = 2_IK, ipph
                lc = ipp2 - l
                ar1h = dcp * a1%re - dsp * a1%im
                a1%im = dcp * a1%im + dsp * a1%re
                a1%re = ar1h
                do jk = 1_IK, idl1
                    C2(jk, l) = Ch2(jk, 1) + a1%re * Ch2(jk, 2)
                    C2(jk, lc) = a1%im * Ch2(jk, ip)
                end do
                dc2 = a1%re
                ds2 = a1%im
                a2%re = a1%re
                a2%im = a1%im
                do j = 3_IK, ipph
                    jc = ipp2 - j
                    ar2h = dc2 * a2%re - ds2 * a2%im
                    a2%im = dc2 * a2%im + ds2 * a2%re
                    a2%re = ar2h
                    do jk = 1_IK, idl1
                        C2(jk, l) = C2(jk, l) + a2%re * Ch2(jk, j)
                        C2(jk, lc) = C2(jk, lc) + a2%im * Ch2(jk, jc)
                    end do
                end do
            end do
            do j = 2_IK, ipph
                do jk = 1_IK, idl1
                    Ch2(jk, 1) = Ch2(jk, 1) + Ch2(jk, j)
                end do
            end do
            do j = 2_IK, ipph
                jc = ipp2 - j
                do k = 1_IK, l1
                    work(1, k, j) = C1(1, k, j) - C1(1, k, jc)
                    work(1, k, jc) = C1(1, k, j) + C1(1, k, jc)
                end do
            end do
            if (ido /= 1_IK) then
                if (nbd < l1) then
                    do j = 2_IK, ipph
                        jc = ipp2 - j
                        do i = 3_IK, ido, 2_IK
                            do k = 1_IK, l1
                                work(i - 1, k, j) = C1(i - 1, k, j) - C1(i, k, jc)
                                work(i - 1, k, jc) = C1(i - 1, k, j) + C1(i, k, jc)
                                work(i, k, j) = C1(i, k, j) + C1(i - 1, k, jc)
                                work(i, k, jc) = C1(i, k, j) - C1(i - 1, k, jc)
                            end do
                        end do
                    end do
                else
                    do j = 2_IK, ipph
                        jc = ipp2 - j
                        do k = 1_IK, l1
                            do i = 3_IK, ido, 2_IK
                                work(i - 1, k, j) = C1(i - 1, k, j) - C1(i, k, jc)
                                work(i - 1, k, jc) = C1(i - 1, k, j) + C1(i, k, jc)
                                work(i, k, j) = C1(i, k, j) + C1(i - 1, k, jc)
                                work(i, k, jc) = C1(i, k, j) - C1(i - 1, k, jc)
                            end do
                        end do
                    end do
                end if
            end if
            if (ido == 1_IK) return
            C2(1 : idl1, 1) = Ch2(1 : idl1, 1)
            C1(1, 1 : l1, 2 : ip) = work(1, 1 : l1, 2 : ip)
            if (nbd > l1) then
                is = -ido
                do j = 2_IK, ip
                    is = is + ido
                    do k = 1_IK, l1
                        idij = is
                        do i = 3_IK, ido, 2_IK
                            idij = idij + 2_IK
                            C1(i - 1, k, j) = coef(ix1 + idij - 1) * work(i - 1, k, j) - coef(ix1 + idij) * work(i, k, j)
                            C1(i, k, j) = coef(ix1 + idij - 1) * work(i, k, j) + coef(ix1 + idij) * work(i - 1, k, j)
                        end do
                    end do
                end do
            else
                is = -ido
                do j = 2_IK, ip
                    is = is + ido
                    idij = is
                    do i = 3_IK, ido, 2_IK
                        idij = idij + 2_IK
                        do k = 1_IK, l1
                            C1(i - 1, k, j) = coef(ix1 + idij - 1) * work(i - 1, k, j) - coef(ix1 + idij) * work(i, k, j)
                            C1(i, k, j) = coef(ix1 + idij - 1) * work(i, k, j) + coef(ix1 + idij) * work(i - 1, k, j)
                        end do
                    end do
                end do
            end if
        end subroutine

        !%%%%%%%%%%%%%%
#elif   setFFTI_ENABLED
        !%%%%%%%%%%%%%%

        integer(IK) :: i
        real(TKC) :: invLenData
        call setFFTR(factor, coef, data, work, inwork)
        invLenData = 1._TKC / real(size(data, 1, IK), TKC)
        if (inwork) then
            do concurrent(i = 1 : size(data, 1, IK))
                work(i) = work(i) * invLenData
            end do
        else
            do concurrent(i = 1 : size(data, 1, IK))
                data(i) = data(i) * invLenData
            end do
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getFFTF_ENABLED || getFFTR_ENABLED || getFFTI_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
        complex(TKC), dimension(size(data, 1, IK)) :: coef, work
#elif   RK_ENABLED
        real(TKC), dimension(size(data, 1, IK)) :: coef, work
#else
#error  "Unrecognized interface."
#endif
        logical(LK) :: inwork
        integer(IK) , allocatable :: factor(:)
        factor = getFactorFFT(data, coef)
        fft = data
#if     getFFTF_ENABLED
        call setFFTF(factor, coef, fft, work, inwork)
#elif   getFFTR_ENABLED
        call setFFTR(factor, coef, fft, work, inwork)
#elif   getFFTI_ENABLED
        call setFFTI(factor, coef, fft, work, inwork)
#else
#error  "Unrecognized interface."
#endif
        if (inwork) fft = work
        ! Normalize the result.
!#if     getFFTF_ENABLED || getFFTR_ENABLED
!        if (inwork) fft = work
!#elif   getFFTI_ENABLED
!        block
!            integer(IK) :: i
!            real(TKC) :: invLenData
!            invLenData = 1._TKC / real(size(data, 1, IK), TKC)
!            if (inwork) then
!                do concurrent(i = 1 : size(data, 1, IK))
!                    fft(i) = fft(i) * invLenData
!                end do
!            else
!                do concurrent(i = 1 : size(data, 1, IK))
!                    fft(i) = work(i) * invLenData
!                end do
!            end if
!        end block
!#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif