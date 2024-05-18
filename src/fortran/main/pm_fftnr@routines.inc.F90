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

#if     setFFTF_ENABLED
        real(TKC), parameter :: trifac = +2 * acos(-1._TKC)
#elif   setFFTR_ENABLED
        real(TKC), parameter :: trifac = -2 * acos(-1._TKC)
#endif
        !%%%%%%%%%%%%%%
#if     setFFTI_ENABLED
        !%%%%%%%%%%%%%%

        call setFFTR(data)
#if     CK_ENABLED
        data = data / size(data, 1, IK)
#elif   RK_ENABLED
        data = data * (2._TKC / size(data, 1, IK))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getFFTF_ENABLED || getFFTR_ENABLED || getFFTI_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fft(1 : size(data, 1, IK)) = data
        fft(size(data, 1, IK) + 1 :) = 0._TKC
#if     getFFTF_ENABLED
        call setFFTF(fft)
#elif   getFFTI_ENABLED
        call setFFTI(fft)
#elif   getFFTR_ENABLED
        call setFFTR(fft)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (setFFTF_ENABLED || setFFTR_ENABLED) && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lendata, lenDataHalf, i, j, m, step, maxm
        complex(TKC) :: temp, w, wp
        real(TKC) :: theta, wtmp
        lendata = size(data, 1, IK)
        CHECK_ASSERTION(__LINE__, isIntPow(lenData) .and. 1 < lenData, SK_"@setFFTF/setFFTR(): The condition `isIntPow(size(data)) .and. 1 < size(data)` must hold. size(data) = "//getStr(lendata))
        lenDataHalf = lenData / 2_IK
        j = 1_IK
        do i = 1, lendata
            if (i < j) then
                ! swap elements.
                temp = data(j)
                data(j) = data(i)
                data(i) = temp
            endif
            m = lenDataHalf
            do
                if (m < 2_IK .or. j <= m) exit
                j = j - m
                m = m / 2_IK
            end do
            j = j + m
        end do
        maxm = 1_IK
        ! Danielson-Lanczos iteration. Repeat `log2(lenDataHalf)` times.
        do
            if (lendata <= maxm) exit
            step = 2_IK * maxm
            theta = trifac / step
            wp%re = -2._TKC * sin(0.5_TKC * theta)**2
            wp%im = sin(theta)
            w = (1._TKC, 0._TKC)
            do m = 1, maxm
                do i = m, lendata, step
                    j = i + maxm
                    temp%re = w%re * data(j)%re - w%im * data(j)%im
                    temp%im = w%re * data(j)%im + w%im * data(j)%re
                    data(j) = data(i) - temp
                    data(i) = data(i) + temp
                end do
                wtmp = w%re
                w%re = w%re * wp%re - w%im * wp%im + w%re
                w%im = w%im * wp%re + wtmp * wp%im + w%im
            end do
            maxm = step
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (setFFTF_ENABLED || setFFTR_ENABLED) && RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        complex(TKC) :: h1, h2, w, wp
        real(TKC) :: c1, c2, theta, wtmp
        integer(IK) :: i, i1, i2, i3, i4, lenData, lenDataPlus3
        lenData = size(data, 1, IK)
        CHECK_ASSERTION(__LINE__, isIntPow(lenData) .and. 1 < lenData, SK_"@setFFTF/setFFTR(): The condition `isIntPow(size(data)) .and. 1 < size(data)` must hold. size(data) = "//getStr(lendata))
        theta = trifac / real(lenData, TKC)
        c1 = 0.5_TKC
        c2 = 0.5_TKC
#if     setFFTF_ENABLED
        call setFFTC(data)
        c2 = -c2
#endif
        wp%re = -2._TKC * sin(0.5_TKC * theta)**2
        wp%im = sin(theta)
        w%im = wp%im
        w%re = 1._TKC + wp%re
        lenDataPlus3 = lenData + 3_IK
        do i = 2_IK, lenData / 4_IK
            i1 = 2_IK * i - 1_IK
            i2 = i1 + 1_IK
            i3 = lenDataPlus3 - i2
            i4 = i3 + 1_IK
            h1%re = +c1 * (data(i1) + data(i3))
            h1%im = +c1 * (data(i2) - data(i4))
            h2%re = -c2 * (data(i2) + data(i4))
            h2%im = +c2 * (data(i1) - data(i3))
            data(i1) = +h1%re + w%re * h2%re - w%im * h2%im
            data(i2) = +h1%im + w%re * h2%im + w%im * h2%re
            data(i3) = +h1%re - w%re * h2%re + w%im * h2%im
            data(i4) = -h1%im + w%re * h2%im + w%im * h2%re
            wtmp = w%re
            w%re = w%re * wp%re - w%im * wp%im + w%re
            w%im = w%im * wp%re + wtmp * wp%im + w%im
        end do
#if     setFFTF_ENABLED
        h1%re = data(1)
        data(1) = h1%re + data(2)
        data(2) = h1%re - data(2)
#elif   setFFTR_ENABLED
        h1%re = data(1)
        data(1) = c1 * (h1%re + data(2))
        data(2) = c1 * (h1%re - data(2))
        call setFFTC(data)
#else
#error  "Unrecognized interface."
#endif
    contains
        pure subroutine setFFTC(data)
            real(TKC), intent(inout), contiguous :: data(:)
            integer(IK) :: i, j, step, m, maxm, lenData, lenDataHalf
            real(TKC) :: tempi, tempr, theta, wi, wpi, wpr, wr, wtemp
            lenData = size(data, 1, IK)
            lenDataHalf = lenData / 2_IK
            j = 1_IK
            do i = 1_IK, lenData, 2_IK
                if(i < j)then
                    tempr = data(j)
                    tempi = data(j + 1)
                    data(j) = data(i)
                    data(j + 1) = data(i + 1)
                    data(i) = tempr
                    data(i + 1) = tempi
                endif
                m = lenDataHalf
                do
                    if (m < 2_IK .or. j <= m) exit
                    j = j - m
                    m = m / 2_IK
                end do
                j = j + m
            end do
            maxm = 2_IK
            do
                if (lenData <= maxm) exit
                step = 2_IK * maxm
                theta = trifac / maxm
                wpr = -2._TKC * sin(0.5_TKC * theta)**2
                wpi = sin(theta)
                wr = 1._TKC
                wi = 0._TKC
                do m = 1_IK, maxm, 2_IK
                    do i = m, lenData, step
                        j = i + maxm
                        tempr = wr * data(j) - wi * data(j + 1)
                        tempi = wr * data(j + 1) + wi * data(j)
                        data(j) = data(i) - tempr
                        data(j + 1) = data(i + 1) - tempi
                        data(i) = data(i) + tempr
                        data(i + 1) = data(i + 1) + tempi
                    end do
                    wtemp = wr
                    wr = wr * wpr - wi * wpi + wr
                    wi = wi * wpr + wtemp * wpi + wi
                end do
                maxm = step
            end do
        end subroutine
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif