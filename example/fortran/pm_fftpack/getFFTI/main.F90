program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_fftpack, only: getfactorFFT
    use pm_fftpack, only: getFFTF, getFFTI
    use pm_distUnif, only: getUnifRand
    use pm_mathCompare, only: isClose
    use pm_err, only: setAsserted

    implicit none

    integer(IK) :: i
    integer(IK) :: lenData
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    do i = 1, 5
    block
        use pm_kind, only: TKC => CK32
        real(TKC) :: reltol
        complex(TKC), allocatable :: data(:), adat(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKC + getUnifRand((0._TKC, 0._TKC), (1._TKC, 1._TKC), lenData)")
                        data = 1._TKC + getUnifRand((0._TKC, 0._TKC), (1._TKC, 1._TKC), lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("adat = getFFTI(getFFTF(data))")
                        adat = getFFTI(getFFTF(data))
        call disp%show("adat")
        call disp%show( adat )
        call disp%show("reltol = sqrt(epsilon(1._TKC))")
                        reltol = sqrt(epsilon(1._TKC))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, adat, reltol = reltol)")
        call disp%show( isClose(data, adat, reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, adat, reltol = reltol)))")
                        call setAsserted(all(isClose(data, adat, reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => CK64
        real(TKC) :: reltol
        complex(TKC), allocatable :: data(:), adat(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKC + getUnifRand((0._TKC, 0._TKC), (1._TKC, 1._TKC), lenData)")
                        data = 1._TKC + getUnifRand((0._TKC, 0._TKC), (1._TKC, 1._TKC), lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("adat = getFFTI(getFFTF(data))")
                        adat = getFFTI(getFFTF(data))
        call disp%show("adat")
        call disp%show( adat )
        call disp%show("reltol = sqrt(epsilon(1._TKC))")
                        reltol = sqrt(epsilon(1._TKC))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, adat, reltol = reltol)")
        call disp%show( isClose(data, adat, reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, adat, reltol = reltol)))")
                        call setAsserted(all(isClose(data, adat, reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => CK128
        real(TKC) :: reltol
        complex(TKC), allocatable :: data(:), adat(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKC + getUnifRand((0._TKC, 0._TKC), (1._TKC, 1._TKC), lenData)")
                        data = 1._TKC + getUnifRand((0._TKC, 0._TKC), (1._TKC, 1._TKC), lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("adat = getFFTI(getFFTF(data))")
                        adat = getFFTI(getFFTF(data))
        call disp%show("adat")
        call disp%show( adat )
        call disp%show("reltol = sqrt(epsilon(1._TKC))")
                        reltol = sqrt(epsilon(1._TKC))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, adat, reltol = reltol)")
        call disp%show( isClose(data, adat, reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, adat, reltol = reltol)))")
                        call setAsserted(all(isClose(data, adat, reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RK32
        real(TKC) :: reltol
        real(TKC), allocatable :: data(:), adat(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKC + getUnifRand(0._TKC, 1._TKC, lenData)")
                        data = 1._TKC + getUnifRand(0._TKC, 1._TKC, lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("adat = getFFTI(getFFTF(data))")
                        adat = getFFTI(getFFTF(data))
        call disp%show("adat")
        call disp%show( adat )
        call disp%show("reltol = sqrt(epsilon(1._TKC))")
                        reltol = sqrt(epsilon(1._TKC))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, adat, reltol = reltol)")
        call disp%show( isClose(data, adat, reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, adat, reltol = reltol)))")
                        call setAsserted(all(isClose(data, adat, reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RK64
        real(TKC) :: reltol
        real(TKC), allocatable :: data(:), adat(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKC + getUnifRand(0._TKC, 1._TKC, lenData)")
                        data = 1._TKC + getUnifRand(0._TKC, 1._TKC, lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("adat = getFFTI(getFFTF(data))")
                        adat = getFFTI(getFFTF(data))
        call disp%show("adat")
        call disp%show( adat )
        call disp%show("reltol = sqrt(epsilon(1._TKC))")
                        reltol = sqrt(epsilon(1._TKC))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, adat, reltol = reltol)")
        call disp%show( isClose(data, adat, reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, adat, reltol = reltol)))")
                        call setAsserted(all(isClose(data, adat, reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RK128
        real(TKC) :: reltol
        real(TKC), allocatable :: data(:), adat(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKC + getUnifRand(0._TKC, 1._TKC, lenData)")
                        data = 1._TKC + getUnifRand(0._TKC, 1._TKC, lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("adat = getFFTI(getFFTF(data))")
                        adat = getFFTI(getFFTF(data))
        call disp%show("adat")
        call disp%show( adat )
        call disp%show("reltol = sqrt(epsilon(1._TKC))")
                        reltol = sqrt(epsilon(1._TKC))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, adat, reltol = reltol)")
        call disp%show( isClose(data, adat, reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, adat, reltol = reltol)))")
                        call setAsserted(all(isClose(data, adat, reltol = reltol)))
        call disp%skip()
    end block
    end do

end program example