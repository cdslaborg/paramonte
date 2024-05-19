program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_fftnr, only: getFFTF, getFFTI
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
        use pm_kind, only: TKG => CKS
        real(TKG) :: reltol
        complex(TKG), allocatable :: data(:), atad(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKG + getUnifRand((0._TKG, 0._TKG), (1._TKG, 1._TKG), lenData)")
                        data = 1._TKG + getUnifRand((0._TKG, 0._TKG), (1._TKG, 1._TKG), lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("atad = getFFTI(getFFTF(data))")
                        atad = getFFTI(getFFTF(data))
        call disp%show("atad")
        call disp%show( atad )
        call disp%show("reltol = sqrt(epsilon(1._TKG))")
                        reltol = sqrt(epsilon(1._TKG))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, atad(1:size(data)), reltol = reltol)")
        call disp%show( isClose(data, atad(1:size(data)), reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))")
                        call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => CKD
        real(TKG) :: reltol
        complex(TKG), allocatable :: data(:), atad(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKG + getUnifRand((0._TKG, 0._TKG), (1._TKG, 1._TKG), lenData)")
                        data = 1._TKG + getUnifRand((0._TKG, 0._TKG), (1._TKG, 1._TKG), lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("atad = getFFTI(getFFTF(data))")
                        atad = getFFTI(getFFTF(data))
        call disp%show("atad")
        call disp%show( atad )
        call disp%show("reltol = sqrt(epsilon(1._TKG))")
                        reltol = sqrt(epsilon(1._TKG))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, atad(1:size(data)), reltol = reltol)")
        call disp%show( isClose(data, atad(1:size(data)), reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))")
                        call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => CKH
        real(TKG) :: reltol
        complex(TKG), allocatable :: data(:), atad(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKG + getUnifRand((0._TKG, 0._TKG), (1._TKG, 1._TKG), lenData)")
                        data = 1._TKG + getUnifRand((0._TKG, 0._TKG), (1._TKG, 1._TKG), lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("atad = getFFTI(getFFTF(data))")
                        atad = getFFTI(getFFTF(data))
        call disp%show("atad")
        call disp%show( atad )
        call disp%show("reltol = sqrt(epsilon(1._TKG))")
                        reltol = sqrt(epsilon(1._TKG))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, atad(1:size(data)), reltol = reltol)")
        call disp%show( isClose(data, atad(1:size(data)), reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))")
                        call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKS
        real(TKG) :: reltol
        real(TKG), allocatable :: data(:), atad(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKG + getUnifRand(0._TKG, 1._TKG, lenData)")
                        data = 1._TKG + getUnifRand(0._TKG, 1._TKG, lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("atad = getFFTI(getFFTF(data))")
                        atad = getFFTI(getFFTF(data))
        call disp%show("atad")
        call disp%show( atad )
        call disp%show("reltol = sqrt(epsilon(1._TKG))")
                        reltol = sqrt(epsilon(1._TKG))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, atad(1:size(data)), reltol = reltol)")
        call disp%show( isClose(data, atad(1:size(data)), reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))")
                        call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKD
        real(TKG) :: reltol
        real(TKG), allocatable :: data(:), atad(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKG + getUnifRand(0._TKG, 1._TKG, lenData)")
                        data = 1._TKG + getUnifRand(0._TKG, 1._TKG, lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("atad = getFFTI(getFFTF(data))")
                        atad = getFFTI(getFFTF(data))
        call disp%show("atad")
        call disp%show( atad )
        call disp%show("reltol = sqrt(epsilon(1._TKG))")
                        reltol = sqrt(epsilon(1._TKG))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, atad(1:size(data)), reltol = reltol)")
        call disp%show( isClose(data, atad(1:size(data)), reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))")
                        call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKH
        real(TKG) :: reltol
        real(TKG), allocatable :: data(:), atad(:)
        call disp%show("lenData = getUnifRand(5, 100)")
                        lenData = getUnifRand(5, 100)
        call disp%show("lenData")
        call disp%show( lenData )
        call disp%show("data = 1._TKG + getUnifRand(0._TKG, 1._TKG, lenData)")
                        data = 1._TKG + getUnifRand(0._TKG, 1._TKG, lenData)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("atad = getFFTI(getFFTF(data))")
                        atad = getFFTI(getFFTF(data))
        call disp%show("atad")
        call disp%show( atad )
        call disp%show("reltol = sqrt(epsilon(1._TKG))")
                        reltol = sqrt(epsilon(1._TKG))
        call disp%show("reltol")
        call disp%show( reltol )
        call disp%show("isClose(data, atad(1:size(data)), reltol = reltol)")
        call disp%show( isClose(data, atad(1:size(data)), reltol = reltol) )
        call disp%show("call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))")
                        call setAsserted(all(isClose(data, atad(1:size(data)), reltol = reltol)))
        call disp%skip()
    end block
    end do

end program example