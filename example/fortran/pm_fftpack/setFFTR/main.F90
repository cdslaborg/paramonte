program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_fftpack, only: getfactorFFT
    use pm_fftpack, only: setFFTF, setFFTR
    use pm_arrayResize, only: setResized

    implicit none

    logical(LK) :: inwork
    integer(IK), allocatable :: factor(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: CKG => CKS
        complex(CKG), allocatable :: data(:), coef(:), work(:)
        call disp%skip()
        call disp%show("data = [complex(CKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.)]")
                        data = [complex(CKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.)]
        call disp%show("if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)")
                        if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)
        call disp%show("if (allocated(work)) deallocate(work); allocate(work, mold = data)")
                        if (allocated(work)) deallocate(work); allocate(work, mold = data)
        call disp%show("factor = getfactorFFT(data, coef)")
                        factor = getfactorFFT(data, coef)
        call disp%show("coef")
        call disp%show( coef )
        call disp%show("factor")
        call disp%show( factor )
        call disp%show("call setFFTF(factor, coef, data, work, inwork)")
                        call setFFTF(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(factor, coef, data, work, inwork)")
                        call setFFTR(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data / size(data)")
        call disp%show( data / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: CKG => CKD
        complex(CKG), allocatable :: data(:), coef(:), work(:)
        call disp%skip()
        call disp%show("data = [complex(CKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.)]")
                        data = [complex(CKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.)]
        call disp%show("if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)")
                        if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)
        call disp%show("if (allocated(work)) deallocate(work); allocate(work, mold = data)")
                        if (allocated(work)) deallocate(work); allocate(work, mold = data)
        call disp%show("factor = getfactorFFT(data, coef)")
                        factor = getfactorFFT(data, coef)
        call disp%show("coef")
        call disp%show( coef )
        call disp%show("factor")
        call disp%show( factor )
        call disp%show("call setFFTF(factor, coef, data, work, inwork)")
                        call setFFTF(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(factor, coef, data, work, inwork)")
                        call setFFTR(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data / size(data)")
        call disp%show( data / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: CKG => CKH
        complex(CKG), allocatable :: data(:), coef(:), work(:)
        call disp%skip()
        call disp%show("data = [complex(CKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.)]")
                        data = [complex(CKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.)]
        call disp%show("if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)")
                        if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)
        call disp%show("if (allocated(work)) deallocate(work); allocate(work, mold = data)")
                        if (allocated(work)) deallocate(work); allocate(work, mold = data)
        call disp%show("factor = getfactorFFT(data, coef)")
                        factor = getfactorFFT(data, coef)
        call disp%show("coef")
        call disp%show( coef )
        call disp%show("factor")
        call disp%show( factor )
        call disp%show("call setFFTF(factor, coef, data, work, inwork)")
                        call setFFTF(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(factor, coef, data, work, inwork)")
                        call setFFTR(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data / size(data)")
        call disp%show( data / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: RKG => RKS
        real(RKG), allocatable :: data(:), coef(:), work(:)
        call disp%skip()
        call disp%show("data = [real(RKG) :: 1., 2., 3., 4., 5., 6.5]")
                        data = [real(RKG) :: 1., 2., 3., 4., 5., 6.5]
        call disp%show("if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)")
                        if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)
        call disp%show("if (allocated(work)) deallocate(work); allocate(work, mold = data)")
                        if (allocated(work)) deallocate(work); allocate(work, mold = data)
        call disp%show("factor = getfactorFFT(data, coef)")
                        factor = getfactorFFT(data, coef)
        call disp%show("coef")
        call disp%show( coef )
        call disp%show("factor")
        call disp%show( factor )
        call disp%show("call setFFTF(factor, coef, data, work, inwork)")
                        call setFFTF(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(factor, coef, data, work, inwork)")
                        call setFFTR(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data / size(data)")
        call disp%show( data / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: RKG => RKD
        real(RKG), allocatable :: data(:), coef(:), work(:)
        call disp%skip()
        call disp%show("data = [real(RKG) :: 1., 2., 3., 4., 5., 6.5]")
                        data = [real(RKG) :: 1., 2., 3., 4., 5., 6.5]
        call disp%show("if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)")
                        if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)
        call disp%show("if (allocated(work)) deallocate(work); allocate(work, mold = data)")
                        if (allocated(work)) deallocate(work); allocate(work, mold = data)
        call disp%show("factor = getfactorFFT(data, coef)")
                        factor = getfactorFFT(data, coef)
        call disp%show("coef")
        call disp%show( coef )
        call disp%show("factor")
        call disp%show( factor )
        call disp%show("call setFFTF(factor, coef, data, work, inwork)")
                        call setFFTF(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(factor, coef, data, work, inwork)")
                        call setFFTR(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data / size(data)")
        call disp%show( data / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: RKG => RKH
        real(RKG), allocatable :: data(:), coef(:), work(:)
        call disp%skip()
        call disp%show("data = [real(RKG) :: 1., 2., 3., 4., 5., 6.5]")
                        data = [real(RKG) :: 1., 2., 3., 4., 5., 6.5]
        call disp%show("if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)")
                        if (allocated(coef)) deallocate(coef); allocate(coef, mold = data)
        call disp%show("if (allocated(work)) deallocate(work); allocate(work, mold = data)")
                        if (allocated(work)) deallocate(work); allocate(work, mold = data)
        call disp%show("factor = getfactorFFT(data, coef)")
                        factor = getfactorFFT(data, coef)
        call disp%show("coef")
        call disp%show( coef )
        call disp%show("factor")
        call disp%show( factor )
        call disp%show("call setFFTF(factor, coef, data, work, inwork)")
                        call setFFTF(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(factor, coef, data, work, inwork)")
                        call setFFTR(factor, coef, data, work, inwork)
        call disp%show("if (inwork) data = work")
                        if (inwork) data = work
        call disp%show("data / size(data)")
        call disp%show( data / size(data) )
        call disp%skip()
    end block

end program example