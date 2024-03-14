program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_fftnr, only: setFFTF, setFFTR
    use pm_arrayResize, only: setResized

    implicit none

    logical(LK) :: inwork
    integer(IK), allocatable :: factor(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => CK32
        complex(TKC), allocatable :: data(:)
        complex(TKC), parameter :: ZERO = 0._TKC
        call disp%skip()
        call disp%show("data = [complex(TKC) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]")
                        data = [complex(TKC) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(data)")
                        call setFFTR(data)
        call disp%show("data * 2 / size(data)")
        call disp%show( data * 2 / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => CK64
        complex(TKC), allocatable :: data(:)
        complex(TKC), parameter :: ZERO = 0._TKC
        call disp%skip()
        call disp%show("data = [complex(TKC) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]")
                        data = [complex(TKC) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(data)")
                        call setFFTR(data)
        call disp%show("data * 2 / size(data)")
        call disp%show( data * 2 / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => CKH
        complex(TKC), allocatable :: data(:)
        complex(TKC), parameter :: ZERO = 0._TKC
        call disp%skip()
        call disp%show("data = [complex(TKC) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]")
                        data = [complex(TKC) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(data)")
                        call setFFTR(data)
        call disp%show("data * 2 / size(data)")
        call disp%show( data * 2 / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RK32
        real(TKC), allocatable :: data(:)
        real(TKC), parameter :: ZERO = 0._TKC
        call disp%skip()
        call disp%show("data = [real(TKC) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]")
                        data = [real(TKC) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(data)")
                        call setFFTR(data)
        call disp%show("data * 2 / size(data)")
        call disp%show( data * 2 / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RK64
        real(TKC), allocatable :: data(:)
        real(TKC), parameter :: ZERO = 0._TKC
        call disp%skip()
        call disp%show("data = [real(TKC) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]")
                        data = [real(TKC) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(data)")
                        call setFFTR(data)
        call disp%show("data * 2 / size(data)")
        call disp%show( data * 2 / size(data) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKH
        real(TKC), allocatable :: data(:)
        real(TKC), parameter :: ZERO = 0._TKC
        call disp%skip()
        call disp%show("data = [real(TKC) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]")
                        data = [real(TKC) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTR(data)")
                        call setFFTR(data)
        call disp%show("data * 2 / size(data)")
        call disp%show( data * 2 / size(data) )
        call disp%skip()
    end block

end program example