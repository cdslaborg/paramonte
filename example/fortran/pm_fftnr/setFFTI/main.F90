program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_fftnr, only: setFFTF, setFFTI
    use pm_arrayResize, only: setResized

    implicit none

    logical(LK) :: inwork
    integer(IK), allocatable :: factor(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKG => CKS
        complex(TKG), allocatable :: data(:)
        complex(TKG), parameter :: ZERO = 0._TKG
        call disp%skip()
        call disp%show("data = [complex(TKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]")
                        data = [complex(TKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTI(data)")
                        call setFFTI(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => CKD
        complex(TKG), allocatable :: data(:)
        complex(TKG), parameter :: ZERO = 0._TKG
        call disp%skip()
        call disp%show("data = [complex(TKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]")
                        data = [complex(TKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTI(data)")
                        call setFFTI(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => CKH
        complex(TKG), allocatable :: data(:)
        complex(TKG), parameter :: ZERO = 0._TKG
        call disp%skip()
        call disp%show("data = [complex(TKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]")
                        data = [complex(TKG) :: (1., -6.), (2., -5.), (3., -4.), (4., -3.), (5., -2.), (6., -1.), ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTI(data)")
                        call setFFTI(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: data(:)
        real(TKG), parameter :: ZERO = 0._TKG
        call disp%skip()
        call disp%show("data = [real(TKG) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]")
                        data = [real(TKG) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTI(data)")
                        call setFFTI(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKD
        real(TKG), allocatable :: data(:)
        real(TKG), parameter :: ZERO = 0._TKG
        call disp%skip()
        call disp%show("data = [real(TKG) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]")
                        data = [real(TKG) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTI(data)")
                        call setFFTI(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKH
        real(TKG), allocatable :: data(:)
        real(TKG), parameter :: ZERO = 0._TKG
        call disp%skip()
        call disp%show("data = [real(TKG) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]")
                        data = [real(TKG) :: 1., 2., 3., 4., 5., 6.5, ZERO, ZERO]
        call disp%show("call setFFTF(data)")
                        call setFFTF(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
        call disp%show("call setFFTI(data)")
                        call setFFTI(data)
        call disp%show("data")
        call disp%show( data )
        call disp%skip()
    end block

end program example