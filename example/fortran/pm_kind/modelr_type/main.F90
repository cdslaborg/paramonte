program example

    use pm_kind, only: SK, IK
    use pm_kind, only: modelr_type
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => RKS
        type(modelr_type) :: model
        real(TKC) :: mold
        call disp%skip()
        call disp%show("kind(mold)")
        call disp%show( kind(mold) )
        call disp%show("model = modelr_type(mold)")
                        model = modelr_type(mold)
        call disp%show("[model%digits, model%maxexponent, model%minexponent, model%precision, model%radix, model%range, model%storage_size]")
        call disp%show( [model%digits, model%maxexponent, model%minexponent, model%precision, model%radix, model%range, model%storage_size] )
        call disp%show("[model%epsilon, model%huge, model%tiny]")
        call disp%show( [model%epsilon, model%huge, model%tiny] )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKD
        type(modelr_type) :: model
        real(TKC) :: mold
        call disp%skip()
        call disp%show("kind(mold)")
        call disp%show( kind(mold) )
        call disp%show("model = modelr_type(mold)")
                        model = modelr_type(mold)
        call disp%show("[model%digits, model%maxexponent, model%minexponent, model%precision, model%radix, model%range, model%storage_size]")
        call disp%show( [model%digits, model%maxexponent, model%minexponent, model%precision, model%radix, model%range, model%storage_size] )
        call disp%show("[model%epsilon, model%huge, model%tiny]")
        call disp%show( [model%epsilon, model%huge, model%tiny] )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKH
        type(modelr_type) :: model
        real(TKC) :: mold
        call disp%skip()
        call disp%show("kind(mold)")
        call disp%show( kind(mold) )
        call disp%show("model = modelr_type(mold)")
                        model = modelr_type(mold)
        call disp%show("[model%digits, model%maxexponent, model%minexponent, model%precision, model%radix, model%range, model%storage_size]")
        call disp%show( [model%digits, model%maxexponent, model%minexponent, model%precision, model%radix, model%range, model%storage_size] )
        call disp%show("[model%epsilon, model%huge, model%tiny]")
        call disp%show( [model%epsilon, model%huge, model%tiny] )
        call disp%skip()
    end block

end program example