program example

    use pm_kind, only: SK, IK
    use pm_kind, only: modeli_type
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => IKS
        type(modeli_type) :: model
        integer(TKC) :: mold
        call disp%skip()
        call disp%show("kind(mold)")
        call disp%show( kind(mold) )
        call disp%show("model = modeli_type(mold)")
                        model = modeli_type(mold)
        call disp%show("[model%digits, model%radix, model%range, model%storage_size]")
        call disp%show( [model%digits, model%radix, model%range, model%storage_size] )
        call disp%show("model%huge")
        call disp%show( model%huge )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => IKD
        type(modeli_type) :: model
        integer(TKC) :: mold
        call disp%skip()
        call disp%show("kind(mold)")
        call disp%show( kind(mold) )
        call disp%show("model = modeli_type(mold)")
                        model = modeli_type(mold)
        call disp%show("[model%digits, model%radix, model%range, model%storage_size]")
        call disp%show( [model%digits, model%radix, model%range, model%storage_size] )
        call disp%show("model%huge")
        call disp%show( model%huge )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => IKH
        type(modeli_type) :: model
        integer(TKC) :: mold
        call disp%skip()
        call disp%show("kind(mold)")
        call disp%show( kind(mold) )
        call disp%show("model = modeli_type(mold)")
                        model = modeli_type(mold)
        call disp%show("[model%digits, model%radix, model%range, model%storage_size]")
        call disp%show( [model%digits, model%radix, model%range, model%storage_size] )
        call disp%show("model%huge")
        call disp%show( model%huge )
        call disp%skip()
    end block

end program example