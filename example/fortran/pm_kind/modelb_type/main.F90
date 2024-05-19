program example

    use pm_kind, only: SK, IK
    use pm_kind, only: modelb_type
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKG => IKS
        type(modelb_type) :: model
        integer(TKG) :: mold
        call disp%skip()
        call disp%show("kind(mold)")
        call disp%show( kind(mold) )
        call disp%show("model = modelb_type(mold)")
                        model = modelb_type(mold)
        call disp%show("[model%bit_size, model%digits, model%radix, model%range, model%storage_size]")
        call disp%show( [model%bit_size, model%digits, model%radix, model%range, model%storage_size] )
        call disp%show("model%huge")
        call disp%show( model%huge )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => IKD
        type(modelb_type) :: model
        integer(TKG) :: mold
        call disp%skip()
        call disp%show("kind(mold)")
        call disp%show( kind(mold) )
        call disp%show("model = modelb_type(mold)")
                        model = modelb_type(mold)
        call disp%show("[model%bit_size, model%digits, model%radix, model%range, model%storage_size]")
        call disp%show( [model%bit_size, model%digits, model%radix, model%range, model%storage_size] )
        call disp%show("model%huge")
        call disp%show( model%huge )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => IKH
        type(modelb_type) :: model
        integer(TKG) :: mold
        call disp%skip()
        call disp%show("kind(mold)")
        call disp%show( kind(mold) )
        call disp%show("model = modelb_type(mold)")
                        model = modelb_type(mold)
        call disp%show("[model%bit_size, model%digits, model%radix, model%range, model%storage_size]")
        call disp%show( [model%bit_size, model%digits, model%radix, model%range, model%storage_size] )
        call disp%show("model%huge")
        call disp%show( model%huge )
        call disp%skip()
    end block

end program example