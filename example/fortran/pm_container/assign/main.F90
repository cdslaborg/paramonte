program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_container, only: assignment(=)

    implicit none

    block
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type) :: destin(3)

        type(display_type) :: disp
        disp = display_type(file = "main.out.F90")

        call disp%skip()
        call disp%show("destin(1) = css_type(SKG_'ParaMonte')")
                        destin(1) = css_type(SKG_'ParaMonte')
        call disp%show("destin(1)")
        call disp%show( destin(1) , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("destin = css_type(SKG_'ParaMonte')")
                        destin = css_type(SKG_'ParaMonte')
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("destin = [css_type(SKG_'a'), css_type(SKG_'b'), css_type(SKG_'c')]")
                        destin = [css_type(SKG_'a'), css_type(SKG_'b'), css_type(SKG_'c')]
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()
    end block

#if PDT_ENABLED
    block
        use pm_kind, only: SKG => SK
        use pm_container, only: css_pdt
        type(css_pdt(SKG)) :: destin(3)

        type(display_type) :: disp
        disp = display_type(file = "main.out.F90")

        call disp%skip()
        call disp%show("destin(1) = css_pdt(SKG_'ParaMonte')")
                        destin(1) = css_pdt(SKG_'ParaMonte')
        call disp%show("destin(1)")
        call disp%show( destin(1) , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("destin = css_pdt(SKG_'ParaMonte')")
                        destin = css_pdt(SKG_'ParaMonte')
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("destin = [css_pdt(SKG_'a'), css_pdt(SKG_'b'), css_pdt(SKG_'c')]")
                        destin = [css_pdt(SKG_'a'), css_pdt(SKG_'b'), css_pdt(SKG_'c')]
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()
    end block
#endif

end program example