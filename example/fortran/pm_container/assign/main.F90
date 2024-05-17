program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_container, only: assignment(=)

    implicit none

    block
        use pm_kind, only: SKC => SK
        use pm_container, only: css_type
        type(css_type) :: destin(3)

        type(display_type) :: disp
        disp = display_type(file = "main.out.F90")

        call disp%skip()
        call disp%show("destin(1) = css_type(SKC_'ParaMonte')")
                        destin(1) = css_type(SKC_'ParaMonte')
        call disp%show("destin(1)")
        call disp%show( destin(1) , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("destin = css_type(SKC_'ParaMonte')")
                        destin = css_type(SKC_'ParaMonte')
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("destin = [css_type(SKC_'a'), css_type(SKC_'b'), css_type(SKC_'c')]")
                        destin = [css_type(SKC_'a'), css_type(SKC_'b'), css_type(SKC_'c')]
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()
    end block

#if PDT_ENABLED
    block
        use pm_kind, only: SKC => SK
        use pm_container, only: css_pdt
        type(css_pdt(SKC)) :: destin(3)

        type(display_type) :: disp
        disp = display_type(file = "main.out.F90")

        call disp%skip()
        call disp%show("destin(1) = css_pdt(SKC_'ParaMonte')")
                        destin(1) = css_pdt(SKC_'ParaMonte')
        call disp%show("destin(1)")
        call disp%show( destin(1) , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("destin = css_pdt(SKC_'ParaMonte')")
                        destin = css_pdt(SKC_'ParaMonte')
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("destin = [css_pdt(SKC_'a'), css_pdt(SKC_'b'), css_pdt(SKC_'c')]")
                        destin = [css_pdt(SKC_'a'), css_pdt(SKC_'b'), css_pdt(SKC_'c')]
        call disp%show("destin")
        call disp%show( destin , deliml = SK_"""" )
        call disp%skip()
    end block
#endif

end program example