program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_container, only: operator(<)
    use pm_container, only: css_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block

        use pm_kind, only: SKC => SK
        use pm_container, only: css_type
        type(css_type) :: con1, con2

        call disp%skip()
        call disp%show("con1 = css_type(SKC_'PARAMONTE')")
                        con1 = css_type(SKC_'PARAMONTE')
        call disp%show("con2 = css_type(SKC_'PARAMONTE')")
                        con2 = css_type(SKC_'PARAMONTE')
        call disp%show("[con1, con2]")
        call disp%show( [con1, con2] , deliml = SK_"""" )
        call disp%show("con1 < con2")
        call disp%show( con1 < con2 )
        call disp%skip()

        call disp%skip()
        call disp%show("con1 = css_type(SKC_'PARAMONTE')")
                        con1 = css_type(SKC_'PARAMONTE')
        call disp%show("con2 = css_type(SKC_'paramonte')")
                        con2 = css_type(SKC_'paramonte')
        call disp%show("[con1, con2]")
        call disp%show( [con1, con2] , deliml = SK_"""" )
        call disp%show("con1 < con2")
        call disp%show( con1 < con2 )
        call disp%skip()

        call disp%show("con1 = css_type(SKC_'paramonte')")
                        con1 = css_type(SKC_'paramonte')
        call disp%show("con2 = css_type(SKC_'PARAMONTE')")
                        con2 = css_type(SKC_'PARAMONTE')
        call disp%skip()
        call disp%show("[con1, con2]")
        call disp%show( [con1, con2] , deliml = SK_"""" )
        call disp%show("con1 < con2")
        call disp%show( con1 < con2 )
        call disp%skip()

        call disp%skip()
        call disp%show("[css_type(SKC_'a'), css_type(SKC_'b'), css_type(SKC_'c')] < css_type(SKC_'b')")
        call disp%show( [css_type(SKC_'a'), css_type(SKC_'b'), css_type(SKC_'c')] < css_type(SKC_'b') )
        call disp%skip()

        call disp%skip()
        call disp%show("[css_type(SKC_'a'), css_type(SKC_'b'), css_type(SKC_'c')] < [css_type(SKC_'c'), css_type(SKC_'b'), css_type(SKC_'a')]")
        call disp%show( [css_type(SKC_'a'), css_type(SKC_'b'), css_type(SKC_'c')] < [css_type(SKC_'c'), css_type(SKC_'b'), css_type(SKC_'a')] )
        call disp%skip()

    end block

#if PDT_ENABLED
    block

        use pm_kind, only: SKC => SK
        use pm_container, only: css_pdt
        type(css_pdt) :: con1, con2

        call disp%skip()
        call disp%show("con1 = css_pdt(SKC_'PARAMONTE')")
                        con1 = css_pdt(SKC_'PARAMONTE')
        call disp%show("con2 = css_pdt(SKC_'PARAMONTE')")
                        con2 = css_pdt(SKC_'PARAMONTE')
        call disp%show("[con1, con2]")
        call disp%show( [con1, con2] , deliml = SK_"""" )
        call disp%show("con1 < con2")
        call disp%show( con1 < con2 )
        call disp%skip()

        call disp%skip()
        call disp%show("con1 = css_pdt(SKC_'PARAMONTE')")
                        con1 = css_pdt(SKC_'PARAMONTE')
        call disp%show("con2 = css_pdt(SKC_'paramonte')")
                        con2 = css_pdt(SKC_'paramonte')
        call disp%show("[con1, con2]")
        call disp%show( [con1, con2] , deliml = SK_"""" )
        call disp%show("con1 < con2")
        call disp%show( con1 < con2 )
        call disp%skip()

        call disp%show("con1 = css_pdt(SKC_'paramonte')")
                        con1 = css_pdt(SKC_'paramonte')
        call disp%show("con2 = css_pdt(SKC_'PARAMONTE')")
                        con2 = css_pdt(SKC_'PARAMONTE')
        call disp%skip()
        call disp%show("[con1, con2]")
        call disp%show( [con1, con2] , deliml = SK_"""" )
        call disp%show("con1 < con2")
        call disp%show( con1 < con2 )
        call disp%skip()

        call disp%skip()
        call disp%show("[css_pdt(SKC_'a'), css_pdt(SKC_'b'), css_pdt(SKC_'c')] < css_pdt(SKC_'b')")
        call disp%show( [css_pdt(SKC_'a'), css_pdt(SKC_'b'), css_pdt(SKC_'c')] < css_pdt(SKC_'b') )
        call disp%skip()

        call disp%skip()
        call disp%show("[css_pdt(SKC_'a'), css_pdt(SKC_'b'), css_pdt(SKC_'c')] < [css_pdt(SKC_'c'), css_pdt(SKC_'b'), css_pdt(SKC_'a')]")
        call disp%show( [css_pdt(SKC_'a'), css_pdt(SKC_'b'), css_pdt(SKC_'c')] < [css_pdt(SKC_'c'), css_pdt(SKC_'b'), css_pdt(SKC_'a')] )
        call disp%skip()

    end block
#endif

end program example