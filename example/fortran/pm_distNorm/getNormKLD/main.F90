program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RK32, RK64, RK128 ! all real kinds are supported.
    use pm_distNorm, only: getNormKLD
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Kullback–Leibler Divergence (KLD) distance the Normal distribution Q from P.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getNormKLD(meanDiffSq = (0.5_RK32 - 1.5_RK32)**2) ! assuming standard deviations of 1.")
    call disp%show( getNormKLD(meanDiffSq = (0.5_RK32 - 1.5_RK32)**2) )
    call disp%show("getNormKLD(meanDiffSq = (0.5_RK64 - 1.5_RK64)**2)")
    call disp%show( getNormKLD(meanDiffSq = (0.5_RK64 - 1.5_RK64)**2) )
    call disp%show("getNormKLD(meanDiffSq = (0.5_RK128 - 1.5_RK128)**2)")
    call disp%show( getNormKLD(meanDiffSq = (0.5_RK128 - 1.5_RK128)**2) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormKLD(varP = 2._RK32 , varQ = 3._RK32) ! assuming similar means.")
    call disp%show( getNormKLD(varP = 2._RK32 , varQ = 3._RK32) )
    call disp%show("getNormKLD(varP = 2._RK64 , varQ = 3._RK64)")
    call disp%show( getNormKLD(varP = 2._RK64 , varQ = 3._RK64) )
    call disp%show("getNormKLD(varP = 2._RK128, varQ = 3._RK128)")
    call disp%show( getNormKLD(varP = 2._RK128, varQ = 3._RK128) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormKLD(varP = 3._RK32 , varQ = 2._RK32 ) ! assuming similar means.")
    call disp%show( getNormKLD(varP = 3._RK32 , varQ = 2._RK32 ) )
    call disp%show("getNormKLD(varP = 3._RK64 , varQ = 2._RK64 )")
    call disp%show( getNormKLD(varP = 3._RK64 , varQ = 2._RK64 ) )
    call disp%show("getNormKLD(varP = 3._RK128, varQ = 2._RK128)")
    call disp%show( getNormKLD(varP = 3._RK128, varQ = 2._RK128) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormKLD(meanDiffSq = 0.2_RK32 , varP = 3._RK32 , varQ = 2._RK32 ! assuming similar means.")
    call disp%show( getNormKLD(meanDiffSq = 0.2_RK32 , varP = 3._RK32 , varQ = 2._RK32 ) )
    call disp%show("getNormKLD(meanDiffSq = 0.2_RK64 , varP = 3._RK64 , varQ = 2._RK64 )")
    call disp%show( getNormKLD(meanDiffSq = 0.2_RK64 , varP = 3._RK64 , varQ = 2._RK64 ) )
    call disp%show("getNormKLD(meanDiffSq = 0.2_RK128, varP = 3._RK128, varQ = 2._RK128)")
    call disp%show( getNormKLD(meanDiffSq = 0.2_RK128, varP = 3._RK128, varQ = 2._RK128) )
    call disp%skip()

end program example