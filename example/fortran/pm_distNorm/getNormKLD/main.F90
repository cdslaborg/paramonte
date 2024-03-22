program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RKS, RKD, RKH ! all real kinds are supported.
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
    call disp%show("getNormKLD(meanDiffSq = (0.5_RKS - 1.5_RKS)**2) ! assuming standard deviations of 1.")
    call disp%show( getNormKLD(meanDiffSq = (0.5_RKS - 1.5_RKS)**2) )
    call disp%show("getNormKLD(meanDiffSq = (0.5_RKD - 1.5_RKD)**2)")
    call disp%show( getNormKLD(meanDiffSq = (0.5_RKD - 1.5_RKD)**2) )
    call disp%show("getNormKLD(meanDiffSq = (0.5_RKH - 1.5_RKH)**2)")
    call disp%show( getNormKLD(meanDiffSq = (0.5_RKH - 1.5_RKH)**2) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormKLD(varP = 2._RKS , varQ = 3._RKS) ! assuming similar means.")
    call disp%show( getNormKLD(varP = 2._RKS , varQ = 3._RKS) )
    call disp%show("getNormKLD(varP = 2._RKD , varQ = 3._RKD)")
    call disp%show( getNormKLD(varP = 2._RKD , varQ = 3._RKD) )
    call disp%show("getNormKLD(varP = 2._RKH, varQ = 3._RKH)")
    call disp%show( getNormKLD(varP = 2._RKH, varQ = 3._RKH) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormKLD(varP = 3._RKS , varQ = 2._RKS ) ! assuming similar means.")
    call disp%show( getNormKLD(varP = 3._RKS , varQ = 2._RKS ) )
    call disp%show("getNormKLD(varP = 3._RKD , varQ = 2._RKD )")
    call disp%show( getNormKLD(varP = 3._RKD , varQ = 2._RKD ) )
    call disp%show("getNormKLD(varP = 3._RKH, varQ = 2._RKH)")
    call disp%show( getNormKLD(varP = 3._RKH, varQ = 2._RKH) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormKLD(meanDiffSq = 0.2_RKS , varP = 3._RKS , varQ = 2._RKS ! assuming similar means.")
    call disp%show( getNormKLD(meanDiffSq = 0.2_RKS , varP = 3._RKS , varQ = 2._RKS ) )
    call disp%show("getNormKLD(meanDiffSq = 0.2_RKD , varP = 3._RKD , varQ = 2._RKD )")
    call disp%show( getNormKLD(meanDiffSq = 0.2_RKD , varP = 3._RKD , varQ = 2._RKD ) )
    call disp%show("getNormKLD(meanDiffSq = 0.2_RKH, varP = 3._RKH, varQ = 2._RKH)")
    call disp%show( getNormKLD(meanDiffSq = 0.2_RKH, varP = 3._RKH, varQ = 2._RKH) )
    call disp%skip()

end program example