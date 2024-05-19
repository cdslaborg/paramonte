program example

    use pm_kind, only: SK, IK, LK, RKG => RKD
    use pm_parallelism, only: setForkJoinScaling
    use pm_arraySpace, only: getLinSpace
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    real(RKG) :: redshift = 5.5_RKG
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Fisher z transformation of a random bounded variable.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKS
        integer(IK), allocatable :: numproc(:)
        real, allocatable :: scaling(:)
        integer :: maxloc
        real :: maxval
        call disp%skip()
        call disp%show("call setForkJoinScaling(conProb = .23, seqSecTime = 0., parSecTime = 1., comSecTime = .001, scaling = scaling, numproc = numproc, scalingMaxVal = maxval, scalingMaxLoc = maxloc)")
                        call setForkJoinScaling(conProb = .23, seqSecTime = 0., parSecTime = 1., comSecTime = .001, scaling = scaling, numproc = numproc, scalingMaxVal = maxval, scalingMaxLoc = maxloc)
        call disp%show("scaling")
        call disp%show( scaling )
        call disp%show("numproc")
        call disp%show( numproc )
        call disp%show("maxloc")
        call disp%show( maxloc )
        call disp%show("maxval")
        call disp%show( maxval )
        call disp%skip()
    end block

    ! Generate both the cosmic rate and the rate density.

    block
        use pm_io, only: getErrTableWrite
        real :: maxval
        integer :: maxloc
        real, allocatable :: scaling(:)
        integer(IK), allocatable :: numproc(:)
        call setForkJoinScaling(conProb = .23, seqSecTime = 0., parSecTime = 1., comSecTime = .001, scaling = scaling, numproc = numproc, scalingMaxVal = maxval, scalingMaxLoc = maxloc, scalingMinLen = 32_IK)
        if (0 /= getErrTableWrite("setForkJoinScaling.csv", reshape([real(numproc), scaling], [size(numproc), 2]), header = "Processor Count,Scaling")) error stop "Table writing failed."
    end block

end program example