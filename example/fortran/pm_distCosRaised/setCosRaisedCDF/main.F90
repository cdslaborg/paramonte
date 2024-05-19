program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RK ! all real kinds are supported.
    use pm_arrayMembership, only: operator(.inrange.)
    use pm_mathSubAdd, only: operator(.subadd.)
    use pm_distCosRaised, only: setCosRaisedCDF
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 5_IK
    real(RKG), dimension(NP) :: mu, invSigma, CDF

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%show("! The standard CDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call setLinSpace(mu, x1 = -5._RKG, x2 = +5._RKG)
    call setLogSpace(invSigma, logx1 = log(0.1_RKG), logx2 = log(10._RKG))

    call disp%skip()
    call disp%show("call setCosRaisedCDF(CDF(1), 0._RKG)")
                    call setCosRaisedCDF(CDF(1), 0._RKG)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setCosRaisedCDF(CDF(1:5), [real(RKG) :: -1., -.5, 0., .5, 1.])")
                    call setCosRaisedCDF(CDF(1:5), [real(RKG) :: -1., -.5, 0., .5, 1.])
    call disp%show("CDF(1:5)")
    call disp%show( CDF(1:5) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("! CDF with a mean.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("call setCosRaisedCDF(CDF(1), mu(1), mu(1))")
                    call setCosRaisedCDF(CDF(1), mu(1), mu(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! CDF with a mean and a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("invSigma(1)")
    call disp%show( invSigma(1) )
    call disp%show("call setCosRaisedCDF(CDF(1:2), mu(1) .subadd. invSigma(1) / 2, mu(1), invSigma(1))")
                    call setCosRaisedCDF(CDF(1:2), mu(1) .subadd. invSigma(1) / 2, mu(1), invSigma(1))
    call disp%show("CDF(1:2)")
    call disp%show( CDF(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with the same mean and standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("invSigma(1)")
    call disp%show( invSigma(1) )
    call disp%show("call setCosRaisedCDF(CDF(1:NP), getLinSpace(mu(1) - invSigma(1), mu(1) + invSigma(1), count = NP), mu(1), invSigma(1))")
                    call setCosRaisedCDF(CDF(1:NP), getLinSpace(mu(1) - invSigma(1), mu(1) + invSigma(1), count = NP), mu(1), invSigma(1))
    call disp%show("CDF(1:NP)")
    call disp%show( CDF(1:NP) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at the same point but with different means and standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1:NP)")
    call disp%show( mu(1:NP) )
    call disp%show("invSigma(1:NP)")
    call disp%show( invSigma(1:NP) )
    call disp%show("call setCosRaisedCDF(CDF(1:NP), mu(1:NP), mu(1:NP), invSigma(1:NP))")
                    call setCosRaisedCDF(CDF(1:NP), mu(1:NP), mu(1:NP), invSigma(1:NP))
    call disp%show("CDF(1:NP)")
    call disp%show( CDF(1:NP) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with different means and a standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1:NP)")
    call disp%show( mu(1:NP) )
    call disp%show("invSigma(1:NP)")
    call disp%show( invSigma(1:NP) )
    call disp%show("call setCosRaisedCDF(CDF(1:NP), mu(1:NP), mu(1:NP), invSigma(1:NP))")
                    call setCosRaisedCDF(CDF(1:NP), mu(1:NP), mu(1:NP), invSigma(1:NP))
    call disp%show("CDF(1:NP)")
    call disp%show( CDF(1:NP) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i, j
        real(RKG) :: point(1000), CDF(4), mu(4), invSigma(4)
        open(newunit = fileUnit, file = "setCosRaisedCDF.RK.txt")
        call setLinSpace(point, x1 = -4._RKG, x2 = +4._RKG)
        invSigma = 1._RKG / [+3._RKG, +1._RKG, +.3_RKG, 1._RKG]
        mu = [+0._RKG, +0._RKG, +0._RKG, -2._RKG]
        do i = 1, size(point)
            do j = 1, size(CDF)
                if(point(i) .inrange. (mu(j) .subadd. 1._RKG / invSigma(j))) then
                    call setCosRaisedCDF(CDF(j), point(i), mu(j), invSigma(j))
                elseif (point(i) > mu(j) + 1._RKG / invSigma(j)) then
                    CDF(j) = 1._RKG
                else
                    CDF(j) = 0._RKG
                end if
            end do
            write(fileUnit,"(5(g0,:,' '))") point(i), CDF
        end do
        close(fileUnit)
    end block

end program example