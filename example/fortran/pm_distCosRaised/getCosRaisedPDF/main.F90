program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RK ! all real kinds are supported.
    use pm_arrayMembership, only: operator(.inrange.)
    use pm_mathSubAdd, only: operator(.subadd.)
    use pm_distCosRaised, only: getCosRaisedPDF
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 5_IK
    real(RKG), dimension(NP) :: mu, sigma, PDF

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%show("! The standard PDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call setLinSpace(mu, x1 = -5._RKG, x2 = +5._RKG)
    call setLogSpace(sigma, logx1 = log(0.1_RKG), logx2 = log(10._RKG))

    call disp%skip()
    call disp%show("PDF(1) = getCosRaisedPDF(0._RKG)")
                    PDF(1) = getCosRaisedPDF(0._RKG)
    call disp%show("PDF(1)")
    call disp%show( PDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("PDF(1:5) = getCosRaisedPDF([real(RKG) :: -1., -.5, 0., .5, 1.])")
                    PDF(1:5) = getCosRaisedPDF([real(RKG) :: -1., -.5, 0., .5, 1.])
    call disp%show("PDF(1:5)")
    call disp%show( PDF(1:5) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("! PDF with a mean.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("PDF(1) = getCosRaisedPDF(mu(1), mu(1))")
                    PDF(1) = getCosRaisedPDF(mu(1), mu(1))
    call disp%show("PDF(1)")
    call disp%show( PDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! PDF with a scale parameter.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("sigma(1)")
    call disp%show( sigma(1) )
    call disp%show("PDF(1) = getCosRaisedPDF(0._RKG, sigma = sigma(1))")
                    PDF(1) = getCosRaisedPDF(0._RKG, sigma = sigma(1))
    call disp%show("PDF(1)")
    call disp%show( PDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! PDF with a mean and a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("sigma(1)")
    call disp%show( sigma(1) )
    call disp%show("PDF(1:2) = getCosRaisedPDF(mu(1) .subadd. sigma(1) / 2, mu(1), sigma(1))")
                    PDF(1:2) = getCosRaisedPDF(mu(1) .subadd. sigma(1) / 2, mu(1), sigma(1))
    call disp%show("PDF(1:2)")
    call disp%show( PDF(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at different points with the same mean and standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("sigma(1)")
    call disp%show( sigma(1) )
    call disp%show("PDF(1:NP) = getCosRaisedPDF(getLinSpace(mu(1) - sigma(1), mu(1) + sigma(1), count = NP), mu(1), sigma(1))")
                    PDF(1:NP) = getCosRaisedPDF(getLinSpace(mu(1) - sigma(1), mu(1) + sigma(1), count = NP), mu(1), sigma(1))
    call disp%show("PDF(1:NP)")
    call disp%show( PDF(1:NP) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at the same point but with different means and standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1:NP)")
    call disp%show( mu(1:NP) )
    call disp%show("sigma(1:NP)")
    call disp%show( sigma(1:NP) )
    call disp%show("PDF(1:NP) = getCosRaisedPDF(mu(1:NP), mu(1:NP), sigma(1:NP))")
                    PDF(1:NP) = getCosRaisedPDF(mu(1:NP), mu(1:NP), sigma(1:NP))
    call disp%show("PDF(1:NP)")
    call disp%show( PDF(1:NP) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at different points with different means and a standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1:NP)")
    call disp%show( mu(1:NP) )
    call disp%show("sigma(1:NP)")
    call disp%show( sigma(1:NP) )
    call disp%show("PDF(1:NP) = getCosRaisedPDF(mu(1:NP), mu(1:NP), sigma(1:NP))")
                    PDF(1:NP) = getCosRaisedPDF(mu(1:NP), mu(1:NP), sigma(1:NP))
    call disp%show("PDF(1:NP)")
    call disp%show( PDF(1:NP) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i, j
        real(RKG) :: point(1000), PDF(4), mu(4), sigma(4)
        open(newunit = fileUnit, file = "getCosRaisedPDF.RK.txt")
        call setLinSpace(point, x1 = -4._RKG, x2 = +4._RKG)
        sigma = [+3._RKG, +1._RKG, +.3_RKG, 1._RKG]
        mu = [+0._RKG, +0._RKG, +0._RKG, -2._RKG]
        do i = 1, size(point)
            do j = 1, size(PDF)
                if(point(i) .inrange. (mu(j) .subadd. sigma(j))) then
                    PDF(j) = getCosRaisedPDF(point(i), mu(j), sigma(j))
                else
                    PDF(j) = 0._RKG
                end if
            end do
            write(fileUnit,"(5(g0,:,' '))") point(i), PDF
        end do
        close(fileUnit)
    end block

end program example