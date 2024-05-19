program fermi



    use pm_kind, only: SK, IK, LK, RKG => RKD
    use pm_val2str, only: getStr
    use pm_distBand, only: getBandEbreak
    use pm_distBand, only: setBandEnergy
    use pm_cosmology, only: setVolComDiffNormed
    use pm_cosmology, only: LOG_HUBBLE_VOLUME_CM3
    use pm_cosmology, only: LOG_HUBBLE_DISTANCE_CM
    use pm_cosmology, only: getHubbleParamNormedSq
    use pm_cosmology, only: getDisComTransNormedWU10
    use pm_cosmicRate, only: getLogRateDensityB10
    use pm_matrixCopy, only: lowDia, transSymm
    use pm_matrixCopy, only: rdpack, lfpack
    use pm_arraySpace, only: setLinSpace
    use pm_matrixCopy, only: setMatCopy
    use pm_sampleCov, only: getCov
    use pm_io, only: display_type
    use pm_io, only: getErrTableRead
    use pm_arrayResize, only: setResized
    use pm_io, only: getErrTableWrite, trans
    use pm_optimization, only: isFailedMinPowell



    implicit none



    type(display_type) :: disp
    character(:, SK), allocatable :: header
    integer(IK) :: ivar, ired, igrb, ngrb, ngrbz, unit, info
    integer(IK), parameter :: nred = 1000, nvar = 4
    real(RKG), parameter :: ZMIN = 0.1_RKG, ZMAX = 1.01e2_RKG
    real(RKG), parameter :: LOGZPLUS1_MIN = log(ZMIN + 1._RKG), LOGZPLUS1Z_MAX = log(ZMAX + 1._RKG)
    real(RKG), parameter :: DELTA_LOGZPLUS1 = (LOGZPLUS1Z_MAX - LOGZPLUS1_MIN) / (nred - 1)
    real(RKG), parameter :: TIME_DILATION_EXPO = 1._RKG ! kfacor correction. .66 is the other possibility.
    real(RKG), parameter :: LBNEW = .1_RKG, UBNEW = 20000._RKG, LB = 10._RKG, UB = 1000._RKG ! kev unit.
    real(RKG), parameter :: ALPHA = -1.1_RKG, BETA = -2.3_RKG
    real(RKG), allocatable :: grbz(:,:) ! (1:nvar, size(zindex) pbol, epk, sbol, t90
    real(RKG), allocatable :: logz(:) ! log of known z values.
    Logical :: findOptimalShiftEnabled = .false.
    Logical :: shiftEnabled = .true.
    Logical :: cleanSampleEnabled = .false.



    real(RKG), parameter :: PI = acos(-1._RKG)
    real(RKG), parameter :: LN10 = log(10._RKG)
    real(RKG), parameter :: LOG_FULL_SKY_SOLID_ANGLE = log(4._RKG * PI)



    type :: zinfo_type
        real(RKG) :: z, logSFR
        real(RKG) :: logInt2Obs(nvar) ! the redshift mapping: the subtraction of the log GRB attributes in rest frame from observer frame.
    end type
    type(zinfo_type) :: zinfo(nred)



    real(RKG), allocatable :: grb(:,:) ! pbol, epk, sbol, t90



    type :: model_type
        real(RKG) :: avg(nvar)
        real(RKG) :: invcov(nvar, nvar)
    end type
    type(model_type) :: model
    model%avg = [51.95, 2.59, 52.24, 1.] * LN10 ! pbol, epk, sbol, t90
    ! This matrix is computed for log_10 parameter values and must be converted to `log()` scale.
    model%invcov = reshape( [ +21.183484566160494_RKG, +7.095003864339406_RKG, -14.827122427564412_RKG, -0.310607307395483_RKG &
                            , +7.0950038643394060_RKG, +4.573466054775477_RKG, -5.4639334908475360_RKG, -0.193710749126834_RKG &
                            , -14.827122427564412_RKG, -5.463933490847536_RKG, +10.954404830891077_RKG, -0.415177537893988_RKG &
                            ,  -0.310607307395483_RKG, -0.193710749126834_RKG, -0.4151775378939880_RKG, +1.997104942888680_RKG &
                            ], shape = [nvar, nvar])



    ! Compute the model parameters.



    !block
    !    real(RKG) :: std(nvar), rho(nvar*(nvar + 1)/2), cor(nvar, nvar), cov(nvar, nvar)
    !    avg = [51.95, 2.59, 52.24, 1.] * LN10 ! pbol, epk, sbol, t90
    !    std = [0.56, 0.36, 0.87, 0.42] * LN10 ! pbol, epk, sbol, t90
    !    rho = [1., .42, .97, .66, 1., .58, .42, 1., .68, 1.]
    !    call setMatCopy(cor, rdpack, rho, lfpack, lowDia)
    !    call setMatCopy(cor(1:3,2:4), rdpack, cor(2:4,1:3), rdpack, lowDia, transSymm)
    !    cov = getCov(cor, std)
    !    !call disp%show(cor)
    !    !call disp%skip
    !    !call disp%show(cov)
    !end block



    ! Compute redshift properties.



    block
        real(RKG) :: zplus1, logzplus1, logLumDisCM, disComTransNormed
        open(newunit = unit, file = SK_"zgrid.txt", status = "replace")
        write(unit,"(*(g0))") "z,logLumDisCM,logSFR"
        do ired = 1, nred
            logzplus1 = LOGZPLUS1_MIN + real(ired - 1, RKG) * DELTA_LOGZPLUS1
            zplus1 = exp(logzplus1)
            zinfo(ired)%z = zplus1 - 1._RKG
            disComTransNormed = getDisComTransNormedWU10(zplus1)
            logLumDisCM = logzplus1 + log(disComTransNormed) + LOG_HUBBLE_DISTANCE_CM
            ! this is log(4*pi*dl^2) where dl is luminosity distance in units of mpc.
            zinfo(ired)%logInt2Obs(1) = LOG_FULL_SKY_SOLID_ANGLE + logLumDisCM * 2 ! log(liso/pbol)
            zinfo(ired)%logInt2Obs(2) = logzplus1 ! log(epkz/epko)
            zinfo(ired)%logInt2Obs(3) = zinfo(ired)%logInt2Obs(1) - logzplus1 ! log(eiso/sbol)
            zinfo(ired)%logInt2Obs(4) = -logzplus1 * TIME_DILATION_EXPO ! log(t90z/t90)
            call setVolComDiffNormed(zinfo(ired)%logSFR, disComTransNormedSq = disComTransNormed**2, hubbleParamNormed = sqrt(getHubbleParamNormedSq(zplus1)))
            zinfo(ired)%logSFR = getLogRateDensityB10(logzplus1) + log(zinfo(ired)%logSFR) + LOG_FULL_SKY_SOLID_ANGLE + LOG_HUBBLE_VOLUME_CM3
            ! compute redshift probability.
            write(unit,"(*(g0,:,','))") zinfo(ired)%z, logLumDisCM, zinfo(ired)%logSFR
        end do
    end block



    ! Read the GRB properties.


    if (cleanSampleEnabled) THEN
        print*,"using fermiClean.in"
         if (0 /= getErrTableRead(SK_"fermiClean.in", grb, trans, header = header)) error stop "Failed to read the fermi.txt file." ! pbol(10-1000 ph/cm2/s), epk(kev), sbol(10-1000 kev/cm2), t90(s), redshift (-1 implies missing info).
    else
        if (0 /= getErrTableRead(SK_"fermi.in", grb, trans, header = header)) error stop "Failed to read the fermi.txt file." ! pbol(10-1000 ph/cm2/s), epk(kev), sbol(10-1000 kev/cm2), t90(s), redshift (-1 implies missing info).
    end if

    ! Compute the GRB properties.

    ngrb = size(grb, 2, IK)
    !print*,grb(1:nvar,1)
    !print*,grb(1:nvar,2)
    if (size(grb, 1, IK) /= nvar + 1) error stop "The input GRB data table is unrecognized."
    block
        use pm_physUnit, only: KEV2ERGS
        real(RKG) :: ebreak, fluxph
        do igrb = 1, size(grb, 2, IK)
            ebreak = getBandEbreak(ALPHA, BETA, epeak = grb(2, igrb))
            ! Correct the flux energy window.
            fluxph = grb(1, igrb)
            call setBandEnergy(grb(1, igrb), LBNEW, UBNEW, fluxph, LB, UB, ALPHA, BETA, ebreak, info)
            if (info < 0) error stop "Failed to compute the bolometric peak flux."
            grb(1, igrb) = grb(1, igrb) * KEV2ERGS
            ! Correct the fluence energy window.
            call setBandEnergy(grb(3, igrb), LBNEW, UBNEW, LB, UB, ALPHA, BETA, ebreak, info)
            if (info < 0) error stop "Failed to compute the bolometric peak flux."
            grb(1 : 4, igrb) = log(grb(1 : 4, igrb))
        end do
        if (0 /= getErrTableWrite(SK_"fermi.out", grb, trans, header = header)) error stop "Failed to write the redshift prediction file."
    end block

    !shift the GRB attributes, comes after block because the shift was calculated using the calculated fluence above
    if(shiftEnabled) THEN
        do igrb=1,size(grb(1,:))
            grb(1:nvar,igrb) = grb(1:nvar,igrb) + [-0.118033988749895_RKG,3.61803398874989_RKG,1.50000000000000000_RKG,0.881966011250105_RKG]
        end do
        ![0.0000000000000000_RKG, 0.0000000000000000_RKG, 0.0000000000000000_RKG, -0.23606797749978969_RKG]!0.0000000000000000_RKG        0.0000000000000000_RKG        0.0000000000000000_RKG      -0.23606797749978969_RKG
        !clean shift -0.11803398874989479_RKG,1.5000000000000000_RKG,0.50000000000000000_RKG,0.50000000000000000_RKG
    end if

    ! Default predicted fermi redshifts with no data transformations.



    block
        integer(IK) :: maxLocPDF
        real(RKG) :: zsum(2, ngrb)  ! table of (zmode, zreal) for all events.
        real(RKG) :: zpdf(nred, 2)  ! table of (z, zpdf) for individual events.
        do igrb = 1, ngrb
            !call disp%show(SK_"processing fermi grb "//getStr(igrb))
            zpdf(:,1) = zinfo%z
            call setRedshiftPDF(zpdf(:,2), maxLocPDF, zinfo, grb(1 : nvar, igrb))
            zsum(1, igrb) = zinfo(maxLocPDF)%z ! predicted z dist mode.
            zsum(2, igrb) = grb(5, igrb) ! actual z.
            if (0 /= getErrTableWrite(SK_"fermi.z."//getStr(igrb)//SK_".out", zpdf, header = SK_"z,zpdf")) error stop "Failed to write the redshift prediction file."
        end do
        if (0 /= getErrTableWrite(SK_"fermi.z.out", zsum, trans, header = SK_"zmode,zreal")) error stop "Failed to write the redshift prediction file."
    end block



    ! Optimize the predicted fermi redshifts with linear data transformations.



    block
        use pm_arrayRange, only: getRange
        real(RKG) :: optimalShift(nvar) ! pbol, epk, sbol, t90 : corrections to add to GRB attributes to make predictions optimal.
        integer(IK), allocatable :: zindex(:) ! index of all z-known grbs in `grb(:,:)`.
        ! generate `grbz` sample containing only redshift-known GRBs.
        zindex = pack(getRange(1_IK, ngrb), nint(grb(5, :)) /= -1)
        logz = log(grb(nvar + 1, zindex)) ! all known zs.
        grbz = grb(1 : nvar, zindex) ! all z-known grbs.
        !print*,grbz(1:nvar,1)
        ngrbz = size(zindex, 1, IK)
        !disp = display_type(file = "temp.txt")
        optimalShift = [-0.118033988749895_RKG,3.61803398874989_RKG,1.50000000000000000_RKG,0.881966011250105_RKG]!0.5_RKG
        !print *,grbz(1:nvar, :)
        if(findOptimalShiftEnabled) Then
            if (isFailedMinPowell(getFunc = getErrSq, xmin = optimalShift)) error stop "Powell minimization failed." !**************************
            !print *, optimalShift
        end if
    end block



contains

    !function getScaling(grbScale) result(errSq)
    !    implicit none
    !    real(RKG), intent(in) :: grbScale(nvar) ! pbol, epk, sbol, t90 : scaling corrections to add to GRB attributes to make predictions optimal.
    !    real(RKG) :: errSq
    !    integer(IK) :: igrbz, maxLocPDF
    !    real(RKG) :: grbzScaled(nvar, ngrbz), zpdf(nred) ! table of zpdf for individual z-known events.
    !    do igrb =1,size(grbz(1,:))
    !        grbzScaled(1:nvar, igrb) = grbz(1:nvar,igrb) * grbScale
    !    end do
    !    errSq = 0._RKG
    !    do igrbz = 1, ngrbz
    !        call setRedshiftPDF(zpdf, maxLocPDF, zinfo, grbzScaled(1 : nvar, igrbz))
    !        errSq = errSq + (logz(igrbz) - log(zinfo(maxLocPDF)%z))**2
    !    end do
    !    !print*,"grbSscale:",grbScale,"error squared",errSq
    !end function

    function getErrSq(grbShift) result(errSq)
        use pm_sampleShift, only: getShifted
        implicit none
        real(RKG), intent(in) :: grbShift(nvar) ! pbol, epk, sbol, t90 : corrections to add to GRB attributes to make predictions optimal.
        real(RKG) :: errSq, grbzScaled(nvar,ngrbz),grbScale(nvar)
        integer(IK) :: igrbz, maxLocPDF
        real(RKG) :: grbzShifted(nvar, ngrbz), zpdf(nred) ! table of zpdf for individual z-known events.

    !    grbScale = 1.3_RKG
    !    if (isFailedMinPowell(getFunc = getScaling, xmin = grbScale)) error stop "Powell minimization failed." !**************************
    !    !print*,"getScaling:",grbScale
    !    do igrb =1,size(grbz(1,:))
    !        grbzScaled(1:nvar, igrb) = grbz(1:nvar,igrb) * grbScale
    !    end do

        !grbzShifted = getShifted(sample = grbzScaled, amount = grbShift, dim = 2_IK)
        grbzShifted = getShifted(sample = grbz, amount = grbShift, dim = 2_IK)
        errSq = 0._RKG

        do igrbz = 1, ngrbz
            call setRedshiftPDF(zpdf, maxLocPDF, zinfo, grbzShifted(1 : nvar, igrbz))
            errSq = errSq + (logz(igrbz) - log(zinfo(maxLocPDF)%z))**2
        end do
        !print*,"grbSscale:",grbScale,"grbShift:",grbShift,"error squared",errSq
        print*,"grbShift:",grbShift,"error squared",errSq

    end function

    !>  \brief
    !>  Generate and return the vector of (normalized) probabilities of an input `event` having the redshifts specified in the input `zinfo%z`.
    !>
    !>  \param[out] pdf         :   The output vector of type `real` of kind \RKD of size `nred`, containing the redshift grid probabilities for the input `event`.
    !>  \param[out] maxLocPDF   :   The output scalar of type `integer` of default kind \IK, the index of `zinfo` and `pdf` corresponding to the most likely redshift of the `event`.
    !>  \param[in]  zinfo       :   The input vector of size `nred` of type [zinfo_type](@ref zinfo_type) containing redshift grid information.
    !>  \param[in]  event       :   The input vector of size `nvar` of type [zinfo_type](@ref zinfo_type) containing redshift grid information.
    pure subroutine setRedshiftPDF(pdf, maxLocPDF, zinfo, event)
        type(zinfo_type), intent(in)    , contiguous    :: zinfo(:)
        real(RKG)       , intent(in)    , contiguous    :: event(:)
        real(RKG)       , intent(out)   , contiguous    :: pdf(:)
        integer(IK)     , intent(out)                   :: maxLocPDF
        real(RKG)                                       :: grbIntShifted(size(event)) ! shifted w.r.t. the mean of the model.
        integer(IK)                                     :: ired
        if (size(event) /= nvar) error stop
        if (size(zinfo) /= nred) error stop
        if (size(pdf) /= nred) error stop
        do ired = 1, size(zinfo)
            grbIntShifted = event + zinfo(ired)%logInt2Obs - model%avg
            ! compute the unnormalized probability at the current redshift.
            pdf(ired) = zinfo(ired)%logSFR - 0.5_RKG * dot_product(grbIntShifted, matmul(model%invcov, grbIntShifted))
        end do
        maxLocPDF = maxloc(pdf, 1) ! the index of maximum probability.
        pdf = exp(pdf - pdf(maxLocPDF)) ! normalize the UDF to PDF and compute the PDF mode.
        pdf = pdf / sum(pdf)
    end subroutine



end program fermi