!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of the ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module Batse_mod

    use Constants_mod, only: IK, RK, PI, LN10

    implicit none

    character(*), parameter :: MODULE_NAME = "@Batse_mod"

    integer(IK) , parameter :: NLGRB = 1366, NSGRB = 565, NVAR = 4
    real(RK)    , parameter :: MIN_LOG10PH53_4_LOGPBOLZERO  = 4.92_RK
    real(RK)    , parameter :: MAX_LOG10PH53_4_LOGPBOLZERO  = 6.318167895318538_RK
    real(RK)    , parameter :: MIN_LOGPH53_4_LOGPBOLZERO    = MIN_LOG10PH53_4_LOGPBOLZERO * LN10
    real(RK)    , parameter :: MAX_LOGPH53_4_LOGPBOLZERO    = MAX_LOG10PH53_4_LOGPBOLZERO * LN10

    ! Parameters for the effective peak flux of SGRBs (P64ms conversion to P1024):

    ! The scale of the change in BATSE efficiency for different GRB durations.
    real(RK)    , parameter :: ERFC_HEIGHT                  = LN10 * (+0.282313526464596_RK)

    ! Mean duration in the Error function used to model the connection between the peak fluxes in 64 and 1024 ms.
    real(RK)    , parameter :: ERFC_AVG                     = LN10 * (-0.483553339256463_RK)

    ! scale of the duration in the Error function used to model the connection between the peak fluxes in 64 and 1024 ms.
    real(RK)    , parameter :: ERFC_STD                     = LN10 * (+1.051469898469480_RK)

    ! GRB attributes
    type :: Event_type
        real(RK) :: logPbol, logEpk, logSbol, logT90, logPF53
    end type Event_type

    type :: GRB_type
        integer(IK) :: count
        type(Event_type), allocatable :: Event(:)
        !type(Event_type) :: Event(NLGRB)
    end type GRB_type

#if defined CAF_ENABLED
    type(GRB_type) :: GRB[*]
#else
    type(GRB_type) :: GRB
#endif

    real(RK), allocatable :: Trigger(:)
    !integer(IK) :: Trigger(NLGRB)
    !integer(IK) :: TriggerSGRB(NSGRB)

!***********************************************************************************************************************************
!***********************************************************************************************************************************
        
contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! returns all log of data in natural (Neper) base.
    subroutine readDataGRB(inFilePath,outFilePath,isLgrb)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: readDataGRB
#endif

        use Constants_mod, only: IK, RK
        implicit none
        character(*), intent(in)            :: inFilePath, outFilePath
        integer(IK)                         :: inFileUnit, outFileUnit, igrb
        logical     , intent(in)            :: isLgrb

        if (isLgrb) then
            GRB%count = NLGRB
        else
            GRB%count = NSGRB
        end if

        allocate(GRB%Event(GRB%count), Trigger(GRB%count))

#if defined CAF_ENABLED
        if (this_image()==1) then
#endif
            open(newunit=inFileUnit,file=inFilePath,status="old")
            open(newunit=outFileUnit,file=outFilePath,status="replace")
            write(outFileUnit,"(9a30)"  ) "trigger"             &
                                        , "logPbol_1eV_20MeV"   &
                                        , "logSbol_1eV_20MeV"   &
                                        , "logEpk"              &
                                        , "logEPR1024"          &
                                        , "logEFR"              &
                                        , "logFPR1024"          &
                                        , "logT90"              &
                                        , "logEffPF53"

            ! skip the header row in the input file
            read(inFileUnit,*)

            ! read BATSE GRB data
            do igrb = 1, GRB%count

                if (isLgrb) then
                    read(inFileUnit,*   ) Trigger(igrb)             &
                                        , GRB%Event(igrb)%logPF53   &
                                        , GRB%Event(igrb)%logEpk    &
                                        , GRB%Event(igrb)%logSbol   &
                                        , GRB%Event(igrb)%logT90
                else
                    read(inFileUnit,*   ) Trigger(igrb)             &
                                        , GRB%Event(igrb)%logPbol   &
                                        , GRB%Event(igrb)%logSbol   &
                                        , GRB%Event(igrb)%logEpk    &
                                        , GRB%Event(igrb)%logT90    &
                                        , GRB%Event(igrb)%logPF53
                end if

                ! convert all values to logarithm in base Neper
                GRB%Event(igrb)%logPF53 = LN10 * GRB%Event(igrb)%logPF53
                GRB%Event(igrb)%logEpk  = LN10 * GRB%Event(igrb)%logEpk
                GRB%Event(igrb)%logSbol = LN10 * GRB%Event(igrb)%logSbol
                GRB%Event(igrb)%logT90  = LN10 * GRB%Event(igrb)%logT90

                if (isLgrb) then
                    ! convert photon count data to energy in units of ergs
                    GRB%Event(igrb)%logPbol = getLogPbol( GRB%Event(igrb)%logEpk, GRB%Event(igrb)%logPF53 )
                    GRB%Event(igrb)%logSbol = getLogPbol( GRB%Event(igrb)%logEpk, GRB%Event(igrb)%logSbol )
                else
                    GRB%Event(igrb)%logPbol = LN10 * GRB%Event(igrb)%logPbol
                end if

                ! write the converted data to output file
                write(outFileUnit,"(i30,7e30.6)") Trigger(igrb)                                     &
                                                , GRB%Event(igrb)%logPbol                           &
                                                , GRB%Event(igrb)%logSbol                           &
                                                , GRB%Event(igrb)%logEpk                            &
                                                , GRB%Event(igrb)%logEpk-GRB%Event(igrb)%logPbol    &
                                                , GRB%Event(igrb)%logEpk-GRB%Event(igrb)%logSbol    &
                                                , GRB%Event(igrb)%logSbol-GRB%Event(igrb)%logPbol   &
                                                , GRB%Event(igrb)%logT90                            &
                                                , GRB%Event(igrb)%logPF53

            end do
        close(outFileUnit)
        close(inFileUnit)

#if defined CAF_ENABLED
        sync images(*)

    else

        sync images(*)
        do igrb = 1, GRB%count
            GRB%Event(igrb)%logPbol = GRB[1]%Event(igrb)%logPbol
            GRB%Event(igrb)%logSbol = GRB[1]%Event(igrb)%logSbol
            GRB%Event(igrb)%logPF53 = GRB[1]%Event(igrb)%logPF53
            GRB%Event(igrb)%logEpk  = GRB[1]%Event(igrb)%logEpk
            GRB%Event(igrb)%logT90  = GRB[1]%Event(igrb)%logT90
        end do

    end if
#endif

    end subroutine readDataGRB

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getLog10PF53(log10epk,log10pbol) result(log10PF53)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog10PF53
#endif
        ! Given Log10(Epk [KeV]) of an LGRB and its bolometric (0.0001-20000 KeV) peak flux [in units of Ergs/s], Log(Pbol), 
        ! this function calculates the corresponding Log10(peak photon flux) in the BATSE detection energy range [50,300] KeV.
        ! Amir Shahmoradi, Wednesday June 27, 2012, 7:15 PM, IFS, The University of Texas at Austin.
        ! Amir Shahmoradi, Wednesday June 27, 2012, 9:29 PM, IFS, The University of Texas at Austin.
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: log10epk,log10pbol
        real(RK)                :: log10PF53
        if (log10epk<-2.915056638230699_RK) then
            log10PF53 = log10pbol    + 4.9200000000000000_RK
        elseif (log10epk<1.5_RK) then
            log10PF53 = log10pbol   + 5.7361200000000000_RK + log10epk * &
                                    ( 0.3093600000000000_RK + log10epk * &
                                    ( 0.0045600000000000_RK + log10epk * &
                                    ( 0.0015900000000000_RK + log10epk * &
                                    ( 0.0001533360000000_RK - log10epk   &
                                    * 0.0003574800000000_RK ))))
        elseif (log10epk<2.5_RK) then
            log10PF53 = log10pbol   + 1.9112800000000000_RK + log10epk * &
                                    ( 39.710390000000000_RK - log10epk * &
                                    ( 96.606280000000000_RK - log10epk * &
                                    ( 109.24696000000000_RK - log10epk * &
                                    ( 67.271800000000000_RK - log10epk * &
                                    ( 23.402390000000000_RK - log10epk * &
                                    ( 4.3454400000000000_RK - log10epk   &
                                    * 0.3360600000000000_RK ))))))
        elseif (log10epk<4._RK) then
            log10PF53 = log10pbol   + 2.8020600000000000_RK + log10epk * &
                                    ( 4.5690700000000000_RK - log10epk * &
                                    ( 1.9277200000000000_RK - log10epk * &
                                    ( 0.2938100000000000_RK - log10epk   &
                                    * 0.0148900000000000_RK )))
        elseif (log10epk<5.4093868613659435_RK) then
            log10PF53 = log10pbol   - 10.465330000000000_RK + log10epk * &
                                    ( 26.706370000000000_RK - log10epk * &
                                    ( 14.476310000000000_RK - log10epk * &
                                    ( 3.5404100000000000_RK - log10epk * &
                                    ( 0.4095700000000000_RK - log10epk   &
                                    * 0.0183100000000000_RK ))))
        else ! if (log10epk>=5.4093868613659435_RK) then
            log10PF53 = log10pbol   + 4.9200000000000000_RK
        end if
    end function getLog10PF53

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getLogPF53(logEpk,logPbol) result(logPF53)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPF53
#endif
        ! Given Log(Epk [KeV]) of an LGRB and its bolometric (0.0001-20000 KeV) peak flux [in units of Ergs/s], Log(Pbol), 
        ! this function calculates the corresponding Log(peak photon flux) in the BATSE detection energy range [50,300] KeV.
        ! Amir Shahmoradi, Wednesday June 27, 2012, 7:15 PM, IFS, The University of Texas at Austin.
        ! Amir Shahmoradi, Wednesday June 27, 2012, 9:29 PM, IFS, The University of Texas at Austin.
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: logEpk,logPbol
        real(RK)                :: logPF53
        if (logEpk<-6.712165960423344_RK) then
            logPF53 = logPbol    + 11.32871865753070600_RK
        elseif (logEpk<3.453877639491069_RK) then
            logPF53 = logPbol   + 13.20790440362500600_RK + logEpk * &
                                ( 0.309360000000000000_RK + logEpk * &
                                ( 0.001980382837478830_RK + logEpk * &
                                ( 0.000299892598248466_RK + logEpk * &
                                ( 1.25602147173493e-05_RK - logEpk   &
                                * 1.27171265917873e-05_RK ))))
        elseif (logEpk<5.756462732485115_RK) then
            logPF53 = logPbol   + 4.400884836537660000_RK + logEpk * &
                                ( 39.71039000000000000_RK - logEpk * &
                                ( 41.95557432120050000_RK - logEpk * &
                                ( 20.60525451895990000_RK - logEpk * &
                                ( 5.510436247342930000_RK - logEpk * &
                                ( 0.832525333390336000_RK - logEpk * &
                                ( 0.067135977132248900_RK - logEpk   &
                                * 0.002254876138523550_RK ))))))
        elseif (logEpk<9.210340371976184_RK) then
            logPF53 = logPbol   + 6.451981585674900000_RK + logEpk * &
                                ( 4.569070000000000000_RK - logEpk * &
                                ( 0.837198158654537000_RK - logEpk * &
                                ( 0.055416002698982300_RK - logEpk   &
                                * 0.001219684856402480_RK )))
        elseif (logEpk<12.455573549219071_RK) then
            logPF53 = logPbol   - 24.09731285126340000_RK + logEpk * &
                                ( 26.70637000000000000_RK - logEpk * &
                                ( 6.286981551320860000_RK - logEpk * &
                                ( 0.667762738216888000_RK - logEpk * &
                                ( 0.033549115287895400_RK - logEpk   &
                                * 0.000651366755890191_RK ))))
        else
            logPF53 = logPbol   + 11.32871865753070600_RK
        end if
    end function getLogPF53

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getLogPbol(logEpk,logPF53) result(logPbol)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogPbol
#endif
        ! Given Log10(Epk [KeV]) of an LGRB and its bolometric (0.0001-20000 KeV) peak flux [in units of Ergs/s], Log(Pbol), 
        ! this function calculates the corresponding Log10(peak photon flux) in the BATSE detection energy range [50,300] KeV.
        ! Amir Shahmoradi, Wednesday June 27, 2012, 7:15 PM, IFS, The University of Texas at Austin.
        ! Amir Shahmoradi, Wednesday June 27, 2012, 9:29 PM, IFS, The University of Texas at Austin.
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: logEpk, logPF53
        real(RK)                :: logPbol
        logPbol = logPF53 - getLogPF53(logEpk,0._RK)
    end function getLogPbol

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Batse_mod
