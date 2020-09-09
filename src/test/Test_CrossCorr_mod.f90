!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a 
!!!!   copy of this software and associated documentation files (the "Software"), 
!!!!   to deal in the Software without restriction, including without limitation 
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!!!!   and/or sell copies of the Software, and to permit persons to whom the 
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be 
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of 
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
!!!!   work (education/research/industry/development/...) by citing the ParaMonte 
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Test_CrossCorr_mod

    use Test_mod, only: Test_type
    use CrossCorr_mod
    implicit none

    private
    public :: test_CrossCorr

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_CrossCorr()

        use Constants_mod, only: IK, RK
        use Timer_mod, only: Timer_type
        use Err_mod, only: Err_type
        implicit none
        integer(IK)                 :: nd,np,id,ip,fileUnit,nlag,ilag, paddedLen
        integer(IK), allocatable    :: Lag(:), Weight(:)
        real(RK)   , allocatable    :: InputData(:,:), NormedData(:,:), AutoCorr(:,:), InverseSumNormedDataSq(:), Mean(:)
        real(RK)   , allocatable    :: NormedDataFFT1(:,:), NormedDataFFT2(:,:), AutoCorrFFT(:,:)
        type(Timer_type)            :: Timer
        type(Err_type)              :: Err

        Test = Test_type(moduleName=MODULE_NAME)

        blockFirstImage: if (Test%Image%id==1_IK) then

            call Test%testing("AutoCorr")
            Timer = Timer_type(Err)
            call Test%checkForErr(Err)

            ! read input data

            open(newunit=fileUnit,file=Test%inDir//"autoCorrDataWeightedCompact.txt",status="old")
            nd = 1
            np = 9985
            allocate(InputData(nd,np),Weight(np))
            do ip = 1,np
                read(fileUnit,*) Weight(ip), InputData(1:nd,ip)
            end do
            close(fileUnit)

            ! normalize data

            allocate(Mean(nd),NormedData(nd,np))
            Mean = 0._RK
            do id = 1,nd
                do ip = 1,np
                    Mean(id) = Mean(id) + InputData(id,ip)
                end do
                Mean(id) = Mean(id) / np
                NormedData(id,1:np) = InputData(id,1:np) - Mean(id)
            end do

            ! Generate Lags

            nlag = -1
            do
                nlag = nlag + 1
                if (2**nlag<np) cycle
                exit
            end do
            nlag = nlag - 1
            Lag = [(2**ilag,ilag=1,nlag)]
#if defined DBG_ENABLED
            write(*,*) nlag, Lag
#endif

            ! compute AutoCorr using Classic definition

            allocate(AutoCorr(nd,nlag),InverseSumNormedDataSq(nd))
            open(newunit=fileUnit,file=Test%outDir//"autoCorrOutput.txt",status="replace")
            write(*,*) "Computing AutoCorr for lag number using the classical AutoCorr definition..."

            InverseSumNormedDataSq = getInverseSumNormedDataSq(1_IK,np,NormedData)
            call Timer%tic()
            call getAutoCorrSlow(nd,np,NormedData(1:nd,1:np),nlag,Lag,AutoCorr,InverseSumNormedDataSq)
            call Timer%toc()
#if defined DBG_ENABLED
            write(*,*) "Done in ", Timer%Time%delta, " seconds."
#endif
            do ilag = 1,nlag
                write(fileUnit,*) Lag(ilag), AutoCorr(1:nd,ilag)
            end do
            close(fileUnit)

            !***********************************************************************************************************************
            ! compute AutoCorr using FFT
            !***********************************************************************************************************************

            paddedLen = getPaddedLen(np)
            allocate( NormedDataFFT1(paddedLen,nd) &
                    , NormedDataFFT2(paddedLen,nd) &
                    , AutoCorrFFT(paddedLen,nd) )
            do id = 1, nd
                NormedDataFFT1(1:paddedLen,id) = padZero(np,NormedData(id,1:np),paddedLen)
            end do
            NormedDataFFT2 = NormedDataFFT1
            !do id = 1, nd
            !    NormedDataFFT1(1:np,id) = NormedData(id,1:np)
            !    NormedDataFFT1(np+1:paddedLen,id) = 0
            !end do
            open(newunit=fileUnit,file=Test%outDir//"autoCorrOutputFFTCompact.txt",status="replace")
            write(*,*) "Computing AutoCorr for Compact Data number using FFT..."
            call Timer%tic()
            do id = 1,nd
                AutoCorrFFT(1:paddedLen,id) = getCrossCorrFFT( paddedLen                      &
                                                             , NormedDataFFT1(1:paddedLen,id) &
                                                             , NormedDataFFT1(1:paddedLen,id) )
                AutoCorrFFT(1:paddedLen,id) = AutoCorrFFT(1:paddedLen,id) / AutoCorrFFT(1,id)
                ! ensure NormedDataFFT1 does not change upon entering and exiting getCrossCorrFFT()
#if defined DBG_ENABLED
                do ip = 1,paddedLen
                    if (NormedDataFFT2(ip,id)/=NormedDataFFT1(ip,id)) then
                        write(*,*) "nonsense: ", NormedDataFFT1(ip,id), NormedDataFFT2(ip,id)
                    end if
                end do
#endif
            end do
            call Timer%toc()
#if defined DBG_ENABLED
            write(*,*) "Done in ", Timer%Time%delta, " seconds."
#endif
            do ip = 1, paddedLen
                write(fileUnit,*) ip-1, AutoCorrFFT(ip,1)
            end do
            close(fileUnit)

            !***********************************************************************************************************************
            ! compute AutoCorr using FFT with missing Weighted Data
            !***********************************************************************************************************************

            open(newunit=fileUnit,file=Test%outDir//"autoCorrOutputFFTCompactWeightMissing.txt",status="replace")
            write(*,*) "Computing AutoCorr for Compact Data using Weighted FFT Implementation..."
            call Timer%tic()
            do id = 1,nd
                AutoCorrFFT(1:paddedLen,id) = getCrossCorrFFTweighted   ( nData1 = np                   &
                                                                        , nData2 = np                   &
                                                                        , paddedLen = paddedLen         &
                                                                        , Data1 = NormedData(id,1:np)   &
                                                                        , Data2 = NormedData(id,1:np)   &
                                                                        )
                AutoCorrFFT(1:paddedLen,id) = AutoCorrFFT(1:paddedLen,id) / AutoCorrFFT(1,id)
            end do
            call Timer%toc()
#if defined DBG_ENABLED
            write(*,*) "Done in ", Timer%Time%delta, " seconds."
#endif
            do ip = 1, paddedLen
                write(fileUnit,*) ip-1, AutoCorrFFT(ip,1)
            end do
            close(fileUnit)

            !***********************************************************************************************************************
            ! compute AutoCorr using FFT with Weighted Data
            !***********************************************************************************************************************

            ! normalize data

            deallocate(Mean,NormedData,NormedDataFFT1,AutoCorrFFT)
            allocate(Mean(nd),NormedData(np,nd),NormedDataFFT1(np,nd))
            Mean = 0._RK
            do id = 1,nd
                do ip = 1,np
                    Mean(id) = Mean(id) + Weight(ip)*InputData(id,ip)
                end do
                Mean(id) = Mean(id) / sum(Weight)
                NormedData(1:np,id) = InputData(id,1:np) - Mean(id) 
            end do
            NormedDataFFT1 = NormedData

            paddedLen = getPaddedLen(sum(Weight))
            allocate(AutoCorrFFT(paddedLen,nd))
            open(newunit=fileUnit,file=Test%outDir//"autoCorrOutputFFTVerbose.txt",status="replace")
            write(*,*) "Computing AutoCorr for Weighted Compact Data using Weighted FFT implementation..."
            call Timer%tic()
            do id = 1,nd
                AutoCorrFFT(1:paddedLen,id) = getCrossCorrFFTweighted   ( nData1 = np                   &
                                                                        , nData2 = np                   &
                                                                        , paddedLen = paddedLen         &
                                                                        , Data1 = NormedData(1:np,id)   &
                                                                        , Data2 = NormedData(1:np,id)   &
                                                                        , Weight1 = Weight(1:np)        &
                                                                        , Weight2 = Weight(1:np)        &
                                                                        )
                AutoCorrFFT(1:paddedLen,id) = AutoCorrFFT(1:paddedLen,id) / AutoCorrFFT(1,id)
                ! ensure NormedData does not change upon entering and exiting getCrossCorrFFT()
#if defined DBG_ENABLED
                do ip = 1,np
                    if (NormedDataFFT1(ip,id)/=NormedData(ip,id)) then
                        write(*,*) "nonsense: ", NormedData(ip,id), NormedDataFFT1(ip,id)
                    end if
                end do
#endif
            end do
            call Timer%toc()
#if defined DBG_ENABLED
            write(*,*) "Done in ", Timer%Time%delta, " seconds."
#endif
            do ip = 1, paddedLen
                write(fileUnit,*) ip-1, AutoCorrFFT(ip,1)
            end do
            close(fileUnit)

            !***********************************************************************************************************************
            ! compute AutoCorr using FFT for a vector of ones
            !***********************************************************************************************************************

            ! normalize data

            nd = 1
            np = 100
            if (allocated(InputData)) deallocate(InputData)
            if (allocated(Weight)) deallocate(Weight)
            allocate(InputData(nd,np),Weight(np))
            deallocate(Mean,NormedData,NormedDataFFT1,AutoCorrFFT)
            allocate(Mean(nd),NormedData(np,nd),NormedDataFFT1(np,nd))
            InputData = 1._RK
            InputData(1,np) = .999999_RK
            Weight = 1_IK
            Mean = 0._RK
            do id = 1,nd
                do ip = 1,np
                    Mean(id) = Mean(id) + Weight(ip)*InputData(id,ip)
                end do
                Mean(id) = Mean(id) / sum(Weight)
                NormedData(1:np,id) = InputData(id,1:np) - Mean(id) 
            end do
            NormedDataFFT1 = NormedData

            paddedLen = getPaddedLen(sum(Weight))
            allocate(AutoCorrFFT(paddedLen,nd))
            open(newunit=fileUnit,file=Test%outDir//"autoCorrOutputFFTIdentity.txt",status="replace")
            write(*,*) "Computing AutoCorr for Identity vector using Weighted FFT implementation..."
            call Timer%tic()
            do id = 1,nd
                AutoCorrFFT(1:paddedLen,id) = getCrossCorrFFTweighted   ( nData1 = np                   &
                                                                        , nData2 = np                   &
                                                                        , paddedLen = paddedLen         &
                                                                        , Data1 = NormedData(1:np,id)   &
                                                                        , Data2 = NormedData(1:np,id)   &
                                                                        , Weight1 = Weight(1:np)        &
                                                                        , Weight2 = Weight(1:np)        &
                                                                        )
                AutoCorrFFT(1:paddedLen,id) = AutoCorrFFT(1:paddedLen,id) / AutoCorrFFT(1,id)
                ! ensure NormedData does not change upon entering and exiting getCrossCorrFFT()
#if defined DBG_ENABLED
                do ip = 1,np
                    if (NormedDataFFT1(ip,id)/=NormedData(ip,id)) then
                        write(*,*) "nonsense: ", NormedData(ip,id), NormedDataFFT1(ip,id)
                    end if
                end do
#endif
            end do
            call Timer%toc()
#if defined DBG_ENABLED
            write(*,*) "Done in ", Timer%Time%delta, " seconds."
#endif
            do ip = 1, paddedLen
                write(fileUnit,*) ip-1, AutoCorrFFT(ip,1)
            end do
            close(fileUnit)

        end if blockFirstImage

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()
        call Test%finalize()

    end subroutine test_CrossCorr

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_CrossCorr_mod