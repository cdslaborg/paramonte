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
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief
!> This file is part of PartitionOptDen_mod.f90 source file.
!> \author
!> Amir Shahmoradi
!> Saturday 00:37 am, March 13, 2021, Dallas, TX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Partition%neopt = Partition%neopt + 1
Partition%Size(Partition%neopt) = Partition%Try(it)%Size(Partition%Try(it)%ic)
Partition%Center(1:nd,Partition%neopt) = Partition%Try(it)%Center(1:nd,Partition%Try(it)%ic)
Partition%Membership( Partition%Try(it)%CumSumSize(Partition%Try(it)%ic-1)+1 : Partition%Try(it)%CumSumSize(Partition%Try(it)%ic) ) = Partition%neopt

if (Partition%Try(it)%LogVolRatio(Partition%Try(it)%ic) < 0._RK) then ! blockExpansionCheck : enlarge, no questions asked.

    Partition%LogLikeFitness(Partition%neopt) = 0._RK
    Partition%LogVolNormed(Partition%neopt) = Partition%LogVolNormed(Partition%neopt) - Partition%Try(it)%LogVolRatio(Partition%Try(it)%ic)

    scaleFactor = exp( -ndInverse * Partition%Try(it)%LogVolRatio(Partition%Try(it)%ic) )
    scaleFactorSq = scaleFactor**2
    scaleFactorSqInverse = 1._RK / scaleFactorSq
    Partition%Try(it)%ScaleFactorSq(Partition%Try(it)%ic) = Partition%Try(it)%ScaleFactorSq(Partition%Try(it)%ic) * scaleFactorSq

    Partition%EffectiveSize(Partition%neopt) = count(Partition%Try(it)%MahalSq(1:Partition%np,Partition%Try(it)%ic) <= Partition%Try(it)%ScaleFactorSq(Partition%Try(it)%ic))
    Partition%InvCovMat(1:nd,1:nd,Partition%neopt) = Partition%Try(it)%InvCovMat(1:nd,1:nd,Partition%Try(it)%ic) * scaleFactorSqInverse
    Partition%ChoDia(1:nd,Partition%neopt) = Partition%Try(it)%ChoDia(1:nd,Partition%Try(it)%ic) * scaleFactor
    do j = 1, nd
        Partition%ChoLowCovUpp(j+1:nd,j,Partition%neopt) = Partition%Try(it)%ChoLowCovUpp(j+1:nd,j,Partition%Try(it)%ic) * scaleFactor
    end do

else ! blockExpansionCheck

    Partition%ChoDia(1:nd,Partition%neopt) = Partition%Try(it)%ChoDia(1:nd,Partition%Try(it)%ic)
    Partition%LogVolNormed(Partition%neopt) = Partition%Try(it)%LogVolNormed(Partition%Try(it)%ic)
    Partition%InvCovMat(1:nd,1:nd,Partition%neopt) = Partition%Try(it)%InvCovMat(1:nd,1:nd,Partition%Try(it)%ic)
    Partition%ChoLowCovUpp(1:nd,1:nd,Partition%neopt) = Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,Partition%Try(it)%ic)

    Partition%EffectiveSize(Partition%neopt) = Partition%Try(it)%EffectiveSize(Partition%Try(it)%ic)
    if (Partition%Try(it)%IsTooLarge(Partition%Try(it)%ic)) Partition%LogLikeFitness(Partition%neopt) = -Partition%LogLikeFitness(Partition%neopt)

end if ! blockExpansionCheck

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
block
    integer :: ip
    integer :: ipend
    integer :: ipstart
    logical :: isInside
    real(RK) :: mahalSqScalar
    real(RK), allocatable :: NormedPoint(:)
    ipend = Partition%Try(it)%CumSumSize(Partition%Try(it)%ic)
    ipstart = Partition%Try(it)%CumSumSize(Partition%Try(it)%ic-1)+1
    do ip = ipstart, ipend
        NormedPoint = Point(:,ip) - Partition%Center(1:nd,Partition%neopt)
        mahalSqScalar = dot_product(NormedPoint,matmul(Partition%InvCovMat(:,:,Partition%neopt),NormedPoint))
        isInside = mahalSqScalar - 1._RK <= 1.e-6_RK
        if (.not. isInside) then
            ! LCOV_EXCL_START
            write(*,"(*(g0.15,:,' '))") new_line("a"), PROCEDURE_NAME//": FATAL - POINT NOT INSIDE!, MAHAL = ", sqrt(mahalSqScalar), new_line("a")
            write(*,"(60(g0,:,','))") "numRecursiveCall", Partition%numRecursiveCall
            write(*,"(60(g0,:,','))") "Partition%Try(it)%IsTooLarge(Partition%Try(it)%ic)", Partition%Try(it)%IsTooLarge(Partition%Try(it)%ic)
            !write(*,"(60(g0,:,','))") "KmeansNemax(Partition%neopt), Partition%neopt", KmeansNemax(Partition%neopt), Partition%neopt
            write(*,"(60(g0,:,','))") "ip, IpStart, IpEnd, size(Membership)", ip, ipstart, ipend, size(Partition%Membership(ipstart:ipend))
            write(*,"(60(g0,:,','))") "Membership", Partition%Membership(ipstart:ipend)
            write(*,"(60(g0,:,','))") "ChoLowCovUpp", Partition%ChoLowCovUpp(:,:,Partition%neopt)
            write(*,"(60(g0,:,','))") "InvCovMat", Partition%InvCovMat(:,:,Partition%neopt)
            write(*,"(60(g0,:,','))") "ChoDia", Partition%ChoDia(:,Partition%neopt)
            write(*,"(60(g0,:,','))") "NormedPoint", NormedPoint
            error stop
            ! LCOV_EXCL_STOP
        end if
    end do
end block
#endif

!#include "PartitionOptDen_mod.optimize.inc.f90"
