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

module Optimization_mod

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@Optimization_mod"

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: BrentMinimum_type
        integer(IK)     :: niter                        ! the number of iterations to reach the minimum of the function
        real(RK)        :: Bracket(3)                   ! the initial 3 Bracketing points that envelop the minimum
        real(RK)        :: xtol = sqrt(epsilon(1._RK))  ! the stopping rule tolerance
        real(RK)        :: xmin                         ! the x-value at the minimum of the function
        real(RK)        :: fmin                         ! the minimum of the function
        type(Err_type)  :: Err
    end type BrentMinimum_type

    ! overload Brent

    interface BrentMinimum_type
        module procedure :: minimizeBrent
    end interface BrentMinimum_type

    ! interfaces of the objective functions

    abstract interface
        function getFunc1D_proc(x) result(funcVal)
            import :: RK
            real(RK)    , intent(in)    :: x
            real(RK)                    :: funcVal
        end function getFunc1D_proc
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: PowellMinimum_type
        integer(IK)             :: niter                        ! the number of iterations to reach the minimum of the function
        integer(IK)             :: ndim                         ! the number of dimensions of the function
        real(RK)                :: ftol = sqrt(epsilon(1._RK))  ! the stopping rule tolerance for the value of function
        real(RK), allocatable   :: xmin(:)                      ! the x-value at the minimum of the function
        real(RK), allocatable   :: DirMat(:,:)                  ! an initial (ndim,ndim) matrix whose columns contain the initial set of directions (usually the ndim unit vectors)
       !real(RK), allocatable   :: StartVec(:)                  ! an initial (ndim) vector representing the start of the search
        real(RK)                :: fmin                         ! the minimum of the function
        type(Err_type)          :: Err
    end type PowellMinimum_type

    ! overload Powell

    interface PowellMinimum_type
        module procedure :: minimizePowell
    end interface PowellMinimum_type

    abstract interface
        function getFuncMD_proc(ndim,Point) result(funcVal)
            import :: RK, IK
            integer(IK) , intent(in)    :: ndim
            real(RK)    , intent(in)    :: Point(ndim)
            real(RK)                    :: funcVal
        end function getFuncMD_proc
    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface minimize
        module procedure :: minimizeBrent, minimizePowell
    end interface minimize

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function minimizeBrent(getFunc, x0, x1, x2, xtol) result(BrentMinimum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: minimizeBrent
#endif
        !   Compute the minimum of the input 1-dimensional function isolated to 
        !   a fractional precision of about xtol using Brent's method.
        !
        !   Input
        !   =====
        !
        !       getFunc
        !
        !           The 1-dimensional function which will have to be minimized.
        !
        !       x0, x1, x2
        !
        !           A set of optional bracketing triplet of abscissas that bracket the minimum of 
        !           the function such that, x0 < x1 < x2 and getFunc(x0) > getFunc(x1) < getFunc(x2).
        !
        !       xtol
        !
        !           An optional fractional precision within which is the minimum is returned.
        !           The default value is sqrt(epsilon(1._RK)).
        !
        !   Output
        !   ======
        !
        !       BrentMinimum
        !
        !           An object of type BrentMinimum_type that contains the minimum of the function (xmin) 
        !           and the function value at the minimum (fmin) as well as other relevant information.
        !
        use Constants_mod, only: IK, RK
        procedure(getFunc1D_proc)           :: getFunc
        real(RK)    , intent(in), optional  :: x0, x1, x2
        real(RK)    , intent(in), optional  :: xtol
        type(BrentMinimum_type)             :: BrentMinimum

        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME//"@minimizeBrent"
        integer(IK) , parameter             :: ITMAX = 1000_IK                  ! the maximum number of iterations
        real(RK)    , parameter             :: CGOLD = 0.3819660_RK             ! Golden Section switch criterion
        real(RK)    , parameter             :: ZEPS = sqrt(epsilon(1._RK))**3   ! tiny nonzero value

        integer(IK) :: iter
        logical     :: isPresentX0, isPresentX1, isPresentX2
        real(RK)    :: a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm
        real(RK)    :: ax, bx, cx, fa, fb, fc

        BrentMinimum%Err%occurred = .false.

        if (present(xtol)) BrentMinimum%xtol = xtol

        isPresentX0 = present(x0)
        isPresentX1 = present(x1)
        isPresentX2 = present(x2)
        if (isPresentX0 .and. isPresentX1 .and. isPresentX2) then
            BrentMinimum%Bracket = [x0, x1, x2]
        else
            if (isPresentX0) then
                ax = x0
            else
                ax = 0._RK ! assume an initial starting point of zero. This is the worst cae scenario.
            end if
            if (isPresentX1) then
                bx = x1
            else
                bx = ax + 1._RK
            end if
            call getBracket ( ax = ax &
                            , bx = bx &
                            , cx = cx &
                            , fa = fa &
                            , fb = fb &
                            , fc = fc &
                            , getFunc = getFunc &
                            )
            BrentMinimum%Bracket = [ax, bx, cx]
        end if

        a = min( BrentMinimum%Bracket(1), BrentMinimum%Bracket(3) )
        b = max( BrentMinimum%Bracket(1), BrentMinimum%Bracket(3) )
        v = BrentMinimum%Bracket(2)
        w = v
        x = v
        e = 0._RK
        fx = getFunc(x)
        fv = fx
        fw = fx
        do iter = 1, ITMAX
            BrentMinimum%niter = iter
            xm = 0.5_RK * (a+b)
            tol1 = xtol*abs(x) + ZEPS
            tol2 = 2.0_RK*tol1
            if (abs(x-xm) <= (tol2-0.5_RK*(b-a))) then
                BrentMinimum%xmin = x
                BrentMinimum%fmin = fx
                return
            end if
            if (abs(e) > tol1) then
                r = (x-w)*(fx-fv)
                q = (x-v)*(fx-fw)
                p = (x-v)*q-(x-w)*r
                q = 2.0_RK*(q-r)
                if (q > 0._RK) p = -p
                q = abs(q)
                etemp = e
                e = d
                if (abs(p) >= abs(0.5_RK*q*etemp) .or. p <= q*(a-x) .or. p >= q*(b-x)) then
                    e = merge(a-x,b-x, x >= xm )
                    d = cgold*e
                else
                    d = p/q
                    u = x+d
                    if (u-a < tol2 .or. b-u < tol2) d = sign(tol1,xm-x)
                end if
            else
                e = merge(a-x,b-x, x >= xm )
                d = cgold*e
            end if
            u = merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
            fu = getFunc(u)
            if (fu <= fx) then
                if (u >= x) then
                    a = x
                else
                    b = x
                end if
                call shft(v,w,x,u)
                call shft(fv,fw,fx,fu)
            else
                if (u < x) then
                    a = u
                else
                    b = u
                end if
                if (fu <= fw .or. w == x) then
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                else if (fu <= fv .or. v == x .or. v == w) then
                    v = u
                    fv = fu
                end if
            end if
        end do

        BrentMinimum%Err%occurred = .true.
        BrentMinimum%Err%msg = PROCEDURE_NAME//": maximum number of iterations exceeded."
        return

    contains

        subroutine shft(a,b,c,d)
            implicit none
            real(RK), intent(out) :: a
            real(RK), intent(inout) :: b,c
            real(RK), intent(in) :: d
            a = b
            b = c
            c = d
        end subroutine shft

    end function minimizeBrent

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getBracket(ax,bx,cx,fa,fb,fc,getFunc)

        use Constants_mod, only: IK, RK, TINY_RK
        use Misc_mod, only : swap
        implicit none

        real(RK), intent(inout)     :: ax, bx
        real(RK), intent(out)       :: cx, fa, fb, fc
        procedure(getFunc1D_proc)   :: getFunc

        real(RK), parameter         :: GOLD = 1.618034_RK
        real(RK), parameter         :: GLIMIT = 100.0_RK
        real(RK), parameter         :: TINY = TINY_RK ! 1.0e-20_RK
        real(RK)                    :: fu, q, r, u, ulim

        fa = getFunc(ax)
        fb = getFunc(bx)
        if (fb > fa) then
            call swap(ax,bx)
            call swap(fa,fb)
        end if

        cx = bx + GOLD*(bx-ax)
        fc = getFunc(cx)
        do
            if (fb < fc) return
            r = (bx-ax)*(fb-fc)
            q = (bx-cx)*(fb-fa)
            u = bx-((bx-cx)*q-(bx-ax)*r)/(2._RK*sign(max(abs(q-r),TINY),q-r))
            ulim = bx+GLIMIT*(cx-bx)
            if ((bx-u)*(u-cx) > 0._RK) then
                fu = getFunc(u)
                if (fu < fc) then
                    ax = bx
                    fa = fb
                    bx = u
                    fb = fu
                    return
                else if (fu > fb) then
                    cx = u
                    fc = fu
                    return
                end if
                u = cx+GOLD*(cx-bx)
                fu = getFunc(u)
            else if ((cx-u)*(u-ulim) > 0._RK) then
                fu = getFunc(u)
                if (fu < fc) then
                    bx = cx
                    cx = u
                    u = cx+GOLD*(cx-bx)
                    call shft(fb,fc,fu,getFunc(u))
                end if
            else if ((u-ulim)*(ulim-cx) >= 0._RK) then
                u = ulim
                fu = getFunc(u)
            else
                u = cx+GOLD*(cx-bx)
                fu = getFunc(u)
            end if
            call shft(ax,bx,cx,u)
            call shft(fa,fb,fc,fu)
        end do

    contains

        subroutine shft(a,b,c,d)
            real(RK), intent(out) :: a
            real(RK), intent(inout) :: b,c
            real(RK), intent(in) :: d
            a = b
            b = c
            c = d
        end subroutine shft

    end subroutine getBracket

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function minimizePowell(ndim, getFuncMD, StartVec, DirMat, ftol) result(PowellMinimum)

        use Constants_mod, only: IK, RK, TINY_RK ! tiny = 1.0e-25_RK
        implicit none

        procedure(getFuncMD_proc)               :: getFuncMD
        integer(IK) , intent(in)                :: ndim
        real(RK)    , intent(in)                :: StartVec(ndim)
        real(RK)    , intent(in)    , optional  :: DirMat(ndim,ndim)
        real(RK)    , intent(in)    , optional  :: ftol
        type(PowellMinimum_type)                :: PowellMinimum

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@minimizeBrent"
        integer(IK), parameter                  :: ITMAX = 1000

        integer(IK)                             :: i, ibig
        real(RK)                                :: del,fp,fptt,t
        real(RK)                                :: pt(ndim), ptt(ndim), xit(ndim)

        PowellMinimum%ndim = ndim
        PowellMinimum%Err%occurred = .false.

        allocate(PowellMinimum%xmin, source = StartVec)

        if (present(DirMat)) then
            PowellMinimum%DirMat = DirMat
        else
            allocate(PowellMinimum%DirMat(ndim,ndim), source = 0._RK)
            do i = 1, ndim
                PowellMinimum%DirMat(i,i) = 1._RK
            end do
        end if
        if (present(ftol)) PowellMinimum%ftol = ftol

        PowellMinimum%fmin = getFuncMD(ndim,PowellMinimum%xmin)
        pt = PowellMinimum%xmin
        PowellMinimum%niter = 0
        do

            PowellMinimum%niter = PowellMinimum%niter + 1
            fp = PowellMinimum%fmin
            ibig = 0
            del = 0._RK
            do i = 1, ndim
                xit = PowellMinimum%DirMat(1:ndim,i)
                fptt = PowellMinimum%fmin
                call linmin(getFuncMD, ndim, PowellMinimum%xmin, xit, PowellMinimum%fmin, PowellMinimum%Err)
                if (PowellMinimum%Err%occurred) then
                    PowellMinimum%Err%msg = PROCEDURE_NAME//PowellMinimum%Err%msg
                    return
                end if
                if (fptt - PowellMinimum%fmin > del) then
                    del = fptt - PowellMinimum%fmin
                    ibig = i
                end if
            end do

            if ( 2._RK*(fp-PowellMinimum%fmin) <= PowellMinimum%ftol*(abs(fp)+abs(PowellMinimum%fmin)) + TINY_RK ) return

            if (PowellMinimum%niter == ITMAX) then
                PowellMinimum%Err%occurred = .true.
                PowellMinimum%Err%msg = PROCEDURE_NAME//": maximum number of iterations exceeded."
                return
            end if

            ptt = 2._RK * PowellMinimum%xmin - pt
            xit = PowellMinimum%xmin - pt
            pt = PowellMinimum%xmin
            fptt = getFuncMD(ndim,ptt)
            if (fptt >= fp) cycle
            t = 2._RK * (fp-2._RK*PowellMinimum%fmin+fptt) * (fp-PowellMinimum%fmin-del)**2 - del*(fp-fptt)**2
            if (t >= 0.0) cycle
            call linmin(getFuncMD, ndim, PowellMinimum%xmin, xit, PowellMinimum%fmin, PowellMinimum%Err)
            if (PowellMinimum%Err%occurred) then
                PowellMinimum%Err%msg = PROCEDURE_NAME//PowellMinimum%Err%msg
                return
            end if
            PowellMinimum%DirMat(1:ndim,ibig) = PowellMinimum%DirMat(1:ndim,ndim)
            PowellMinimum%DirMat(1:ndim,ndim) = xit

        end do

    end function minimizePowell

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine linmin(getFuncMD, ndim, StartVec, dirVec, fmin, Err)
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        implicit none
        procedure(getFuncMD_proc)       :: getFuncMD
        integer(IK)     , intent(in)    :: ndim
        real(RK)        , intent(inout) :: StartVec(ndim), dirVec(ndim)
        real(RK)        , intent(out)   :: fmin
        type(Err_type)  , intent(out)   :: Err
        real(RK)        , parameter     :: XTOL = 1.0e-8_RK
        real(RK)                        :: ax, bx, fa, fb, fx, xx
        type(BrentMinimum_type)         :: BrentMinimum
        ax = 0.0
        xx = 1.0
        call getBracket(ax,xx,bx,fa,fx,fb,getFunc1D)
        BrentMinimum = minimizeBrent(getFunc1D, ax, xx, bx, XTOL)
        if (BrentMinimum%Err%occurred) then
            Err = BrentMinimum%Err
            return
        else
            Err%occurred = .false.
        end if
        fmin = BrentMinimum%fmin
        dirVec = BrentMinimum%xmin * dirVec
        StartVec = StartVec + dirVec
    contains
        function getFunc1D(x) result(funcVal)
            implicit none
            real(RK), intent(in)    :: x
            real(RK)                :: funcVal
            real(RK), allocatable   :: xt(:)
            xt = StartVec + x * dirVec
            funcVal = getFuncMD(ndim,xt)
        end function getFunc1D
    end subroutine linmin

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Optimization_mod