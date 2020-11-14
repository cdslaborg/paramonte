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

module QuadPackSPR_mod

    !use, intrinsic :: iso_fortran_env, only: IK => int32, RK => real64
    use Constants_mod, only: IK, RK

    implicit none

    character(*), parameter :: MODULE_NAME = "@QuadPackSPR_mod"
    character(*), parameter :: QAG_MSG_ERR_0 =  "ierr = 0: Normal and reliable termination of the routine: &
                                                &It is assumed that the requested accuracy has been achieved."
    character(*), parameter :: QAG_MSG_ERR_1 =  "ierr = 1: maximum number of subdivisions allowed has been achieved. &
                                                &Allow more subdivisions by increasing the value of LIMIT in QAG. &
                                                &However, if this yields no improvement it is advised to analyze the &
                                                &integrand to determine the integration difficulties.  If the position &
                                                &of a local difficulty can be determined, such as a singularity or &
                                                &discontinuity within the interval) one will probably gain from &
                                                &splitting up the interval at this point and calling the integrator &
                                                &on the subranges.  If possible, an appropriate special-purpose &
                                                &integrator should be used which is designed for handling the type &
                                                &of difficulty involved."
    character(*), parameter :: QAG_MSG_ERR_2 =  "ierr = 2: the occurrence of roundoff error is detected, which prevents the &
                                                &requested tolerance from being achieved."
    character(*), parameter :: QAG_MSG_ERR_3 =  "ierr = 3: extremely bad integrand behavior occurs at some points of the integration interval."
    character(*), parameter :: QAG_MSG_ERR_6 =  "ierr = 6: the input is invalid, because EPSABS < 0 and EPSREL < 0."

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine aaaa

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! AAAA is a dummy subroutine with QUADPACK documentation in its comments.
    !
    ! 1. introduction
    !
    !    quadpack is a fortran subroutine package for the numerical
    !    computation of definite one-dimensional integrals. it originated
    !    from a joint project of r. piessens and e. de doncker (appl.
    !    math. and progr. div.- k.u.leuven, belgium), c. ueberhuber (inst.
    !    fuer math.- techn.u.wien, austria), and d. kahaner (nation. bur.
    !    of standards- washington d.c., u.s.a.).
    !
    ! 2. survey
    !
    !    - qags : is an integrator based on globally adaptive interval
    !             subdivision in connection with extrapolation (de doncker,
    !             1978) by the epsilon algorithm (wynn, 1956).
    !
    !    - qagp : serves the same purposes as qags, but also allows
    !             for eventual user-supplied information, i.e. the
    !             abscissae of internal singularities, discontinuities
    !             and other difficulties of the integrand function.
    !             the algorithm is a modification of that in qags.
    !
    !    - qagi : handles integration over infinite intervals. the
    !             infinite range is mapped onto a finite interval and
    !             then the same strategy as in qags is applied.
    !
    !    - qawo : is a routine for the integration of cos(omega*x)*f(x)
    !             or sin(omega*x)*f(x) over a finite interval (a,b).
    !             omega is is specified by the user
    !             the rule evaluation component is based on the
    !             modified clenshaw-curtis technique.
    !             an adaptive subdivision scheme is used connected with
    !             an extrapolation procedure, which is a modification
    !             of that in qags and provides the possibility to deal
    !             even with singularities in f.
    !
    !    - qawf : calculates the fourier cosine or fourier sine
    !             transform of f(x), for user-supplied interval (a,
    !             infinity), omega, and f. the procedure of qawo is
    !             used on successive finite intervals, and convergence
    !             acceleration by means of the epsilon algorithm (wynn,
    !             1956) is applied to the series of the integral
    !             contributions.
    !
    !    - qaws : integrates w(x)*f(x) over (a,b) with a < b finite,
    !             and   w(x) = ((x-a)**alfa)*((b-x)**beta)*v(x)
    !             where v(x) = 1 or log(x-a) or log(b-x)
    !                            or log(x-a)*log(b-x)
    !             and   -1 < alfa, -1 < beta.
    !             the user specifies a, b, alfa, beta and the type of
    !             the function v.
    !             a globally adaptive subdivision strategy is applied,
    !             with modified clenshaw-curtis integration on the
    !             subintervals which contain a or b.
    !
    !    - qawc : computes the cauchy principal value of f(x)/(x-c)
    !             over a finite interval (a,b) and for
    !             user-determined c.
    !             the strategy is globally adaptive, and modified
    !             clenshaw-curtis integration is used on the subranges
    !             which contain the point x = c.
    !
    !  each of the routines above also has a "more detailed" version
    !    with a name ending in e, as qage.  these provide more
    !    information and control than the easier versions.
    !
    !
    !   the preceeding routines are all automatic.  that is, the user
    !      inputs his problem and an error tolerance.  the routine
    !      attempts to perform the integration to within the requested
    !      absolute or relative error.
    !   there are, in addition, a number of non-automatic integrators.
    !      these are most useful when the problem is such that the
    !      user knows that a fixed rule will provide the accuracy
    !      required.  typically they return an error estimate but make
    !      no attempt to satisfy any particular input error request.
    !
    !      qk15
    !      qk21
    !      qk31
    !      qk41
    !      qk51
    !      qk61
    !           estimate the integral on [a,b] using 15, 21,..., 61
    !           point rule and return an error estimate.
    !      qk15i 15 point rule for (semi)infinite interval.
    !      qk15w 15 point rule for special singular weight functions.
    !      qc25c 25 point rule for cauchy principal values
    !      qc25o 25 point rule for sin/cos integrand.
    !      qmomo integrates k-th degree chebychev polynomial times
    !            function with various explicit singularities.
    !
    ! 3. guidelines for the use of quadpack
    !
    !    here it is not our purpose to investigate the question when
    !    automatic quadrature should be used. we shall rather attempt
    !    to help the user who already made the decision to use quadpack,
    !    with selecting an appropriate routine or a combination of
    !    several routines for handling his problem.
    !
    !    for both quadrature over finite and over infinite intervals,
    !    one of the first questions to be answered by the user is
    !    related to the amount of computer time he wants to spend,
    !    versus his -own- time which would be needed, for example, for
    !    manual subdivision of the interval or other analytic
    !    manipulations.
    !
    !    (1) the user may not care about computer time, or not be
    !        willing to do any analysis of the problem. especially when
    !        only one or a few integrals must be calculated, this attitude
    !        can be perfectly reasonable. in this case it is clear that
    !        either the most sophisticated of the routines for finite
    !        intervals, qags, must be used, or its analogue for infinite
    !        intervals, qagi. these routines are able to cope with
    !        rather difficult, even with improper integrals.
    !        this way of proceeding may be expensive. but the integrator
    !        is supposed to give you an answer in return, with additional
    !        information in the case of a failure, through its error
    !        estimate and flag. yet it must be stressed that the programs
    !        cannot be totally reliable.
    !
    !    (2) the user may want to examine the integrand function.
    !        if bad local difficulties occur, such as a discontinuity, a
    !        singularity, derivative singularity or high peak at one or
    !        more points within the interval, the first advice is to
    !        split up the interval at these points. the integrand must
    !        then be examinated over each of the subintervals separately,
    !        so that a suitable integrator can be selected for each of
    !        them. if this yields problems involving relative accuracies
    !        to be imposed on -finite- subintervals, one can make use of
    !        qagp, which must be provided with the positions of the local
    !        difficulties. however, if strong singularities are present
    !        and a high accuracy is requested, application of qags on the
    !        subintervals may yield a better result.
    !
    !        for quadrature over finite intervals we thus dispose of qags
    !        and
    !        - qng for well-behaved integrands,
    !        - qag for functions with an oscillating behavior of a non
    !          specific type,
    !        - qawo for functions, eventually singular, containing a
    !          factor cos(omega*x) or sin(omega*x) where omega is known,
    !        - qaws for integrands with algebraico-logarithmic end point
    !          singularities of known type,
    !        - qawc for cauchy principal values.
    !
    !        remark
    !
    !        on return, the work arrays in the argument lists of the
    !        adaptive integrators contain information about the interval
    !        subdivision process and hence about the integrand behavior:
    !        the end points of the subintervals, the local integral
    !        contributions and error estimates, and eventually other
    !        characteristics. for this reason, and because of its simple
    !        globally adaptive nature, the routine qag in particular is
    !        well-suited for integrand examination. difficult spots can
    !        be located by investigating the error estimates on the
    !        subintervals.
    !
    !        for infinite intervals we provide only one general-purpose
    !        routine, qagi. it is based on the qags algorithm applied
    !        after a transformation of the original interval into (0,1).
    !        yet it may eventuate that another type of transformation is
    !        more appropriate, or one might prefer to break up the
    !        original interval and use qagi only on the infinite part
    !        and so on. these kinds of actions suggest a combined use of
    !        different quadpack integrators. note that, when the only
    !        difficulty is an integrand singularity at the finite
    !        integration limit, it will in general not be necessary to
    !        break up the interval, as qagi deals with several types of
    !        singularity at the boundary point of the integration range.
    !        it also handles slowly convergent improper integrals, on
    !        the condition that the integrand does not oscillate over
    !        the entire infinite interval. if it does we would advise
    !        to sum succeeding positive and negative contributions to
    !        the integral -e.g. integrate between the zeros- with one
    !        or more of the finite-range integrators, and apply
    !        convergence acceleration eventually by means of quadpack
    !        subroutine qelg which implements the epsilon algorithm.
    !        such quadrature problems include the fourier transform as
    !        a special case. yet for the latter we have an automatic
    !        integrator available, qawf.
    !
      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qag ( f, a, b, epsabs, epsrel, key, result, abserr, neval, ier )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAG approximates an integral over a finite interval.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral   
    !      I = integral of F over (A,B),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !    QAG is a simple globally adaptive integrator using the strategy of 
    !    Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
    !    Gauss-Kronrod quadrature formulae for the rule evaluation component. 
    !    The pairs of high degree of precision are suitable for handling
    !    integration difficulties due to a strongly oscillating integrand.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Input, integer(IK) KEY, chooses the order of the local integration rule:
    !    1,  7 Gauss points, 15 Gauss-Kronrod points,
    !    2, 10 Gauss points, 21 Gauss-Kronrod points,
    !    3, 15 Gauss points, 31 Gauss-Kronrod points,
    !    4, 20 Gauss points, 41 Gauss-Kronrod points,
    !    5, 25 Gauss points, 51 Gauss-Kronrod points,
    !    6, 30 Gauss points, 61 Gauss-Kronrod points.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !    Output, integer(IK) IER, return code.
    !    0, normal and reliable termination of the routine.  It is assumed that the 
    !      requested accuracy has been achieved.
    !    1, maximum number of subdivisions allowed has been achieved.  One can 
    !      allow more subdivisions by increasing the value of LIMIT in QAG. 
    !      However, if this yields no improvement it is advised to analyze the
    !      integrand to determine the integration difficulties.  If the position
    !      of a local difficulty can be determined, such as a singularity or
    !      discontinuity within the interval) one will probably gain from 
    !      splitting up the interval at this point and calling the integrator 
    !      on the subranges.  If possible, an appropriate special-purpose 
    !      integrator should be used which is designed for handling the type 
    !      of difficulty involved.
    !    2, the occurrence of roundoff error is detected, which prevents the
    !      requested tolerance from being achieved.
    !    3, extremely bad integrand behavior occurs at some points of the
    !      integration interval.
    !    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
    !
    !  Local parameters:
    !
    !    LIMIT is the maximum number of subintervals allowed in
    !    the subdivision process of QAGE.
    !
      implicit none

      integer(IK), parameter :: limit = 500

      real(RK) a
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) b
      real(RK) blist(limit)
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK), external :: f
      integer(IK) ier
      integer(IK) iord(limit)
      integer(IK) key
      integer(IK) last
      integer(IK) neval
      real(RK) result
      real(RK) rlist(limit)

      call qage ( f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
        ier, alist, blist, rlist, elist, iord, last )

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qage ( f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
      ier, alist, blist, rlist, elist, iord, last )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAGE estimates a definite integral.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral   
    !      I = integral of F over (A,B),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Input, integer(IK) KEY, chooses the order of the local integration rule:
    !    1,  7 Gauss points, 15 Gauss-Kronrod points,
    !    2, 10 Gauss points, 21 Gauss-Kronrod points,
    !    3, 15 Gauss points, 31 Gauss-Kronrod points,
    !    4, 20 Gauss points, 41 Gauss-Kronrod points,
    !    5, 25 Gauss points, 51 Gauss-Kronrod points,
    !    6, 30 Gauss points, 61 Gauss-Kronrod points.
    !
    !    Input, integer(IK) LIMIT, the maximum number of subintervals that
    !    can be used.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !    Output, integer(IK) IER, return code.
    !    0, normal and reliable termination of the routine.  It is assumed that the 
    !      requested accuracy has been achieved.
    !    1, maximum number of subdivisions allowed has been achieved.  One can 
    !      allow more subdivisions by increasing the value of LIMIT in QAG. 
    !      However, if this yields no improvement it is advised to analyze the
    !      integrand to determine the integration difficulties.  If the position
    !      of a local difficulty can be determined, such as a singularity or
    !      discontinuity within the interval) one will probably gain from 
    !      splitting up the interval at this point and calling the integrator 
    !      on the subranges.  If possible, an appropriate special-purpose 
    !      integrator should be used which is designed for handling the type 
    !      of difficulty involved.
    !    2, the occurrence of roundoff error is detected, which prevents the
    !      requested tolerance from being achieved.
    !    3, extremely bad integrand behavior occurs at some points of the
    !      integration interval.
    !    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
    !
    !    Workspace, real(RK) ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
    !    through LAST the left and right ends of the partition subintervals.
    !
    !    Workspace, real(RK) RLIST(LIMIT), contains in entries 1 through LAST
    !    the integral approximations on the subintervals.
    !
    !    Workspace, real(RK) ELIST(LIMIT), contains in entries 1 through LAST
    !    the absolute error estimates on the subintervals.
    !
    !    Output, integer(IK) IORD(LIMIT), the first K elements of which are pointers 
    !    to the error estimates over the subintervals, such that
    !    elist(iord(1)), ..., elist(iord(k)) form a decreasing sequence, with
    !    k = last if last <= (limit/2+2), and k = limit+1-last otherwise.
    !
    !    Output, integer(IK) LAST, the number of subintervals actually produced 
    !    in the subdivision process.
    !
    !  Local parameters:
    !
    !    alist     - list of left end points of all subintervals
    !                       considered up to now
    !    blist     - list of right end points of all subintervals
    !                       considered up to now
    !    elist(i)  - error estimate applying to rlist(i)
    !    maxerr    - pointer to the interval with largest error estimate
    !    errmax    - elist(maxerr)
    !    area      - sum of the integrals over the subintervals
    !    errsum    - sum of the errors over the subintervals
    !    errbnd    - requested accuracy max(epsabs,epsrel*abs(result))
    !    *****1    - variable for the left subinterval
    !    *****2    - variable for the right subinterval
    !    last      - index for subdivision
    !
      implicit none

      integer(IK) limit

      real(RK) a
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) area
      real(RK) area1
      real(RK) area12
      real(RK) area2
      real(RK) a1
      real(RK) a2
      real(RK) b
      real(RK) blist(limit)
      real(RK) b1
      real(RK) b2
      real(RK) c
      real(RK) defabs
      real(RK) defab1
      real(RK) defab2
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK) errbnd
      real(RK) errmax
      real(RK) error1
      real(RK) error2
      real(RK) erro12
      real(RK) errsum
      real(RK), external :: f
      integer(IK) ier
      integer(IK) iord(limit)
      integer(IK) iroff1
      integer(IK) iroff2
      integer(IK) key
      integer(IK) keyf
      integer(IK) last
      integer(IK) maxerr
      integer(IK) neval
      integer(IK) nrmax
      real(RK) resabs
      real(RK) result
      real(RK) rlist(limit)
    !
    !  Test on validity of parameters.
    !
      ier = 0
      neval = 0
      last = 0
      result = 0._RK
      abserr = 0._RK
      alist(1) = a
      blist(1) = b
      rlist(1) = 0._RK
      elist(1) = 0._RK
      iord(1) = 0

      if ( epsabs < 0._RK .and. epsrel < 0._RK ) then
        ier = 6
        return
      end if
    !
    !  First approximation to the integral.
    !
      keyf = key
      keyf = max ( keyf, 1 )
      keyf = min ( keyf, 6 )

      c = keyf
      neval = 0

      if ( keyf == 1 ) then
        call qk15 ( f, a, b, result, abserr, defabs, resabs )
      else if ( keyf == 2 ) then
        call qk21 ( f, a, b, result, abserr, defabs, resabs )
      else if ( keyf == 3 ) then
        call qk31 ( f, a, b, result, abserr, defabs, resabs )
      else if ( keyf == 4 ) then
        call qk41 ( f, a, b, result, abserr, defabs, resabs )
      else if ( keyf == 5 ) then
        call qk51 ( f, a, b, result, abserr, defabs, resabs )
      else if ( keyf == 6 ) then
        call qk61 ( f, a, b, result, abserr, defabs, resabs )
      end if

      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
    !
    !  Test on accuracy.
    !
      errbnd = max ( epsabs, epsrel * abs ( result ) )

      if ( abserr <= 5.0E+01 * epsilon ( defabs ) * defabs .and. errbnd < abserr ) then
        ier = 2
      end if

      if ( limit == 1 ) then
        ier = 1
      end if

      if ( ier /= 0 .or. &
        ( abserr <= errbnd .and. abserr /= resabs ) .or. &
        abserr == 0._RK ) then

        if ( keyf /= 1 ) then
          neval = (10*keyf+1) * (2*neval+1)
        else
          neval = 30 * neval + 15
        end if

        return

      end if
    !
    !  Initialization.
    !
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0

      do last = 2, limit
    !
    !  Bisect the subinterval with the largest error estimate.
    !
        a1 = alist(maxerr)
        b1 = 0.5E+00 * ( alist(maxerr) + blist(maxerr) )
        a2 = b1
        b2 = blist(maxerr)

        if ( keyf == 1 ) then
          call qk15 ( f, a1, b1, area1, error1, resabs, defab1 )
        else if ( keyf == 2 ) then
          call qk21 ( f, a1, b1, area1, error1, resabs, defab1 )
        else if ( keyf == 3 ) then
          call qk31 ( f, a1, b1, area1, error1, resabs, defab1 )
        else if ( keyf == 4 ) then
          call qk41 ( f, a1, b1, area1, error1, resabs, defab1)
        else if ( keyf == 5 ) then
          call qk51 ( f, a1, b1, area1, error1, resabs, defab1 )
        else if ( keyf == 6 ) then
          call qk61 ( f, a1, b1, area1, error1, resabs, defab1 )
        end if

        if ( keyf == 1 ) then
          call qk15 ( f, a2, b2, area2, error2, resabs, defab2 )
        else if ( keyf == 2 ) then
          call qk21 ( f, a2, b2, area2, error2, resabs, defab2 )
        else if ( keyf == 3 ) then
          call qk31 ( f, a2, b2, area2, error2, resabs, defab2 )
        else if ( keyf == 4 ) then
          call qk41 ( f, a2, b2, area2, error2, resabs, defab2 )
        else if ( keyf == 5 ) then
          call qk51 ( f, a2, b2, area2, error2, resabs, defab2 )
        else if ( keyf == 6 ) then
          call qk61 ( f, a2, b2, area2, error2, resabs, defab2 )
        end if
    !
    !  Improve previous approximations to integral and error and
    !  test for accuracy.
    !
        neval = neval + 1
        area12 = area1 + area2
        erro12 = error1 + error2
        errsum = errsum + erro12 - errmax
        area = area + area12 - rlist(maxerr)

        if ( defab1 /= error1 .and. defab2 /= error2 ) then

          if ( abs ( rlist(maxerr) - area12 ) <= 1.0E-05 * abs ( area12 ) .and. 9.9E-01 * errmax <= erro12 ) then
            iroff1 = iroff1 + 1
          end if

          if ( 10 < last .and. errmax < erro12 ) then
            iroff2 = iroff2 + 1
          end if

        end if

        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max ( epsabs, epsrel * abs ( area ) )
    !
    !  Test for roundoff error and eventually set error flag.
    !
        if ( errbnd < errsum ) then

          if ( 6 <= iroff1 .or. 20 <= iroff2 ) then
            ier = 2
          end if
    !
    !  Set error flag in the case that the number of subintervals
    !  equals limit.
    !
          if ( last == limit ) then
            ier = 1
          end if
    !
    !  Set error flag in the case of bad integrand behavior
    !  at a point of the integration range.
    !
          if ( max ( abs ( a1 ), abs ( b2 ) ) <= ( 1._RK + c * 1.0E+03 * &
            epsilon ( a1 ) ) * ( abs ( a2 ) + 1.0E+04 * tiny ( a2 ) ) ) then
            ier = 3
          end if

        end if
    !
    !  Append the newly-created intervals to the list.
    !
        if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
        else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
        end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with the largest error estimate (to be bisected next).
    !
        call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )
     
        if ( ier /= 0 .or. errsum <= errbnd ) then
          exit
        end if

      end do
    !
    !  Compute final result.
    !
      result = sum ( rlist(1:last) )

      abserr = errsum

      if ( keyf /= 1 ) then
        neval = ( 10 * keyf + 1 ) * ( 2 * neval + 1 )
      else
        neval = 30 * neval + 15
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qagi ( f, bound, inf, epsabs, epsrel, result, abserr, neval, ier )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAGI estimates an integral over a semi-infinite or infinite interval.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral   
    !      I = integral of F over (A, +Infinity), 
    !    or 
    !      I = integral of F over (-Infinity,A)
    !    or 
    !      I = integral of F over (-Infinity,+Infinity),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) BOUND, the value of the finite endpoint of the integration
    !    range, if any, that is, if INF is 1 or -1.
    !
    !    Input, integer(IK) INF, indicates the type of integration range.
    !    1:  (  BOUND,    +Infinity),
    !    -1: ( -Infinity,  BOUND),
    !    2:  ( -Infinity, +Infinity).
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !    Output, integer(IK) IER, error indicator.
    !    0, normal and reliable termination of the routine.  It is assumed that 
    !      the requested accuracy has been achieved.
    !    > 0,  abnormal termination of the routine.  The estimates for result
    !      and error are less reliable.  It is assumed that the requested
    !      accuracy has not been achieved.
    !    1, maximum number of subdivisions allowed has been achieved.  One can 
    !      allow more subdivisions by increasing the data value of LIMIT in QAGI
    !      (and taking the according dimension adjustments into account).
    !      However, if this yields no improvement it is advised to analyze the
    !      integrand in order to determine the integration difficulties.  If the
    !      position of a local difficulty can be determined (e.g. singularity,
    !      discontinuity within the interval) one will probably gain from
    !      splitting up the interval at this point and calling the integrator 
    !      on the subranges.  If possible, an appropriate special-purpose 
    !      integrator should be used, which is designed for handling the type
    !      of difficulty involved.
    !    2, the occurrence of roundoff error is detected, which prevents the
    !      requested tolerance from being achieved.  The error may be
    !      under-estimated.
    !    3, extremely bad integrand behavior occurs at some points of the
    !      integration interval.
    !    4, the algorithm does not converge.  Roundoff error is detected in the
    !      extrapolation table.  It is assumed that the requested tolerance
    !      cannot be achieved, and that the returned result is the best which 
    !      can be obtained.
    !    5, the integral is probably divergent, or slowly convergent.  It must 
    !      be noted that divergence can occur with any other value of IER.
    !    6, the input is invalid, because INF /= 1 and INF /= -1 and INF /= 2, or
    !      epsabs < 0 and epsrel < 0.  result, abserr, neval are set to zero.
    !
    !  Local parameters:
    !
    !            the dimension of rlist2 is determined by the value of
    !            limexp in QEXTR.
    !
    !           alist     - list of left end points of all subintervals
    !                       considered up to now
    !           blist     - list of right end points of all subintervals
    !                       considered up to now
    !           rlist(i)  - approximation to the integral over
    !                       (alist(i),blist(i))
    !           rlist2    - array of dimension at least (limexp+2),
    !                       containing the part of the epsilon table
    !                       which is still needed for further computations
    !           elist(i)  - error estimate applying to rlist(i)
    !           maxerr    - pointer to the interval with largest error
    !                       estimate
    !           errmax    - elist(maxerr)
    !           erlast    - error on the interval currently subdivided
    !                       (before that subdivision has taken place)
    !           area      - sum of the integrals over the subintervals
    !           errsum    - sum of the errors over the subintervals
    !           errbnd    - requested accuracy max(epsabs,epsrel*
    !                       abs(result))
    !           *****1    - variable for the left subinterval
    !           *****2    - variable for the right subinterval
    !           last      - index for subdivision
    !           nres      - number of calls to the extrapolation routine
    !           numrl2    - number of elements currently in rlist2. if an
    !                       appropriate approximation to the compounded
    !                       integral has been obtained, it is put in
    !                       rlist2(numrl2) after numrl2 has been increased
    !                       by one.
    !           small     - length of the smallest interval considered up
    !                       to now, multiplied by 1.5
    !           erlarg    - sum of the errors over the intervals larger
    !                       than the smallest interval considered up to now
    !           extrap    - logical variable denoting that the routine
    !                       is attempting to perform extrapolation. i.e.
    !                       before subdividing the smallest interval we
    !                       try to decrease the value of erlarg.
    !           noext     - logical variable denoting that extrapolation
    !                       is no longer allowed (true-value)
    !
      implicit none

      integer(IK), parameter :: limit = 500

      real(RK) abseps
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) area
      real(RK) area1
      real(RK) area12
      real(RK) area2
      real(RK) a1
      real(RK) a2
      real(RK) blist(limit)
      real(RK) boun
      real(RK) bound
      real(RK) b1
      real(RK) b2
      real(RK) correc
      real(RK) defabs
      real(RK) defab1
      real(RK) defab2
      real(RK) dres
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK) erlarg
      real(RK) erlast
      real(RK) errbnd
      real(RK) errmax
      real(RK) error1
      real(RK) error2
      real(RK) erro12
      real(RK) errsum
      real(RK) ertest
      logical extrap
      real(RK), external :: f
      integer(IK) id
      integer(IK) ier
      integer(IK) ierro
      integer(IK) inf
      integer(IK) iord(limit)
      integer(IK) iroff1
      integer(IK) iroff2
      integer(IK) iroff3
      integer(IK) jupbnd
      integer(IK) k
      integer(IK) ksgn
      integer(IK) ktmin
      integer(IK) last
      integer(IK) maxerr
      integer(IK) neval
      logical noext
      integer(IK) nres
      integer(IK) nrmax
      integer(IK) numrl2
      real(RK) resabs
      real(RK) reseps
      real(RK) result
      real(RK) res3la(3)
      real(RK) rlist(limit)
      real(RK) rlist2(52)
      real(RK) small
    !
    !  Test on validity of parameters.
    !
      ier = 0
      neval = 0
      last = 0
      result = 0._RK
      abserr = 0._RK
      alist(1) = 0._RK
      blist(1) = 1._RK
      rlist(1) = 0._RK
      elist(1) = 0._RK
      iord(1) = 0

      if ( epsabs < 0._RK .and. epsrel < 0._RK ) then
        ier = 6
        return
      end if
    !
    !  First approximation to the integral.
    !
    !  Determine the interval to be mapped onto (0,1).
    !  If INF = 2 the integral is computed as i = i1+i2, where
    !  i1 = integral of f over (-infinity,0),
    !  i2 = integral of f over (0,+infinity).
    !
      if ( inf == 2 ) then
        boun = 0._RK
      else
        boun = bound
      end if

      !call qk15i ( f, boun, inf, 0._RK, 1._RK, result, abserr, defabs, resabs )
      call qk15i ( f, boun, inf, 0._RK, 1._RK, result, abserr, defabs, resabs )
    !
    !  Test on accuracy.
    !
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      dres = abs ( result )
      errbnd = max ( epsabs, epsrel * dres )

      if ( abserr <= 100._RK * epsilon ( defabs ) * defabs .and. &
        errbnd < abserr ) then
        ier = 2
      end if

      if ( limit == 1 ) then
        ier = 1
      end if

      if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
        abserr == 0._RK ) go to 130
    !
    !  Initialization.
    !
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = huge ( abserr )
      nrmax = 1
      nres = 0
      ktmin = 0
      numrl2 = 2
      extrap = .false.
      noext = .false.
      ierro = 0
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0

      if ( ( 1._RK - 5.0E+01 * epsilon ( defabs ) ) * defabs <= dres ) then
        ksgn = 1
      else
        ksgn = -1
      end if

      do last = 2, limit
    !
    !  Bisect the subinterval with nrmax-th largest error estimate.
    !
        a1 = alist(maxerr)
        b1 = 5.0E-01_RK * ( alist(maxerr) + blist(maxerr) )
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call qk15i ( f, boun, inf, a1, b1, area1, error1, resabs, defab1 )
        call qk15i ( f, boun, inf, a2, b2, area2, error2, resabs, defab2 )
    !
    !  Improve previous approximations to integral and error
    !  and test for accuracy.
    !
        area12 = area1 + area2
        erro12 = error1 + error2
        errsum = errsum + erro12 - errmax
        area = area + area12 - rlist(maxerr)

        if ( defab1 /= error1 .and. defab2 /= error2 ) then

          if ( abs ( rlist(maxerr) - area12 ) <= 1.0E-05_RK * abs ( area12 ) &
            .and. 9.9E-01_RK * errmax <= erro12 ) then

            if ( extrap ) then
              iroff2 = iroff2 + 1
            end if

            if ( .not. extrap ) then
              iroff1 = iroff1 + 1
            end if

          end if

          if ( 10 < last .and. errmax < erro12 ) then
            iroff3 = iroff3 + 1
          end if

        end if

        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max ( epsabs, epsrel * abs ( area ) )
    !
    !  Test for roundoff error and eventually set error flag.
    !
        if ( 10 <= iroff1 + iroff2 .or. 20 <= iroff3 ) then
          ier = 2
        end if

        if ( 5 <= iroff2 ) then
          ierro = 3
        end if
    !
    !  Set error flag in the case that the number of subintervals equals LIMIT.
    !
        if ( last == limit ) then
          ier = 1
        end if
    !
    !  Set error flag in the case of bad integrand behavior
    !  at some points of the integration range.
    !
        if ( max ( abs(a1), abs(b2) ) <= (1._RK + 1.0E+03_RK * epsilon ( a1 ) ) * &
        ( abs(a2) + 1.0E+03_RK * tiny ( a2 ) )) then
          ier = 4
        end if
    !
    !  Append the newly-created intervals to the list.
    !
        if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
        else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
        end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with NRMAX-th largest error estimate (to be bisected next).
    !
        call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

        if ( errsum <= errbnd ) go to 115

        if ( ier /= 0 ) then
          exit
        end if

        if ( last == 2 ) then
          small = 3.75E-01_RK
          erlarg = errsum
          ertest = errbnd
          rlist2(2) = area
          cycle
        end if

        if ( noext ) then
          cycle
        end if

        erlarg = erlarg - erlast

        if ( small < abs ( b1 - a1 ) ) then
          erlarg = erlarg + erro12
        end if
    !
    !  Test whether the interval to be bisected next is the
    !  smallest interval.
    !
        if ( .not. extrap ) then

          if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
            cycle
          end if

          extrap = .true.
          nrmax = 2

        end if

        if ( ierro == 3 .or. erlarg <= ertest ) then
          go to 60
        end if
    !
    !  The smallest interval has the largest error.
    !  before bisecting decrease the sum of the errors over the
    !  larger intervals (erlarg) and perform extrapolation.
    !
        id = nrmax
        jupbnd = last

        if ( (2+limit/2) < last ) then
          jupbnd = limit + 3 - last
        end if

        do k = id, jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
            go to 90
          end if
          nrmax = nrmax + 1
        end do
    !
    !  Extrapolate.
    !
    60  continue

        numrl2 = numrl2 + 1
        rlist2(numrl2) = area
        call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres ) 
        ktmin = ktmin+1

        if ( 5 < ktmin .and. abserr < 1.0E-03_RK * errsum ) then
          ier = 5
        end if

        if ( abseps < abserr ) then

          ktmin = 0
          abserr = abseps
          result = reseps
          correc = erlarg
          ertest = max ( epsabs, epsrel * abs(reseps) )

          if ( abserr <= ertest ) then
            exit
          end if

        end if
    !
    !  Prepare bisection of the smallest interval.
    !
        if ( numrl2 == 1 ) then
          noext = .true.
        end if

        if ( ier == 5 ) then
          exit
        end if

        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small * 5.0E-01_RK
        erlarg = errsum

    90  continue

      end do
    !
    !  Set final result and error estimate.
    !
      if ( abserr == huge ( abserr ) ) then
        go to 115
      end if

      if ( ( ier + ierro ) == 0 ) then
        go to 110
      end if

      if ( ierro == 3 ) then
        abserr = abserr + correc
      end if

      if ( ier == 0 ) then
        ier = 3
      end if

      if ( result /= 0._RK .and. area /= 0._RK) then
        go to 105
      end if

      if ( errsum < abserr ) then
        go to 115
      end if

      if ( area == 0._RK ) then
        go to 130
      end if

      go to 110

    105 continue

      if ( errsum / abs ( area ) < abserr / abs ( result )  ) then
        go to 115
      end if
    !
    !  Test on divergence
    !
    110 continue

      if ( ksgn == (-1) .and. max ( abs(result), abs(area) ) <=  defabs * 1.0E-02_RK) go to 130

      if ( 1.0E-02_RK > (result/area) .or. (result/area) > 1.0E+02_RK .or. errsum > abs(area)) then
        ier = 6
      end if

      go to 130
    !
    !  Compute global integral sum.
    !
      115 continue

      result = sum ( rlist(1:last) )

      abserr = errsum
      130 continue

      neval = 30 * last - 15
      if ( inf == 2 ) then
        neval = 2 * neval
      end if

      if ( 2 < ier ) then
        ier = ier - 1
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qagp ( f, a, b, npts2, points, epsabs, epsrel, result, abserr, neval, ier )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAGP computes a definite integral.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral   
    !      I = integral of F over (A,B),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !    Interior break points of the integration interval,
    !    where local difficulties of the integrand may occur, such as
    !    singularities or discontinuities, are provided by the user.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, 
    !    of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, integer(IK) NPTS2, the number of user-supplied break points within 
    !    the integration range, plus 2.  NPTS2 must be at least 2.
    !
    !    Input/output, real(RK) POINTS(NPTS2), contains the user provided interior
    !    breakpoints in entries 1 through NPTS2-2.  If these points are not
    !    in ascending order on input, they will be sorted.
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !    Output, integer(IK) IER, return flag.
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the requested
    !                             accuracy has been achieved.
    !                     ier > 0 abnormal termination of the routine.
    !                             the estimates for integral and error are
    !                             less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                     ier = 1 maximum number of subdivisions allowed
    !                             has been achieved. one can allow more
    !                             subdivisions by increasing the data value
    !                             of limit in qagp(and taking the according
    !                             dimension adjustments into account).
    !                             however, if this yields no improvement
    !                             it is advised to analyze the integrand
    !                             in order to determine the integration
    !                             difficulties. if the position of a local
    !                             difficulty can be determined (i.e.
    !                             singularity, discontinuity within the
    !                             interval), it should be supplied to the
    !                             routine as an element of the vector
    !                             points. if necessary, an appropriate
    !                             special-purpose integrator must be used,
    !                             which is designed for handling the type
    !                             of difficulty involved.
    !                         = 2 the occurrence of roundoff error is
    !                             detected, which prevents the requested
    !                             tolerance from being achieved.
    !                             the error may be under-estimated.
    !                         = 3 extremely bad integrand behavior occurs
    !                             at some points of the integration
    !                             interval.
    !                         = 4 the algorithm does not converge. roundoff
    !                             error is detected in the extrapolation
    !                             table. it is presumed that the requested
    !                             tolerance cannot be achieved, and that
    !                             the returned result is the best which
    !                             can be obtained.
    !                         = 5 the integral is probably divergent, or
    !                             slowly convergent. it must be noted that
    !                             divergence can occur with any other value
    !                             of ier > 0.
    !                         = 6 the input is invalid because
    !                             npts2 < 2 or
    !                             break points are specified outside
    !                             the integration range or
    !                             epsabs < 0 and epsrel < 0,
    !                             or limit < npts2.
    !                             result, abserr, neval are set to zero.
    !
    !  Local parameters:
    !
    !            the dimension of rlist2 is determined by the value of
    !            limexp in QEXTR (rlist2 should be of dimension
    !            (limexp+2) at least).
    !
    !           alist     - list of left end points of all subintervals
    !                       considered up to now
    !           blist     - list of right end points of all subintervals
    !                       considered up to now
    !           rlist(i)  - approximation to the integral over
    !                       (alist(i),blist(i))
    !           rlist2    - array of dimension at least limexp+2
    !                       containing the part of the epsilon table which
    !                       is still needed for further computations
    !           elist(i)  - error estimate applying to rlist(i)
    !           maxerr    - pointer to the interval with largest error
    !                       estimate
    !           errmax    - elist(maxerr)
    !           erlast    - error on the interval currently subdivided
    !                       (before that subdivision has taken place)
    !           area      - sum of the integrals over the subintervals
    !           errsum    - sum of the errors over the subintervals
    !           errbnd    - requested accuracy max(epsabs,epsrel*
    !                       abs(result))
    !           *****1    - variable for the left subinterval
    !           *****2    - variable for the right subinterval
    !           last      - index for subdivision
    !           nres      - number of calls to the extrapolation routine
    !           numrl2    - number of elements in rlist2. if an appropriate
    !                       approximation to the compounded integral has
    !                       obtained, it is put in rlist2(numrl2) after
    !                       numrl2 has been increased by one.
    !           erlarg    - sum of the errors over the intervals larger
    !                       than the smallest interval considered up to now
    !           extrap    - logical variable denoting that the routine
    !                       is attempting to perform extrapolation. i.e.
    !                       before subdividing the smallest interval we
    !                       try to decrease the value of erlarg.
    !           noext     - logical variable denoting that extrapolation is
    !                       no longer allowed (true-value)
    !
      implicit none

      integer(IK), parameter :: limit = 500

      real(RK) a
      real(RK) abseps
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) area
      real(RK) area1
      real(RK) area12
      real(RK) area2
      real(RK) a1
      real(RK) a2
      real(RK) b
      real(RK) blist(limit)
      real(RK) b1
      real(RK) b2
      real(RK) correc
      real(RK) defabs
      real(RK) defab1
      real(RK) defab2
      real(RK) dres
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK) erlarg
      real(RK) erlast
      real(RK) errbnd
      real(RK) errmax
      real(RK) error1
      real(RK) erro12
      real(RK) error2
      real(RK) errsum
      real(RK) ertest
      logical extrap
      real(RK), external :: f
      integer(IK) i
      integer(IK) id
      integer(IK) ier
      integer(IK) ierro
      integer(IK) ind1
      integer(IK) ind2
      integer(IK) iord(limit)
      integer(IK) iroff1
      integer(IK) iroff2
      integer(IK) iroff3
      integer(IK) j
      integer(IK) jlow
      integer(IK) jupbnd
      integer(IK) k
      integer(IK) ksgn
      integer(IK) ktmin
      integer(IK) last
      integer(IK) levcur
      integer(IK) level(limit)
      integer(IK) levmax
      integer(IK) maxerr
      integer(IK) ndin(40)
      integer(IK) neval
      integer(IK) nint
      logical noext
      integer(IK) npts
      integer(IK) npts2
      integer(IK) nres
      integer(IK) nrmax
      integer(IK) numrl2
      real(RK) points(40)
      real(RK) pts(40)
      real(RK) resa
      real(RK) resabs
      real(RK) reseps
      real(RK) result
      real(RK) res3la(3)
      real(RK) rlist(limit)
      real(RK) rlist2(52)
      real(RK) sign
      real(RK) temp
    !
    !  Test on validity of parameters.
    !
      ier = 0
      neval = 0
      last = 0
      result = 0._RK
      abserr = 0._RK
      alist(1) = a
      blist(1) = b
      rlist(1) = 0._RK
      elist(1) = 0._RK
      iord(1) = 0
      level(1) = 0
      npts = npts2 - 2

      if ( npts2 < 2 ) then
        ier = 6
        return
      else if ( limit <= npts .or. ( epsabs < 0._RK .and. &
        epsrel < 0._RK) ) then
        ier = 6
        return
      end if
    !
    !  If any break points are provided, sort them into an
    !  ascending sequence.
    !
      if ( b < a ) then
        sign = -1._RK
      else
        sign = +1._RK
      end if

      pts(1) = min ( a, b )

      do i = 1, npts
        pts(i+1) = points(i)
      end do

      pts(npts+2) = max ( a, b )
      nint = npts+1
      a1 = pts(1)

      if ( npts /= 0 ) then

        do i = 1, nint
          do j = i+1, nint+1
            if ( pts(j) < pts(i) ) then
              temp   = pts(i)
              pts(i) = pts(j)
              pts(j) = temp
            end if
          end do
        end do

        if ( pts(1) /= min ( a, b ) .or. pts(nint+1) /= max ( a, b ) ) then
          ier = 6
          return
        end if

      end if
    !
    !  Compute first integral and error approximations.
    !
      resabs = 0._RK

      do i = 1, nint

        b1 = pts(i+1)
        call qk21 ( f, a1, b1, area1, error1, defabs, resa )
        abserr = abserr + error1
        result = result + area1
        ndin(i) = 0

        if ( error1 == resa .and. error1 /= 0._RK ) then
          ndin(i) = 1
        end if

        resabs = resabs + defabs
        level(i) = 0
        elist(i) = error1
        alist(i) = a1
        blist(i) = b1
        rlist(i) = area1
        iord(i) = i
        a1 = b1

      end do

      errsum = 0._RK

      do i = 1, nint
        if ( ndin(i) == 1 ) then
          elist(i) = abserr
        end if
        errsum = errsum + elist(i)
      end do
    !
    !  Test on accuracy.
    !
      last = nint
      neval = 21 * nint
      dres = abs ( result )
      errbnd = max ( epsabs, epsrel * dres )

      if ( abserr <= 1.0E+02_RK * epsilon ( resabs ) * resabs .and. &
        abserr > errbnd ) then
        ier = 2
      end if

      if ( nint /= 1 ) then

        do i = 1, npts

          jlow = i+1
          ind1 = iord(i)

          do j = jlow, nint
            ind2 = iord(j)
            if ( elist(ind1) <= elist(ind2) ) then
              ind1 = ind2
              k = j
            end if
          end do

          if ( ind1 /= iord(i) ) then
            iord(k) = iord(i)
            iord(i) = ind1
          end if

        end do

        if ( limit < npts2 ) then
          ier = 1
        end if

      end if

      if ( ier /= 0 .or. abserr <= errbnd ) then
        return
      end if
    !
    !  Initialization
    !
      rlist2(1) = result
      maxerr = iord(1)
      errmax = elist(maxerr)
      area = result
      nrmax = 1
      nres = 0
      numrl2 = 1
      ktmin = 0
      extrap = .false.
      noext = .false.
      erlarg = errsum
      ertest = errbnd
      levmax = 1
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ierro = 0
      abserr = huge ( abserr )

      if ( dres >= ( 1._RK - 0.5E+00_RK * epsilon ( resabs ) ) * resabs ) then
        ksgn = 1
      else
        ksgn = -1
      end if

      do last = npts2, limit
    !
    !  Bisect the subinterval with the NRMAX-th largest error estimate.
    !
        levcur = level(maxerr) + 1
        a1 = alist(maxerr)
        b1 = 0.5E+00_RK * ( alist(maxerr) + blist(maxerr) )
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call qk21 ( f, a1, b1, area1, error1, resa, defab1 )
        call qk21 ( f, a2, b2, area2, error2, resa, defab2 )
    !
    !  Improve previous approximations to integral and error
    !  and test for accuracy.
    !
        neval = neval + 42
        area12 = area1 + area2
        erro12 = error1 + error2
        errsum = errsum + erro12 -errmax
        area = area + area12 - rlist(maxerr)

        if ( defab1 /= error1 .and. defab2 /= error2 ) then

          if ( abs ( rlist ( maxerr ) - area12 ) <= 1.0E-05 * abs(area12) .and. &
            erro12 >= 9.9E-01_RK * errmax ) then

            if ( extrap ) then
              iroff2 = iroff2+1
            else
              iroff1 = iroff1+1
            end if

          end if

          if ( last > 10 .and. erro12 > errmax ) then
            iroff3 = iroff3 + 1
          end if

        end if

        level(maxerr) = levcur
        level(last) = levcur
        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max ( epsabs, epsrel * abs ( area ) )
    !
    !  Test for roundoff error and eventually set error flag.
    !
        if ( 10 <= iroff1 + iroff2 .or. 20 <= iroff3 ) then
          ier = 2
        end if

        if ( 5 <= iroff2 ) then
          ierro = 3
        end if
    !
    !  Set error flag in the case that the number of subintervals
    !  equals limit.
    !
        if ( last == limit ) then
          ier = 1
        end if
    !
    !  Set error flag in the case of bad integrand behavior
    !  at a point of the integration range
    !
        if ( max ( abs(a1), abs(b2)) <= ( 1._RK + 1.0E+03_RK * epsilon ( a1 ) )*( abs ( a2 ) + 1.0E+03 * tiny ( a2 ) ) ) then
          ier = 4
        end if
    !
    !  Append the newly-created intervals to the list.
    !
        if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
        else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
        end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with nrmax-th largest error estimate (to be bisected next).
    !
        call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

        if ( errsum <= errbnd ) then
          go to 190
        end if

        if ( ier /= 0 ) then
          exit
        end if

        if ( noext ) then
          cycle
        end if

        erlarg = erlarg - erlast

        if ( levcur+1 <= levmax ) then
          erlarg = erlarg + erro12
        end if
    !
    !  Test whether the interval to be bisected next is the
    !  smallest interval.
    !
        if ( .not. extrap ) then

          if ( level(maxerr)+1 <= levmax ) then
            cycle
          end if

          extrap = .true.
          nrmax = 2

        end if
    !
    !  The smallest interval has the largest error.
    !  Before bisecting decrease the sum of the errors over the
    !  larger intervals (erlarg) and perform extrapolation.
    !
        if ( ierro /= 3 .and. erlarg > ertest ) then

          id = nrmax
          jupbnd = last
          if ( last > (2+limit/2) ) then
            jupbnd = limit+3-last
          end if

          do k = id, jupbnd
            maxerr = iord(nrmax)
            errmax = elist(maxerr)
            if ( level(maxerr)+1 <= levmax ) then
              go to 160
            end if
            nrmax = nrmax + 1
          end do

        end if
    !
    !  Perform extrapolation.
    !
        numrl2 = numrl2 + 1
        rlist2(numrl2) = area

        if ( numrl2 <= 2 ) then
          go to 155
        end if

        call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
        ktmin = ktmin+1

        if ( 5 < ktmin .and. abserr < 1.0E-03_RK * errsum ) then
          ier = 5
        end if

        if ( abseps < abserr ) then

          ktmin = 0
          abserr = abseps
          result = reseps
          correc = erlarg
          ertest = max ( epsabs, epsrel * abs(reseps) )

          if ( abserr < ertest ) then
            exit
          end if

        end if
    !
    !  Prepare bisection of the smallest interval.
    !
        if ( numrl2 == 1 ) then
          noext = .true.
        end if

        if ( 5 <= ier ) then
          exit
        end if

    155 continue

        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        levmax = levmax + 1
        erlarg = errsum

    160 continue

      end do
    !
    !  Set the final result.
    !
      if ( abserr == huge ( abserr ) ) then
        go to 190
      end if

      if ( ( ier + ierro ) == 0 ) then
        go to 180
      end if

      if ( ierro == 3 ) then
        abserr = abserr + correc
      end if

      if ( ier == 0 ) then
        ier = 3
      end if

      if ( result /= 0._RK .and. area /= 0._RK ) then
        go to 175
      end if

      if ( errsum < abserr ) then
        go to 190
      end if

      if ( area == 0._RK ) then
        go to 210
      end if

      go to 180

    175 continue

      if ( abserr / abs(result) > errsum / abs(area) ) then
        go to 190
      end if
    !
    !  Test on divergence.
    !
      180 continue

      if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <= resabs*1.0E-02_RK ) go to 210

      if ( 1.0E-02_RK > (result/area) .or. (result/area) > 1.0E+02 .or. errsum > abs(area) ) then
        ier = 6
      end if

      go to 210
    !
    !  Compute global integral sum.
    !
    190 continue

      result = sum ( rlist(1:last) )

      abserr = errsum

    210 continue

      if ( 2 < ier ) then
        ier = ier - 1
      end if

      result = result * sign

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qags ( f, a, b, epsabs, epsrel, result, abserr, neval, ier )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAGS estimates the integral of a function.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral   
    !      I = integral of F over (A,B),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !    Output, integer(IK) IER, error flag.
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the requested
    !                             accuracy has been achieved.
    !                     ier > 0 abnormal termination of the routine
    !                             the estimates for integral and error are
    !                             less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                         = 1 maximum number of subdivisions allowed
    !                             has been achieved. one can allow more sub-
    !                             divisions by increasing the data value of
    !                             limit in qags (and taking the according
    !                             dimension adjustments into account).
    !                             however, if this yields no improvement
    !                             it is advised to analyze the integrand
    !                             in order to determine the integration
    !                             difficulties. if the position of a
    !                             local difficulty can be determined (e.g.
    !                             singularity, discontinuity within the
    !                             interval) one will probably gain from
    !                             splitting up the interval at this point
    !                             and calling the integrator on the sub-
    !                             ranges. if possible, an appropriate
    !                             special-purpose integrator should be used,
    !                             which is designed for handling the type
    !                             of difficulty involved.
    !                         = 2 the occurrence of roundoff error is detec-
    !                             ted, which prevents the requested
    !                             tolerance from being achieved.
    !                             the error may be under-estimated.
    !                         = 3 extremely bad integrand behavior occurs
    !                             at some  points of the integration
    !                             interval.
    !                         = 4 the algorithm does not converge. roundoff
    !                             error is detected in the extrapolation
    !                             table. it is presumed that the requested
    !                             tolerance cannot be achieved, and that the
    !                             returned result is the best which can be
    !                             obtained.
    !                         = 5 the integral is probably divergent, or
    !                             slowly convergent. it must be noted that
    !                             divergence can occur with any other value
    !                             of ier.
    !                         = 6 the input is invalid, because
    !                             epsabs < 0 and epsrel < 0,
    !                             result, abserr and neval are set to zero.
    !
    !  Local Parameters:
    !
    !           alist     - list of left end points of all subintervals
    !                       considered up to now
    !           blist     - list of right end points of all subintervals
    !                       considered up to now
    !           rlist(i)  - approximation to the integral over
    !                       (alist(i),blist(i))
    !           rlist2    - array of dimension at least limexp+2 containing
    !                       the part of the epsilon table which is still
    !                       needed for further computations
    !           elist(i)  - error estimate applying to rlist(i)
    !           maxerr    - pointer to the interval with largest error
    !                       estimate
    !           errmax    - elist(maxerr)
    !           erlast    - error on the interval currently subdivided
    !                       (before that subdivision has taken place)
    !           area      - sum of the integrals over the subintervals
    !           errsum    - sum of the errors over the subintervals
    !           errbnd    - requested accuracy max(epsabs,epsrel*
    !                       abs(result))
    !           *****1    - variable for the left interval
    !           *****2    - variable for the right interval
    !           last      - index for subdivision
    !           nres      - number of calls to the extrapolation routine
    !           numrl2    - number of elements currently in rlist2. if an
    !                       appropriate approximation to the compounded
    !                       integral has been obtained it is put in
    !                       rlist2(numrl2) after numrl2 has been increased
    !                       by one.
    !           small     - length of the smallest interval considered
    !                       up to now, multiplied by 1.5
    !           erlarg    - sum of the errors over the intervals larger
    !                       than the smallest interval considered up to now
    !           extrap    - logical variable denoting that the routine is
    !                       attempting to perform extrapolation i.e. before
    !                       subdividing the smallest interval we try to
    !                       decrease the value of erlarg.
    !           noext     - logical variable denoting that extrapolation
    !                       is no longer allowed (true value)
    !
      implicit none

      integer(IK), parameter :: limit = 500

      real(RK) a
      real(RK) abseps
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) area
      real(RK) area1
      real(RK) area12
      real(RK) area2
      real(RK) a1
      real(RK) a2
      real(RK) b
      real(RK) blist(limit)
      real(RK) b1
      real(RK) b2
      real(RK) correc
      real(RK) defabs
      real(RK) defab1
      real(RK) defab2
      real(RK) dres
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK) erlarg
      real(RK) erlast
      real(RK) errbnd
      real(RK) errmax
      real(RK) error1
      real(RK) error2
      real(RK) erro12
      real(RK) errsum
      real(RK) ertest
      logical extrap
      real(RK), external :: f
      integer(IK) id
      integer(IK) ier
      integer(IK) ierro
      integer(IK) iord(limit)
      integer(IK) iroff1
      integer(IK) iroff2
      integer(IK) iroff3
      integer(IK) jupbnd
      integer(IK) k
      integer(IK) ksgn
      integer(IK) ktmin
      integer(IK) last
      logical noext
      integer(IK) maxerr
      integer(IK) neval
      integer(IK) nres
      integer(IK) nrmax
      integer(IK) numrl2
      real(RK) resabs
      real(RK) reseps
      real(RK) result
      real(RK) res3la(3)
      real(RK) rlist(limit)
      real(RK) rlist2(52)
      real(RK) small
    !
    !  The dimension of rlist2 is determined by the value of
    !  limexp in QEXTR (rlist2 should be of dimension
    !  (limexp+2) at least).
    !
    !  Test on validity of parameters.
    !
      ier = 0
      neval = 0
      last = 0
      result = 0._RK
      abserr = 0._RK
      alist(1) = a
      blist(1) = b
      rlist(1) = 0._RK
      elist(1) = 0._RK

      if ( epsabs < 0._RK .and. epsrel < 0._RK ) then
        ier = 6
        return
      end if
    !
    !  First approximation to the integral.
    !
      ierro = 0
      call qk21 ( f, a, b, result, abserr, defabs, resabs )
    !
    !  Test on accuracy.
    !
      dres = abs ( result )
      errbnd = max ( epsabs, epsrel * dres )
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1

      if ( abserr <= 1.0E+02 * epsilon ( defabs ) * defabs .and. &
        abserr > errbnd ) then
        ier = 2
      end if

      if ( limit == 1 ) then
        ier = 1
      end if

      if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
        abserr == 0._RK ) go to 140
    !
    !  Initialization.
    !
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = huge ( abserr )
      nrmax = 1
      nres = 0
      numrl2 = 2
      ktmin = 0
      extrap = .false.
      noext = .false.
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0

      if ( dres >= (1._RK - 5.0E+01_RK* epsilon ( defabs ) ) * defabs ) then
        ksgn = 1
      else
        ksgn = -1
      end if

      do last = 2, limit
    !
    !  Bisect the subinterval with the nrmax-th largest error estimate.
    !
        a1 = alist(maxerr)
        b1 = 5.0E-01_RK * ( alist(maxerr) + blist(maxerr) )
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call qk21 ( f, a1, b1, area1, error1, resabs, defab1 )
        call qk21 ( f, a2, b2, area2, error2, resabs, defab2 )
    !
    !  Improve previous approximations to integral and error
    !  and test for accuracy.
    !
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)

        if ( defab1 == error1 .or. defab2 == error2 ) go to 15

        if ( abs ( rlist(maxerr) - area12) > 1.0E-05_RK * abs(area12) &
          .or. erro12 < 9.9E-01_RK * errmax ) go to 10

        if ( extrap ) then
          iroff2 = iroff2+1
        else
          iroff1 = iroff1+1
        end if

    10  continue

        if ( last > 10 .and. erro12 > errmax ) then
          iroff3 = iroff3+1
        end if

    15  continue

        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max ( epsabs, epsrel*abs(area) )
    !
    !  Test for roundoff error and eventually set error flag.
    !
        if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
          ier = 2
        end if

        if ( iroff2 >= 5 ) then
          ierro = 3
        end if
    !
    !  Set error flag in the case that the number of subintervals
    !  equals limit.
    !
        if ( last == limit ) then
          ier = 1
        end if
    !
    !  Set error flag in the case of bad integrand behavior
    !  at a point of the integration range.
    !
        if ( max ( abs(a1),abs(b2)) <= (1._RK + 1.0E+03_RK * epsilon ( a1 ) )* &
          (abs(a2)+1.0E+03_RK* tiny ( a2 ) ) ) then
          ier = 4
        end if
    !
    !  Append the newly-created intervals to the list.
    !
        if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
        else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
        end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with nrmax-th largest error estimate (to be bisected next).
    !
        call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

        if ( errsum <= errbnd ) go to 115

        if ( ier /= 0 ) then
          exit
        end if

        if ( last == 2 ) go to 80
        if ( noext ) go to 90

        erlarg = erlarg-erlast

        if ( abs(b1-a1) > small ) then
          erlarg = erlarg+erro12
        end if
    !
    !  Test whether the interval to be bisected next is the
    !  smallest interval.
    !
        if ( .not. extrap ) then
          if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
          extrap = .true.
          nrmax = 2
        end if

    !40  continue
    !
    !  The smallest interval has the largest error.
    !  Before bisecting decrease the sum of the errors over the
    !  larger intervals (erlarg) and perform extrapolation.
    !
        if ( ierro /= 3 .and. erlarg > ertest ) then

          id = nrmax
          jupbnd = last

          if ( last > (2+limit/2) ) then
            jupbnd = limit+3-last
          end if

          do k = id, jupbnd
            maxerr = iord(nrmax)
            errmax = elist(maxerr)
            if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
              go to 90
            end if
            nrmax = nrmax+1
          end do

        end if
    !
    !  Perform extrapolation.
    !
    !60  continue

        numrl2 = numrl2+1
        rlist2(numrl2) = area
        call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
        ktmin = ktmin+1

        if ( ktmin > 5 .and. abserr < 1.0E-03_RK * errsum ) then
          ier = 5
        end if

        if ( abseps < abserr ) then

          ktmin = 0
          abserr = abseps
          result = reseps
          correc = erlarg
          ertest = max ( epsabs,epsrel*abs(reseps))

          if ( abserr <= ertest ) then
            exit
          end if

        end if
    !
    !  Prepare bisection of the smallest interval.
    !
        if ( numrl2 == 1 ) then
          noext = .true.
        end if

        if ( ier == 5 ) then
          exit
        end if

        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small * 5.0E-01_RK
        erlarg = errsum
        go to 90

    80  continue

        small = abs ( b - a ) * 3.75E-01_RK
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area

    90  continue

      end do
    !
    !  Set final result and error estimate.
    !
      if ( abserr == huge ( abserr ) ) then
        go to 115
      end if

      if ( ier + ierro == 0 ) then
        go to 110
      end if

      if ( ierro == 3 ) then
        abserr = abserr + correc
      end if

      if ( ier == 0 ) then
        ier = 3
      end if

      if ( result /= 0._RK .and. area /= 0._RK ) then
        go to 105
      end if

      if ( abserr > errsum ) go to 115
      if ( area == 0._RK ) go to 130
      go to 110

    105 continue

      if ( abserr/abs(result) > errsum/abs(area) ) go to 115
    !
    !  Test on divergence.
    !
    110 continue

      if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <= defabs*1.0E-02_RK ) go to 130

      if ( 1.0E-02_RK > (result/area) .or. (result/area) > 1.0E+02_RK .or. errsum > abs(area) ) then
        ier = 6
      end if

      go to 130
    !
    !  Compute global integral sum.
    !
    115 continue

      result = sum ( rlist(1:last) )

      abserr = errsum

    130 continue
     
      if ( 2 < ier ) then
        ier = ier - 1
      end if

    140 continue

      neval = 42*last-21

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qawc ( f, a, b, c, epsabs, epsrel, result, abserr, neval, ier )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAWC computes a Cauchy principal value.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a Cauchy principal
    !    value 
    !      I = integral of F*W over (A,B),
    !    with
    !      W(X) = 1 / (X-C),
    !    with C distinct from A and B, hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) C, a parameter in the weight function, which must
    !    not be equal to A or B.
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !            ier    - integer(IK)
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the requested
    !                             accuracy has been achieved.
    !                     ier > 0 abnormal termination of the routine
    !                             the estimates for integral and error are
    !                             less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                     ier = 1 maximum number of subdivisions allowed
    !                             has been achieved. one can allow more sub-
    !                             divisions by increasing the data value of
    !                             limit in qawc (and taking the according
    !                             dimension adjustments into account).
    !                             however, if this yields no improvement it
    !                             is advised to analyze the integrand in
    !                             order to determine the integration
    !                             difficulties. if the position of a local
    !                             difficulty can be determined (e.g.
    !                             singularity, discontinuity within the
    !                             interval one will probably gain from
    !                             splitting up the interval at this point
    !                             and calling appropriate integrators on the
    !                             subranges.
    !                         = 2 the occurrence of roundoff error is detec-
    !                             ted, which prevents the requested
    !                             tolerance from being achieved.
    !                         = 3 extremely bad integrand behavior occurs
    !                             at some points of the integration
    !                             interval.
    !                         = 6 the input is invalid, because
    !                             c = a or c = b or
    !                             epsabs < 0 and epsrel < 0,
    !                             result, abserr, neval are set to zero.
    !
    !  Local parameters:
    !
    !    LIMIT is the maximum number of subintervals allowed in the
    !    subdivision process of qawce. take care that limit >= 1.
    !
      implicit none

      integer(IK), parameter :: limit = 500

      real(RK) a
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) b
      real(RK) blist(limit)
      real(RK) elist(limit)
      real(RK) c
      real(RK) epsabs
      real(RK) epsrel
      real(RK), external :: f
      integer(IK) ier
      integer(IK) iord(limit)
      integer(IK) last
      integer(IK) neval
      real(RK) result
      real(RK) rlist(limit)

      call qawce ( f, a, b, c, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iord, last )

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qawce ( f, a, b, c, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iord, last )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAWCE computes a Cauchy principal value.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a Cauchy principal
    !    value   
    !      I = integral of F*W over (A,B),
    !    with
    !      W(X) = 1 / ( X - C ),
    !    with C distinct from A and B, hopefully satisfying
    !      | I - RESULT | <= max ( EPSABS, EPSREL * |I| ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) C, a parameter in the weight function, which cannot be
    !    equal to A or B.
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Input, integer(IK) LIMIT, the upper bound on the number of subintervals that
    !    will be used in the partition of [A,B].  LIMIT is typically 500.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !            ier    - integer(IK)
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the requested
    !                             accuracy has been achieved.
    !                     ier > 0 abnormal termination of the routine
    !                             the estimates for integral and error are
    !                             less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                     ier = 1 maximum number of subdivisions allowed
    !                             has been achieved. one can allow more sub-
    !                             divisions by increasing the value of
    !                             limit. however, if this yields no
    !                             improvement it is advised to analyze the
    !                             integrand, in order to determine the
    !                             integration difficulties.  if the position
    !                             of a local difficulty can be determined
    !                             (e.g. singularity, discontinuity within
    !                             the interval) one will probably gain
    !                             from splitting up the interval at this
    !                             point and calling appropriate integrators
    !                             on the subranges.
    !                         = 2 the occurrence of roundoff error is detec-
    !                             ted, which prevents the requested
    !                             tolerance from being achieved.
    !                         = 3 extremely bad integrand behavior occurs
    !                             at some interior points of the integration
    !                             interval.
    !                         = 6 the input is invalid, because
    !                             c = a or c = b or
    !                             epsabs < 0 and epsrel < 0,
    !                             or limit < 1.
    !                             result, abserr, neval, rlist(1), elist(1),
    !                             iord(1) and last are set to zero.
    !                             alist(1) and blist(1) are set to a and b
    !                             respectively.
    !
    !    Workspace, real(RK) ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
    !    through LAST the left and right ends of the partition subintervals.
    !
    !    Workspace, real(RK) RLIST(LIMIT), contains in entries 1 through LAST
    !    the integral approximations on the subintervals.
    !
    !    Workspace, real(RK) ELIST(LIMIT), contains in entries 1 through LAST
    !    the absolute error estimates on the subintervals.
    !
    !            iord    - integer(IK)
    !                      vector of dimension at least limit, the first k
    !                      elements of which are pointers to the error
    !                      estimates over the subintervals, so that
    !                      elist(iord(1)), ...,  elist(iord(k)) with
    !                      k = last if last <= (limit/2+2), and
    !                      k = limit+1-last otherwise, form a decreasing
    !                      sequence.
    !
    !            last    - integer(IK)
    !                      number of subintervals actually produced in
    !                      the subdivision process
    !
    !  Local parameters:
    !
    !           alist     - list of left end points of all subintervals
    !                       considered up to now
    !           blist     - list of right end points of all subintervals
    !                       considered up to now
    !           rlist(i)  - approximation to the integral over
    !                       (alist(i),blist(i))
    !           elist(i)  - error estimate applying to rlist(i)
    !           maxerr    - pointer to the interval with largest error
    !                       estimate
    !           errmax    - elist(maxerr)
    !           area      - sum of the integrals over the subintervals
    !           errsum    - sum of the errors over the subintervals
    !           errbnd    - requested accuracy max(epsabs,epsrel*
    !                       abs(result))
    !           *****1    - variable for the left subinterval
    !           *****2    - variable for the right subinterval
    !           last      - index for subdivision
    !
      implicit none

      integer(IK) limit

      real(RK) a
      real(RK) aa
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) area
      real(RK) area1
      real(RK) area12
      real(RK) area2
      real(RK) a1
      real(RK) a2
      real(RK) b
      real(RK) bb
      real(RK) blist(limit)
      real(RK) b1
      real(RK) b2
      real(RK) c
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK) errbnd
      real(RK) errmax
      real(RK) error1
      real(RK) error2
      real(RK) erro12
      real(RK) errsum
      real(RK), external :: f
      integer(IK) ier
      integer(IK) iord(limit)
      integer(IK) iroff1
      integer(IK) iroff2
      integer(IK) krule
      integer(IK) last
      integer(IK) maxerr
      integer(IK) nev
      integer(IK) neval
      integer(IK) nrmax
      real(RK) result
      real(RK) rlist(limit)
    !
    !  Test on validity of parameters.
    !
      ier = 0
      neval = 0
      last = 0
      alist(1) = a
      blist(1) = b
      rlist(1) = 0._RK
      elist(1) = 0._RK
      iord(1) = 0
      result = 0._RK
      abserr = 0._RK

      if ( c == a  ) then
        ier = 6
        return
      else if ( c == b ) then
        ier = 6
        return
      else if ( epsabs < 0._RK .and. epsrel < 0._RK ) then
        ier = 6
        return
      end if
    !
    !  First approximation to the integral.
    !
      if ( a <= b ) then
        aa = a
        bb = b
      else
        aa = b
        bb = a
      end if

      krule = 1
      call qc25c ( f, aa, bb, c, result, abserr, krule, neval )
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      alist(1) = a
      blist(1) = b
    !
    !  Test on accuracy.
    !
      errbnd = max ( epsabs, epsrel * abs(result) )

      if ( limit == 1 ) then
        ier = 1
        go to 70
      end if

      if ( abserr < min ( 1.0E-02_RK * abs(result), errbnd)  ) then
        go to 70
      end if
    !
    !  Initialization
    !
      alist(1) = aa
      blist(1) = bb
      rlist(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0

      do last = 2, limit
    !
    !  Bisect the subinterval with nrmax-th largest error estimate.
    !
        a1 = alist(maxerr)
        b1 = 5.0E-01_RK*(alist(maxerr)+blist(maxerr))
        b2 = blist(maxerr)

        if ( c <= b1 .and. a1 < c ) then
          b1 = 5.0E-01_RK*(c+b2)
        end if

        if ( b1 < c .and. c < b2 ) then
          b1 = 5.0E-01_RK * ( a1 + c )
        end if

        a2 = b1
        krule = 2

        call qc25c ( f, a1, b1, c, area1, error1, krule, nev )
        neval = neval+nev

        call qc25c ( f, a2, b2, c, area2, error2, krule, nev )
        neval = neval+nev
    !
    !  Improve previous approximations to integral and error
    !  and test for accuracy.
    !
        area12 = area1 + area2
        erro12 = error1 + error2
        errsum = errsum + erro12 - errmax
        area = area + area12 - rlist(maxerr)

        if ( abs ( rlist(maxerr)-area12) < 1.0E-05_RK * abs(area12) .and. erro12 >= 9.9E-01_RK * errmax .and. krule == 0 ) iroff1 = iroff1+1

        if ( last > 10 .and. erro12 > errmax .and. krule == 0 ) then
          iroff2 = iroff2+1
        end if

        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max ( epsabs, epsrel * abs(area) )

        if ( errsum > errbnd ) then
    !
    !  Test for roundoff error and eventually set error flag.
    !
          if ( iroff1 >= 6 .and. iroff2 > 20 ) then
            ier = 2
          end if
    !
    !  Set error flag in the case that number of interval
    !  bisections exceeds limit.
    !
          if ( last == limit ) then
            ier = 1
          end if
    !
    !  Set error flag in the case of bad integrand behavior at
    !  a point of the integration range.
    !
          if ( max ( abs(a1), abs(b2) ) <= ( 1._RK + 1.0E+03_RK * epsilon ( a1 ) ) &
            *( abs(a2)+1.0E+03_RK * tiny ( a2 ) )) then
            ier = 3
          end if

        end if
    !
    !  Append the newly-created intervals to the list.
    !
        if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
        else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
        end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with NRMAX-th largest error estimate (to be bisected next).
    !
        call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

        if ( ier /= 0 .or. errsum <= errbnd ) then
          exit
        end if

      end do
    !
    !  Compute final result.
    !
      result = sum ( rlist(1:last) )

      abserr = errsum

    70 continue 

      if ( aa == b ) then
        result = - result
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qawf ( f, a, omega, integr, epsabs, result, abserr, neval, ier )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAWF computes Fourier integrals over the interval [ A, +Infinity ).
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral  
    ! 
    !      I = integral of F*COS(OMEGA*X) 
    !    or 
    !      I = integral of F*SIN(OMEGA*X) 
    !
    !    over the interval [A,+Infinity), hopefully satisfying
    !
    !      || I - RESULT || <= EPSABS.
    !
    !    If OMEGA = 0 and INTEGR = 1, the integral is calculated by means 
    !    of QAGI, and IER has the meaning as described in the comments of QAGI.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, the lower limit of integration.
    !
    !    Input, real(RK) OMEGA, the parameter in the weight function.
    !
    !    Input, integer(IK) INTEGR, indicates which weight functions is used
    !    = 1, w(x) = cos(omega*x)
    !    = 2, w(x) = sin(omega*x)
    !
    !    Input, real(RK) EPSABS, the absolute accuracy requested.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !            ier    - integer(IK)
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the
    !                             requested accuracy has been achieved.
    !                     ier > 0 abnormal termination of the routine.
    !                             the estimates for integral and error are
    !                             less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                    if omega /= 0
    !                     ier = 6 the input is invalid because
    !                             (integr /= 1 and integr /= 2) or
    !                              epsabs <= 0
    !                              result, abserr, neval, lst are set to
    !                              zero.
    !                         = 7 abnormal termination of the computation
    !                             of one or more subintegrals
    !                         = 8 maximum number of cycles allowed
    !                             has been achieved, i.e. of subintervals
    !                             (a+(k-1)c,a+kc) where
    !                             c = (2*int(abs(omega))+1)*pi/abs(omega),
    !                             for k = 1, 2, ...
    !                         = 9 the extrapolation table constructed for
    !                             convergence acceleration of the series
    !                             formed by the integral contributions
    !                             over the cycles, does not converge to
    !                             within the requested accuracy.
    !
    !  Local parameters:
    !
    !    integer(IK) LIMLST, gives an upper bound on the number of cycles, LIMLST >= 3.
    !    if limlst < 3, the routine will end with ier = 6.
    !
    !    integer(IK) MAXP1, an upper bound on the number of Chebyshev moments which 
    !    can be stored, i.e. for the intervals of lengths abs(b-a)*2**(-l), 
    !    l = 0,1, ..., maxp1-2, maxp1 >= 1.  if maxp1 < 1, the routine will end
    !    with ier = 6.
    !
      implicit none

      integer(IK), parameter :: limit = 500
      integer(IK), parameter :: limlst = 50
      integer(IK), parameter :: maxp1 = 21

      real(RK) a
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) blist(limit)
      real(RK) chebmo(maxp1,25)
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) erlst(limlst)
      real(RK), external :: f
      integer(IK) ier
      integer(IK) integr
      integer(IK) iord(limit)
      integer(IK) ierlst(limlst)
      integer(IK) lst
      integer(IK) neval
      integer(IK) nnlog(limit)
      real(RK) omega
      real(RK) result
      real(RK) rlist(limit)
      real(RK) rslst(limlst)

      ier = 6
      neval = 0
      result = 0._RK
      abserr = 0._RK

      if ( limlst < 3 .or. maxp1 < 1 ) then
        return
      end if

      call qawfe ( f, a, omega, integr, epsabs, limlst, limit, maxp1, &
        result, abserr, neval, ier, rslst, erlst, ierlst, lst, alist, blist, &
        rlist, elist, iord, nnlog, chebmo )

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qawfe ( f, a, omega, integr, epsabs, limlst, limit, maxp1, &
      result, abserr, neval, ier, rslst, erlst, ierlst, lst, alist, blist, &
      rlist, elist, iord, nnlog, chebmo )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAWFE computes Fourier integrals.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral   
    !      I = integral of F*COS(OMEGA*X) or F*SIN(OMEGA*X) over (A,+Infinity),
    !    hopefully satisfying
    !      || I - RESULT || <= EPSABS.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, the lower limit of integration.
    !
    !    Input, real(RK) OMEGA, the parameter in the weight function.
    !
    !    Input, integer(IK) INTEGR, indicates which weight function is used
    !    = 1      w(x) = cos(omega*x)
    !    = 2      w(x) = sin(omega*x)
    !
    !    Input, real(RK) EPSABS, the absolute accuracy requested.
    !
    !    Input, integer(IK) LIMLST, an upper bound on the number of cycles.
    !    LIMLST must be at least 1.  In fact, if LIMLST < 3, the routine 
    !    will end with IER= 6.
    !
    !    Input, integer(IK) LIMIT, an upper bound on the number of subintervals 
    !    allowed in the partition of each cycle, limit >= 1.
    !
    !            maxp1  - integer(IK)
    !                     gives an upper bound on the number of
    !                     Chebyshev moments which can be stored, i.e.
    !                     for the intervals of lengths abs(b-a)*2**(-l),
    !                     l=0,1, ..., maxp1-2, maxp1 >= 1
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !            ier    - ier = 0 normal and reliable termination of
    !                             the routine. it is assumed that the
    !                             requested accuracy has been achieved.
    !                     ier > 0 abnormal termination of the routine
    !                             the estimates for integral and error
    !                             are less reliable. it is assumed that
    !                             the requested accuracy has not been
    !                             achieved.
    !                    if omega /= 0
    !                     ier = 6 the input is invalid because
    !                             (integr /= 1 and integr /= 2) or
    !                              epsabs <= 0 or limlst < 3.
    !                              result, abserr, neval, lst are set
    !                              to zero.
    !                         = 7 bad integrand behavior occurs within
    !                             one or more of the cycles. location
    !                             and type of the difficulty involved
    !                             can be determined from the vector ierlst.
    !                             here lst is the number of cycles actually
    !                             needed (see below).
    !                             ierlst(k) = 1 the maximum number of
    !                                           subdivisions (= limit)
    !                                           has been achieved on the
    !                                           k th cycle.
    !                                       = 2 occurence of roundoff
    !                                           error is detected and
    !                                           prevents the tolerance
    !                                           imposed on the k th cycle
    !                                           from being acheived.
    !                                       = 3 extremely bad integrand
    !                                           behavior occurs at some
    !                                           points of the k th cycle.
    !                                       = 4 the integration procedure
    !                                           over the k th cycle does
    !                                           not converge (to within the
    !                                           required accuracy) due to
    !                                           roundoff in the
    !                                           extrapolation procedure
    !                                           invoked on this cycle. it
    !                                           is assumed that the result
    !                                           on this interval is the
    !                                           best which can be obtained.
    !                                       = 5 the integral over the k th
    !                                           cycle is probably divergent
    !                                           or slowly convergent. it
    !                                           must be noted that
    !                                           divergence can occur with
    !                                           any other value of
    !                                           ierlst(k).
    !                         = 8 maximum number of  cycles  allowed
    !                             has been achieved, i.e. of subintervals
    !                             (a+(k-1)c,a+kc) where
    !                             c = (2*int(abs(omega))+1)*pi/abs(omega),
    !                             for k = 1, 2, ..., lst.
    !                             one can allow more cycles by increasing
    !                             the value of limlst (and taking the
    !                             according dimension adjustments into
    !                             account).
    !                             examine the array iwork which contains
    !                             the error flags over the cycles, in order
    !                             to eventual look for local integration
    !                             difficulties.
    !                             if the position of a local difficulty can
    !                             be determined (e.g. singularity,
    !                             discontinuity within the interval)
    !                             one will probably gain from splitting
    !                             up the interval at this point and
    !                             calling appopriate integrators on the
    !                             subranges.
    !                         = 9 the extrapolation table constructed for
    !                             convergence acceleration of the series
    !                             formed by the integral contributions
    !                             over the cycles, does not converge to
    !                             within the required accuracy.
    !                             as in the case of ier = 8, it is advised
    !                             to examine the array iwork which contains
    !                             the error flags on the cycles.
    !                    if omega = 0 and integr = 1,
    !                    the integral is calculated by means of qagi
    !                    and ier = ierlst(1) (with meaning as described
    !                    for ierlst(k), k = 1).
    !
    !            rslst  - real
    !                     vector of dimension at least limlst
    !                     rslst(k) contains the integral contribution
    !                     over the interval (a+(k-1)c,a+kc) where
    !                     c = (2*int(abs(omega))+1)*pi/abs(omega),
    !                     k = 1, 2, ..., lst.
    !                     note that, if omega = 0, rslst(1) contains
    !                     the value of the integral over (a,infinity).
    !
    !            erlst  - real
    !                     vector of dimension at least limlst
    !                     erlst(k) contains the error estimate
    !                     corresponding with rslst(k).
    !
    !            ierlst - integer(IK)
    !                     vector of dimension at least limlst
    !                     ierlst(k) contains the error flag corresponding
    !                     with rslst(k). for the meaning of the local error
    !                     flags see description of output parameter ier.
    !
    !            lst    - integer(IK)
    !                     number of subintervals needed for the integration
    !                     if omega = 0 then lst is set to 1.
    !
    !            alist, blist, rlist, elist - real
    !                     vector of dimension at least limit,
    !
    !            iord, nnlog - integer(IK)
    !                     vector of dimension at least limit, providing
    !                     space for the quantities needed in the
    !                     subdivision process of each cycle
    !
    !            chebmo - real
    !                     array of dimension at least (maxp1,25),
    !                     providing space for the Chebyshev moments
    !                     needed within the cycles
    !
    !  Local parameters:
    !
    !           c1, c2    - end points of subinterval (of length
    !                       cycle)
    !           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
    !           psum      - vector of dimension at least (limexp+2)
    !                       (see routine qextr)
    !                       psum contains the part of the epsilon table
    !                       which is still needed for further computations.
    !                       each element of psum is a partial sum of
    !                       the series which should sum to the value of
    !                       the integral.
    !           errsum    - sum of error estimates over the
    !                       subintervals, calculated cumulatively
    !           epsa      - absolute tolerance requested over current
    !                       subinterval
    !           chebmo    - array containing the modified Chebyshev
    !                       moments (see also routine qc25o)
    !
      implicit none

      integer(IK) limit
      integer(IK) limlst
      integer(IK) maxp1

      real(RK) a
      real(RK) abseps
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) blist(limit)
      real(RK) chebmo(maxp1,25)
      real(RK) correc
      real(RK) cycle
      real(RK) c1
      real(RK) c2
      real(RK) dl
    ! real(RK) dla
      real(RK) drl
      real(RK) elist(limit)
      real(RK) ep
      real(RK) eps
      real(RK) epsa
      real(RK) epsabs
      real(RK) erlst(limlst)
      real(RK) errsum
      real(RK), external :: f
      real(RK) fact
      integer(IK) ier
      integer(IK) ierlst(limlst)
      integer(IK) integr
      integer(IK) iord(limit)
      integer(IK) ktmin
      integer(IK) l
      integer(IK) ll
      integer(IK) lst
      integer(IK) momcom
      integer(IK) nev
      integer(IK) neval
      integer(IK) nnlog(limit)
      integer(IK) nres
      integer(IK) numrl2
      real(RK) omega
      real(RK), parameter :: p = 0.9E+00_RK
      real(RK), parameter :: pi = 3.1415926535897932E+00_RK
      real(RK) p1
      real(RK) psum(52)
      real(RK) reseps
      real(RK) result
      real(RK) res3la(3)
      real(RK) rlist(limit)
      real(RK) rslst(limlst)
    !
    !  The dimension of  psum  is determined by the value of
    !  limexp in QEXTR (psum must be
    !  of dimension (limexp+2) at least).
    !
    !  Test on validity of parameters.
    !
      result = 0._RK
      abserr = 0._RK
      neval = 0
      lst = 0
      ier = 0

      if ( (integr /= 1 .and. integr /= 2 ) .or. epsabs <= 0._RK .or. limlst < 3 ) then
        ier = 6
        return
      end if

      if ( omega == 0._RK ) then

        if ( integr == 1 ) then
          call qagi ( f, 0._RK, 1, epsabs, 0._RK, result, abserr, neval, ier )
        else
          result = 0._RK
          abserr = 0._RK
          neval = 0
          ier = 0
        end if

        rslst(1) = result
        erlst(1) = abserr
        ierlst(1) = ier
        lst = 1

        return
      end if
    !
    !  Initializations.
    !
      l = int ( abs ( omega ) )
      dl = 2 * l + 1
      cycle = dl * pi / abs ( omega )
      ier = 0
      ktmin = 0
      neval = 0
      numrl2 = 0
      nres = 0
      c1 = a
      c2 = cycle+a
      p1 = 1._RK - p
      eps = epsabs

      if ( epsabs > tiny ( epsabs ) / p1 ) then
        eps = epsabs * p1
      end if

      ep = eps
      fact = 1._RK
      correc = 0._RK
      abserr = 0._RK
      errsum = 0._RK

      do lst = 1, limlst
    !
    !  Integrate over current subinterval.
    !
    !   dla = lst
        epsa = eps * fact

        call qfour ( f, c1, c2, omega, integr, epsa, 0._RK, limit, lst, maxp1, &
          rslst(lst), erlst(lst), nev, ierlst(lst), alist, blist, rlist, elist, &
          iord, nnlog, momcom, chebmo )

        neval = neval + nev
        fact = fact * p
        errsum = errsum + erlst(lst)
        drl = 5.0E+01 * abs(rslst(lst))
    !
    !  Test on accuracy with partial sum.
    !
        if ((errsum+drl) <= epsabs.and.lst >= 6) then
          go to 80
        end if

        correc = max ( correc,erlst(lst))

        if ( ierlst(lst) /= 0 ) then
          eps = max ( ep,correc*p1)
          ier = 7
        end if

        if ( ier == 7 .and. (errsum+drl) <= correc*1.0E+01_RK .and. lst > 5) go to 80

        numrl2 = numrl2+1

        if ( lst <= 1 ) then
          psum(1) = rslst(1)
          go to 40
        end if

        psum(numrl2) = psum(ll) + rslst(lst)

        if ( lst == 2 ) then
          go to 40
        end if
    !
    !  Test on maximum number of subintervals
    !
        if ( lst == limlst ) then
          ier = 8
        end if
    !
    !  Perform new extrapolation
    !
        call qextr ( numrl2, psum, reseps, abseps, res3la, nres )
    !
    !  Test whether extrapolated result is influenced by roundoff
    !
        ktmin = ktmin + 1

        if ( ktmin >= 15 .and. abserr <= 1.0E-03_RK * (errsum+drl) ) then
          ier = 9
        end  if

        if ( abseps <= abserr .or. lst == 3 ) then

          abserr = abseps
          result = reseps
          ktmin = 0
    !
    !  If IER is not 0, check whether direct result (partial
    !  sum) or extrapolated result yields the best integral
    !  approximation
    !
          if ( ( abserr + 1.0E+01_RK * correc ) <= epsabs ) then
            exit
          end if

          if ( abserr <= epsabs .and. 1.0E+01_RK * correc >= epsabs ) then
            exit
          end if

        end if

        if ( ier /= 0 .and. ier /= 7 ) then
          exit
        end if

    40  continue

        ll = numrl2
        c1 = c2
        c2 = c2+cycle

      end do
    !
    !  Set final result and error estimate.
    !
    !60 continue

      abserr = abserr + 1.0E+01_RK * correc

      if ( ier == 0 ) then
        return
      end if

      if ( result /= 0._RK .and. psum(numrl2) /= 0._RK) go to 70

      if ( abserr > errsum ) then
        go to 80
      end if

      if ( psum(numrl2) == 0._RK ) then
        return
      end if

    70 continue

      if ( abserr / abs(result) <= (errsum+drl)/abs(psum(numrl2)) ) then

        if ( ier >= 1 .and. ier /= 7 ) then
          abserr = abserr + drl
        end if

        return

      end if

    80 continue

      result = psum(numrl2)
      abserr = errsum + drl

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qawo ( f, a, b, omega, integr, epsabs, epsrel, result, abserr, neval, ier )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAWO computes the integrals of oscillatory integrands.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a given
    !    definite integral
    !      I = Integral ( A <= X <= B ) F(X) * cos ( OMEGA * X ) dx
    !    or 
    !      I = Integral ( A <= X <= B ) F(X) * sin ( OMEGA * X ) dx
    !    hopefully satisfying following claim for accuracy
    !      | I - RESULT | <= max ( epsabs, epsrel * |I| ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) OMEGA, the parameter in the weight function.
    !
    !    Input, integer(IK) INTEGR, specifies the weight function:
    !    1, W(X) = cos ( OMEGA * X )
    !    2, W(X) = sin ( OMEGA * X )
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !            ier    - integer(IK)
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the
    !                             requested accuracy has been achieved.
    !                   - ier > 0 abnormal termination of the routine.
    !                             the estimates for integral and error are
    !                             less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                     ier = 1 maximum number of subdivisions allowed
    !                             (= leniw/2) has been achieved. one can
    !                             allow more subdivisions by increasing the
    !                             value of leniw (and taking the according
    !                             dimension adjustments into account).
    !                             however, if this yields no improvement it
    !                             is advised to analyze the integrand in
    !                             order to determine the integration
    !                             difficulties. if the position of a local
    !                             difficulty can be determined (e.g.
    !                             singularity, discontinuity within the
    !                             interval) one will probably gain from
    !                             splitting up the interval at this point
    !                             and calling the integrator on the
    !                             subranges. if possible, an appropriate
    !                             special-purpose integrator should
    !                             be used which is designed for handling
    !                             the type of difficulty involved.
    !                         = 2 the occurrence of roundoff error is
    !                             detected, which prevents the requested
    !                             tolerance from being achieved.
    !                             the error may be under-estimated.
    !                         = 3 extremely bad integrand behavior occurs
    !                             at some interior points of the integration
    !                             interval.
    !                         = 4 the algorithm does not converge. roundoff
    !                             error is detected in the extrapolation
    !                             table. it is presumed that the requested
    !                             tolerance cannot be achieved due to
    !                             roundoff in the extrapolation table,
    !                             and that the returned result is the best
    !                             which can be obtained.
    !                         = 5 the integral is probably divergent, or
    !                             slowly convergent. it must be noted that
    !                             divergence can occur with any other value
    !                             of ier.
    !                         = 6 the input is invalid, because
    !                             epsabs < 0 and epsrel < 0,
    !                             result, abserr, neval are set to zero.
    !
    !  Local parameters:
    !
    !    limit is the maximum number of subintervals allowed in the
    !    subdivision process of QFOUR. take care that limit >= 1.
    !
    !    maxp1 gives an upper bound on the number of Chebyshev moments
    !    which can be stored, i.e. for the intervals of lengths
    !    abs(b-a)*2**(-l), l = 0, 1, ... , maxp1-2. take care that
    !    maxp1 >= 1.

      implicit none

      integer(IK), parameter :: limit = 500
      integer(IK), parameter :: maxp1 = 21

      real(RK) a
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) b
      real(RK) blist(limit)
      real(RK) chebmo(maxp1,25)
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK), external :: f
      integer(IK) ier
      integer(IK) integr
      integer(IK) iord(limit)
      integer(IK) momcom
      integer(IK) neval
      integer(IK) nnlog(limit)
      real(RK) omega
      real(RK) result
      real(RK) rlist(limit)

      call qfour ( f, a, b, omega, integr, epsabs, epsrel, limit, 1, maxp1, &
        result, abserr, neval, ier, alist, blist, rlist, elist, iord, nnlog, &
        momcom, chebmo )

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qaws ( f, a, b, alfa, beta, integr, epsabs, epsrel, result, &
      abserr, neval, ier )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAWS estimates integrals with algebraico-logarithmic endpoint singularities.
    !
    !  Discussion:
    !
    !    This routine calculates an approximation RESULT to a given
    !    definite integral   
    !      I = integral of f*w over (a,b) 
    !    where w shows a singular behavior at the end points, see parameter
    !    integr, hopefully satisfying following claim for accuracy
    !      abs(i-result) <= max(epsabs,epsrel*abs(i)).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) ALFA, BETA, parameters used in the weight function.
    !    ALFA and BETA should be greater than -1.
    !
    !    Input, integer(IK) INTEGR, indicates which weight function is to be used
    !    = 1  (x-a)**alfa*(b-x)**beta
    !    = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
    !    = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
    !    = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !            ier    - integer(IK)
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the requested
    !                             accuracy has been achieved.
    !                     ier > 0 abnormal termination of the routine
    !                             the estimates for the integral and error
    !                             are less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                     ier = 1 maximum number of subdivisions allowed
    !                             has been achieved. one can allow more
    !                             subdivisions by increasing the data value
    !                             of limit in qaws (and taking the according
    !                             dimension adjustments into account).
    !                             however, if this yields no improvement it
    !                             is advised to analyze the integrand, in
    !                             order to determine the integration
    !                             difficulties which prevent the requested
    !                             tolerance from being achieved. in case of
    !                             a jump discontinuity or a local
    !                             singularity of algebraico-logarithmic type
    !                             at one or more interior points of the
    !                             integration range, one should proceed by
    !                             splitting up the interval at these points
    !                             and calling the integrator on the
    !                             subranges.
    !                         = 2 the occurrence of roundoff error is
    !                             detected, which prevents the requested
    !                             tolerance from being achieved.
    !                         = 3 extremely bad integrand behavior occurs
    !                             at some points of the integration
    !                             interval.
    !                         = 6 the input is invalid, because
    !                             b <= a or alfa <= (-1) or beta <= (-1) or
    !                             integr < 1 or integr > 4 or
    !                             epsabs < 0 and epsrel < 0,
    !                             result, abserr, neval are set to zero.
    !
    !  Local parameters:
    !
    !    LIMIT is the maximum number of subintervals allowed in the
    !    subdivision process of qawse. take care that limit >= 2.
    !
      implicit none

      integer(IK), parameter :: limit = 500

      real(RK) a
      real(RK) abserr
      real(RK) alfa
      real(RK) alist(limit)
      real(RK) b
      real(RK) blist(limit)
      real(RK) beta
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK), external :: f
      integer(IK) ier
      integer(IK) integr
      integer(IK) iord(limit)
      integer(IK) last
      integer(IK) neval
      real(RK) result
      real(RK) rlist(limit)

      call qawse ( f, a, b, alfa, beta, integr, epsabs, epsrel, limit, result, &
        abserr, neval, ier, alist, blist, rlist, elist, iord, last )

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qawse ( f, a, b, alfa, beta, integr, epsabs, epsrel, limit, &
      result, abserr, neval, ier, alist, blist, rlist, elist, iord, last )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QAWSE estimates integrals with algebraico-logarithmic endpoint singularities.
    !
    !  Discussion:
    !
    !    This routine calculates an approximation RESULT to an integral
    !      I = integral of F(X) * W(X) over (a,b), 
    !    where W(X) shows a singular behavior at the endpoints, hopefully 
    !    satisfying:
    !      | I - RESULT | <= max ( epsabs, epsrel * |I| ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) ALFA, BETA, parameters used in the weight function.
    !    ALFA and BETA should be greater than -1.
    !
    !    Input, integer(IK) INTEGR, indicates which weight function is used:
    !    = 1  (x-a)**alfa*(b-x)**beta
    !    = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
    !    = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
    !    = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Input, integer(IK) LIMIT, an upper bound on the number of subintervals
    !    in the partition of (A,B), LIMIT >= 2.  If LIMIT < 2, the routine 
    !     will end with IER = 6.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !            ier    - integer(IK)
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the requested
    !                             accuracy has been achieved.
    !                     ier > 0 abnormal termination of the routine
    !                             the estimates for the integral and error
    !                             are less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                         = 1 maximum number of subdivisions allowed
    !                             has been achieved. one can allow more
    !                             subdivisions by increasing the value of
    !                             limit. however, if this yields no
    !                             improvement it is advised to analyze the
    !                             integrand, in order to determine the
    !                             integration difficulties which prevent
    !                             the requested tolerance from being
    !                             achieved. in case of a jump discontinuity
    !                             or a local singularity of algebraico-
    !                             logarithmic type at one or more interior
    !                             points of the integration range, one
    !                             should proceed by splitting up the
    !                             interval at these points and calling the
    !                             integrator on the subranges.
    !                         = 2 the occurrence of roundoff error is
    !                             detected, which prevents the requested
    !                             tolerance from being achieved.
    !                         = 3 extremely bad integrand behavior occurs
    !                             at some points of the integration
    !                             interval.
    !                         = 6 the input is invalid, because
    !                             b <= a or alfa <= (-1) or beta <= (-1) or
    !                             integr < 1 or integr > 4, or
    !                             epsabs < 0 and epsrel < 0,
    !                             or limit < 2.
    !                             result, abserr, neval, rlist(1), elist(1),
    !                             iord(1) and last are set to zero.
    !                             alist(1) and blist(1) are set to a and b
    !                             respectively.
    !
    !    Workspace, real(RK) ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
    !    through LAST the left and right ends of the partition subintervals.
    !
    !    Workspace, real(RK) RLIST(LIMIT), contains in entries 1 through LAST
    !    the integral approximations on the subintervals.
    !
    !    Workspace, real(RK) ELIST(LIMIT), contains in entries 1 through LAST
    !    the absolute error estimates on the subintervals.
    !
    !            iord   - integer(IK)
    !                     vector of dimension at least limit, the first k
    !                     elements of which are pointers to the error
    !                     estimates over the subintervals, so that
    !                     elist(iord(1)), ..., elist(iord(k)) with k = last
    !                     if last <= (limit/2+2), and k = limit+1-last
    !                     otherwise, form a decreasing sequence.
    !
    !    Output, integer(IK) LAST, the number of subintervals actually produced in 
    !    the subdivision process.
    !
    !  Local parameters:
    !
    !           alist     - list of left end points of all subintervals
    !                       considered up to now
    !           blist     - list of right end points of all subintervals
    !                       considered up to now
    !           rlist(i)  - approximation to the integral over
    !                       (alist(i),blist(i))
    !           elist(i)  - error estimate applying to rlist(i)
    !           maxerr    - pointer to the interval with largest error
    !                       estimate
    !           errmax    - elist(maxerr)
    !           area      - sum of the integrals over the subintervals
    !           errsum    - sum of the errors over the subintervals
    !           errbnd    - requested accuracy max(epsabs,epsrel*
    !                       abs(result))
    !           *****1    - variable for the left subinterval
    !           *****2    - variable for the right subinterval
    !           last      - index for subdivision
    !
      implicit none

      integer(IK) limit

      real(RK) a
      real(RK) abserr
      real(RK) alfa
      real(RK) alist(limit)
      real(RK) area
      real(RK) area1
      real(RK) area12
      real(RK) area2
      real(RK) a1
      real(RK) a2
      real(RK) b
      real(RK) beta
      real(RK) blist(limit)
      real(RK) b1
      real(RK) b2
      real(RK) centre
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK) errbnd
      real(RK) errmax
      real(RK) error1
      real(RK) erro12
      real(RK) error2
      real(RK) errsum
      real(RK), external :: f
      integer(IK) ier
      integer(IK) integr
      integer(IK) iord(limit)
      integer(IK) iroff1
      integer(IK) iroff2
      integer(IK) last
      integer(IK) maxerr
      integer(IK) nev
      integer(IK) neval
      integer(IK) nrmax
      real(RK) resas1
      real(RK) resas2
      real(RK) result
      real(RK) rg(25)
      real(RK) rh(25)
      real(RK) ri(25)
      real(RK) rj(25)
      real(RK) rlist(limit)
    !
    !  Test on validity of parameters.
    !
      ier = 0
      neval = 0
      last = 0
      rlist(1) = 0._RK
      elist(1) = 0._RK
      iord(1) = 0
      result = 0._RK
      abserr = 0._RK

      if ( b <= a .or. &
        (epsabs < 0._RK .and. epsrel < 0._RK) .or. &
        alfa <= (-1._RK) .or. &
        beta <= (-1._RK) .or. &
        integr < 1 .or. &
        integr > 4 .or. &
        limit < 2 ) then
        ier = 6
        return
      end if
    !
    !  Compute the modified Chebyshev moments.
    !
      call qmomo ( alfa, beta, ri, rj, rg, rh, integr )
    !
    !  Integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
    !
      centre = 5.0E-01 * ( b + a )

      call qc25s ( f, a, b, a, centre, alfa, beta, ri, rj, rg, rh, area1, &
        error1, resas1, integr, nev )

      neval = nev

      call qc25s ( f, a, b, centre, b, alfa, beta, ri, rj, rg, rh, area2, &
        error2, resas2, integr, nev )

      last = 2
      neval = neval+nev
      result = area1+area2
      abserr = error1+error2
    !
    !  Test on accuracy.
    !
      errbnd = max ( epsabs, epsrel * abs ( result ) )
    !
    !  Initialization.
    !
      if ( error2 <= error1 ) then
        alist(1) = a
        alist(2) = centre
        blist(1) = centre
        blist(2) = b
        rlist(1) = area1
        rlist(2) = area2
        elist(1) = error1
        elist(2) = error2
      else
        alist(1) = centre
        alist(2) = a
        blist(1) = b
        blist(2) = centre
        rlist(1) = area2
        rlist(2) = area1
        elist(1) = error2
        elist(2) = error1
      end if

      iord(1) = 1
      iord(2) = 2

      if ( limit == 2 ) then
        ier = 1
        return
      end if

      if ( abserr <= errbnd ) then
        return
      end if

      errmax = elist(1)
      maxerr = 1
      nrmax = 1
      area = result
      errsum = abserr
      iroff1 = 0
      iroff2 = 0

      do last = 3, limit
    !
    !  Bisect the subinterval with largest error estimate.
    !
        a1 = alist(maxerr)
        b1 = 5.0E-01_RK * ( alist(maxerr) + blist(maxerr) )
        a2 = b1
        b2 = blist(maxerr)

        call qc25s ( f, a, b, a1, b1, alfa, beta, ri, rj, rg, rh, area1, &
          error1, resas1, integr, nev )

        neval = neval + nev

        call qc25s ( f, a, b, a2, b2, alfa, beta, ri, rj, rg, rh, area2, &
          error2, resas2, integr, nev )

        neval = neval + nev
    !
    !  Improve previous approximations integral and error and
    !  test for accuracy.
    !
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
    !
    !  Test for roundoff error.
    !
        if ( a /= a1 .and. b /= b2 ) then

          if ( resas1 /= error1 .and. resas2 /= error2 ) then

            if ( abs ( rlist(maxerr) - area12 ) < 1.0E-05_RK * abs ( area12 ) &
              .and.erro12 >= 9.9E-01_RK*errmax ) then
              iroff1 = iroff1 + 1
            end if

            if ( last > 10 .and. erro12 > errmax ) then
              iroff2 = iroff2 + 1
            end if

          end if

        end if

        rlist(maxerr) = area1
        rlist(last) = area2
    !
    !  Test on accuracy.
    !
        errbnd = max ( epsabs, epsrel * abs ( area ) )

        if ( errsum > errbnd ) then
    !
    !  Set error flag in the case that the number of interval
    !  bisections exceeds limit.
    !
          if ( last == limit ) then
            ier = 1
          end if
    !
    !  Set error flag in the case of roundoff error.
    !
          if ( iroff1 >= 6 .or. iroff2 >= 20 ) then
            ier = 2
         end if
    !
    !  Set error flag in the case of bad integrand behavior
    !  at interior points of integration range.
    !
          if ( max ( abs(a1),abs(b2)) <= (1._RK+1.0E+03_RK* epsilon ( a1 ) )*( abs(a2) + 1.0E+03_RK * tiny ( a2) )) then
            ier = 3
          end if

        end if
    !
    !  Append the newly-created intervals to the list.
    !
        if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
        else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
        end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with largest error estimate (to be bisected next).
    !
        call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

        if ( ier /= 0 .or. errsum <= errbnd ) then
          exit
        end if

      end do
    !
    !  Compute final result.
    !
      result = sum ( rlist(1:last) )

      abserr = errsum

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qc25c ( f, a, b, c, result, abserr, krul, neval )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QC25C returns integration rules for Cauchy Principal Value integrals.
    !
    !  Discussion:
    !
    !    This routine estimates 
    !      I = integral of F(X) * W(X) over (a,b) 
    !    with error estimate, where 
    !      w(x) = 1/(x-c)
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) C, the parameter in the weight function.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !    RESULT is computed by using a generalized Clenshaw-Curtis method if
    !    C lies within ten percent of the integration interval.  In the 
    !    other case the 15-point Kronrod rule obtained by optimal addition
    !    of abscissae to the 7-point Gauss rule, is applied.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !           krul   - integer(IK)
    !                    key which is decreased by 1 if the 15-point
    !                    Gauss-Kronrod scheme has been used
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !  Local parameters:
    !
    !           fval   - value of the function f at the points
    !                    cos(k*pi/24),  k = 0, ..., 24
    !           cheb12 - Chebyshev series expansion coefficients, for the
    !                    function f, of degree 12
    !           cheb24 - Chebyshev series expansion coefficients, for the
    !                    function f, of degree 24
    !           res12  - approximation to the integral corresponding to the
    !                    use of cheb12
    !           res24  - approximation to the integral corresponding to the
    !                    use of cheb24
    !           qwgtc  - external function subprogram defining the weight
    !                    function
    !           hlgth  - half-length of the interval
    !           centr  - mid point of the interval
    !
      implicit none

      real(RK) a
      real(RK) abserr
      real(RK) ak22
      real(RK) amom0
      real(RK) amom1
      real(RK) amom2
      real(RK) b
      real(RK) c
      real(RK) cc
      real(RK) centr
      real(RK) cheb12(13)
      real(RK) cheb24(25)
      real(RK), external :: f
      real(RK) fval(25)
      real(RK) hlgth
      integer(IK) i
      integer(IK) isym
      integer(IK) k 
      integer(IK) kp
      integer(IK) krul
      integer(IK) neval
      real(RK) p2
      real(RK) p3
      real(RK) p4
      !real(RK), external :: qwgtc
      real(RK) resabs
      real(RK) resasc
      real(RK) result
      real(RK) res12
      real(RK) res24
      real(RK) u
      real(RK), parameter, dimension ( 11 ) :: x = (/ &
        9.914448613738104E-01_RK, 9.659258262890683E-01_RK, &
        9.238795325112868E-01_RK, 8.660254037844386E-01_RK, &
        7.933533402912352E-01_RK, 7.071067811865475E-01_RK, &
        6.087614290087206E-01_RK, 5.000000000000000E-01_RK, &
        3.826834323650898E-01_RK, 2.588190451025208E-01_RK, &
        1.305261922200516E-01_RK /)
    !
    !  Check the position of C.
    !
      cc = ( 2._RK * c - b - a ) / ( b - a )
    !
    !  Apply the 15-point Gauss-Kronrod scheme.
    !
      if ( abs ( cc ) >= 1.1E+00_RK ) then
        krul = krul - 1
        call qk15w ( f, qwgtc, c, p2, p3, p4, kp, a, b, result, abserr, &
          resabs, resasc )
        neval = 15
        if ( resasc == abserr ) then
          krul = krul+1
        end if
        return
      end if
    !
    !  Use the generalized Clenshaw-Curtis method.
    !
      hlgth = 5.0E-01_RK * ( b - a )
      centr = 5.0E-01_RK * ( b + a )
      neval = 25
      fval(1) = 5.0E-01_RK * f(hlgth+centr)
      fval(13) = f(centr)
      fval(25) = 5.0E-01_RK * f(centr-hlgth)

      do i = 2, 12
        u = hlgth * x(i-1)
        isym = 26 - i
        fval(i) = f(u+centr)
        fval(isym) = f(centr-u)
      end do
    !
    !  Compute the Chebyshev series expansion.
    !
      call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  The modified Chebyshev moments are computed by forward
    !  recursion, using AMOM0 and AMOM1 as starting values.
    !
      amom0 = log ( abs ( ( 1._RK - cc ) / ( 1._RK + cc ) ) )
      amom1 = 2._RK + cc * amom0
      res12 = cheb12(1) * amom0 + cheb12(2) * amom1
      res24 = cheb24(1) * amom0 + cheb24(2) * amom1

      do k = 3, 13
        amom2 = 2._RK * cc * amom1 - amom0
        ak22 = ( k - 2 ) * ( k - 2 )
        if ( ( k / 2 ) * 2 == k ) then
          amom2 = amom2 - 4._RK / ( ak22 - 1._RK )
        end if
        res12 = res12 + cheb12(k) * amom2
        res24 = res24 + cheb24(k) * amom2
        amom0 = amom1
        amom1 = amom2
      end do

      do k = 14, 25
        amom2 = 2._RK * cc * amom1 - amom0
        ak22 = ( k - 2 ) * ( k - 2 )
        if ( ( k / 2 ) * 2 == k ) then
          amom2 = amom2 - 4._RK / ( ak22 - 1._RK )
        end if
        res24 = res24 + cheb24(k) * amom2
        amom0 = amom1
        amom1 = amom2
      end do

      result = res24
      abserr = abs ( res24 - res12 )

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qc25o ( f, a, b, omega, integr, nrmom, maxp1, ksave, result, &
      abserr, neval, resabs, resasc, momcom, chebmo )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QC25O returns integration rules for integrands with a COS or SIN factor.
    !
    !  Discussion:
    !
    !    This routine estimates the integral
    !      I = integral of f(x) * w(x) over (a,b)
    !    where
    !      w(x) = cos(omega*x)
    !    or 
    !      w(x) = sin(omega*x),
    !    and estimates
    !      J = integral ( A <= X <= B ) |F(X)| dx.
    !
    !    For small values of OMEGA or small intervals (a,b) the 15-point
    !    Gauss-Kronrod rule is used.  In all other cases a generalized
    !    Clenshaw-Curtis method is used, that is, a truncated Chebyshev 
    !    expansion of the function F is computed on (a,b), so that the 
    !    integrand can be written as a sum of terms of the form W(X)*T(K,X), 
    !    where T(K,X) is the Chebyshev polynomial of degree K.  The Chebyshev
    !    moments are computed with use of a linear recurrence relation.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) OMEGA, the parameter in the weight function.
    !
    !    Input, integer(IK) INTEGR, indicates which weight function is to be used
    !    = 1, w(x) = cos(omega*x)
    !    = 2, w(x) = sin(omega*x)
    !
    !    ?, integer(IK) NRMOM, the length of interval (a,b) is equal to the length
    !    of the original integration interval divided by
    !    2**nrmom (we suppose that the routine is used in an
    !    adaptive integration process, otherwise set
    !    nrmom = 0).  nrmom must be zero at the first call.
    !
    !           maxp1  - integer(IK)
    !                    gives an upper bound on the number of Chebyshev
    !                    moments which can be stored, i.e. for the intervals
    !                    of lengths abs(bb-aa)*2**(-l), l = 0,1,2, ...,
    !                    maxp1-2.
    !
    !           ksave  - integer(IK)
    !                    key which is one when the moments for the
    !                    current interval have been computed
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !           abserr - real
    !                    estimate of the modulus of the absolute
    !                    error, which should equal or exceed abs(i-result)
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !    Output, real(RK) RESABS, approximation to the integral J.
    !
    !    Output, real(RK) RESASC, approximation to the integral of abs(F-I/(B-A)).
    !
    !         on entry and return
    !           momcom - integer(IK)
    !                    for each interval length we need to compute
    !                    the Chebyshev moments. momcom counts the number
    !                    of intervals for which these moments have already
    !                    been computed. if nrmom < momcom or ksave = 1,
    !                    the Chebyshev moments for the interval (a,b)
    !                    have already been computed and stored, otherwise
    !                    we compute them and we increase momcom.
    !
    !           chebmo - real
    !                    array of dimension at least (maxp1,25) containing
    !                    the modified Chebyshev moments for the first momcom
    !                    interval lengths
    !
    !  Local parameters:
    !
    !    maxp1 gives an upper bound
    !           on the number of Chebyshev moments which can be
    !           computed, i.e. for the interval (bb-aa), ...,
    !           (bb-aa)/2**(maxp1-2).
    !           should this number be altered, the first dimension of
    !           chebmo needs to be adapted.
    !
    !    x contains the values cos(k*pi/24)
    !           k = 1, ...,11, to be used for the Chebyshev expansion of f
    !
    !           centr  - mid point of the integration interval
    !           hlgth  - half length of the integration interval
    !           fval   - value of the function f at the points
    !                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5
    !                    k = 0, ...,24
    !           cheb12 - coefficients of the Chebyshev series expansion
    !                    of degree 12, for the function f, in the
    !                    interval (a,b)
    !           cheb24 - coefficients of the Chebyshev series expansion
    !                    of degree 24, for the function f, in the
    !                    interval (a,b)
    !           resc12 - approximation to the integral of
    !                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
    !                    over (-1,+1), using the Chebyshev series
    !                    expansion of degree 12
    !           resc24 - approximation to the same integral, using the
    !                    Chebyshev series expansion of degree 24
    !           ress12 - the analogue of resc12 for the sine
    !           ress24 - the analogue of resc24 for the sine
    !
      implicit none

      integer(IK) maxp1

      real(RK) a
      real(RK) abserr
      real(RK) ac
      real(RK) an
      real(RK) an2
      real(RK) as
      real(RK) asap
      real(RK) ass
      real(RK) b
      real(RK) centr
      real(RK) chebmo(maxp1,25)
      real(RK) cheb12(13)
      real(RK) cheb24(25)
      real(RK) conc
      real(RK) cons
      real(RK) cospar
      real(RK) d(28)
      real(RK) d1(28)
      real(RK) d2(28)
      real(RK) d3(28)
      real(RK) estc
      real(RK) ests
      real(RK), external :: f
      real(RK) fval(25)
      real(RK) hlgth
      integer(IK) i
      integer(IK) integr
      integer(IK) isym
      integer(IK) j
      integer(IK) k
      integer(IK) ksave
      integer(IK) m
      integer(IK) momcom
      integer(IK) neval
      integer(IK), parameter :: nmac = 28
      integer(IK) noeq1
      integer(IK) noequ
      integer(IK) nrmom
      real(RK) omega
      real(RK) parint
      real(RK) par2
      real(RK) par22
      real(RK) p2
      real(RK) p3
      real(RK) p4
      !real(RK), external :: qwgto
      real(RK) resabs
      real(RK) resasc
      real(RK) resc12
      real(RK) resc24
      real(RK) ress12
      real(RK) ress24
      real(RK) result
      real(RK) sinpar
      real(RK) v(28)
      real(RK), dimension ( 11 ) :: x = (/ &
        9.914448613738104E-01_RK,     9.659258262890683E-01_RK, &
        9.238795325112868E-01_RK,     8.660254037844386E-01_RK, &
        7.933533402912352E-01_RK,     7.071067811865475E-01_RK, &
        6.087614290087206E-01_RK,     5.000000000000000E-01_RK, &
        3.826834323650898E-01_RK,     2.588190451025208E-01_RK, &
        1.305261922200516E-01_RK /)

      centr = 5.0E-01 * ( b + a )
      hlgth = 5.0E-01 * ( b - a )
      parint = omega * hlgth
    !
    !  Compute the integral using the 15-point Gauss-Kronrod
    !  formula if the value of the parameter in the integrand
    !  is small or if the length of the integration interval
    !  is less than (bb-aa)/2**(maxp1-2), where (aa,bb) is the
    !  original integration interval.
    !
      if ( abs ( parint ) <= 2._RK ) then

        call qk15w ( f, qwgto, omega, p2, p3, p4, integr, a, b, result, &
          abserr, resabs, resasc )

        neval = 15
        return

      end if
    !
    !  Compute the integral using the generalized clenshaw-curtis method.
    !
      conc = hlgth * cos(centr*omega)
      cons = hlgth * sin(centr*omega)
      resasc = huge ( resasc )
      neval = 25
    !
    !  Check whether the Chebyshev moments for this interval
    !  have already been computed.
    !
      if ( nrmom < momcom .or. ksave == 1 ) then
        go to 140
      end if
    !
    !  Compute a new set of Chebyshev moments.
    !
      m = momcom + 1
      par2 = parint * parint
      par22 = par2 + 2._RK
      sinpar = sin(parint)
      cospar = cos(parint)
    !
    !  Compute the Chebyshev moments with respect to cosine.
    !
      v(1) = 2._RK * sinpar / parint
      v(2) = (8._RK*cospar+(par2+par2-8._RK)*sinpar/ parint)/par2
      v(3) = (3.2E+01_RK*(par2-1.2E+01_RK)*cospar+(2._RK*((par2-8.0E+01_RK)*par2+1.92E+02_RK)*sinpar)/parint)/(par2*par2)
      ac = 8._RK*cospar
      as = 2.4E+01*parint*sinpar

      if ( abs ( parint ) > 2.4E+01_RK ) then
        go to 70
      end if
    !
    !  Compute the Chebyshev moments as the solutions of a boundary value 
    !  problem with one initial value (v(3)) and one end value computed
    !  using an asymptotic formula.
    !
      noequ = nmac-3
      noeq1 = noequ-1
      an = 6._RK

      do k = 1, noeq1
        an2 = an*an
        d(k) = -2._RK*(an2-4._RK) * (par22-an2-an2)
        d2(k) = (an-1._RK)*(an-2._RK) * par2
        d1(k) = (an+3._RK)*(an+4._RK) * par2
        v(k+3) = as-(an2-4._RK) * ac
        an = an+2._RK
      end do

      an2 = an*an
      d(noequ) = -2._RK*(an2-4._RK) * (par22-an2-an2)
      v(noequ+3) = as - ( an2 - 4._RK ) * ac
      v(4) = v(4) - 5.6E+01_RK * par2 * v(3)
      ass = parint * sinpar
      asap  = (((((2.10E+02_RK*par2-1._RK)*cospar-(1.05E+02_RK*par2-6.3E+01_RK)*ass)/an2-(1._RK-1.5E+01_RK*par2)*cospar &
            + 1.5E+01_RK*ass)/an2-cospar+3._RK*ass)/an2-cospar)/an2
      v(noequ+3) = v(noequ+3)-2._RK*asap*par2*(an-1._RK)* &
          (an-2._RK)
    !
    !  Solve the tridiagonal system by means of Gaussian
    !  elimination with partial pivoting.
    !
      d3(1:noequ) = 0._RK

      d2(noequ) = 0._RK

      do i = 1, noeq1

        if ( abs(d1(i)) > abs(d(i)) ) then
          an = d1(i)
          d1(i) = d(i)
          d(i) = an
          an = d2(i)
          d2(i) = d(i+1)
          d(i+1) = an
          d3(i) = d2(i+1)
          d2(i+1) = 0._RK
          an = v(i+4)
          v(i+4) = v(i+3)
          v(i+3) = an
        end if

        d(i+1) = d(i+1)-d2(i)*d1(i)/d(i)
        d2(i+1) = d2(i+1)-d3(i)*d1(i)/d(i)
        v(i+4) = v(i+4)-v(i+3)*d1(i)/d(i)

      end do

      v(noequ+3) = v(noequ+3) / d(noequ)
      v(noequ+2) = (v(noequ+2)-d2(noeq1)*v(noequ+3))/d(noeq1)

      do i = 2, noeq1
        k = noequ-i
        v(k+3) = (v(k+3)-d3(k)*v(k+5)-d2(k)*v(k+4))/d(k)
      end do

      go to 90
    !
    !  Compute the Chebyshev moments by means of forward recursion
    !
    70 continue

      an = 4._RK

      do i = 4, 13
        an2 = an*an
        v(i) = ((an2-4._RK)*(2._RK*(par22-an2-an2)*v(i-1)-ac) &
         +as-par2*(an+1._RK)*(an+2._RK)*v(i-2))/ &
         (par2*(an-1._RK)*(an-2._RK))
        an = an+2._RK
      end do

    90 continue

      do j = 1, 13
        chebmo(m,2*j-1) = v(j)
      end do
    !
    !  Compute the Chebyshev moments with respect to sine.
    !
      v(1) = 2._RK*(sinpar-parint*cospar)/par2
      v(2) = (1.8E+01_RK-4.8E+01_RK/par2)*sinpar/par2 &
         +(-2._RK+4.8E+01_RK/par2)*cospar/parint
      ac = -2.4E+01_RK*parint*cospar
      as = -8._RK*sinpar
      chebmo(m,2) = v(1)
      chebmo(m,4) = v(2)

      if ( abs(parint) <= 2.4E+01_RK ) then

        do k = 3, 12
          an = k
          chebmo(m,2*k) = -sinpar/(an*(2._RK*an-2._RK)) &
                     -2.5E-01_RK*parint*(v(k+1)/an-v(k)/(an-1._RK))
        end do
    !
    !  Compute the Chebyshev moments by means of forward recursion.
    !
      else

        an = 3._RK

        do i = 3, 12
          an2 = an*an
          v(i) = ((an2-4._RK)*(2._RK*(par22-an2-an2)*v(i-1)+as) &
           +ac-par2*(an+1._RK)*(an+2._RK)*v(i-2)) &
           /(par2*(an-1._RK)*(an-2._RK))
          an = an+2._RK
          chebmo(m,2*i) = v(i)
        end do

      end if

    140 continue

      if ( nrmom < momcom ) then
        m = nrmom + 1
      end if

      if ( momcom < maxp1 - 1 .and. nrmom >= momcom ) then
        momcom = momcom + 1
      end if
    !
    !  Compute the coefficients of the Chebyshev expansions
    !  of degrees 12 and 24 of the function F.
    !
      fval(1) = 5.0E-01_RK * f(centr+hlgth)
      fval(13) = f(centr)
      fval(25) = 5.0E-01_RK * f(centr-hlgth)

      do i = 2, 12
        isym = 26-i
        fval(i) = f(hlgth*x(i-1)+centr)
        fval(isym) = f(centr-hlgth*x(i-1))
      end do

      call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  Compute the integral and error estimates.
    !
      resc12 = cheb12(13) * chebmo(m,13)
      ress12 = 0._RK
      estc = abs ( cheb24(25)*chebmo(m,25))+abs((cheb12(13)- &
        cheb24(13))*chebmo(m,13) )
      ests = 0._RK
      k = 11

      do j = 1, 6
        resc12 = resc12+cheb12(k)*chebmo(m,k)
        ress12 = ress12+cheb12(k+1)*chebmo(m,k+1)
        estc = estc+abs((cheb12(k)-cheb24(k))*chebmo(m,k))
        ests = ests+abs((cheb12(k+1)-cheb24(k+1))*chebmo(m,k+1))
        k = k-2
      end do

      resc24 = cheb24(25)*chebmo(m,25)
      ress24 = 0._RK
      resabs = abs(cheb24(25))
      k = 23

      do j = 1, 12

        resc24 = resc24+cheb24(k)*chebmo(m,k)
        ress24 = ress24+cheb24(k+1)*chebmo(m,k+1)
        resabs = resabs+abs(cheb24(k))+abs(cheb24(k+1))

        if ( j <= 5 ) then
          estc = estc+abs(cheb24(k)*chebmo(m,k))
          ests = ests+abs(cheb24(k+1)*chebmo(m,k+1))
        end if

        k = k-2

      end do

      resabs = resabs * abs ( hlgth )

      if ( integr == 1 ) then
        result = conc * resc24-cons*ress24
        abserr = abs ( conc * estc ) + abs ( cons * ests )
      else
        result = conc*ress24+cons*resc24
        abserr = abs(conc*ests)+abs(cons*estc)
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qc25s ( f, a, b, bl, br, alfa, beta, ri, rj, rg, rh, result, abserr, resasc, integr, neval )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QC25S returns rules for algebraico-logarithmic end point singularities.
    !
    !  Discussion:
    !
    !    This routine computes 
    !      i = integral of F(X) * W(X) over (bl,br), 
    !    with error estimate, where the weight function W(X) has a singular
    !    behavior of algebraico-logarithmic type at the points
    !    a and/or b. 
    !
    !    The interval (bl,br) is a subinterval of (a,b).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) BL, BR, the lower and upper limits of integration.
    !    A <= BL < BR <= B.
    !
    !    Input, real(RK) ALFA, BETA, parameters in the weight function.
    !
    !    Input, real(RK) RI(25), RJ(25), RG(25), RH(25), modified Chebyshev moments 
    !    for the application of the generalized Clenshaw-Curtis method,
    !    computed in QMOMO.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral, computed by 
    !    using a generalized clenshaw-curtis method if b1 = a or br = b.
    !    In all other cases the 15-point Kronrod rule is applied, obtained by
    !    optimal addition of abscissae to the 7-point Gauss rule.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, real(RK) RESASC, approximation to the integral of abs(F*W-I/(B-A)).
    !
    !    Input, integer(IK) INTEGR,  determines the weight function
    !    1, w(x) = (x-a)**alfa*(b-x)**beta
    !    2, w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
    !    3, w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
    !    4, w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !  Local Parameters:
    !
    !           fval   - value of the function f at the points
    !                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
    !                    k = 0, ..., 24
    !           cheb12 - coefficients of the Chebyshev series expansion
    !                    of degree 12, for the function f, in the interval
    !                    (bl,br)
    !           cheb24 - coefficients of the Chebyshev series expansion
    !                    of degree 24, for the function f, in the interval
    !                    (bl,br)
    !           res12  - approximation to the integral obtained from cheb12
    !           res24  - approximation to the integral obtained from cheb24
    !           qwgts  - external function subprogram defining the four
    !                    possible weight functions
    !           hlgth  - half-length of the interval (bl,br)
    !           centr  - mid point of the interval (bl,br)
    !
    !           the vector x contains the values cos(k*pi/24)
    !           k = 1, ..., 11, to be used for the computation of the
    !           Chebyshev series expansion of f.
    !
      implicit none

      real(RK) a
      real(RK) abserr
      real(RK) alfa
      real(RK) b
      real(RK) beta
      real(RK) bl
      real(RK) br
      real(RK) centr
      real(RK) cheb12(13)
      real(RK) cheb24(25)
      real(RK) dc
      real(RK), external :: f
      real(RK) factor
      real(RK) fix
      real(RK) fval(25)
      real(RK) hlgth
      integer(IK) i
      integer(IK) integr
      integer(IK) isym
      integer(IK) neval
      !real(RK), external :: qwgts
      real(RK) resabs
      real(RK) resasc
      real(RK) result
      real(RK) res12
      real(RK) res24
      real(RK) rg(25)
      real(RK) rh(25)
      real(RK) ri(25)
      real(RK) rj(25)
      real(RK) u
      real(RK), dimension ( 11 ) :: x = [ &
        9.914448613738104E-01_RK,     9.659258262890683E-01_RK, &
        9.238795325112868E-01_RK,     8.660254037844386E-01_RK, &
        7.933533402912352E-01_RK,     7.071067811865475E-01_RK, &
        6.087614290087206E-01_RK,     5.000000000000000E-01_RK, &
        3.826834323650898E-01_RK,     2.588190451025208E-01_RK, &
        1.305261922200516E-01_RK]

      neval = 25

      if ( bl == a .and. (alfa /= 0._RK .or. integr == 2 .or. integr == 4)) then
        go to 10
      end if

      if ( br == b .and. (beta /= 0._RK .or. integr == 3 .or. integr == 4)) &
        go to 140
    !
    !  If a > bl and b < br, apply the 15-point Gauss-Kronrod scheme.
    !
      call qk15w ( f, qwgts, a, b, alfa, beta, integr, bl, br, result, abserr, resabs, resasc )

      neval = 15
      return
    !
    !  This part of the program is executed only if a = bl.
    !
    !  Compute the Chebyshev series expansion of the function
    !  f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta*f(0.5*(br-a)*x+0.5*(br+a))
    !
    10 continue

      hlgth = 5.0E-01_RK*(br-bl)
      centr = 5.0E-01_RK*(br+bl)
      fix = b-centr
      fval(1) = 5.0E-01_RK*f(hlgth+centr)*(fix-hlgth)**beta
      fval(13) = f(centr)*(fix**beta)
      fval(25) = 5.0E-01_RK*f(centr-hlgth)*(fix+hlgth)**beta

      do i = 2, 12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix-u)**beta
        fval(isym) = f(centr-u)*(fix+u)**beta
      end do

      factor = hlgth**(alfa+1._RK)
      result = 0._RK
      abserr = 0._RK
      res12 = 0._RK
      res24 = 0._RK

      if ( integr > 2 ) go to 70

      call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  integr = 1  (or 2)
    !
      do i = 1, 13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
      end do

      do i = 14, 25
        res24 = res24 + cheb24(i) * ri(i)
      end do

      if ( integr == 1 ) go to 130
    !
    !  integr = 2
    !
      dc = log ( br - bl )
      result = res24 * dc
      abserr = abs((res24-res12)*dc)
      res12 = 0._RK
      res24 = 0._RK

      do i = 1, 13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res24+cheb24(i)*rg(i)
      end do

      do i = 14, 25
        res24 = res24+cheb24(i)*rg(i)
      end do

      go to 130
    !
    !  Compute the Chebyshev series expansion of the function
    !  F4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
    !
    70 continue

      fval(1) = fval(1) * log ( fix - hlgth )
      fval(13) = fval(13) * log ( fix )
      fval(25) = fval(25) * log ( fix + hlgth )

      do i = 2, 12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i) * log ( fix - u )
        fval(isym) = fval(isym) * log ( fix + u )
      end do

      call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  integr = 3  (or 4)
    !
      do i = 1, 13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
      end do

      do i = 14, 25
        res24 = res24+cheb24(i)*ri(i)
      end do

      if ( integr == 3 ) then
        go to 130
      end if
    !
    !  integr = 4
    !
      dc = log ( br - bl )
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0._RK
      res24 = 0._RK

      do i = 1, 13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res24+cheb24(i)*rg(i)
      end do

      do i = 14, 25
        res24 = res24+cheb24(i)*rg(i)
      end do

    130 continue

      result = (result+res24)*factor
      abserr = (abserr+abs(res24-res12))*factor
      go to 270
    !
    !  This part of the program is executed only if b = br.
    !
    !  Compute the Chebyshev series expansion of the function
    !  f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa*f(0.5*(b-bl)*x+0.5*(b+bl))
    !
    140 continue

      hlgth = 5.0E-01_RK*(br-bl)
      centr = 5.0E-01_RK*(br+bl)
      fix = centr-a
      fval(1) = 5.0E-01_RK*f(hlgth+centr)*(fix+hlgth)**alfa
      fval(13) = f(centr)*(fix**alfa)
      fval(25) = 5.0E-01_RK*f(centr-hlgth)*(fix-hlgth)**alfa

      do i = 2, 12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix+u)**alfa
        fval(isym) = f(centr-u)*(fix-u)**alfa
      end do

      factor = hlgth**(beta+1._RK)
      result = 0._RK
      abserr = 0._RK
      res12 = 0._RK
      res24 = 0._RK

      if ( integr == 2 .or. integr == 4 ) then
        go to 200
      end if
    !
    !  integr = 1  (or 3)
    !
      call qcheb ( x, fval, cheb12, cheb24 )

      do i = 1, 13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
      end do

      do i = 14, 25
        res24 = res24+cheb24(i)*rj(i)
      end do

      if ( integr == 1 ) go to 260
    !
    !  integr = 3
    !
      dc = log ( br - bl )
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0._RK
      res24 = 0._RK

      do i = 1, 13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
      end do

      do i = 14, 25
        res24 = res24+cheb24(i)*rh(i)
      end do

      go to 260
    !
    !  Compute the Chebyshev series expansion of the function
    !  f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
    !
    200 continue

      fval(1) = fval(1) * log ( hlgth + fix )
      fval(13) = fval(13) * log ( fix )
      fval(25) = fval(25) * log ( fix - hlgth )

      do i = 2, 12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i) * log(u+fix)
        fval(isym) = fval(isym) * log(fix-u)
      end do

      call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  integr = 2  (or 4)
    !
      do i = 1, 13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
      end do

      do i = 14, 25
        res24 = res24+cheb24(i)*rj(i)
      end do

      if ( integr == 2 ) go to 260

      dc = log(br-bl)
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0._RK
      res24 = 0._RK
    !
    !  integr = 4
    !
      do i = 1, 13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
      end do

      do i = 14, 25
        res24 = res24+cheb24(i)*rh(i)
      end do

    260 continue

      result = (result+res24)*factor
      abserr = (abserr+abs(res24-res12))*factor

    270 continue

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qcheb ( x, fval, cheb12, cheb24 )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QCHEB computes the Chebyshev series expansion.
    !
    !  Discussion:
    !
    !    This routine computes the Chebyshev series expansion
    !    of degrees 12 and 24 of a function using a fast Fourier transform method
    !
    !      f(x) = sum(k=1, ...,13) (cheb12(k)*t(k-1,x)),
    !      f(x) = sum(k=1, ...,25) (cheb24(k)*t(k-1,x)),
    !
    !    where T(K,X) is the Chebyshev polynomial of degree K.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(RK) X(11), contains the values of COS(K*PI/24), for K = 1 to 11.
    !
    !    Input/output, real(RK) FVAL(25), the function values at the points
    !    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24, where (a,b) is the 
    !    approximation interval.  FVAL(1) and FVAL(25) are divided by two
    !    These values are destroyed at output.
    !
    !    Output, real(RK) CHEB12(13), the Chebyshev coefficients for degree 12.
    !
    !    Output, real(RK) CHEB24(25), the Chebyshev coefficients for degree 24.
    !
      implicit none

      real(RK) alam
      real(RK) alam1
      real(RK) alam2
      real(RK) cheb12(13)
      real(RK) cheb24(25)
      real(RK) fval(25)
      integer(IK) i
      integer(IK) j
      real(RK) part1
      real(RK) part2
      real(RK) part3
      real(RK) v(12)
      real(RK) x(11)

      do i = 1, 12
        j = 26-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
      end do

      alam1 = v(1)-v(9)
      alam2 = x(6)*(v(3)-v(7)-v(11))
      cheb12(4) = alam1+alam2
      cheb12(10) = alam1-alam2
      alam1 = v(2)-v(8)-v(10)
      alam2 = v(4)-v(6)-v(12)
      alam = x(3)*alam1+x(9)*alam2
      cheb24(4) = cheb12(4)+alam
      cheb24(22) = cheb12(4)-alam
      alam = x(9)*alam1-x(3)*alam2
      cheb24(10) = cheb12(10)+alam
      cheb24(16) = cheb12(10)-alam
      part1 = x(4)*v(5)
      part2 = x(8)*v(9)
      part3 = x(6)*v(7)
      alam1 = v(1)+part1+part2
      alam2 = x(2)*v(3)+part3+x(10)*v(11)
      cheb12(2) = alam1+alam2
      cheb12(12) = alam1-alam2
      alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8)+x(9)*v(10)+x(11)*v(12)
      cheb24(2) = cheb12(2)+alam
      cheb24(24) = cheb12(2)-alam
      alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8)+x(3)*v(10)-x(1)*v(12)
      cheb24(12) = cheb12(12)+alam
      cheb24(14) = cheb12(12)-alam
      alam1 = v(1)-part1+part2
      alam2 = x(10)*v(3)-part3+x(2)*v(11)
      cheb12(6) = alam1+alam2
      cheb12(8) = alam1-alam2
      alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6)-x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
      cheb24(6) = cheb12(6)+alam
      cheb24(20) = cheb12(6)-alam
      alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8)-x(9)*v(10)-x(5)*v(12)
      cheb24(8) = cheb12(8)+alam
      cheb24(18) = cheb12(8)-alam

      do i = 1, 6
        j = 14-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
      end do

      alam1 = v(1)+x(8)*v(5)
      alam2 = x(4)*v(3)
      cheb12(3) = alam1+alam2
      cheb12(11) = alam1-alam2
      cheb12(7) = v(1)-v(5)
      alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
      cheb24(3) = cheb12(3)+alam
      cheb24(23) = cheb12(3)-alam
      alam = x(6)*(v(2)-v(4)-v(6))
      cheb24(7) = cheb12(7)+alam
      cheb24(19) = cheb12(7)-alam
      alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
      cheb24(11) = cheb12(11)+alam
      cheb24(15) = cheb12(11)-alam

      do i = 1, 3
        j = 8-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
      end do

      cheb12(5) = v(1)+x(8)*v(3)
      cheb12(9) = fval(1)-x(8)*fval(3)
      alam = x(4)*v(2)
      cheb24(5) = cheb12(5)+alam
      cheb24(21) = cheb12(5)-alam
      alam = x(8)*fval(2)-fval(4)
      cheb24(9) = cheb12(9)+alam
      cheb24(17) = cheb12(9)-alam
      cheb12(1) = fval(1)+fval(3)
      alam = fval(2)+fval(4)
      cheb24(1) = cheb12(1)+alam
      cheb24(25) = cheb12(1)-alam
      cheb12(13) = v(1)-v(3)
      cheb24(13) = cheb12(13)
      alam = 1._RK/6._RK

      do i = 2, 12
        cheb12(i) = cheb12(i)*alam
      end do

      alam = 5.0E-01*alam
      cheb12(1) = cheb12(1)*alam
      cheb12(13) = cheb12(13)*alam

      do i = 2, 24
        cheb24(i) = cheb24(i)*alam
      end do

      cheb24(1) = 0.5E+00_RK * alam*cheb24(1)
      cheb24(25) = 0.5E+00_RK * alam*cheb24(25)

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qextr ( n, epstab, result, abserr, res3la, nres )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QEXTR carries out the Epsilon extrapolation algorithm.
    !
    !  Discussion:
    !
    !    The routine determines the limit of a given sequence of approximations, 
    !    by means of the epsilon algorithm of P. Wynn.  An estimate of the 
    !    absolute error is also given.  The condensed epsilon table is computed.
    !    Only those elements needed for the computation of the next diagonal
    !    are preserved.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, integer(IK) N, indicates the entry of EPSTAB which contains
    !    the new element in the first column of the epsilon table.
    !
    !    Input/output, real(RK) EPSTAB(52), the two lower diagonals of the triangular
    !    epsilon table.  The elements are numbered starting at the right-hand 
    !    corner of the triangle.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, estimate of the absolute error computed from
    !    RESULT and the 3 previous results.
    !
    !    ?, real(RK) RES3LA(3), the last 3 results.
    !
    !    Input/output, integer(IK) NRES, the number of calls to the routine.  This
    !    should be zero on the first call, and is automatically updated
    !    before return.
    !
    !  Local Parameters:
    !
    !           e0     - the 4 elements on which the
    !           e1       computation of a new element in
    !           e2       the epsilon table is based
    !           e3                 e0
    !                        e3    e1    new
    !                              e2
    !           newelm - number of elements to be computed in the new
    !                    diagonal
    !           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
    !           result - the element in the new diagonal with least value
    !                    of error
    !           limexp is the maximum number of elements the epsilon table
    !           can contain. if this number is reached, the upper diagonal
    !           of the epsilon table is deleted.
    !
      implicit none

      real(RK) abserr
      real(RK) delta1
      real(RK) delta2
      real(RK) delta3
      real(RK) epsinf
      real(RK) epstab(52)
      real(RK) error
      real(RK) err1
      real(RK) err2
      real(RK) err3
      real(RK) e0
      real(RK) e1
      real(RK) e1abs
      real(RK) e2
      real(RK) e3
      integer(IK) i
      integer(IK) ib
      integer(IK) ib2
      integer(IK) ie
      integer(IK) indx
      integer(IK) k1
      integer(IK) k2
      integer(IK) k3
      integer(IK) limexp
      integer(IK) n
      integer(IK) newelm
      integer(IK) nres
      integer(IK) num
      real(RK) res
      real(RK) result
      real(RK) res3la(3)
      real(RK) ss
      real(RK) tol1
      real(RK) tol2
      real(RK) tol3

      nres = nres+1
      abserr = huge ( abserr )
      result = epstab(n)

      if ( n < 3 ) then
        abserr = max ( abserr,0.5E+00_RK * epsilon ( result ) *abs(result))
        return
      end if

      limexp = 50
      epstab(n+2) = epstab(n)
      newelm = (n-1)/2
      epstab(n) = huge ( epstab(n) )
      num = n
      k1 = n

      do i = 1, newelm

        k2 = k1-1
        k3 = k1-2
        res = epstab(k1+2)
        e0 = epstab(k3)
        e1 = epstab(k2)
        e2 = res
        e1abs = abs(e1)
        delta2 = e2-e1
        err2 = abs(delta2)
        tol2 = max ( abs(e2),e1abs)* epsilon ( e2 )
        delta3 = e1-e0
        err3 = abs(delta3)
        tol3 = max ( e1abs,abs(e0))* epsilon ( e0 )
    !
    !  If e0, e1 and e2 are equal to within machine accuracy, convergence 
    !  is assumed.
    !
        if ( err2 <= tol2 .and. err3 <= tol3 ) then
          result = res
          abserr = err2+err3
          abserr = max ( abserr,0.5E+00_RK* epsilon ( result ) *abs(result))
          return
        end if

        e3 = epstab(k1)
        epstab(k1) = e1
        delta1 = e1-e3
        err1 = abs(delta1)
        tol1 = max ( e1abs,abs(e3))* epsilon ( e3 )
    !
    !  If two elements are very close to each other, omit a part
    !  of the table by adjusting the value of N.
    !
        if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) go to 20

        ss = 1._RK/delta1+1._RK/delta2-1._RK/delta3
        epsinf = abs ( ss*e1 )
    !
    !  Test to detect irregular behavior in the table, and
    !  eventually omit a part of the table adjusting the value of N.
    !
        if ( epsinf > 1.0E-04_RK ) go to 30

    20  continue

        n = i+i-1
        exit
    !
    !  Compute a new element and eventually adjust the value of RESULT.
    !
    30  continue

        res = e1+1._RK/ss
        epstab(k1) = res
        k1 = k1-2
        error = err2+abs(res-e2)+err3

        if ( error <= abserr ) then
          abserr = error
          result = res
        end if

      end do
    !
    !  Shift the table.
    !
      if ( n == limexp ) then
        n = 2*(limexp/2)-1
      end if

      if ( (num/2)*2 == num ) then
        ib = 2
      else
        ib = 1
      end if

      ie = newelm+1

      do i = 1, ie
        ib2 = ib+2
        epstab(ib) = epstab(ib2)
        ib = ib2
      end do

      if ( num /= n ) then

        indx = num-n+1

        do i = 1, n
          epstab(i)= epstab(indx)
          indx = indx+1
        end do

      end if

      if ( nres < 4 ) then
        res3la(nres) = result
        abserr = huge ( abserr )
      else
        abserr = abs(result-res3la(3))+abs(result-res3la(2)) &
          +abs(result-res3la(1))
        res3la(1) = res3la(2)
        res3la(2) = res3la(3)
        res3la(3) = result
      end if

      abserr = max ( abserr,0.5E+00_RK* epsilon ( result ) *abs(result))

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qfour ( f, a, b, omega, integr, epsabs, epsrel, limit, icall, &
      maxp1, result, abserr, neval, ier, alist, blist, rlist, elist, iord, &
      nnlog, momcom, chebmo )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QFOUR estimates the integrals of oscillatory functions.
    !
    !  Discussion:
    !
    !    This routine calculates an approximation RESULT to a definite integral
    !      I = integral of F(X) * COS(OMEGA*X) 
    !    or
    !      I = integral of F(X) * SIN(OMEGA*X) 
    !    over (A,B), hopefully satisfying:
    !      | I - RESULT | <= max ( epsabs, epsrel * |I| ) ).
    !
    !    QFOUR is called by QAWO and QAWF.  It can also be called directly in 
    !    a user-written program.  In the latter case it is possible for the 
    !    user to determine the first dimension of array CHEBMO(MAXP1,25).
    !    See also parameter description of MAXP1.  Additionally see
    !    parameter description of ICALL for eventually re-using
    !    Chebyshev moments computed during former call on subinterval
    !    of equal length abs(B-A).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) OMEGA, the multiplier of X in the weight function.
    !
    !    Input, integer(IK) INTEGR, indicates the weight functions to be used.
    !    = 1, w(x) = cos(omega*x)
    !    = 2, w(x) = sin(omega*x)
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Input, integer(IK) LIMIT, the maximum number of subintervals of [A,B]
    !    that can be generated.
    !
    !    icall  - integer(IK)
    !                     if qfour is to be used only once, ICALL must
    !                     be set to 1.  assume that during this call, the
    !                     Chebyshev moments (for clenshaw-curtis integration
    !                     of degree 24) have been computed for intervals of
    !                     lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
    !                     the Chebyshev moments already computed can be
    !                     re-used in subsequent calls, if qfour must be
    !                     called twice or more times on intervals of the
    !                     same length abs(b-a). from the second call on, one
    !                     has to put then ICALL > 1.
    !                     if ICALL < 1, the routine will end with ier = 6.
    !
    !            maxp1  - integer(IK)
    !                     gives an upper bound on the number of
    !                     Chebyshev moments which can be stored, i.e.
    !                     for the intervals of lenghts abs(b-a)*2**(-l),
    !                     l=0,1, ..., maxp1-2, maxp1 >= 1.
    !                     if maxp1 < 1, the routine will end with ier = 6.
    !                     increasing (decreasing) the value of maxp1
    !                     decreases (increases) the computational time but
    !                     increases (decreases) the required memory space.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !            ier    - integer(IK)
    !                     ier = 0 normal and reliable termination of the
    !                             routine. it is assumed that the
    !                             requested accuracy has been achieved.
    !                   - ier > 0 abnormal termination of the routine.
    !                             the estimates for integral and error are
    !                             less reliable. it is assumed that the
    !                             requested accuracy has not been achieved.
    !                     ier = 1 maximum number of subdivisions allowed
    !                             has been achieved. one can allow more
    !                             subdivisions by increasing the value of
    !                             limit (and taking according dimension
    !                             adjustments into account). however, if
    !                             this yields no improvement it is advised
    !                             to analyze the integrand, in order to
    !                             determine the integration difficulties.
    !                             if the position of a local difficulty can
    !                             be determined (e.g. singularity,
    !                             discontinuity within the interval) one
    !                             will probably gain from splitting up the
    !                             interval at this point and calling the
    !                             integrator on the subranges. if possible,
    !                             an appropriate special-purpose integrator
    !                             should be used which is designed for
    !                             handling the type of difficulty involved.
    !                         = 2 the occurrence of roundoff error is
    !                             detected, which prevents the requested
    !                             tolerance from being achieved.
    !                             the error may be under-estimated.
    !                         = 3 extremely bad integrand behavior occurs
    !                             at some points of the integration
    !                             interval.
    !                         = 4 the algorithm does not converge. roundoff
    !                             error is detected in the extrapolation
    !                             table. it is presumed that the requested
    !                             tolerance cannot be achieved due to
    !                             roundoff in the extrapolation table, and
    !                             that the returned result is the best which
    !                             can be obtained.
    !                         = 5 the integral is probably divergent, or
    !                             slowly convergent. it must be noted that
    !                             divergence can occur with any other value
    !                             of ier > 0.
    !                         = 6 the input is invalid, because
    !                             epsabs < 0 and epsrel < 0,
    !                             or (integr /= 1 and integr /= 2) or
    !                             ICALL < 1 or maxp1 < 1.
    !                             result, abserr, neval, last, rlist(1),
    !                             elist(1), iord(1) and nnlog(1) are set to
    !                             zero. alist(1) and blist(1) are set to a
    !                             and b respectively.
    !
    !    Workspace, real(RK) ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
    !    through LAST the left and right ends of the partition subintervals.
    !
    !    Workspace, real(RK) RLIST(LIMIT), contains in entries 1 through LAST
    !    the integral approximations on the subintervals.
    !
    !    Workspace, real(RK) ELIST(LIMIT), contains in entries 1 through LAST
    !    the absolute error estimates on the subintervals.
    !
    !            iord   - integer(IK)
    !                     vector of dimension at least limit, the first k
    !                     elements of which are pointers to the error
    !                     estimates over the subintervals, such that
    !                     elist(iord(1)), ..., elist(iord(k)), form
    !                     a decreasing sequence, with k = last
    !                     if last <= (limit/2+2), and
    !                     k = limit+1-last otherwise.
    !
    !            nnlog  - integer(IK)
    !                     vector of dimension at least limit, indicating the
    !                     subdivision levels of the subintervals, i.e.
    !                     iwork(i) = l means that the subinterval numbered
    !                     i is of length abs(b-a)*2**(1-l)
    !
    !         on entry and return
    !            momcom - integer(IK)
    !                     indicating that the Chebyshev moments have been
    !                     computed for intervals of lengths
    !                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
    !                     momcom < maxp1
    !
    !            chebmo - real
    !                     array of dimension (maxp1,25) containing the
    !                     Chebyshev moments
    !
    !  Local Parameters:
    !
    !           alist     - list of left end points of all subintervals
    !                       considered up to now
    !           blist     - list of right end points of all subintervals
    !                       considered up to now
    !           rlist(i)  - approximation to the integral over
    !                       (alist(i),blist(i))
    !           rlist2    - array of dimension at least limexp+2 containing
    !                       the part of the epsilon table which is still
    !                       needed for further computations
    !           elist(i)  - error estimate applying to rlist(i)
    !           maxerr    - pointer to the interval with largest error
    !                       estimate
    !           errmax    - elist(maxerr)
    !           erlast    - error on the interval currently subdivided
    !           area      - sum of the integrals over the subintervals
    !           errsum    - sum of the errors over the subintervals
    !           errbnd    - requested accuracy max(epsabs,epsrel*
    !                       abs(result))
    !           *****1    - variable for the left subinterval
    !           *****2    - variable for the right subinterval
    !           last      - index for subdivision
    !           nres      - number of calls to the extrapolation routine
    !           numrl2    - number of elements in rlist2. if an appropriate
    !                       approximation to the compounded integral has
    !                       been obtained it is put in rlist2(numrl2) after
    !                       numrl2 has been increased by one
    !           small     - length of the smallest interval considered
    !                       up to now, multiplied by 1.5
    !           erlarg    - sum of the errors over the intervals larger
    !                       than the smallest interval considered up to now
    !           extrap    - logical variable denoting that the routine is
    !                       attempting to perform extrapolation, i.e. before
    !                       subdividing the smallest interval we try to
    !                       decrease the value of erlarg
    !           noext     - logical variable denoting that extrapolation
    !                       is no longer allowed (true value)
    !
      implicit none

      integer(IK) limit
      integer(IK) maxp1

      real(RK) a
      real(RK) abseps
      real(RK) abserr
      real(RK) alist(limit)
      real(RK) area
      real(RK) area1
      real(RK) area12
      real(RK) area2
      real(RK) a1
      real(RK) a2
      real(RK) b
      real(RK) blist(limit)
      real(RK) b1
      real(RK) b2
      real(RK) chebmo(maxp1,25)
      real(RK) correc
      real(RK) defab1
      real(RK) defab2
      real(RK) defabs
      real(RK) domega
      real(RK) dres
      real(RK) elist(limit)
      real(RK) epsabs
      real(RK) epsrel
      real(RK) erlarg
      real(RK) erlast
      real(RK) errbnd
      real(RK) errmax
      real(RK) error1
      real(RK) erro12
      real(RK) error2
      real(RK) errsum
      real(RK) ertest
      logical extall
      logical extrap
      real(RK), external :: f
      integer(IK) icall
      integer(IK) id
      integer(IK) ier
      integer(IK) ierro
      integer(IK) integr
      integer(IK) iord(limit)
      integer(IK) iroff1
      integer(IK) iroff2
      integer(IK) iroff3
      integer(IK) jupbnd
      integer(IK) k
      integer(IK) ksgn
      integer(IK) ktmin
      integer(IK) last
      integer(IK) maxerr
      integer(IK) momcom
      integer(IK) nev
      integer(IK) neval
      integer(IK) nnlog(limit)
      logical noext
      integer(IK) nres
      integer(IK) nrmax
      integer(IK) nrmom
      integer(IK) numrl2
      real(RK) omega
      real(RK) resabs
      real(RK) reseps
      real(RK) result
      real(RK) res3la(3)
      real(RK) rlist(limit)
      real(RK) rlist2(52)
      real(RK) small
      real(RK) width
    !
    !  the dimension of rlist2 is determined by  the value of
    !  limexp in QEXTR (rlist2 should be of dimension
    !  (limexp+2) at least).
    !
    !  Test on validity of parameters.
    !
      ier = 0
      neval = 0
      last = 0
      result = 0._RK
      abserr = 0._RK
      alist(1) = a
      blist(1) = b
      rlist(1) = 0._RK
      elist(1) = 0._RK
      iord(1) = 0
      nnlog(1) = 0

      if ( (integr /= 1 .and. integr /= 2) .or. (epsabs < 0._RK .and. epsrel < 0._RK) .or. icall < 1 .or. maxp1 < 1 ) then
        ier = 6
        return
      end if
    !
    !  First approximation to the integral.
    !
      domega = abs ( omega )
      nrmom = 0

      if ( icall <= 1 ) then
        momcom = 0
      end if

      call qc25o ( f, a, b, domega, integr, nrmom, maxp1, 0, result, abserr, neval, defabs, resabs, momcom, chebmo )
    !
    !  Test on accuracy.
    !
      dres = abs(result)
      errbnd = max ( epsabs,epsrel*dres)
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if ( abserr <= 1.0E+02_RK * epsilon ( defabs ) *defabs .and. &
        abserr > errbnd ) ier = 2

      if ( limit == 1 ) then
        ier = 1
      end if

      if ( ier /= 0 .or. abserr <= errbnd ) then
        go to 200
      end if
    !
    !  Initializations
    !
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = huge ( abserr )
      nrmax = 1
      extrap = .false.
      noext = .false.
      ierro = 0
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ktmin = 0
      small = abs(b-a)*7.5E-01_RK
      nres = 0
      numrl2 = 0
      extall = .false.

      if ( 5.0E-01*abs(b-a)*domega <= 2._RK) then
        numrl2 = 1
        extall = .true.
        rlist2(1) = result
      end if

      if ( 2.5E-01_RK * abs(b-a) * domega <= 2._RK ) then
        extall = .true.
      end if

      if ( dres >= (1._RK - 5.0E+01_RK * epsilon ( defabs ) )*defabs ) then
        ksgn = 1
      else
        ksgn = -1
      end if
    !
    !  main do-loop
    !
      do last = 2, limit
    !
    !  Bisect the subinterval with the nrmax-th largest error estimate.
    !
        nrmom = nnlog(maxerr)+1
        a1 = alist(maxerr)
        b1 = 5.0E-01*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax

        call qc25o ( f, a1, b1, domega, integr, nrmom, maxp1, 0, area1, error1, nev, resabs, defab1, momcom, chebmo )

        neval = neval+nev

        call qc25o ( f, a2, b2, domega, integr, nrmom, maxp1, 1, area2, error2, nev, resabs, defab2, momcom, chebmo )

        neval = neval+nev
    !
    !  Improve previous approximations to integral and error and
    !  test for accuracy.
    !
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if ( defab1 == error1 .or. defab2 == error2 ) go to 25
        if ( abs(rlist(maxerr)-area12) > 1.0E-05_RK*abs(area12) &
        .or. erro12 < 9.9E-01_RK*errmax ) go to 20
        if ( extrap ) iroff2 = iroff2+1

        if ( .not.extrap ) then
          iroff1 = iroff1+1
        end if

    20  continue

        if ( last > 10 .and. erro12 > errmax ) iroff3 = iroff3 + 1

    25  continue

        rlist(maxerr) = area1
        rlist(last) = area2
        nnlog(maxerr) = nrmom
        nnlog(last) = nrmom
        errbnd = max ( epsabs,epsrel*abs(area))
    !
    !  Test for roundoff error and eventually set error flag
    !
        if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) ier = 2

        if ( iroff2 >= 5) ierro = 3
    !
    !  Set error flag in the case that the number of subintervals
    !  equals limit.
    !
        if ( last == limit ) then
          ier = 1
        end if
    !
    !  Set error flag in the case of bad integrand behavior at
    !  a point of the integration range.
    !
        if ( max ( abs(a1),abs(b2)) <= (1._RK + 1.0E+03_RK * epsilon ( a1 ) )*(abs(a2) + 1.0E+03_RK * tiny ( a2 ) )) then
          ier = 4
        end if
    !
    !  Append the newly-created intervals to the list.
    !
        if ( error2 <= error1 ) then
          alist(last) = a2
          blist(maxerr) = b1
          blist(last) = b2
          elist(maxerr) = error1
          elist(last) = error2
        else
          alist(maxerr) = a2
          alist(last) = a1
          blist(last) = b1
          rlist(maxerr) = area2
          rlist(last) = area1
          elist(maxerr) = error2
          elist(last) = error1
        end if
    !
    !  Call QSORT to maintain the descending ordering
    !  in the list of error estimates and select the subinterval
    !  with nrmax-th largest error estimate (to be bisected next).
    !

        call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

        if ( errsum <= errbnd ) then
          go to 170
        end if

        if ( ier /= 0 ) then
          exit
        end if

        if ( last == 2 .and. extall ) go to 120

        if ( noext ) then
          cycle
        end if

        if ( .not. extall ) go to 50
        erlarg = erlarg-erlast
        if ( abs(b1-a1) > small ) erlarg = erlarg+erro12
        if ( extrap ) go to 70
    !
    !  Test whether the interval to be bisected next is the
    !  smallest interval.
    !
    50  continue

        width = abs(blist(maxerr)-alist(maxerr))

        if ( width > small ) then
          cycle
        end if

        if ( extall ) go to 60
    !
    !  Test whether we can start with the extrapolation procedure
    !  (we do this if we integrate over the next interval with
    !  use of a Gauss-Kronrod rule - see QC25O).
    !
        small = small*5.0E-01_RK

        if ( 2.5E-01_RK*width*domega > 2._RK ) then
          cycle
        end if

        extall = .true.
        go to 130

    60  continue

        extrap = .true.
        nrmax = 2

    70  continue

        if ( ierro == 3 .or. erlarg <= ertest ) go to 90
    !
    !  The smallest interval has the largest error.
    !  Before bisecting decrease the sum of the errors over the
    !  larger intervals (ERLARG) and perform extrapolation.
    !
        jupbnd = last

        if ( last > (limit/2+2) ) then
          jupbnd = limit+3-last
        end if

        id = nrmax

        do k = id, jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 140
          nrmax = nrmax+1
        end do
    !
    !  Perform extrapolation.
    !
    90  continue

        numrl2 = numrl2+1
        rlist2(numrl2) = area

        if ( numrl2 < 3 ) go to 110

        call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
        ktmin = ktmin+1

        if ( ktmin > 5 .and. abserr < 1.0E-03_RK*errsum ) then
          ier = 5
        end if

        if ( abseps >= abserr ) go to 100

        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = max ( epsabs, epsrel*abs(reseps))

        if ( abserr <= ertest ) then
          exit
        end if
    !
    !  Prepare bisection of the smallest interval.
    !
    100 continue

        if ( numrl2 == 1 ) then
          noext = .true.
        end if

        if ( ier == 5 ) then
          exit
        end if

    110 continue

        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*5.0E-01_RK
        erlarg = errsum
        cycle

    120 continue

        small = small * 5.0E-01_RK
        numrl2 = numrl2 + 1
        rlist2(numrl2) = area

    130 continue

        ertest = errbnd
        erlarg = errsum

    140 continue

      end do
    !
    !  set the final result.
    !
      if ( abserr == huge ( abserr ) .or. nres == 0 ) then
        go to 170
      end if

      if ( ier+ierro == 0 ) go to 165
      if ( ierro == 3 ) abserr = abserr+correc
      if ( ier == 0 ) ier = 3
      if ( result /= 0._RK .and. area /= 0._RK ) go to 160
      if ( abserr > errsum ) go to 170
      if ( area == 0._RK ) go to 190
      go to 165

    160 continue

      if ( abserr/abs(result) > errsum/abs(area) ) go to 170
    !
    !  Test on divergence.
    !
      165 continue

      if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <= defabs*1.0E-02_RK ) goto 190

      if ( 1.0E-02_RK > (result/area) .or. (result/area) > 1.0E+02_RK .or. errsum >= abs(area) ) ier = 6

      go to 190
    !
    !  Compute global integral sum.
    !
    170 continue

      result = sum ( rlist(1:last) )

      abserr = errsum

    190 continue

      if (ier > 2) ier=ier-1

    200 continue

      if ( integr == 2 .and. omega < 0._RK ) then
        result = -result
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qk15 ( f, a, b, result, abserr, resabs, resasc )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 15-point Kronrod rule (RESK) 
    !    obtained by optimal addition of abscissae to the 7-point Gauss rule 
    !    (RESG).
    !
    !    Output, real(RK) ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(RK) RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(RK) RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 15-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 7-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 7-point Gauss rule
    !
    !           wgk    - weights of the 15-point Kronrod rule
    !
    !           wg     - weights of the 7-point Gauss rule
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 7-point Gauss formula
    !           resk   - result of the 15-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
      implicit none

      real(RK) a
      real(RK) absc
      real(RK) abserr
      real(RK) b
      real(RK) centr
      real(RK) dhlgth
      real(RK), external :: f
      real(RK) fc
      real(RK) fsum
      real(RK) fval1
      real(RK) fval2
      real(RK) fv1(7)
      real(RK) fv2(7)
      real(RK) hlgth
      integer(IK) j
      integer(IK) jtw
      integer(IK) jtwm1
      real(RK) resabs
      real(RK) resasc
      real(RK) resg
      real(RK) resk
      real(RK) reskh
      real(RK) result
      real(RK) wg(4)
      real(RK) wgk(8)
      real(RK) xgk(8)

      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
           9.914553711208126E-01_RK,   9.491079123427585E-01_RK, &
           8.648644233597691E-01_RK,   7.415311855993944E-01_RK, &
           5.860872354676911E-01_RK,   4.058451513773972E-01_RK, &
           2.077849550078985E-01_RK,   0._RK              /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
           2.293532201052922E-02_RK,   6.309209262997855E-02_RK, &
           1.047900103222502E-01_RK,   1.406532597155259E-01_RK, &
           1.690047266392679E-01_RK,   1.903505780647854E-01_RK, &
           2.044329400752989E-01_RK,   2.094821410847278E-01_RK/
      data wg(1),wg(2),wg(3),wg(4)/ &
           1.294849661688697E-01_RK,   2.797053914892767E-01_RK, &
           3.818300505051189E-01_RK,   4.179591836734694E-01_RK/
    !
      centr = 5.0E-01_RK*(a+b)
      hlgth = 5.0E-01_RK*(b-a)
      dhlgth = abs(hlgth)
    !
    !  Compute the 15-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
      fc = f(centr)
      resg = fc*wg(4)
      resk = fc*wgk(8)
      resabs = abs(resk)

      do j = 1, 3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
      end do

      do j = 1, 4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
      end do

      reskh = resk * 5.0E-01_RK
      resasc = wgk(8)*abs(fc-reskh)

      do j = 1, 7
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      end do

      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)

      if ( resasc /= 0._RK .and. abserr /= 0._RK ) then
        abserr = resasc*min ( 1._RK,(2.0E+02_RK*abserr/resasc)**1.5E+00_RK)
      end if

      if ( resabs > tiny ( resabs ) / (5.0E+01_RK* epsilon ( resabs ) ) ) then
        abserr = max (( epsilon ( resabs ) *5.0E+01_RK)*resabs,abserr)
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qk15i ( f, boun, inf, a, b, result, abserr, resabs, resasc )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QK15I applies a 15 point Gauss-Kronrod quadrature on an infinite interval.
    !
    !  Discussion:
    !
    !    The original infinite integration range is mapped onto the interval 
    !    (0,1) and (a,b) is a part of (0,1).  The routine then computes:
    !
    !    i = integral of transformed integrand over (a,b),
    !    j = integral of abs(transformed integrand) over (a,b).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) BOUN, the finite bound of the original integration range,
    !    or zero if INF is 2.
    !
    !    Input, integer(IK) INF, indicates the type of the interval.
    !    -1: the original interval is (-infinity,BOUN),
    !    +1, the original interval is (BOUN,+infinity),
    !    +2, the original interval is (-infinity,+infinity) and
    !    the integral is computed as the sum of two integrals, one 
    !    over (-infinity,0) and one over (0,+infinity).
    !
    !    Input, real(RK) A, B, the limits of integration, over a subrange of [0,1].
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained 
    !    by optimal addition of abscissae to the 7-point Gauss rule (RESG).
    !
    !    Output, real(RK) ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(RK) RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(RK) RESASC, approximation to the integral of the
    !    transformated integrand | F-I/(B-A) | over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc*  - abscissa
    !           tabsc* - transformed abscissa
    !           fval*  - function value
    !           resg   - result of the 7-point Gauss formula
    !           resk   - result of the 15-point Kronrod formula
    !           reskh  - approximation to the mean value of the transformed
    !                    integrand over (a,b), i.e. to i/(b-a)
    !
      implicit none

      real(RK) a
      real(RK) absc
      real(RK) absc1
      real(RK) absc2
      real(RK) abserr
      real(RK) b
      real(RK) boun
      real(RK) centr
      real(RK) dinf
      real(RK), external :: f
      real(RK) fc
      real(RK) fsum
      real(RK) fval1
      real(RK) fval2
      real(RK) fv1(7)
      real(RK) fv2(7)
      real(RK) hlgth
      integer(IK) inf
      integer(IK) j
      real(RK) resabs
      real(RK) resasc
      real(RK) resg
      real(RK) resk
      real(RK) reskh
      real(RK) result
      real(RK) tabsc1
      real(RK) tabsc2
      real(RK) wg(8)
      real(RK) wgk(8)
      real(RK) xgk(8)
    !
    !  the abscissae and weights are supplied for the interval
    !  (-1,1).  because of symmetry only the positive abscissae and
    !  their corresponding weights are given.
    !
    !           xgk    - abscissae of the 15-point Kronrod rule
    !                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
    !                    rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 7-point Gauss rule
    !
    !           wgk    - weights of the 15-point Kronrod rule
    !
    !           wg     - weights of the 7-point Gauss rule, corresponding
    !                    to the abscissae xgk(2), xgk(4), ...
    !                    wg(1), wg(3), ... are set to zero.
    !
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
           9.914553711208126E-01_RK,     9.491079123427585E-01_RK, &
           8.648644233597691E-01_RK,     7.415311855993944E-01_RK, &
           5.860872354676911E-01_RK,     4.058451513773972E-01_RK, &
           2.077849550078985E-01_RK,     0.000000000000000E+00_RK/

      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
           2.293532201052922E-02_RK,     6.309209262997855E-02_RK, &
           1.047900103222502E-01_RK,     1.406532597155259E-01_RK, &
           1.690047266392679E-01_RK,     1.903505780647854E-01_RK, &
           2.044329400752989E-01_RK,     2.094821410847278E-01_RK/

      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
           0.0000000000000000E+00_RK,     1.294849661688697E-01_RK, &
           0.0000000000000000E+00_RK,     2.797053914892767E-01_RK, &
           0.0000000000000000E+00_RK,     3.818300505051189E-01_RK, &
           0.0000000000000000E+00_RK,     4.179591836734694E-01_RK/

      dinf = min ( 1, inf )

      centr = 5.0E-01_RK*(a+b)
      hlgth = 5.0E-01_RK*(b-a)
      tabsc1 = boun+dinf*(1._RK-centr)/centr
      fval1 = f(tabsc1)
      if ( inf == 2 ) fval1 = fval1+f(-tabsc1)
      fc = (fval1/centr)/centr
    !
    !  Compute the 15-point Kronrod approximation to the integral,
    !  and estimate the error.
    !
      resg = wg(8)*fc
      resk = wgk(8)*fc
      resabs = abs(resk)

      do j = 1, 7

        absc = hlgth*xgk(j)
        absc1 = centr-absc
        absc2 = centr+absc
        tabsc1 = boun+dinf*(1._RK-absc1)/absc1
        tabsc2 = boun+dinf*(1._RK-absc2)/absc2
        fval1 = f(tabsc1)
        fval2 = f(tabsc2)

        if ( inf == 2 ) then
          fval1 = fval1+f(-tabsc1)
          fval2 = fval2+f(-tabsc2)
        end if

        fval1 = (fval1/absc1)/absc1
        fval2 = (fval2/absc2)/absc2
        fv1(j) = fval1
        fv2(j) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(j)*fsum
        resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
      end do

      reskh = resk * 5.0E-01
      resasc = wgk(8) * abs(fc-reskh)

      do j = 1, 7
        resasc = resasc + wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      end do

      result = resk * hlgth
      resasc = resasc * hlgth
      resabs = resabs * hlgth
      abserr = abs ( ( resk - resg ) * hlgth )

      if ( resasc /= 0._RK .and. abserr /= 0._RK) then
        abserr = resasc* min ( 1._RK,(2.0E+02_RK*abserr/resasc)**1.5E+00_RK)
      end if

      if ( resabs > tiny ( resabs ) / ( 5.0E+01_RK * epsilon ( resabs ) ) ) then
        abserr = max (( epsilon ( resabs ) *5.0E+01_RK)*resabs,abserr)
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qk15w ( f, w, p1, p2, p3, p4, kp, a, b, result, abserr, resabs, &
      resasc )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QK15W applies a 15 point Gauss-Kronrod rule for a weighted integrand.
    !
    !  Discussion:
    !
    !    This routine approximates 
    !      i = integral of f*w over (a,b), 
    !    with error estimate, and
    !      j = integral of abs(f*w) over (a,b)
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !              w      - real
    !                       function subprogram defining the integrand
    !                       weight function w(x). the actual name for w
    !                       needs to be declared e x t e r n a l in the
    !                       calling program.
    !
    !    ?, real(RK) P1, P2, P3, P4, parameters in the weight function
    !
    !    Input, integer(IK) KP, key for indicating the type of weight function
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained by
    !    optimal addition of abscissae to the 7-point Gauss rule (RESG).
    !
    !    Output, real(RK) ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(RK) RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(RK) RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc*  - abscissa
    !           fval*  - function value
    !           resg   - result of the 7-point Gauss formula
    !           resk   - result of the 15-point Kronrod formula
    !           reskh  - approximation to the mean value of f*w over (a,b),
    !                    i.e. to i/(b-a)
    !
      implicit none

      real(RK) a
      real(RK) absc
      real(RK) absc1
      real(RK) absc2
      real(RK) abserr
      real(RK) b
      real(RK) centr
      real(RK) dhlgth
      real(RK), external :: f
      real(RK) fc
      real(RK) fsum
      real(RK) fval1
      real(RK) fval2
      real(RK) fv1(7)
      real(RK) fv2(7)
      real(RK) hlgth
      integer(IK) j
      integer(IK) jtw
      integer(IK) jtwm1
      integer(IK) kp
      real(RK) p1
      real(RK) p2
      real(RK) p3
      real(RK) p4
      real(RK) resabs
      real(RK) resasc
      real(RK) resg
      real(RK) resk
      real(RK) reskh
      real(RK) result
      real(RK), external :: w
      real(RK), dimension ( 4 ) :: wg = (/ &
        1.294849661688697E-01_RK,     2.797053914892767E-01_RK, &
        3.818300505051889E-01_RK,     4.179591836734694E-01_RK /)
      real(RK) wgk(8)
      real(RK) xgk(8)
    !
    !  the abscissae and weights are given for the interval (-1,1).
    !  because of symmetry only the positive abscissae and their
    !  corresponding weights are given.
    !
    !           xgk    - abscissae of the 15-point Gauss-Kronrod rule
    !                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
    !                    rule
    !                    xgk(1), xgk(3), ... abscissae which are optimally
    !                    added to the 7-point Gauss rule
    !
    !           wgk    - weights of the 15-point Gauss-Kronrod rule
    !
    !           wg     - weights of the 7-point Gauss rule
    !
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
           9.914553711208126E-01_RK,     9.491079123427585E-01_RK, &
           8.648644233597691E-01_RK,     7.415311855993944E-01_RK, &
           5.860872354676911E-01_RK,     4.058451513773972E-01_RK, &
           2.077849550789850E-01_RK,     0.000000000000000E+00_RK/

      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
           2.293532201052922E-02_RK,     6.309209262997855E-02_RK, &
           1.047900103222502E-01_RK,     1.406532597155259E-01_RK, &
           1.690047266392679E-01_RK,     1.903505780647854E-01_RK, &
           2.044329400752989E-01_RK,     2.094821410847278E-01_RK/
    !
      centr = 5.0E-01_RK*(a+b)
      hlgth = 5.0E-01_RK*(b-a)
      dhlgth = abs(hlgth)
    !
    !  Compute the 15-point Kronrod approximation to the integral,
    !  and estimate the error.
    !
      fc = f(centr)*w(centr,p1,p2,p3,p4,kp)
      resg = wg(4)*fc
      resk = wgk(8)*fc
      resabs = abs(resk)

      do j = 1, 3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
      end do

      do j = 1, 4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
      end do

      reskh = resk*5.0E-01_RK
      resasc = wgk(8)*abs(fc-reskh)

      do j = 1, 7
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      end do

      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)

      if ( resasc /= 0._RK .and. abserr /= 0._RK) then
        abserr = resasc*min ( 1._RK,(2.0E+02_RK*abserr/resasc)**1.5E+00_RK)
      end if

      if ( resabs > tiny ( resabs ) /(5.0E+01_RK* epsilon ( resabs ) ) ) then
        abserr = max ( ( epsilon ( resabs ) * 5.0E+01_RK)*resabs,abserr)
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qk21 ( f, a, b, result, abserr, resabs, resasc )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 21-point Kronrod rule (resk) 
    !    obtained by optimal addition of abscissae to the 10-point Gauss 
    !    rule (resg).
    !
    !    Output, real(RK) ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(RK) RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(RK) RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
      implicit none

      real(RK) a
      real(RK) absc
      real(RK) abserr
      real(RK) b
      real(RK) centr
      real(RK) dhlgth
      real(RK), external :: f
      real(RK) fc
      real(RK) fsum
      real(RK) fval1
      real(RK) fval2
      real(RK) fv1(10)
      real(RK) fv2(10)
      real(RK) hlgth
      integer(IK) j
      integer(IK) jtw
      integer(IK) jtwm1
      real(RK) resabs
      real(RK) resasc
      real(RK) resg
      real(RK) resk
      real(RK) reskh
      real(RK) result
      real(RK) wg(5)
      real(RK) wgk(11)
      real(RK) xgk(11)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 21-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 10-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 10-point Gauss rule
    !
    !           wgk    - weights of the 21-point Kronrod rule
    !
    !           wg     - weights of the 10-point Gauss rule
    !
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
        xgk(9),xgk(10),xgk(11)/ &
           9.956571630258081E-01_RK,     9.739065285171717E-01_RK, &
           9.301574913557082E-01_RK,     8.650633666889845E-01_RK, &
           7.808177265864169E-01_RK,     6.794095682990244E-01_RK, &
           5.627571346686047E-01_RK,     4.333953941292472E-01_RK, &
           2.943928627014602E-01_RK,     1.488743389816312E-01_RK, &
           0.000000000000000E+00_RK/
    !
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
        wgk(9),wgk(10),wgk(11)/ &
           1.169463886737187E-02_RK,     3.255816230796473E-02_RK, &
           5.475589657435200E-02_RK,     7.503967481091995E-02_RK, &
           9.312545458369761E-02_RK,     1.093871588022976E-01_RK, &
           1.234919762620659E-01_RK,     1.347092173114733E-01_RK, &
           1.427759385770601E-01_RK,     1.477391049013385E-01_RK, &
           1.494455540029169E-01_RK/
    !
      data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
           6.667134430868814E-02_RK,     1.494513491505806E-01_RK, &
           2.190863625159820E-01_RK,     2.692667193099964E-01_RK, &
           2.955242247147529E-01_RK/
    !
    !
    !           list of major variables
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 10-point Gauss formula
    !           resk   - result of the 21-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
      centr = 5.0E-01_RK*(a+b)
      hlgth = 5.0E-01_RK*(b-a)
      dhlgth = abs(hlgth)
    !
    !  Compute the 21-point Kronrod approximation to the
    !  integral, and estimate the absolute error.
    !
      resg = 0._RK
      fc = f(centr)
      resk = wgk(11)*fc
      resabs = abs(resk)

      do j = 1, 5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
      end do

      do j = 1, 5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
      end do

      reskh = resk*5.0E-01_RK
      resasc = wgk(11)*abs(fc-reskh)

      do j = 1, 10
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      end do

      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)

      if ( resasc /= 0._RK .and. abserr /= 0._RK) then
        abserr = resasc*min ( 1._RK,(2.0E+02_RK*abserr/resasc)**1.5E+00_RK)
      end if

      if ( resabs > tiny ( resabs ) /(5.0E+01_RK * epsilon ( resabs ) )) then
        abserr = max (( epsilon ( resabs ) *5.0E+01_RK)*resabs,abserr)
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qk31 ( f, a, b, result, abserr, resabs, resasc )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QK31 carries out a 31 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !                       result is computed by applying the 31-point
    !                       Gauss-Kronrod rule (resk), obtained by optimal
    !                       addition of abscissae to the 15-point Gauss
    !                       rule (resg).
    !
    !    Output, real(RK) ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(RK) RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(RK) RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
      implicit none

      real(RK) a
      real(RK) absc
      real(RK) abserr
      real(RK) b
      real(RK) centr
      real(RK) dhlgth
      real(RK), external :: f
      real(RK) fc
      real(RK) fsum
      real(RK) fval1
      real(RK) fval2
      real(RK) fv1(15)
      real(RK) fv2(15)
      real(RK) hlgth
      integer(IK) j
      integer(IK) jtw
      integer(IK) jtwm1
      real(RK) resabs
      real(RK) resasc
      real(RK) resg
      real(RK) resk
      real(RK) reskh
      real(RK) result
      real(RK) wg(8)
      real(RK) wgk(16)
      real(RK) xgk(16)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 31-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 15-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 15-point Gauss rule
    !
    !           wgk    - weights of the 31-point Kronrod rule
    !
    !           wg     - weights of the 15-point Gauss rule
    !
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
        xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16)/ &
           9.980022986933971E-01_RK,   9.879925180204854E-01_RK, &
           9.677390756791391E-01_RK,   9.372733924007059E-01_RK, &
           8.972645323440819E-01_RK,   8.482065834104272E-01_RK, &
           7.904185014424659E-01_RK,   7.244177313601700E-01_RK, &
           6.509967412974170E-01_RK,   5.709721726085388E-01_RK, &
           4.850818636402397E-01_RK,   3.941513470775634E-01_RK, &
           2.991800071531688E-01_RK,   2.011940939974345E-01_RK, &
           1.011420669187175E-01_RK,   0._RK               /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
        wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16)/ &
           5.377479872923349E-03_RK,   1.500794732931612E-02_RK, &
           2.546084732671532E-02_RK,   3.534636079137585E-02_RK, &
           4.458975132476488E-02_RK,   5.348152469092809E-02_RK, &
           6.200956780067064E-02_RK,   6.985412131872826E-02_RK, &
           7.684968075772038E-02_RK,   8.308050282313302E-02_RK, &
           8.856444305621177E-02_RK,   9.312659817082532E-02_RK, &
           9.664272698362368E-02_RK,   9.917359872179196E-02_RK, &
           1.007698455238756E-01_RK,   1.013300070147915E-01_RK/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
           3.075324199611727E-02_RK,   7.036604748810812E-02_RK, &
           1.071592204671719E-01_RK,   1.395706779261543E-01_RK, &
           1.662692058169939E-01_RK,   1.861610000155622E-01_RK, &
           1.984314853271116E-01_RK,   2.025782419255613E-01_RK/
    !
    !
    !           list of major variables
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 15-point Gauss formula
    !           resk   - result of the 31-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
      centr = 5.0E-01_RK*(a+b)
      hlgth = 5.0E-01_RK*(b-a)
      dhlgth = abs(hlgth)
    !
    !  Compute the 31-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
      fc = f(centr)
      resg = wg(8)*fc
      resk = wgk(16)*fc
      resabs = abs(resk)

      do j = 1, 7
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
      end do

      do j = 1, 8
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
      end do

      reskh = resk*5.0E-01
      resasc = wgk(16)*abs(fc-reskh)

      do j = 1, 15
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
     end do

      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)

      if ( resasc /= 0._RK .and. abserr /= 0._RK) &
        abserr = resasc*min ( 1._RK,(2.0E+02_RK*abserr/resasc)**1.5E+00_RK)

      if ( resabs > tiny ( resabs ) /(5.0E+01_RK* epsilon ( resabs ) )) then
        abserr = max (( epsilon ( resabs ) *5.0E+01_RK)*resabs,abserr)
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qk41 ( f, a, b, result, abserr, resabs, resasc )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QK41 carries out a 41 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !                       result is computed by applying the 41-point
    !                       Gauss-Kronrod rule (resk) obtained by optimal
    !                       addition of abscissae to the 20-point Gauss
    !                       rule (resg).
    !
    !    Output, real(RK) ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(RK) RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(RK) RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 20-point Gauss formula
    !           resk   - result of the 41-point Kronrod formula
    !           reskh  - approximation to mean value of f over (a,b), i.e.
    !                    to i/(b-a)
    !
      implicit none

      real(RK) a
      real(RK) absc
      real(RK) abserr
      real(RK) b
      real(RK) centr
      real(RK) dhlgth
      real(RK), external :: f
      real(RK) fc
      real(RK) fsum
      real(RK) fval1
      real(RK) fval2
      real(RK) fv1(20)
      real(RK) fv2(20)
      real(RK) hlgth
      integer(IK) j
      integer(IK) jtw
      integer(IK) jtwm1
      real(RK) resabs
      real(RK) resasc
      real(RK) resg
      real(RK) resk
      real(RK) reskh
      real(RK) result
      real(RK) wg(10)
      real(RK) wgk(21)
      real(RK) xgk(21)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 41-point Gauss-Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 20-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 20-point Gauss rule
    !
    !           wgk    - weights of the 41-point Gauss-Kronrod rule
    !
    !           wg     - weights of the 20-point Gauss rule
    !
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
        xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16), &
        xgk(17),xgk(18),xgk(19),xgk(20),xgk(21)/ &
           9.988590315882777E-01_RK,   9.931285991850949E-01_RK, &
           9.815078774502503E-01_RK,   9.639719272779138E-01_RK, &
           9.408226338317548E-01_RK,   9.122344282513259E-01_RK, &
           8.782768112522820E-01_RK,   8.391169718222188E-01_RK, &
           7.950414288375512E-01_RK,   7.463319064601508E-01_RK, &
           6.932376563347514E-01_RK,   6.360536807265150E-01_RK, &
           5.751404468197103E-01_RK,   5.108670019508271E-01_RK, &
           4.435931752387251E-01_RK,   3.737060887154196E-01_RK, &
           3.016278681149130E-01_RK,   2.277858511416451E-01_RK, &
           1.526054652409227E-01_RK,   7.652652113349733E-02_RK, &
           0._RK               /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
        wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16), &
        wgk(17),wgk(18),wgk(19),wgk(20),wgk(21)/ &
           3.073583718520532E-03_RK,   8.600269855642942E-03_RK, &
           1.462616925697125E-02_RK,   2.038837346126652E-02_RK, &
           2.588213360495116E-02_RK,   3.128730677703280E-02_RK, &
           3.660016975820080E-02_RK,   4.166887332797369E-02_RK, &
           4.643482186749767E-02_RK,   5.094457392372869E-02_RK, &
           5.519510534828599E-02_RK,   5.911140088063957E-02_RK, &
           6.265323755478117E-02_RK,   6.583459713361842E-02_RK, &
           6.864867292852162E-02_RK,   7.105442355344407E-02_RK, &
           7.303069033278667E-02_RK,   7.458287540049919E-02_RK, &
           7.570449768455667E-02_RK,   7.637786767208074E-02_RK, &
           7.660071191799966E-02_RK/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10)/ &
           1.761400713915212E-02_RK,   4.060142980038694E-02_RK, &
           6.267204833410906E-02_RK,   8.327674157670475E-02_RK, &
           1.019301198172404E-01_RK,   1.181945319615184E-01_RK, &
           1.316886384491766E-01_RK,   1.420961093183821E-01_RK, &
           1.491729864726037E-01_RK,   1.527533871307259E-01_RK/
    !
      centr = 5.0E-01_RK*(a+b)
      hlgth = 5.0E-01_RK*(b-a)
      dhlgth = abs(hlgth)
    !
    !  Compute 41-point Gauss-Kronrod approximation to the
    !  the integral, and estimate the absolute error.
    !
      resg = 0._RK
      fc = f(centr)
      resk = wgk(21)*fc
      resabs = abs(resk)

      do j = 1, 10
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
      end do

      do j = 1, 10
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
      end do

      reskh = resk*5.0E-01_RK
      resasc = wgk(21)*abs(fc-reskh)

      do j = 1, 20
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      end do

      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)

      if ( resasc /= 0._RK .and. abserr /= 0._RK) abserr = resasc*min ( 1._RK,(2.0E+02_RK*abserr/resasc)**1.5E+00_RK)

      if ( resabs > tiny ( resabs ) /(5.0E+01_RK* epsilon ( resabs ) )) then
        abserr = max (( epsilon ( resabs ) *5.0E+01_RK)*resabs,abserr)
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qk51 ( f, a, b, result, abserr, resabs, resasc )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QK51 carries out a 51 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !                       result is computed by applying the 51-point
    !                       Kronrod rule (resk) obtained by optimal addition
    !                       of abscissae to the 25-point Gauss rule (resg).
    !
    !    Output, real(RK) ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(RK) RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(RK) RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 25-point Gauss formula
    !           resk   - result of the 51-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
      implicit none

      real(RK) a
      real(RK) absc
      real(RK) abserr
      real(RK) b
      real(RK) centr
      real(RK) dhlgth
      real(RK), external :: f
      real(RK) fc
      real(RK) fsum
      real(RK) fval1
      real(RK) fval2
      real(RK) fv1(25)
      real(RK) fv2(25)
      real(RK) hlgth
      integer(IK) j
      integer(IK) jtw
      integer(IK) jtwm1
      real(RK) resabs
      real(RK) resasc
      real(RK) resg
      real(RK) resk
      real(RK) reskh
      real(RK) result
      real(RK) wg(13)
      real(RK) wgk(26)
      real(RK) xgk(26)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 51-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 25-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 25-point Gauss rule
    !
    !           wgk    - weights of the 51-point Kronrod rule
    !
    !           wg     - weights of the 25-point Gauss rule
    !
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
        xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14)/ &
           9.992621049926098E-01_RK,   9.955569697904981E-01_RK, &
           9.880357945340772E-01_RK,   9.766639214595175E-01_RK, &
           9.616149864258425E-01_RK,   9.429745712289743E-01_RK, &
           9.207471152817016E-01_RK,   8.949919978782754E-01_RK, &
           8.658470652932756E-01_RK,   8.334426287608340E-01_RK, &
           7.978737979985001E-01_RK,   7.592592630373576E-01_RK, &
           7.177664068130844E-01_RK,   6.735663684734684E-01_RK/
       data xgk(15),xgk(16),xgk(17),xgk(18),xgk(19),xgk(20),xgk(21), &
        xgk(22),xgk(23),xgk(24),xgk(25),xgk(26)/ &
           6.268100990103174E-01_RK,   5.776629302412230E-01_RK, &
           5.263252843347192E-01_RK,   4.730027314457150E-01_RK, &
           4.178853821930377E-01_RK,   3.611723058093878E-01_RK, &
           3.030895389311078E-01_RK,   2.438668837209884E-01_RK, &
           1.837189394210489E-01_RK,   1.228646926107104E-01_RK, &
           6.154448300568508E-02_RK,   0._RK               /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
        wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14)/ &
           1.987383892330316E-03_RK,   5.561932135356714E-03_RK, &
           9.473973386174152E-03_RK,   1.323622919557167E-02_RK, &
           1.684781770912830E-02_RK,   2.043537114588284E-02_RK, &
           2.400994560695322E-02_RK,   2.747531758785174E-02_RK, &
           3.079230016738749E-02_RK,   3.400213027432934E-02_RK, &
           3.711627148341554E-02_RK,   4.008382550403238E-02_RK, &
           4.287284502017005E-02_RK,   4.550291304992179E-02_RK/
       data wgk(15),wgk(16),wgk(17),wgk(18),wgk(19),wgk(20),wgk(21), &
        wgk(22),wgk(23),wgk(24),wgk(25),wgk(26)/ &
           4.798253713883671E-02_RK,   5.027767908071567E-02_RK, &
           5.236288580640748E-02_RK,   5.425112988854549E-02_RK, &
           5.595081122041232E-02_RK,   5.743711636156783E-02_RK, &
           5.868968002239421E-02_RK,   5.972034032417406E-02_RK, &
           6.053945537604586E-02_RK,   6.112850971705305E-02_RK, &
           6.147118987142532E-02_RK,   6.158081806783294E-02_RK/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10), &
        wg(11),wg(12),wg(13)/ &
           1.139379850102629E-02_RK,   2.635498661503214E-02_RK, &
           4.093915670130631E-02_RK,   5.490469597583519E-02_RK, &
           6.803833381235692E-02_RK,   8.014070033500102E-02_RK, &
           9.102826198296365E-02_RK,   1.005359490670506E-01_RK, &
           1.085196244742637E-01_RK,   1.148582591457116E-01_RK, &
           1.194557635357848E-01_RK,   1.222424429903100E-01_RK, &
           1.231760537267155E-01_RK/
    !
      centr = 5.0E-01_RK*(a+b)
      hlgth = 5.0E-01_RK*(b-a)
      dhlgth = abs(hlgth)
    !
    !  Compute the 51-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
      fc = f(centr)
      resg = wg(13)*fc
      resk = wgk(26)*fc
      resabs = abs(resk)

      do j = 1, 12
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
      end do

      do j = 1, 13
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
      end do

      reskh = resk*5.0E-01_RK
      resasc = wgk(26)*abs(fc-reskh)

      do j = 1, 25
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      end do

      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)

      if ( resasc /= 0._RK .and. abserr /= 0._RK) then
        abserr = resasc*min ( 1._RK,(2.0E+02_RK*abserr/resasc)**1.5E+00_RK)
      end if

      if ( resabs > tiny ( resabs ) / (5.0E+01_RK* epsilon ( resabs ) ) ) then
        abserr = max (( epsilon ( resabs ) *5.0E+01_RK)*resabs,abserr)
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qk61 ( f, a, b, result, abserr, resabs, resasc ) 

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QK61 carries out a 61 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !                    result is computed by applying the 61-point
    !                    Kronrod rule (resk) obtained by optimal addition of
    !                    abscissae to the 30-point Gauss rule (resg).
    !
    !    Output, real(RK) ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(RK) RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(RK) RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 30-point Gauss rule
    !           resk   - result of the 61-point Kronrod rule
    !           reskh  - approximation to the mean value of f
    !                    over (a,b), i.e. to i/(b-a)
    !
      implicit none

      real(RK) a
      real(RK) absc
      real(RK) abserr
      real(RK) b
      real(RK) centr
      real(RK) dhlgth
      real(RK), external :: f
      real(RK) fc
      real(RK) fsum
      real(RK) fval1
      real(RK) fval2
      real(RK) fv1(30)
      real(RK) fv2(30)
      real(RK) hlgth
      integer(IK) j
      integer(IK) jtw
      integer(IK) jtwm1
      real(RK) resabs
      real(RK) resasc
      real(RK) resg
      real(RK) resk
      real(RK) reskh
      real(RK) result
      real(RK) wg(15)
      real(RK) wgk(31)
      real(RK) xgk(31)
    !
    !           the abscissae and weights are given for the
    !           interval (-1,1). because of symmetry only the positive
    !           abscissae and their corresponding weights are given.
    !
    !           xgk   - abscissae of the 61-point Kronrod rule
    !                   xgk(2), xgk(4)  ... abscissae of the 30-point
    !                   Gauss rule
    !                   xgk(1), xgk(3)  ... optimally added abscissae
    !                   to the 30-point Gauss rule
    !
    !           wgk   - weights of the 61-point Kronrod rule
    !
    !           wg    - weigths of the 30-point Gauss rule
    !
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10)/ &
           9.994844100504906E-01_RK,     9.968934840746495E-01_RK, &
           9.916309968704046E-01_RK,     9.836681232797472E-01_RK, &
           9.731163225011263E-01_RK,     9.600218649683075E-01_RK, &
           9.443744447485600E-01_RK,     9.262000474292743E-01_RK, &
           9.055733076999078E-01_RK,     8.825605357920527E-01_RK/
      data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),xgk(17), &
        xgk(18),xgk(19),xgk(20)/ &
           8.572052335460611E-01_RK,     8.295657623827684E-01_RK, &
           7.997278358218391E-01_RK,     7.677774321048262E-01_RK, &
           7.337900624532268E-01_RK,     6.978504947933158E-01_RK, &
           6.600610641266270E-01_RK,     6.205261829892429E-01_RK, &
           5.793452358263617E-01_RK,     5.366241481420199E-01_RK/
      data xgk(21),xgk(22),xgk(23),xgk(24),xgk(25),xgk(26),xgk(27), &
        xgk(28),xgk(29),xgk(30),xgk(31)/ &
           4.924804678617786E-01_RK,     4.470337695380892E-01_RK, &
           4.004012548303944E-01_RK,     3.527047255308781E-01_RK, &
           3.040732022736251E-01_RK,     2.546369261678898E-01_RK, &
           2.045251166823099E-01_RK,     1.538699136085835E-01_RK, &
           1.028069379667370E-01_RK,     5.147184255531770E-02_RK, &
           0._RK                   /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
        wgk(9),wgk(10)/ &
           1.389013698677008E-03_RK,     3.890461127099884E-03_RK, &
           6.630703915931292E-03_RK,     9.273279659517763E-03_RK, &
           1.182301525349634E-02_RK,     1.436972950704580E-02_RK, &
           1.692088918905327E-02_RK,     1.941414119394238E-02_RK, &
           2.182803582160919E-02_RK,     2.419116207808060E-02_RK/
      data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),wgk(17), &
        wgk(18),wgk(19),wgk(20)/ &
           2.650995488233310E-02_RK,     2.875404876504129E-02_RK, &
           3.090725756238776E-02_RK,     3.298144705748373E-02_RK, &
           3.497933802806002E-02_RK,     3.688236465182123E-02_RK, &
           3.867894562472759E-02_RK,     4.037453895153596E-02_RK, &
           4.196981021516425E-02_RK,     4.345253970135607E-02_RK/
      data wgk(21),wgk(22),wgk(23),wgk(24),wgk(25),wgk(26),wgk(27), &
        wgk(28),wgk(29),wgk(30),wgk(31)/ &
           4.481480013316266E-02_RK,     4.605923827100699E-02_RK, &
           4.718554656929915E-02_RK,     4.818586175708713E-02_RK, &
           4.905543455502978E-02_RK,     4.979568342707421E-02_RK, &
           5.040592140278235E-02_RK,     5.088179589874961E-02_RK, &
           5.122154784925877E-02_RK,     5.142612853745903E-02_RK, &
           5.149472942945157E-02_RK/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
           7.968192496166606E-03_RK,     1.846646831109096E-02_RK, &
           2.878470788332337E-02_RK,     3.879919256962705E-02_RK, &
           4.840267283059405E-02_RK,     5.749315621761907E-02_RK, &
           6.597422988218050E-02_RK,     7.375597473770521E-02_RK/
      data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/ &
           8.075589522942022E-02_RK,     8.689978720108298E-02_RK, &
           9.212252223778613E-02_RK,     9.636873717464426E-02_RK, &
           9.959342058679527E-02_RK,     1.017623897484055E-01_RK, &
           1.028526528935588E-01_RK/

      centr = 5.0E-01_RK*(b+a)
      hlgth = 5.0E-01_RK*(b-a)
      dhlgth = abs(hlgth)
    !
    !  Compute the 61-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
      resg = 0._RK
      fc = f(centr)
      resk = wgk(31)*fc
      resabs = abs(resk)

      do j = 1, 15
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
      end do

      do j = 1, 15
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
      end do

      reskh = resk * 5.0E-01_RK
      resasc = wgk(31)*abs(fc-reskh)

      do j = 1, 30
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
      end do

      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)

      if ( resasc /= 0._RK .and. abserr /= 0._RK) then
        abserr = resasc*min ( 1._RK,(2.0E+02_RK*abserr/resasc)**1.5E+00_RK)
      end if

      if ( resabs > tiny ( resabs ) / (5.0E+01_RK* epsilon ( resabs ) )) then
        abserr = max ( ( epsilon ( resabs ) *5.0E+01_RK)*resabs, abserr )
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qmomo ( alfa, beta, ri, rj, rg, rh, integr )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QMOMO computes modified Chebyshev moments.
    !
    !  Discussion:
    !
    !    This routine computes modified Chebyshev moments.
    !    The K-th modified Chebyshev moment is defined as the
    !    integral over (-1,1) of W(X)*T(K,X), where T(K,X) is the
    !    Chebyshev polynomial of degree K.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(RK) ALFA, a parameter in the weight function w(x), ALFA > -1.
    !
    !    Input, real(RK) BETA, a parameter in the weight function w(x), BETA > -1.
    !
    !           ri     - real
    !                    vector of dimension 25
    !                    ri(k) is the integral over (-1,1) of
    !                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
    !
    !           rj     - real
    !                    vector of dimension 25
    !                    rj(k) is the integral over (-1,1) of
    !                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
    !
    !           rg     - real
    !                    vector of dimension 25
    !                    rg(k) is the integral over (-1,1) of
    !                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ...,25.
    !
    !           rh     - real
    !                    vector of dimension 25
    !                    rh(k) is the integral over (-1,1) of
    !                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
    !
    !           integr - integer(IK)
    !                    input parameter indicating the modified moments
    !                    to be computed
    !                    integr = 1 compute ri, rj
    !                           = 2 compute ri, rj, rg
    !                           = 3 compute ri, rj, rh
    !                           = 4 compute ri, rj, rg, rh
    !
      implicit none

      real(RK) alfa
      real(RK) alfp1
      real(RK) alfp2
      real(RK) an
      real(RK) anm1
      real(RK) beta
      real(RK) betp1
      real(RK) betp2
      integer(IK) i
      integer(IK) im1
      integer(IK) integr
      real(RK) ralf
      real(RK) rbet
      real(RK) rg(25)
      real(RK) rh(25)
      real(RK) ri(25)
      real(RK) rj(25)
    !
      alfp1 = alfa+1._RK
      betp1 = beta+1._RK
      alfp2 = alfa+2._RK
      betp2 = beta+2._RK
      ralf = 2._RK**alfp1
      rbet = 2._RK**betp1
    !
    !  Compute RI, RJ using a forward recurrence relation.
    !
      ri(1) = ralf/alfp1
      rj(1) = rbet/betp1
      ri(2) = ri(1)*alfa/alfp2
      rj(2) = rj(1)*beta/betp2
      an = 2._RK
      anm1 = 1._RK

      do i = 3, 25
        ri(i) = -(ralf+an*(an-alfp2)*ri(i-1))/(anm1*(an+alfp1))
        rj(i) = -(rbet+an*(an-betp2)*rj(i-1))/(anm1*(an+betp1))
        anm1 = an
        an = an+1._RK
      end do

      if ( integr == 1 ) go to 70
      if ( integr == 3 ) go to 40
    !
    !  Compute RG using a forward recurrence relation.
    !
      rg(1) = -ri(1)/alfp1
      rg(2) = -(ralf+ralf)/(alfp2*alfp2)-rg(1)
      an = 2._RK
      anm1 = 1._RK
      im1 = 2

      do i = 3, 25
        rg(i) = -(an*(an-alfp2)*rg(im1)-an*ri(im1)+anm1*ri(i))/ &
        (anm1*(an+alfp1))
        anm1 = an
        an = an+1._RK
        im1 = i
      end do

      if ( integr == 2 ) go to 70
    !
    !  Compute RH using a forward recurrence relation.
    !
    40 continue

      rh(1) = -rj(1) / betp1
      rh(2) = -(rbet+rbet)/(betp2*betp2)-rh(1)
      an = 2._RK
      anm1 = 1._RK
      im1 = 2

      do i = 3, 25
        rh(i) = -(an*(an-betp2)*rh(im1)-an*rj(im1)+ &
        anm1*rj(i))/(anm1*(an+betp1))
        anm1 = an
        an = an+1._RK
        im1 = i
      end do

      do i = 2, 25, 2
        rh(i) = -rh(i)
      end do

       70 continue

      do i = 2, 25, 2
        rj(i) = -rj(i)
      end do

    !  90 continue

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qng ( f, a, b, epsabs, epsrel, result, abserr, neval, ier )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QNG estimates an integral, using non-adaptive integration.
    !
    !  Discussion:
    !
    !    The routine calculates an approximation RESULT to a definite integral   
    !      I = integral of F over (A,B),
    !    hopefully satisfying
    !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
    !
    !    The routine is a simple non-adaptive automatic integrator, based on
    !    a sequence of rules with increasing degree of algebraic
    !    precision (Patterson, 1968).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(RK) F, the name of the function routine, of the form
    !      function f ( x )
    !      real(RK) f
    !      real(RK) x
    !    which evaluates the integrand function.
    !
    !    Input, real(RK) A, B, the limits of integration.
    !
    !    Input, real(RK) EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Output, real(RK) RESULT, the estimated value of the integral.
    !    RESULT is obtained by applying the 21-point Gauss-Kronrod rule (RES21)
    !    obtained  by optimal addition of abscissae to the 10-point Gauss rule
    !    (RES10), or by applying the 43-point rule (RES43) obtained by optimal
    !    addition of abscissae to the 21-point Gauss-Kronrod rule, or by 
    !    applying the 87-point rule (RES87) obtained by optimal addition of
    !    abscissae to the 43-point rule.
    !
    !    Output, real(RK) ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, integer(IK) NEVAL, the number of times the integral was evaluated.
    !
    !           ier    - ier = 0 normal and reliable termination of the
    !                            routine. it is assumed that the requested
    !                            accuracy has been achieved.
    !                    ier > 0 abnormal termination of the routine. it is
    !                            assumed that the requested accuracy has
    !                            not been achieved.
    !                    ier = 1 the maximum number of steps has been
    !                            executed. the integral is probably too
    !                            difficult to be calculated by qng.
    !                        = 6 the input is invalid, because
    !                            epsabs < 0 and epsrel < 0,
    !                            result, abserr and neval are set to zero.
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the integration interval
    !           hlgth  - half-length of the integration interval
    !           fcentr - function value at mid point
    !           absc   - abscissa
    !           fval   - function value
    !           savfun - array of function values which have already
    !                    been computed
    !           res10  - 10-point Gauss result
    !           res21  - 21-point Kronrod result
    !           res43  - 43-point result
    !           res87  - 87-point result
    !           resabs - approximation to the integral of abs(f)
    !           resasc - approximation to the integral of abs(f-i/(b-a))
    !
      implicit none

      real(RK) a
      real(RK) absc
      real(RK) abserr
      real(RK) b
      real(RK) centr
      real(RK) dhlgth
      real(RK) epsabs
      real(RK) epsrel
      real(RK), external :: f
      real(RK) fcentr
      real(RK) fval
      real(RK) fval1
      real(RK) fval2
      real(RK) fv1(5)
      real(RK) fv2(5)
      real(RK) fv3(5)
      real(RK) fv4(5)
      real(RK) hlgth
      integer(IK) ier
      integer(IK) ipx
      integer(IK) k
      integer(IK) l
      integer(IK) neval
      real(RK) result
      real(RK) res10
      real(RK) res21
      real(RK) res43
      real(RK) res87
      real(RK) resabs
      real(RK) resasc
      real(RK) reskh
      real(RK) savfun(21)
      real(RK) w10(5)
      real(RK) w21a(5)
      real(RK) w21b(6)
      real(RK) w43a(10)
      real(RK) w43b(12)
      real(RK) w87a(21)
      real(RK) w87b(23)
      real(RK) x1(5)
      real(RK) x2(5)
      real(RK) x3(11)
      real(RK) x4(22)
    !
    !           the following data statements contain the abscissae
    !           and weights of the integration rules used.
    !
    !           x1      abscissae common to the 10-, 21-, 43- and 87-point
    !                   rule
    !           x2      abscissae common to the 21-, 43- and 87-point rule
    !           x3      abscissae common to the 43- and 87-point rule
    !           x4      abscissae of the 87-point rule
    !           w10     weights of the 10-point formula
    !           w21a    weights of the 21-point formula for abscissae x1
    !           w21b    weights of the 21-point formula for abscissae x2
    !           w43a    weights of the 43-point formula for absissae x1, x3
    !           w43b    weights of the 43-point formula for abscissae x3
    !           w87a    weights of the 87-point formula for abscissae x1,
    !                   x2 and x3
    !           w87b    weights of the 87-point formula for abscissae x4
    !
      data x1(1),x1(2),x1(3),x1(4),x1(5)/ &
           9.739065285171717E-01_RK,     8.650633666889845E-01_RK, &
           6.794095682990244E-01_RK,     4.333953941292472E-01_RK, &
           1.488743389816312E-01_RK/
      data x2(1),x2(2),x2(3),x2(4),x2(5)/ &
           9.956571630258081E-01_RK,     9.301574913557082E-01_RK, &
           7.808177265864169E-01_RK,     5.627571346686047E-01_RK, &
           2.943928627014602E-01_RK/
      data x3(1),x3(2),x3(3),x3(4),x3(5),x3(6),x3(7),x3(8),x3(9),x3(10), &
        x3(11)/ &
           9.993333609019321E-01_RK,     9.874334029080889E-01_RK, &
           9.548079348142663E-01_RK,     9.001486957483283E-01_RK, &
           8.251983149831142E-01_RK,     7.321483889893050E-01_RK, &
           6.228479705377252E-01_RK,     4.994795740710565E-01_RK, &
           3.649016613465808E-01_RK,     2.222549197766013E-01_RK, &
           7.465061746138332E-02_RK/
      data x4(1),x4(2),x4(3),x4(4),x4(5),x4(6),x4(7),x4(8),x4(9),x4(10), &
        x4(11),x4(12),x4(13),x4(14),x4(15),x4(16),x4(17),x4(18),x4(19), &
        x4(20),x4(21),x4(22)/            9.999029772627292E-01_RK, &
           9.979898959866787E-01_RK,     9.921754978606872E-01_RK, &
           9.813581635727128E-01_RK,     9.650576238583846E-01_RK, &
           9.431676131336706E-01_RK,     9.158064146855072E-01_RK, &
           8.832216577713165E-01_RK,     8.457107484624157E-01_RK, &
           8.035576580352310E-01_RK,     7.570057306854956E-01_RK, &
           7.062732097873218E-01_RK,     6.515894665011779E-01_RK, &
           5.932233740579611E-01_RK,     5.314936059708319E-01_RK, &
           4.667636230420228E-01_RK,     3.994248478592188E-01_RK, &
           3.298748771061883E-01_RK,     2.585035592021616E-01_RK, &
           1.856953965683467E-01_RK,     1.118422131799075E-01_RK, &
           3.735212339461987E-02_RK/
      data w10(1),w10(2),w10(3),w10(4),w10(5)/ &
           6.667134430868814E-02_RK,     1.494513491505806E-01_RK, &
           2.190863625159820E-01_RK,     2.692667193099964E-01_RK, &
           2.955242247147529E-01_RK/
      data w21a(1),w21a(2),w21a(3),w21a(4),w21a(5)/ &
           3.255816230796473E-02_RK,     7.503967481091995E-02_RK, &
           1.093871588022976E-01_RK,     1.347092173114733E-01_RK, &
           1.477391049013385E-01_RK/
      data w21b(1),w21b(2),w21b(3),w21b(4),w21b(5),w21b(6)/ &
           1.169463886737187E-02_RK,     5.475589657435200E-02_RK, &
           9.312545458369761E-02_RK,     1.234919762620659E-01_RK, &
           1.427759385770601E-01_RK,     1.494455540029169E-01_RK/
      data w43a(1),w43a(2),w43a(3),w43a(4),w43a(5),w43a(6),w43a(7), &
        w43a(8),w43a(9),w43a(10)/     1.629673428966656E-02_RK, &
           3.752287612086950E-02_RK,     5.469490205825544E-02_RK, &
           6.735541460947809E-02_RK,     7.387019963239395E-02_RK, &
           5.768556059769796E-03_RK,     2.737189059324884E-02_RK, &
           4.656082691042883E-02_RK,     6.174499520144256E-02_RK, &
           7.138726726869340E-02_RK/
      data w43b(1),w43b(2),w43b(3),w43b(4),w43b(5),w43b(6),w43b(7), &
        w43b(8),w43b(9),w43b(10),w43b(11),w43b(12)/ &
           1.844477640212414E-03_RK,     1.079868958589165E-02_RK, &
           2.189536386779543E-02_RK,     3.259746397534569E-02_RK, &
           4.216313793519181E-02_RK,     5.074193960018458E-02_RK, &
           5.837939554261925E-02_RK,     6.474640495144589E-02_RK, &
           6.956619791235648E-02_RK,     7.282444147183321E-02_RK, &
           7.450775101417512E-02_RK,     7.472214751740301E-02_RK/
      data w87a(1),w87a(2),w87a(3),w87a(4),w87a(5),w87a(6),w87a(7), &
        w87a(8),w87a(9),w87a(10),w87a(11),w87a(12),w87a(13),w87a(14), &
        w87a(15),w87a(16),w87a(17),w87a(18),w87a(19),w87a(20),w87a(21)/ &
           8.148377384149173E-03_RK,     1.876143820156282E-02_RK, &
           2.734745105005229E-02_RK,     3.367770731163793E-02_RK, &
           3.693509982042791E-02_RK,     2.884872430211531E-03_RK, &
           1.368594602271270E-02_RK,     2.328041350288831E-02_RK, &
           3.087249761171336E-02_RK,     3.569363363941877E-02_RK, &
           9.152833452022414E-04_RK,     5.399280219300471E-03_RK, &
           1.094767960111893E-02_RK,     1.629873169678734E-02_RK, &
           2.108156888920384E-02_RK,     2.537096976925383E-02_RK, &
           2.918969775647575E-02_RK,     3.237320246720279E-02_RK, &
           3.478309895036514E-02_RK,     3.641222073135179E-02_RK, &
           3.725387550304771E-02_RK/
      data w87b(1),w87b(2),w87b(3),w87b(4),w87b(5),w87b(6),w87b(7), &
        w87b(8),w87b(9),w87b(10),w87b(11),w87b(12),w87b(13),w87b(14), &
        w87b(15),w87b(16),w87b(17),w87b(18),w87b(19),w87b(20),w87b(21), &
        w87b(22),w87b(23)/               2.741455637620724E-04_RK, &
           1.807124155057943E-03_RK,     4.096869282759165E-03_RK, &
           6.758290051847379E-03_RK,     9.549957672201647E-03_RK, &
           1.232944765224485E-02_RK,     1.501044734638895E-02_RK, &
           1.754896798624319E-02_RK,     1.993803778644089E-02_RK, &
           2.219493596101229E-02_RK,     2.433914712600081E-02_RK, &
           2.637450541483921E-02_RK,     2.828691078877120E-02_RK, &
           3.005258112809270E-02_RK,     3.164675137143993E-02_RK, &
           3.305041341997850E-02_RK,     3.425509970422606E-02_RK, &
           3.526241266015668E-02_RK,     3.607698962288870E-02_RK, &
           3.669860449845609E-02_RK,     3.712054926983258E-02_RK, &
           3.733422875193504E-02_RK,     3.736107376267902E-02_RK/
    !
    !  Test on validity of parameters.
    !
      result = 0._RK
      abserr = 0._RK
      neval = 0

      if ( epsabs < 0._RK .and. epsrel < 0._RK ) then
        ier = 6
        return
      end if

      hlgth = 5.0E-01_RK * ( b - a )
      dhlgth = abs ( hlgth )
      centr = 5.0E-01_RK * ( b + a )
      fcentr = f(centr)
      neval = 21
      ier = 1
    !
    !  Compute the integral using the 10- and 21-point formula.
    !
      do l = 1, 3

        if ( l == 1 ) then

          res10 = 0._RK
          res21 = w21b(6) * fcentr
          resabs = w21b(6) * abs(fcentr)

          do k = 1, 5
            absc = hlgth * x1(k)
            fval1 = f(centr+absc)
            fval2 = f(centr-absc)
            fval = fval1 + fval2
            res10 = res10 + w10(k)*fval
            res21 = res21 + w21a(k)*fval
            resabs = resabs + w21a(k)*(abs(fval1)+abs(fval2))
            savfun(k) = fval
            fv1(k) = fval1
            fv2(k) = fval2
          end do

          ipx = 5

          do k = 1, 5
            ipx = ipx + 1
            absc = hlgth * x2(k)
            fval1 = f(centr+absc)
            fval2 = f(centr-absc)
            fval = fval1 + fval2
            res21 = res21 + w21b(k) * fval
            resabs = resabs + w21b(k) * ( abs ( fval1 ) + abs ( fval2 ) )
            savfun(ipx) = fval
            fv3(k) = fval1
            fv4(k) = fval2
          end do
    !
    !  Test for convergence.
    !
          result = res21 * hlgth
          resabs = resabs * dhlgth
          reskh = 5.0E-01_RK * res21
          resasc = w21b(6) * abs ( fcentr - reskh )

          do k = 1, 5
            resasc = resasc+w21a(k)*(abs(fv1(k)-reskh)+abs(fv2(k)-reskh))+w21b(k)*(abs(fv3(k)-reskh)+abs(fv4(k)-reskh))
          end do

          abserr = abs ( ( res21 - res10 ) * hlgth )
          resasc = resasc * dhlgth
    !
    !  Compute the integral using the 43-point formula.
    !
        else if ( l == 2 ) then

          res43 = w43b(12)*fcentr
          neval = 43

          do k = 1, 10
            res43 = res43 + savfun(k) * w43a(k)
          end do

          do k = 1, 11
            ipx = ipx + 1
            absc = hlgth * x3(k)
            fval = f(absc+centr) + f(centr-absc)
            res43 = res43 + fval * w43b(k)
            savfun(ipx) = fval
          end do
    !
    !  Test for convergence.
    !
          result = res43 * hlgth
          abserr = abs((res43-res21)*hlgth)
    !
    !  Compute the integral using the 87-point formula.
    !
        else if ( l == 3 ) then

          res87 = w87b(23) * fcentr
          neval = 87

          do k = 1, 21
            res87 = res87 + savfun(k) * w87a(k)
          end do

          do k = 1, 22
            absc = hlgth * x4(k)
            res87 = res87 + w87b(k) * ( f(absc+centr) + f(centr-absc) )
          end do

          result = res87 * hlgth
          abserr = abs ( ( res87 - res43) * hlgth )

        end if

        if ( resasc /= 0._RK .and. abserr /= 0._RK ) then
          abserr = resasc * min ( 1._RK,(2.0E+02_RK*abserr/resasc)**1.5E+00_RK)
        end if

        if ( resabs > tiny ( resabs ) / ( 5.0E+01_RK * epsilon ( resabs ) ) ) then
          abserr = max (( epsilon ( resabs ) *5.0E+01_RK) * resabs, abserr )
        end if

        if ( abserr <= max ( epsabs, epsrel*abs(result))) then
          ier = 0
        end if

        if ( ier == 0 ) then
          exit
        end if

      end do

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine qsort ( limit, last, maxerr, ermax, elist, iord, nrmax )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QSORT maintains the order of a list of local error estimates.
    !
    !  Discussion:
    !
    !    This routine maintains the descending ordering in the list of the 
    !    local error estimates resulting from the interval subdivision process. 
    !    At each call two error estimates are inserted using the sequential 
    !    search top-down for the largest error estimate and bottom-up for the
    !    smallest error estimate.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, integer(IK) LIMIT, the maximum number of error estimates the list can
    !    contain.
    !
    !    Input, integer(IK) LAST, the current number of error estimates.
    !
    !    Input/output, integer(IK) MAXERR, the index in the list of the NRMAX-th 
    !    largest error.
    !
    !    Output, real(RK) ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
    !
    !    Input, real(RK) ELIST(LIMIT), contains the error estimates.
    !
    !    Input/output, integer(IK) IORD(LAST).  The first K elements contain 
    !    pointers to the error estimates such that ELIST(IORD(1)) through
    !    ELIST(IORD(K)) form a decreasing sequence, with
    !      K = LAST 
    !    if 
    !      LAST <= (LIMIT/2+2), 
    !    and otherwise
    !      K = LIMIT+1-LAST.
    !
    !    Input/output, integer(IK) NRMAX.
    !
      implicit none

      integer(IK) last

      real(RK) elist(last)
      real(RK) ermax
      real(RK) errmax
      real(RK) errmin
      integer(IK) i
      integer(IK) ibeg
      integer(IK) iord(last)
      integer(IK) isucc
      integer(IK) j
      integer(IK) jbnd
      integer(IK) jupbn
      integer(IK) k
      integer(IK) limit
      integer(IK) maxerr
      integer(IK) nrmax
    !
    !  Check whether the list contains more than two error estimates.
    !
      if ( last <= 2 ) then
        iord(1) = 1
        iord(2) = 2
        go to 90
      end if
    !
    !  This part of the routine is only executed if, due to a
    !  difficult integrand, subdivision increased the error
    !  estimate. in the normal case the insert procedure should
    !  start after the nrmax-th largest error estimate.
    !
      errmax = elist(maxerr)

      do i = 1, nrmax-1

        isucc = iord(nrmax-1)

        if ( errmax <= elist(isucc) ) then
          exit
        end if

        iord(nrmax) = isucc
        nrmax = nrmax-1

      end do
    !
    !  Compute the number of elements in the list to be maintained
    !  in descending order.  This number depends on the number of
    !  subdivisions still allowed.
    !
      jupbn = last

      if ( (limit/2+2) < last ) then
        jupbn = limit+3-last
      end if

      errmin = elist(last)
    !
    !  Insert errmax by traversing the list top-down, starting
    !  comparison from the element elist(iord(nrmax+1)).
    !
      jbnd = jupbn-1
      ibeg = nrmax+1

      do i = ibeg, jbnd
        isucc = iord(i)
        if ( elist(isucc) <= errmax ) then
          go to 60
        end if
        iord(i-1) = isucc
      end do

      iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
    !
    !  Insert errmin by traversing the list bottom-up.
    !
    60 continue

      iord(i-1) = maxerr
      k = jbnd

      do j = i, jbnd
        isucc = iord(k)
        if ( errmin < elist(isucc) ) then
          go to 80
        end if
        iord(k+1) = isucc
        k = k-1
      end do

      iord(i) = last
      go to 90

    80 continue

      iord(k+1) = last
    !
    !  Set maxerr and ermax.
    !
    90 continue

      maxerr = iord(nrmax)
      ermax = elist(maxerr)

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function qwgtc ( x, c, p2, p3, p4, kp )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QWGTC defines the weight function used by QC25C.
    !
    !  Discussion:
    !
    !    The weight function has the form 1 / ( X - C ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(RK) X, the point at which the weight function is evaluated.
    !
    !    Input, real(RK) C, the location of the singularity.
    !
    !    Input, real(RK) P2, P3, P4, parameters that are not used.
    !
    !    Input, integer(IK) KP, a parameter that is not used.
    !
    !    Output, real(RK) QWGTC, the value of the weight function at X.
    !
      implicit none

      real(RK) c
      integer(IK) kp
      real(RK) p2
      real(RK) p3
      real(RK) p4
      real(RK) qwgtc
      real(RK) x

      qwgtc = 1._RK / ( x - c )

      return
    end
    function qwgto ( x, omega, p2, p3, p4, integr )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QWGTO defines the weight functions used by QC25O.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(RK) X, the point at which the weight function is evaluated.
    !
    !    Input, real(RK) OMEGA, the factor multiplying X.
    !
    !    Input, real(RK) P2, P3, P4, parameters that are not used.
    !
    !    Input, integer(IK) INTEGR, specifies which weight function is used:
    !    1. W(X) = cos ( OMEGA * X )
    !    2, W(X) = sin ( OMEGA * X )
    !
    !    Output, real(RK) QWGTO, the value of the weight function at X.
    !
      implicit none

      integer(IK) integr
      real(RK) omega
      real(RK) p2
      real(RK) p3
      real(RK) p4
      real(RK) qwgto
      real(RK) x

      if ( integr == 1 ) then
        qwgto = cos ( omega * x )
      else if ( integr == 2 ) then
        qwgto = sin ( omega * x )
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function qwgts ( x, a, b, alfa, beta, integr )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! QWGTS defines the weight functions used by QC25S.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(RK) X, the point at which the weight function is evaluated.
    !
    !    Input, real(RK) A, B, the endpoints of the integration interval.
    !
    !    Input, real(RK) ALFA, BETA, exponents that occur in the weight function.
    !
    !    Input, integer(IK) INTEGR, specifies which weight function is used:
    !    1. W(X) = (X-A)**ALFA * (B-X)**BETA
    !    2, W(X) = (X-A)**ALFA * (B-X)**BETA * log (X-A)
    !    3, W(X) = (X-A)**ALFA * (B-X)**BETA * log (B-X)
    !    4, W(X) = (X-A)**ALFA * (B-X)**BETA * log (X-A) * log(B-X)
    !
    !    Output, real(RK) QWGTS, the value of the weight function at X.
    !
      implicit none

      real(RK) a
      real(RK) alfa
      real(RK) b
      real(RK) beta
      integer(IK) integr
      real(RK) qwgts
      real(RK) x

      if ( integr == 1 ) then
        qwgts = ( x - a )**alfa * ( b - x )**beta
      else if ( integr == 2 ) then
        qwgts = ( x - a )**alfa * ( b - x )**beta * log ( x - a )
      else if ( integr == 3 ) then
        qwgts = ( x - a )**alfa * ( b - x )**beta * log ( b - x )
      else if ( integr == 4 ) then
        qwgts = ( x - a )**alfa * ( b - x )**beta * log ( x - a ) * log ( b - x )
      end if

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine timestamp ( )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
      implicit none

      character ( len = 8 ) ampm
      integer(IK) d
      integer(IK) h
      integer(IK) m
      integer(IK) mm
      character ( len = 9 ), parameter, dimension(12) :: month = [ &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December '  ]
      integer(IK) n
      integer(IK) s
      integer(IK) values(8)
      integer(IK) y

      call date_and_time ( values = values )

      y = values(1)
      m = values(2)
      d = values(3)
      h = values(5)
      n = values(6)
      s = values(7)
      mm = values(8)

      if ( h < 12 ) then
        ampm = 'AM'
      else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h < 12 ) then
          ampm = 'PM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
        d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

      return
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module QuadPackSPR_mod