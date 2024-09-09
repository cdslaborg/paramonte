program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKS, RKD, RKH ! all processor real/complex kinds are supported.
    use pm_io, only: display_type
    use pm_polynomial, only: setPolyRoot
    use pm_polynomial, only: eigen, jenkins, laguerre, sgl
    use pm_polynomial, only: eigen_type, jenkins_type, laguerre_type, sgl_type
    use pm_arrayResize, only: setResized
    use pm_polynomial, only: getPolyVal

    implicit none

    integer(IK) :: count

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

! Template to avoid code duplication in the
! example for various types and kind parameters.
! See the output below for the actual code.
! developer guideline:
! TYPE - type of coefficient vector: `real`, `complex`
! RKG  - kind of coefficient vector: any supported by the processor
! CREP - Number of coefficient elements to display per line: use `1` for `real` and `2` for `complex` coefficients.
#define GET_ROOT(CREP, TYPE, RKG, METHOD) \
block; \
use pm_val2str, only: getStr; \
TYPE(RKG), allocatable :: coef(:); \
complex(RKG), allocatable :: root(:); \
type(METHOD) :: method; \
call disp%skip(); \
coef = COEF; \
call disp%show("coef"); \
call disp%show( coef , format = "(sp,"//getStr(CREP)//"(g0,:,', '))"); \
call disp%show("call setResized(root, size(coef, 1, IK) - 1_IK)"); \
                call setResized(root, size(coef, 1, IK) - 1_IK); \
call disp%show("[same_type_as(method, sgl), same_type_as(method, eigen), same_type_as(method, jenkins), same_type_as(method, laguerre)]"); \
call disp%show( [same_type_as(method, sgl), same_type_as(method, eigen), same_type_as(method, jenkins), same_type_as(method, laguerre)] ); \
call disp%show("call setPolyRoot(root, count, coef, eigen)"); \
                call setPolyRoot(root, count, coef, eigen); \
call disp%show("count"); \
call disp%show( count ); \
call disp%show("root(1:count)"); \
call disp%show( root(1:count) , format = "(sp,2(g0,:,', '))"); \
call disp%show("getPolyVal(coef, root(1:count))"); \
call disp%show( getPolyVal(coef, root(1:count)) , format = "(sp,2(g0,:,', '))"); \
call disp%skip(); \
end block;

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the roots of polynomials with real coefficients.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define COEF \
    [8, -8, 16, -16, 8, -8]
    GET_ROOT(1, real, RKS, sgl_type)
    GET_ROOT(1, real, RKD, sgl_type)
    GET_ROOT(1, real, RKH, sgl_type)
    GET_ROOT(1, real, RKS, eigen_type)
    GET_ROOT(1, real, RKD, eigen_type)
    GET_ROOT(1, real, RKH, eigen_type)
    GET_ROOT(1, real, RKS, jenkins_type)
    GET_ROOT(1, real, RKD, jenkins_type)
    GET_ROOT(1, real, RKH, jenkins_type)
    GET_ROOT(1, real, RKS, laguerre_type)
    GET_ROOT(1, real, RKD, laguerre_type)
    GET_ROOT(1, real, RKH, laguerre_type)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the roots of polynomials with real coefficients.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define COEF \
    [3628800, -10628640, 12753576, -8409500, 3416930, -902055, 157773, -18150, 1320, -55, 1]
    GET_ROOT(1, real, RKS, sgl_type)
    GET_ROOT(1, real, RKD, sgl_type)
    GET_ROOT(1, real, RKH, sgl_type)
    GET_ROOT(1, real, RKS, eigen_type)
    GET_ROOT(1, real, RKD, eigen_type)
    GET_ROOT(1, real, RKH, eigen_type)
    GET_ROOT(1, real, RKS, jenkins_type)
    GET_ROOT(1, real, RKD, jenkins_type)
    GET_ROOT(1, real, RKH, jenkins_type)
    GET_ROOT(1, real, RKS, laguerre_type)
    GET_ROOT(1, real, RKD, laguerre_type)
    GET_ROOT(1, real, RKH, laguerre_type)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the roots of polynomials with complex coefficients with zeros 1,2,...,10.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define COEF \
    [3628800, -10628640, 12753576, -8409500, 3416930, -902055, 157773, -18150, 1320, -55, 1]
    GET_ROOT(2, complex, RKS, sgl_type)
    GET_ROOT(2, complex, RKD, sgl_type)
    GET_ROOT(2, complex, RKH, sgl_type)
    GET_ROOT(2, complex, RKS, eigen_type)
    GET_ROOT(2, complex, RKD, eigen_type)
    GET_ROOT(2, complex, RKH, eigen_type)
    GET_ROOT(2, complex, RKS, jenkins_type)
    GET_ROOT(2, complex, RKD, jenkins_type)
    GET_ROOT(2, complex, RKH, jenkins_type)
    GET_ROOT(2, complex, RKS, laguerre_type)
    GET_ROOT(2, complex, RKD, laguerre_type)
    GET_ROOT(2, complex, RKH, laguerre_type)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the roots of polynomials with complex coefficients with zeros on imaginary axis.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define COEF \
    [(0._RKH, 1._RKH), (-10001.0001_RKH, 0._RKH), (0._RKH, -10001.0001_RKH), (1._RKH, 0._RKH)]
    GET_ROOT(2, complex, RKD, sgl_type)
    GET_ROOT(2, complex, RKH, sgl_type)
    GET_ROOT(2, complex, RKD, eigen_type)
    GET_ROOT(2, complex, RKH, eigen_type)
    GET_ROOT(2, complex, RKD, jenkins_type)
    GET_ROOT(2, complex, RKH, jenkins_type)
    GET_ROOT(2, complex, RKD, laguerre_type)
    GET_ROOT(2, complex, RKH, laguerre_type)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the roots of polynomials with complex coefficients with zeros at 1+i,1/2*(1+i)....1/(2**-9)*(1+i).")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define COEF \
[ (0._RKH, 9.094947017729282e-13_RKH) \
, (-4.652065399568528e-10_RKH, -4.652065399568528e-10_RKH) \
, (1.584803612786345e-7_RKH, 0._RKH) \
, (-1.154642632172909e-5_RKH, 1.154642632172909e-5_RKH) \
, (0._RKH, -7.820779428584501e-4_RKH) \
, (1.271507365163416e-2_RKH, 1.271507365163416e-2_RKH) \
, (-.2002119533717632e0_RKH, 0._RKH) \
, (.7567065954208374e0_RKH, -7.567065954208374e-1_RKH) \
, (0._RKH, 2.658859252929688e0_RKH) \
, (-1.998046875_RKH, -1.998046875_RKH) \
, (1._RKH, 0._RKH) \
]
GET_ROOT(2, complex, RKD, sgl_type)
GET_ROOT(2, complex, RKH, sgl_type)
GET_ROOT(2, complex, RKD, eigen_type)
GET_ROOT(2, complex, RKH, eigen_type)
GET_ROOT(2, complex, RKD, jenkins_type)
GET_ROOT(2, complex, RKH, jenkins_type)
GET_ROOT(2, complex, RKD, laguerre_type)
GET_ROOT(2, complex, RKH, laguerre_type)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the roots of polynomials with complex coefficients with multiple zeros.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define COEF \
[ (288._RKH, 0._RKH) \
, (-1344._RKH, 504._RKH) \
, (2204._RKH, -2352._RKH) \
, (-920._RKH, 4334._RKH) \
, (-1587._RKH, -3836._RKH) \
, (2374._RKH, 1394._RKH) \
, (-1293._RKH, 200._RKH) \
, (284._RKH, -334._RKH) \
, (3._RKH, 100._RKH) \
, (-10._RKH, -10._RKH) \
, (1._RKH, 0._RKH) \
]
GET_ROOT(2, complex, RKS, sgl_type)
GET_ROOT(2, complex, RKD, sgl_type)
GET_ROOT(2, complex, RKH, sgl_type)
GET_ROOT(2, complex, RKS, eigen_type)
GET_ROOT(2, complex, RKD, eigen_type)
GET_ROOT(2, complex, RKH, eigen_type)
GET_ROOT(2, complex, RKS, jenkins_type)
GET_ROOT(2, complex, RKD, jenkins_type)
GET_ROOT(2, complex, RKH, jenkins_type)
GET_ROOT(2, complex, RKS, laguerre_type)
GET_ROOT(2, complex, RKD, laguerre_type)
GET_ROOT(2, complex, RKH, laguerre_type)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the roots of polynomials with complex coefficients with 12 zeros evenly distribute on a circle of radius 1. centered at 0+2i.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#define COEF \
[ (4095._RKH, 0._RKH) \
, (0._RKH, 24576._RKH) \
, (-67584._RKH, 0._RKH) \
, (0._RKH, -112640._RKH) \
, (126720._RKH, 0._RKH) \
, (0._RKH, 101376._RKH) \
, (-59136._RKH, 0._RKH) \
, (0._RKH, -25344._RKH) \
, (7920._RKH, 0._RKH) \
, (0._RKH, 1760._RKH) \
, (-264._RKH, 0._RKH) \
, (0._RKH, -24._RKH) \
, (1._RKH, 0._RKH) \
]
GET_ROOT(2, complex, RKS, sgl_type)
GET_ROOT(2, complex, RKD, sgl_type)
GET_ROOT(2, complex, RKH, sgl_type)
GET_ROOT(2, complex, RKS, eigen_type)
GET_ROOT(2, complex, RKD, eigen_type)
GET_ROOT(2, complex, RKH, eigen_type)
GET_ROOT(2, complex, RKS, jenkins_type)
GET_ROOT(2, complex, RKD, jenkins_type)
GET_ROOT(2, complex, RKH, jenkins_type)
GET_ROOT(2, complex, RKS, laguerre_type)
GET_ROOT(2, complex, RKD, laguerre_type)
GET_ROOT(2, complex, RKH, laguerre_type)

end program example