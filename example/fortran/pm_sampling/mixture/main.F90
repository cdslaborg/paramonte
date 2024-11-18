module auxil_mod

    use pm_kind, only: RKG => RKD ! All real kinds are supported.

    implicit none

    real(RKG), parameter :: EPS = sqrt(epsilon(0._RKG))
    real(RKG), parameter :: LARGE = sqrt(huge(0._RKG))

end module auxil_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module data_mod

    use auxil_mod, only: RKG
    use pm_kind, only: IK
    implicit none

    private
    public :: obs_type, posuniv_type, stat_type, data_type

    type, abstract :: obs_type
    end type

    type, extends(obs_type) :: univ_type
        real(RKG) :: val
    end type

    type, extends(univ_type) :: posuniv_type
        real(RKG) :: log
    end type

    interface posuniv_type
        module procedure :: posuniv_typer
    end interface

    type stat_type
        real(RKG) :: lim(2)
    end type

    type :: data_type
        integer(IK) :: nobs
        type(stat_type) :: stat
        class(obs_type), allocatable :: obs(:)
    end type

    interface data_type
        module procedure :: data_typer
    end interface

contains

    function posuniv_typer(val) result(self)
        real(RKG), intent(in), contiguous :: val(:)
        type(posuniv_type), allocatable :: self(:)
        integer(IK) :: iobs
        allocate(self(size(val)))
        do concurrent(iobs = 1 : size(val))
            self(iobs)%log = log(val(iobs))
            self(iobs)%val = val(iobs)
        end do
    end function

    function data_typer(obs, stat) result(self)
        class(obs_type), intent(in), contiguous :: obs(:)
        type(stat_type), intent(in) :: stat
        type(data_type) :: self
        self%stat = stat
        self%nobs = size(obs, 1, IK)
        allocate(self%obs, source = obs)
    end function

end module data_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module domain_mod

    !use pm_container, only: css_type
    use auxil_mod, only: RKG
    use pm_kind, only: SK
    implicit none

    type :: domain_type
        real(RKG), allocatable :: lower(:)
        real(RKG), allocatable :: upper(:)
        !type(css_type), allocatable :: names(:)
        character(63, SK), allocatable :: names(:)
    end type

end module domain_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module model_mod

    use pm_kind, only: SK, IK
    use data_mod, only: obs_type
    use domain_mod, only: domain_type
    use auxil_mod, only: RKG
    implicit none

    type, abstract :: model_type
        type(domain_type) :: domain
        !real(RKG), allocatable :: param(:)
        integer(IK) :: npar
    contains
        procedure(getParam_proc), deferred :: getParam
        procedure(setParam_proc), deferred :: setParam
        procedure(getLogPDF_proc), deferred :: getLogPDF
    end type

    type :: con_type
        class(model_type), allocatable :: model
    end type

    interface con_type
        module procedure :: con_typer
    end interface

    abstract interface
        function getParam_proc(self) result(param)
            import :: RKG, model_type
            class(model_type), intent(in) :: self
            real(RKG) :: param(self%npar)
        end function
        subroutine setParam_proc(self, param)
            import :: RKG, model_type
            class(model_type), intent(inout) :: self
            real(RKG), intent(in), contiguous :: param(:)
        end subroutine
        impure elemental function getLogPDF_proc(self, obs) result(logPDF)
            import :: RKG, model_type, obs_type
            class(model_type), intent(in) :: self
            class(obs_type), intent(in) :: obs
            real(RKG) :: logPDF
        end function
    end interface

contains

    pure elemental function con_typer(model) result(self)
        class(model_type), intent(in) :: model
        type(con_type) :: self
        self%model = model
    end function

end module model_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module lognorm_mod

    use model_mod, only: domain_type, model_type, IK
    use auxil_mod, only: RKG, LARGE
    implicit none

    type, extends(model_type) :: lognorm_type
        real(RKG), private :: avg, invstd, loginvstd
    contains
        procedure :: getParam, setParam, getLogPDF
    end type

    interface lognorm_type
        module procedure :: lognorm_typer
    end interface

contains

    function lognorm_typer(param, domain) result(self)
        real(RKG), intent(in), contiguous, optional :: param(:)
        type(domain_type), intent(in), optional :: domain
        type(lognorm_type) :: self
        self%npar = 2
        if (present(domain)) then
            self%domain = domain
        else
            self%domain = domain_type([-LARGE, -LARGE], [+LARGE, +LARGE], [character(63) :: "lognorm_avg", "lognorm_logstd"])
        end if
        if (present(param)) call self%setParam(param)
    end function

    subroutine setParam(self, param)
        use pm_kind, only: SK, IK
        class(lognorm_type), intent(inout) :: self
        real(RKG), intent(in), contiguous :: param(:)
        if (size(param) /= self%npar) error stop "`lognorm_type%setParam()` takes a vector of two parameters representing `avg` and `logstd`."
        self%invstd = exp(-param(2))
        self%loginvstd = -param(2)
        self%avg = param(1)
    end subroutine

    function getParam(self) result(param)
        class(lognorm_type), intent(in) :: self
        real(RKG) :: param(self%npar)
        param = [self%avg, -self%loginvstd]
    end function

    impure elemental function getLogPDF(self, obs) result(logPDF)
        use pm_distNorm, only: setNormLogPDF
        use data_mod, only: obs_type, posuniv_type
        class(lognorm_type), intent(in) :: self
        class(obs_type), intent(in):: obs
        real(RKG) :: logPDF
        select type (obs)
        type is (posuniv_type)
            call setNormLogPDF(logPDF, obs%log, mu = self%avg, invSigma = self%invstd, logInvSigma = self%loginvstd)
        class default
            error stop "Unrecognized data."
        end select
    end function

end module lognorm_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module flatPoweto_mod

    use auxil_mod, only: RKG, LARGE
    use model_mod, only: domain_type, model_type
    implicit none

    type, extends(model_type) :: flatPoweto_type
        real(RKG), private :: logbreak, alpha, break, alphap1, logPDFNF, limx(2), loglimx(2)
    contains
        procedure, pass :: getParam, setParam, getLogPDF
    end type

    interface flatPoweto_type
        module procedure :: flatPoweto_typer
    end interface

contains

    function flatPoweto_typer(limx, param, domain) result(self)
        real(RKG), intent(in), contiguous, optional :: param(:)
        type(domain_type), intent(in), optional :: domain
        real(RKG), intent(in), contiguous :: limx(:)
        character(:), allocatable :: msg
        type(flatPoweto_type) :: self
        self%npar = 2
        self%limx = limx
        self%loglimx = log(limx)
        if (present(domain)) then
            self%domain = domain
        else
            self%domain = domain_type([-LARGE, -LARGE], [+LARGE, +LARGE], [character(63) :: "flatPoweto_logbreak", "flatPoweto_alpha"])
        end if
        if (present(param)) call self%setParam(param)
    end function

    function getParam(self) result(param)
        class(flatPoweto_type), intent(in) :: self
        real(RKG) :: param(self%npar)
        param = [self%logbreak, self%alpha]
    end function

    subroutine setParam(self, param)
        use pm_kind, only: SK, IK
        use pm_val2str, only: getStr
        use pm_mathMinMax, only: getMinMax
        use pm_mathLogAddExp, only: getLogAddExp
        use pm_mathLogSubExp, only: getLogSubExp
        real(RKG), intent(in), contiguous :: param(:)
        class(flatPoweto_type), intent(inout) :: self
        real(RKG) :: logModelInt1, logModelInt2, small, large
        if (size(param) /= self%npar) error stop "`lognorm_type%setParam()` takes a vector of two parameters representing `avg` and `logstd`."
        self%alpha = param(2)
        self%logbreak = param(1)
        self%alphap1 = self%alpha + 1
        self%break = exp(self%logbreak)
        logModelInt1 = log(self%break - self%limx(1))
        if (self%alphap1 < 0._RKG) then
            logModelInt2 = -log(-self%alphap1)
            large = self%alphap1 * self%logbreak
            small = self%alphap1 * self%loglimx(2)
        elseif (0._RKG < self%alphap1) then
            logModelInt2 = -log(self%alphap1)
            small = self%alphap1 * self%logbreak
            large = self%alphap1 * self%loglimx(2)
        else
            logModelInt2 = 0._RKG
            small = self%logbreak
            large = self%loglimx(2)
        end if
        logModelInt2 = logModelInt2 - self%alpha * self%logbreak + getLogSubExp(smaller = small, larger = large)
        self%logPDFNF = -getLogAddExp(getMinMax(logModelInt1, logModelInt2))
    end subroutine

    pure elemental function getLogPDF(self, obs) result(logPDF)
        use data_mod, only: obs_type, posuniv_type
        class(flatPoweto_type), intent(in) :: self
        class(obs_type), intent(in) :: obs
        real(RKG) :: logPDF
        select type (obs)
        type is (posuniv_type)
            if (obs%log < self%logbreak) then
                logPDF = obs%log
            else
                logPDF = self%alphap1 * obs%log - self%alpha * self%logbreak
            end if
            logPDF = logPDF + self%logPDFNF
        class default
            error stop "Unrecognized data."
        end select
    end function

end module flatPoweto_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module flatPowetoTapered_mod

    use auxil_mod, only: RKG, LARGE
    use model_mod, only: domain_type, model_type
    implicit none

    type, extends(model_type) :: flatPowetoTapered_type
        real(RKG), private  :: logbreak, alpha, beta, break, alphap1, logPDFNF, limx(2), loglimx(2)
    contains
        procedure, private :: getLogUDF
        procedure :: getParam, setParam, getLogPDF
    end type

    interface flatPowetoTapered_type
        module procedure :: flatPowetoTapered_typer
    end interface

contains

    function flatPowetoTapered_typer(limx, param, domain) result(self)
        real(RKG), intent(in), contiguous, optional :: param(:)
        type(domain_type), intent(in), optional :: domain
        real(RKG), intent(in), contiguous :: limx(:)
        type(flatPowetoTapered_type) :: self
        self%npar = 3
        self%limx = limx
        self%loglimx = log(limx)
        if (present(domain)) then
            self%domain = domain
        else
            self%domain = domain_type([-LARGE, -LARGE, -LARGE], [+LARGE, +LARGE, +LARGE], [character(63) :: "flatPowetoTapered_logbreak", "flatPowetoTapered_alpha", "flatPowetoTapered_beta"])
        end if
        if (present(param)) call self%setParam(param)
    end function

    function getParam(self) result(param)
        class(flatPowetoTapered_type), intent(in) :: self
        real(RKG) :: param(self%npar)
        param = [self%logbreak, self%alpha, self%beta]
    end function

    subroutine setParam(self, param)
        use pm_val2str, only: getStr
        use pm_kind, only: SK, IK, LK
        use pm_mathMinMax, only: getMinMax
        use pm_mathLogAddExp, only: getLogAddExp
        use pm_quadPack, only: isFailedQuad, getQuadErr, weps
        class(flatPowetoTapered_type), intent(inout) :: self
        real(RKG), intent(in), contiguous :: param(:)
        real(RKG) :: logModelInt1, logModelInt2, logNF
        character(255, SK) :: msg
        logical(LK) :: failed
        if (size(param) /= self%npar) error stop "`flatPowetoFlatPowetoTapered_type` takes a vector of three parameters"&
        //"representing `logbreak`, `alpha`, `beta`. size(param), npar = "//getStr([size(param, 1, IK), self%npar])
        self%logbreak = param(1)
        self%alpha = param(2)
        self%beta = param(3)
        self%alphap1 = self%alpha + 1
        self%break = exp(self%logbreak)
        logModelInt1 = log(self%break - self%limx(1))
        logNF = max(self%logbreak, self%alphap1 * self%loglimx(2) - self%alpha * self%logbreak - self%beta * (self%limx(2) - self%break))
        failed = isFailedQuad(getDensity, lb = self%logbreak, ub = self%loglimx(2), help = weps, integral = logModelInt2, msg = msg)
        if (failed) then
            !error stop "@getLogModelInt(): "//trim(msg)
            self%logPDFNF = sqrt(huge(self%logPDFNF))
        else
            logModelInt2 = log(logModelInt2) + logNF
            self%logPDFNF = -getLogAddExp(getMinMax(logModelInt1, logModelInt2))
        end if
    contains
        function getDensity(logx) result(density)
            real(RKG), intent(in) :: logx
            real(RKG) :: density
            density = exp(self%getLogUDF(logx, exp(logx)) - logNF)
        end function
    end subroutine

    pure elemental function getLogPDF(self, obs) result(logPDF)
        use data_mod, only: obs_type, posuniv_type
        class(flatPowetoTapered_type), intent(in) :: self
        class(obs_type), intent(in):: obs
        real(RKG) :: logPDF
        select type (obs)
        type is (posuniv_type)
            if (obs%log < self%logbreak) then
                logPDF = obs%log
            else
                logPDF = self%getLogUDF(obs%log, obs%val)
            end if
            logPDF = logPDF + self%logPDFNF
        class default
            error stop "Unrecognized data."
        end select
    end function

    pure elemental function getLogUDF(self, logx, x) result(logUDF)
        class(flatPowetoTapered_type), intent(in) :: self
        real(RKG), intent(in) :: logx, x
        real(RKG) :: logUDF
        logUDF = self%alphap1 * logx - self%alpha * self%logbreak - self%beta * (x - self%break)
    end function

end module flatPowetoTapered_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module mixture_mod

    use pm_val2str, only: getStr
    use pm_arrayResize, only: setResized
    use pm_mathLogAddExp, only: getLogAddExp
    use pm_mathLogSumExp, only: getLogSumExp
    use model_mod, only: domain_type, model_type, con_type
    use auxil_mod, only: RKG, EPS, LARGE
    use data_mod, only: obs_type
    use pm_kind, only: SK, IK
    implicit none

    type, extends(model_type) :: mixture_type
        type(con_type), allocatable :: comp(:)
        real(RKG), allocatable :: logfrac(:)
        integer(IK), private :: ncomp
    contains
        procedure :: getParam, setParam, getLogPDF
    end type

    interface mixture_type
        module procedure :: mixture_typer
    end interface

contains

    function mixture_typer(comp, param) result(self)
        type(con_type), intent(in), contiguous :: comp(:)
        real(RKG), intent(in), contiguous, optional :: param(:)
        character(:, SK), allocatable :: prefix
        integer(IK) :: icomp, lpar, upar, ipar
        type(mixture_type) :: self
        self%ncomp = size(comp, 1, IK)
        allocate(self%comp, source = comp)
        call setResized(self%logfrac, self%ncomp)
        self%npar = sum([(self%comp(icomp)%model%npar, icomp = 1, self%ncomp)]) + self%ncomp - 1
        call setResized(self%domain%lower, self%npar)
        call setResized(self%domain%upper, self%npar)
        call setResized(self%domain%names, self%npar)
        lpar = 0
        do icomp = 1, self%ncomp - 1
            upar = lpar + self%comp(icomp)%model%npar
            prefix = SK_"comp"//getStr(icomp)//SK_"_"
            self%domain%lower(lpar + 1 : upar + 1) = [self%comp(icomp)%model%domain%lower, -LARGE]
            self%domain%upper(lpar + 1 : upar + 1) = [self%comp(icomp)%model%domain%upper, +LARGE]
            self%domain%names(lpar + 1 : upar + 1) = [character(63, SK) :: prefix//self%comp(icomp)%model%domain%names, prefix//SK_"fisherz(frac)"]
            lpar = upar + 1
        end do
        prefix = SK_"comp"//getStr(icomp)//SK_"_"
        self%domain%lower(lpar + 1 :) = self%comp(icomp)%model%domain%lower
        self%domain%upper(lpar + 1 :) = self%comp(icomp)%model%domain%upper
        self%domain%names(lpar + 1 :) = prefix//self%comp(icomp)%model%domain%names
        if (present(param)) call self%setParam(param)
    end function

    subroutine setParam(self, param)
        use pm_mathFisher, only: getFisherInv
        class(mixture_type), intent(inout) :: self
        real(RKG), intent(in), contiguous :: param(:)
        integer(IK) :: icomp, lpar, upar
        real(RKG) :: sumfrac, mixfrac
        if (size(param) /= self%npar) error stop "@mixture_type%setParam(): The condition `size(param) == self%npar` must hold. size(param), self%npar = "//getStr([size(param), self%npar])
        lpar = 0
        sumfrac = 0
        do icomp = 1, self%ncomp - 1
            upar = lpar + self%comp(icomp)%model%npar
            call self%comp(icomp)%model%setParam(param(lpar + 1 : upar))
            mixfrac = getFisherInv(param(upar + 1), 0._RKG, 1._RKG)
            self%logfrac(icomp) = log(mixfrac)
            sumfrac = sumfrac + mixfrac
            lpar = upar + 1
        end do
        self%logfrac(icomp) = log(1 - sumfrac)
        call self%comp(icomp)%model%setParam(param(lpar + 1 :))
        if (.not. abs(1 - sum(exp(self%logfrac))) < EPS) error stop "The condition `sum(exp(self%logfrac)) == 1` must hold. self%logfrac = "//getStr(self%logfrac)
    end subroutine

    function getParam(self) result(param)
        use pm_mathFisher, only: getFisher
        class(mixture_type), intent(in) :: self
        integer(IK) :: icomp, lpar, upar
        real(RKG) :: param(self%npar)
        param = [(self%comp(icomp)%model%getParam(), getFisher(exp(self%logfrac(icomp)), 0._RKG, 1._RKG), icomp = 1, self%ncomp - 1), self%comp(self%ncomp)%model%getParam()]
    end function

    impure elemental function getLogPDF(self, obs) result(logPDF)
        class(mixture_type), intent(in) :: self
        real(RKG) :: logPDFS(self%ncomp), maxLogPDF
        class(obs_type), intent(in) :: obs
        real(RKG) :: logPDF
        integer(IK) :: icomp
        maxLogPDF = -huge(0._RKG)
        do icomp = 1, self%ncomp
            logPDFS(icomp) = self%logfrac(icomp) + self%comp(icomp)%model%getLogPDF(obs)
            maxLogPDF = max(maxLogPDF, logPDFS(icomp))
        end do
        logPDF = getLogSumExp(logPDFS, maxLogPDF)
    end function

end module mixture_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module fit_mod

    use auxil_mod, only: RKG
    use pm_sampling, only: sampler_type
    use pm_sampling, only: paradram_type
    use model_mod, only: model_type
    use data_mod, only: data_type
    implicit none

    private
    public :: fit_type, paradram_type

    type :: best_type
        real(RKG) :: loglike
        real(RKG), allocatable :: param(:) ! best fit model parameters.
    end type

    type :: fit_type
        real(RKG) :: bic
        type(best_type) :: best
        type(data_type) :: data
        class(model_type), allocatable :: model
        class(sampler_type), allocatable :: sampler
    contains
        procedure, pass :: run
        !procedure, pass :: setHist
    end type

    interface fit_type
        module procedure :: fit_typer
    end interface

contains

    function fit_typer(data, model, sampler) result(self)
        class(sampler_type), intent(in), optional :: sampler
        class(model_type), intent(in) :: model
        type(data_type), intent(in) :: data
        type(fit_type) :: self
        if (present(sampler)) then
            self%sampler = sampler
        else
            self%sampler = paradram_type()
        end if
        self%model = model
        self%data = data
    end function

    subroutine run(self)

        use pm_kind, only: SK, IK
        use pm_err, only: err_type
        use pm_sampling, only: getErrSampling, paradram_type
        class(fit_type), intent(inout) :: self
        type(err_type) :: err

        if (.not. allocated(self%sampler%outputStatus)) self%sampler%outputStatus = "retry"
        if (.not. allocated(self%sampler%domainAxisName)) self%sampler%domainAxisName = self%model%domain%names
        if (.not. allocated(self%sampler%domainCubeLimitLower)) self%sampler%domainCubeLimitLower = self%model%domain%lower
        if (.not. allocated(self%sampler%domainCubeLimitUpper)) self%sampler%domainCubeLimitUpper = self%model%domain%upper

        select type (sampler => self%sampler)
        type is (paradram_type)
            if (.not. allocated(sampler%outputChainSize)) sampler%outputChainSize = 30000
            if (.not. allocated(sampler%proposalScale)) sampler%proposalScale = "gelman"
            if (.not. allocated(sampler%proposalStart)) sampler%proposalStart = self%model%getParam()
            err = getErrSampling(sampler, getLogLike, self%model%npar)
        end select
        if (err%occurred) error stop "sampler failed: "//err%msg

        ! Find the mean best fit parameters and write post-prediction data to output file for visualization.

        block

            use pm_io, only: getErrTableRead
            use pm_sysPath, only: glob, css_type
            use pm_parallelism, only: getImageID

            integer(IK) :: stat, ibest!, offset = 1
            type(css_type), allocatable :: path(:)
            real(RKG), allocatable :: logx(:), logPDF(:)
            real(RKG), allocatable :: table(:,:), state(:)

            if (getImageID() == 1) write(*, "(A)") "Searching for files: "//self%sampler%outputFileName//SK_"*"
            path = glob(self%sampler%outputFileName//SK_"*")
            if (size(path) == 0) error stop "There is no sample file in the output folder."

            stat = getErrTableRead(path(size(path))%val, table, roff = 1_IK)
            if (stat /= 0) error stop "Failed to read output sample."
            ibest = maxloc(table(:, 1), dim = 1_IK)
            self%best%loglike = table(ibest, 1)
            self%best%param = table(ibest, 2:)
            self%bic = self%model%npar * log(real(self%data%nobs, RKG)) - 2 * self%best%loglike

        end block

    contains

        function getLogLike(param) result(logLike)
            real(RKG), intent(in), contiguous :: param(:)
            real(RKG) :: logLike
            call self%model%setParam(param)
            loglike = sum(self%model%getLogPDF(self%data%obs))
        end function

    end subroutine

end module fit_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program example

    use pm_kind, only: SK, IK
    use fit_mod, only: fit_type
    use model_mod, only: con_type
    use auxil_mod, only: RKG, LARGE
    use lognorm_mod, only: lognorm_type
    use mixture_mod, only: mixture_type
    use pm_sampling, only: paradram_type
    use pm_arraySpace, only: getLogSpace
    use flatPoweto_mod, only: flatPoweto_type
    use data_mod, only: posuniv_type, data_type, stat_type
    use flatPowetoTapered_mod, only: flatPowetoTapered_type
    use pm_io, only: display_type, css_type, getErrTableWrite

    integer(IK) :: imodel, stat
    type(paradram_type) :: sampler
    type(data_type) :: data, pred
    type(display_type) :: disp
    type(fit_type) :: fit

    disp = display_type(file = SK_"main.out.F90")

    ! read data.

    block
        use pm_io, only: getErrTableRead
        real(RKG), allocatable :: table(:,:)
        character(255, SK) :: iomsg
        stat = getErrTableRead("data.csv", table, roff = 1_IK, iomsg = iomsg)
        if (stat /= 0) error stop "Failed to read data: "//trim(iomsg)
        associate(val => table(:, 4))
            data = data_type(obs = posuniv_type(val), stat = stat_type([minval(val), maxval(val)]))
        end associate
    end block

    ! perform fitting.

    do imodel = 1, 4

        sampler = paradram_type()
        sampler%parallelismMpiFinalizeEnabled = .false.

        if (imodel == 1) then
            sampler%outputFileName = "./mixLogNormLogNorm"
            sampler%proposalStart = [-.2, 0.4, .4, 3.6, -0.2]
            sampler%domainCubeLimitLower = [real(RKG) :: -LARGE, -LARGE, -LARGE, log(3.), -LARGE]
            sampler%domainCubeLimitUpper = [real(RKG) :: log(3.), +LARGE, +LARGE, +LARGE, +LARGE]
            !sampler%proposalScale = "0.1 * gelman"
            sampler%proposalCov = reshape   ( &
                                            [  1.6E-2, 5.5E-3, 4.1E-3, 3.9E-3,-2.5E-3 &
                                            ,  5.5E-3, 3.1E-3, 1.7E-3, 1.4E-3,-9.5E-4 &
                                            ,  4.1E-3, 1.7E-3, 1.9E-3, 1.2E-3,-8.8E-4 &
                                            ,  3.9E-3, 1.4E-3, 1.2E-3, 2.0E-3,-8.9E-4 &
                                            , -2.5E-3,-9.5E-4,-8.8E-4,-8.9E-4, 1.0E-3 &
                                            ], shape = [size(sampler%proposalStart), size(sampler%proposalStart)])
            fit = fit_type(data, mixture_type([con_type(lognorm_type()), con_type(lognorm_type())]), sampler)
        elseif (imodel == 2) then
            sampler%proposalScale = "0.1 * gelman"
            sampler%outputFileName = "./mixFlatPowetoFlatPowetoTapered"
            sampler%proposalStart = [log(.37), -1.4, .3, log(21.3), -.5, 0.0156]
            sampler%proposalCov = reshape   ( &
                                            [  1.E-2, -3.E-3,  1.E-2, -7.E-3, -1.E-4, -6.E-4 &
                                            , -3.E-3,  3.E-3, -4.E-3,  6.E-3,  1.E-4,  1.E-3 &
                                            ,  1.E-2, -4.E-3,  1.E-1, -7.E-2, -8.E-4,  3.E-4 &
                                            , -7.E-3,  6.E-3, -7.E-2,  7.E-2,  8.E-4,  3.E-3 &
                                            , -1.E-4,  1.E-4, -8.E-4,  8.E-4,  1.E-5,  6.E-5 &
                                            , -6.E-4,  1.E-3,  3.E-4,  3.E-3,  6.E-5,  2.E-3 &
                                            ], shape = [size(sampler%proposalStart), size(sampler%proposalStart)])
            fit = fit_type(data, mixture_type([con_type(flatPoweto_type(data%stat%lim)), con_type(flatPowetoTapered_type(data%stat%lim))]), sampler)
        elseif (imodel == 3) then
            sampler%proposalScale = "0.1 * gelman"
            sampler%outputFileName = "./mixLogNormFlatPowetoTapered"
            sampler%proposalStart = [-.2, 0.4, .3, log(21.3), -.5, 0.0156]
            fit = fit_type(data, mixture_type([con_type(lognorm_type()), con_type(flatPowetoTapered_type(data%stat%lim))]), sampler)
        elseif (imodel == 4) then
            sampler%proposalScale = "0.1 * gelman"
            sampler%outputFileName = "./mixFlatPowetoLogNorm"
            sampler%proposalStart = [log(.37), -1.4, .3, 3.6, -0.2]
            fit = fit_type(data, mixture_type([con_type(flatPoweto_type(data%stat%lim)), con_type(lognorm_type())]), sampler)
        else
            sampler%outputFileName = "./logNorm"
            fit = fit_type(data, lognorm_type([real(RKG) :: 1, 1]), sampler)
        end if

        call fit%run()

        call disp%show(css_type(fit%sampler%domainAxisName))
        call disp%show(fit%best%param)
        call disp%skip()
        call disp%show("[fit%best%loglike, fit%bic]")
        call disp%show( [fit%best%loglike, fit%bic] )
        call disp%skip()

        ! generate synthetic data for histogram reconstruction.

        associate(val => getLogSpace(logx1 = log(data%stat%lim(1)), logx2 = log(data%stat%lim(2)), count = 1000_IK))
            pred = data_type(obs = posuniv_type(val), stat = stat_type([minval(val), maxval(val)]))
            call fit%model%setParam(fit%best%param)
            stat = getErrTableWrite(fit%sampler%outputFileName//SK_".hist", reshape([log(val), fit%model%getLogPDF(pred%obs)], [size(val), 2]), header = SK_"logx,logPDF")
            if (stat /= 0) error stop "Failed to write the histogram visualization data."
        end associate

    end do

end