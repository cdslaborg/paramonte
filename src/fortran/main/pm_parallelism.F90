!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This module contains procedures and generic interfaces for facilitating parallel computations or computing the performance of the parallel Coarray/MPI/OpenMP algorithms.
!>
!>  \details
!>  A primary goal for the design of this type and the associated procedures is to facilitate parallelism-agnostic coding within a codebase.<br>
!>  As such, the attributes and procedures within this module are guaranteed to also work in serial applications.<br>
!>  However, the procedures of this module typically do nothing when the application is serial or the parallelism is not recognized.<br>
!>
!>  \test
!>  [test_pm_parallelism](@ref test_pm_parallelism)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_parallelism

    use pm_kind, only: SK, IK, LK
    use pm_val2str, only: getStr

    implicit none

    character(*, SK), parameter :: MODULE_NAME = SK_"@pm_parallelism"

#if CAF_ENABLED || MPI_ENABLED
    !>  \brief
    !>  The scalar constant of type `character` of default kind \SK, whose value is set to `parallel`
    !>  if the ParaMonte library is built with the preprocessor macro `-DCAF_ENABLED=1` `-DMPI_ENABLED=1`,
    !>  otherwise, it is set to `"serial"`.<br>
    character(*, SK), parameter :: PARALLELIZATION_MODE = SK_"parallel"
#else
    character(*, SK), parameter :: PARALLELIZATION_MODE = SK_"serial"
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the [imageis_type](@ref pm_parallelism::imageis_type) type for generating objects with components
    !>  of type `logical` of default kind \LK that contain information about the current image/processor/thread.<br>
    !>
    !>  \details
    !>  Objects of this type are not meant to be used directly by the end user.<br>
    !>  This type merely exists to create the `is` component of the [image_type](@ref pm_parallelism::image_type) class.<br>
    !>
    !>  \interface{imageis_type}
    !>  \code{.F90}
    !>
    !>      use pm_parallelism, only: imageis_type
    !>      type(imageis_type) :: imageis
    !>
    !>      imageis = imageis_type()
    !>
    !>  \endcode
    !>
    !>  \final{imageis_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type                                :: imageis_type
        logical(LK)                     :: first    = .false._LK    !<  \public The scalar `logical` of default kind \LK indicating whether the current process is ID #1.
        logical(LK)                     :: extra    = .false._LK    !<  \public The scalar `logical` of default kind \LK indicating whether the current process is NOT ID #1.
        logical(LK)                     :: leader   = .false._LK    !<  \public The scalar `logical` of default kind \LK indicating whether the current process is a leader.<br>
                                                                    !!          The default value is `.true.` <b>if and only if</b> the corresponding image ID is `1`.<br>
                                                                    !!          **Otherwise, it must be set by the user after calling the type constructor depending on the parallelism type**.<br>
        logical(LK)                     :: rooter   = .false._LK    !<  \public The scalar `logical` of default kind \LK indicating whether the current process is a follower.<br>
                                                                    !!          The default value is `.true.` <b>if and only if</b> the corresponding image ID is **not** `1`.<br>
                                                                    !!          **Otherwise, it must be set by the user after calling the type constructor depending on the parallelism type**.<br>
    end type

    !>  \brief
    !>  This is the [image_type](@ref pm_parallelism::image_type) type for generating objects that contain
    !>  information about the current image/processor/thread and facilitate its synchronization with other processes,
    !>  or the global finalization of all inter-process parallel communications (e.g., as is done in MPI applications).
    !>
    !>  \details
    !>  <ol>
    !>      <li>    Within the Fortran terminology, an **image** (or **process**) is equivalent to a replication of the program having its own set of data objects.<br>
    !>      <li>    The number of images could be the same as or more than or less than the available number of physical processors.
    !>      <li>    A particular implementation may permit the number of images to be chosen at compile time, at link time, or at execution time.
    !>      <li>    Each image executes asynchronously, and the normal rules of Fortran apply within each image.
    !>      <li>    The execution sequence can differ from image to image as specified by the programmer who, with the help
    !>              of a unique image index, determines the actual path using normal Fortran control constructs and explicit synchronizations.
    !>  </ol>
    !>
    !>  \return
    !>  `image` :   The output scalar of type [image_type](@ref pm_parallelism::image_type).
    !>
    !>  \interface{image_type}
    !>  \code{.F90}
    !>
    !>      use pm_parallelism, only: image_type
    !>      type(image_type) :: image
    !>
    !>      image = image_type()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Even when the application is MPI/OpenMP-parallel, the `id` assigned to the zeroth process is `1` in this type.<br>
    !>  This convention follows the rules of the Fortran Coarray parallel programming language.<br>
    !>
    !>  \warning
    !>  For **OpenMP-enabled** ParaMonte library builds, this routine return the value of the OpenMP-intrinsic `omp_get_num_threads()` routine,
    !>  that is, the number of threads in the team that is executing the parallel region to which the routine region binds.<br>
    !>  If this constructor is called from the sequential part of an OpenMP-enabled program, this routine returns `1`.<br>
    !>  To return a value larger than `1`, this routine must be called from within an OpenMP-enabled parallel region.<br>
    !>  The desired number of OpenMP threads can be set at runtime via [setImageCount()](@ref pm_parallelism::setImageCount),
    !>  or directly via the OpenMP intrinsic subroutine `omp_set_num_threads()`.<br>
    !>
    !>  \final{image_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type                                :: image_type
        integer(IK)                     :: count        = -huge(1_IK)       !<  \public The scalar `integer` of default kind \IK representing the total count of runtime parallel processes available within the current communication.
        integer(IK)                     :: id           = -huge(1_IK)       !<  \public The scalar `integer` of default kind \IK representing the ID of the runtime parallel process starting with `1`: `1`, `2`, `3`, ...
        type(imageis_type)              :: is           = imageis_type()    !<  \public The scalar of type [imageis_type](@ref pm_parallelism::imageis_type) containing `logical` components that signify the current image role.
        character(:, SK), allocatable   :: label                            !<  \public The `allocatable` scalar `character` of default kind \SK containing the ID of the current process in the format `@process(ID)`.
        contains
        procedure, nopass               :: sync => setImageSynced
        procedure, nopass               :: finalize => setImageFinalized
    end type

    !>  \cond excluded
    interface image_type
        module procedure :: image_typer
    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the predicted parallel Fork-Join speedup scaling behavior for simulations whose image
    !>  contributions are stochastically accepted only from the first successful image starting from image ID `1`.
    !>
    !>  \details
    !>  This generic interface can be used to compute theoretical speedup gained by a
    !>  Fork-Join parallel simulation where the probability of contribution of an image to the simulation
    !>  is less than `1`, meaning that image contributions are stochastic at each step of the Fork-Join cycle.<br>
    !>  If the image contributions are inspected sequentially from the first image to last and only the first successful image contribution
    !>  is kept in the simulation, then it can be shown that the contribution of the images to the parallel Fork-Join simulation follows a [Cyclic Geometric Distribution](@ref pm_distGeomCyclic).<br>
    !>  See [Amir Shahmoradi, Fatemeh Bagheri (2020). ParaDRAM: A Cross-Language Toolbox for Parallel High-Performance Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo Simulations.](https://www.cdslab.org/pubs/2020_Shahmoradi_I.pdf)
    !>  for theoretical details.<br>
    !>
    !>  \param[in]      conProb         :   The input scalar of,
    !>                                      <ol>
    !>                                          <li>    type `real` of kind \RKALL,
    !>                                      </ol>
    !>                                      containing the average probability of contribution of each image to the Fork-Join simulation.<br>
    !>                                      For example, `conProb` can be the effective acceptance rate in parallel Fork-Join Monte Carlo or MCMC sampling simulations
    !>                                      where each image is tasked with returning a single proposal state for subsequent evaluation and if is first to be accepted accepted,
    !>                                      for usage in the next simulation step.<br>
    !>  \param[in]      seqSecTime      :   The input scalar of the same type and kind as `conProb` containing the time (in seconds) **per image (or process or thread)** of the inherently-sequential sections of the entire Fork-Join simulation.<br>
    !>  \param[in]      parSecTime      :   The input scalar of the same type and kind as `conProb` containing the time (in seconds) **per image (or process or thread)** of the inherently-parallel sections of the entire Fork-Join simulation.<br>
    !>  \param[in]      comSecTime      :   The input scalar of the same type and kind as `conProb` containing the time (in seconds) **per image (or process or thread)** of the **pairwise** inter-process communication sections of the entire Fork-Join simulation.<br>
    !>                                      Assuming all rooter images communicate only with the main image and the entire communication takes `overhead` seconds, then `comSecTime = overhead / (nproc - 1)` with `nproc` representing the total number of images.<br>
    !>  \param[out]     scaling         :   The output `allocatable` vector of the same type and kind as the input `conProb` containing the predicted speedup scaling of
    !>                                      the stochastic Fork-Join simulation under a range of possible image (or thread or process) counts given in the output `numproc`.<br>
    !>                                      On output, the vector `scaling` will be resized until the maximum speedup and the corresponding image count is identified in the scaling.<br>
    !>                                      The `i`th element of `scaling` represents theoretically-predicted speedup if `numproc(i)` images (threads or processes) were used in the parallel Fork-Join simulation.<br>
    !>  \param[out]     numproc         :   The output `allocatable` vector of type `integer` of default kind \IK containing the
    !>                                      number of images for the speedup reported in the corresponding element of `scaling`.<br>
    !>                                      If the input argument `comSecTime` is positive, then `numproc` is simply a linear range of integers starting with `1` with a jump size of `1`.<br>
    !>                                      If the input argument `comSecTime` is positive, then `numproc` is simply a log(2)-linear range of integers starting with `1` with a jump size of `1`.<br>
    !>  \param[in]      scalingMaxVal   :   The output scalar of the same type and kind as the input `conProb` containing the maximum achievable speedup (corresponding to `maxval(scaling)`).<br>
    !>  \param[in]      scalingMaxLoc   :   The output scalar of type `integer` of default kind \IK containing the index of element of the output `scaling` containing `scalingMaxVal` (i.e., `maxloc(scaling)`).<br>
    !>  \param[out]     scalingMinLen   :   The input positive scalar of type `integer` of default kind \IK containing the minimum size of the output `scaling`.<br>
    !>                                      (**optional**, default = `int(2 / max(conProb, 0.001))`)
    !>
    !>  \interface{setForkJoinScaling}
    !>  \code{.F90}
    !>
    !>      use pm_parallelism, only: setForkJoinScaling
    !>
    !>      call setForkJoinScaling(conProb, seqSecTime, parSecTime, comSecTime, scaling, numproc, scalingMaxVal, scalingMaxLoc, scalingMinLen = scalingMinLen)
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= seqSecTime` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= parSecTime` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= comSecTime` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= scalingMinLen` must hold for the corresponding input arguments.<br>
    !>  The condition `0 <= conProb .and. conProb <= 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \see
    !>  [pm_sampling](@ref pm_sampling)<br>
    !>  [pm_distGeomCyclic](@ref pm_distGeomCyclic)<br>
    !>  [Amir Shahmoradi, Fatemeh Bagheri (2020). ParaDRAM: A Cross-Language Toolbox for Parallel High-Performance Delayed-Rejection Adaptive Metropolis Markov Chain Monte Carlo Simulations.](https://www.cdslab.org/pubs/2020_Shahmoradi_I.pdf)<br>
    !>
    !>  \example{setForkJoinScaling}
    !>  \include{lineno} example/pm_parallelism/setForkJoinScaling/main.F90
    !>  \compilef{setForkJoinScaling}
    !>  \output{setForkJoinScaling}
    !>  \include{lineno} example/pm_parallelism/setForkJoinScaling/main.out.F90
    !>  \postproc{setForkJoinScaling}
    !>  \include{lineno} example/pm_parallelism/setForkJoinScaling/main.py
    !>  \vis{setForkJoinScaling}
    !>  \image html pm_parallelism/setForkJoinScaling/setForkJoinScaling.png width=700
    !>
    !>  \test
    !>  [test_pm_parallelism](@ref test_pm_parallelism)
    !>
    !>  \final{setForkJoinScaling}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 4:13 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface setForkJoinScaling

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setForkJoinScaling_RK5(conProb, seqSecTime, parSecTime, comSecTime, scaling, numproc, scalingMaxVal, scalingMaxLoc, scalingMinLen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setForkJoinScaling_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: scalingMaxLoc
        integer(IK)         , intent(in)    , optional      :: scalingMinLen
        integer(IK)         , intent(out)   , allocatable   :: numproc(:)
        real(RKG)           , intent(out)   , allocatable   :: scaling(:)
        real(RKG)           , intent(in)                    :: conProb, seqSecTime, parSecTime, comSecTime
        real(RKG)           , intent(out)                   :: scalingMaxVal
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setForkJoinScaling_RK4(conProb, seqSecTime, parSecTime, comSecTime, scaling, numproc, scalingMaxVal, scalingMaxLoc, scalingMinLen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setForkJoinScaling_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: scalingMaxLoc
        integer(IK)         , intent(in)    , optional      :: scalingMinLen
        integer(IK)         , intent(out)   , allocatable   :: numproc(:)
        real(RKG)           , intent(out)   , allocatable   :: scaling(:)
        real(RKG)           , intent(in)                    :: conProb, seqSecTime, parSecTime, comSecTime
        real(RKG)           , intent(out)                   :: scalingMaxVal
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setForkJoinScaling_RK3(conProb, seqSecTime, parSecTime, comSecTime, scaling, numproc, scalingMaxVal, scalingMaxLoc, scalingMinLen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setForkJoinScaling_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: scalingMaxLoc
        integer(IK)         , intent(in)    , optional      :: scalingMinLen
        integer(IK)         , intent(out)   , allocatable   :: numproc(:)
        real(RKG)           , intent(out)   , allocatable   :: scaling(:)
        real(RKG)           , intent(in)                    :: conProb, seqSecTime, parSecTime, comSecTime
        real(RKG)           , intent(out)                   :: scalingMaxVal
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setForkJoinScaling_RK2(conProb, seqSecTime, parSecTime, comSecTime, scaling, numproc, scalingMaxVal, scalingMaxLoc, scalingMinLen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setForkJoinScaling_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: scalingMaxLoc
        integer(IK)         , intent(in)    , optional      :: scalingMinLen
        integer(IK)         , intent(out)   , allocatable   :: numproc(:)
        real(RKG)           , intent(out)   , allocatable   :: scaling(:)
        real(RKG)           , intent(in)                    :: conProb, seqSecTime, parSecTime, comSecTime
        real(RKG)           , intent(out)                   :: scalingMaxVal
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setForkJoinScaling_RK1(conProb, seqSecTime, parSecTime, comSecTime, scaling, numproc, scalingMaxVal, scalingMaxLoc, scalingMinLen)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setForkJoinScaling_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: scalingMaxLoc
        integer(IK)         , intent(in)    , optional      :: scalingMinLen
        integer(IK)         , intent(out)   , allocatable   :: numproc(:)
        real(RKG)           , intent(out)   , allocatable   :: scaling(:)
        real(RKG)           , intent(in)                    :: conProb, seqSecTime, parSecTime, comSecTime
        real(RKG)           , intent(out)                   :: scalingMaxVal
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \cond excluded
#if OMP_ENABLED
    logical(LK) , save :: mv_failed = .false._LK
#endif
!>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an object of class [image_type](@ref pm_parallelism::image_type)
    !>  containing information and statistics of the parallel images/processes/threads available, depending on the type of parallelism requested.<br>
    !>
    !>  \details
    !>  This procedure is the constructor of the type [image_type](@ref pm_parallelism::image_type).<br>
    !>  See the documentation of [image_type](@ref pm_parallelism::image_type) for further details and example usage.<br>
    !>
    !>  \return
    !>  `image`             :   The output scalar of type [image_type](@ref pm_parallelism::image_type).
    !>
    !>  \interface{image_typer}
    !>  \code{.F90}
    !>
    !>      use pm_parallelism, only: image_type
    !>      type(image_type) :: image
    !>
    !>      image = image_type()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Even when the application is MPI/OpenMP-parallel, the `id` assigned to the zeroth process is `1` in this type.<br>
    !>  This convention follows the rules of the Fortran Coarray parallel programming language.<br>
    !>
    !>  \warning
    !>  For **OpenMP-enabled** ParaMonte library builds, this routine return the value of the OpenMP-intrinsic `omp_get_num_threads()` routine,
    !>  that is, the number of threads in the team that is executing the parallel region to which the routine region binds.<br>
    !>  If this constructor is called from the sequential part of an OpenMP-enabled program, this routine returns `1`.<br>
    !>  To return a value larger than `1`, this routine must be called from within an OpenMP-enabled parallel region.<br>
    !>  The desired number of OpenMP threads can be set at runtime via [setImageCount()](@ref pm_parallelism::setImageCount),
    !>  or directly via the OpenMP intrinsic subroutine `omp_set_num_threads()`.<br>
    !>
    !>  \final{image_typer}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    function image_typer() result(image)
        type(image_type) :: image
        !integer(IK), intent(in), optional :: nthread
        ! setup general processor / coarray image variables
        image%id            = getImageID()
        image%count         = getImageCount()
        image%label         = SK_"@process("//getStr(image%id)//SK_")"
        image%is%first       = image%id == 1_IK
        image%is%extra    = image%id /= 1_IK
       !image%is%leader      = .false._LK ! ATTN: this is to be set by the user at runtime, depending on the parallelism type.
       !image%is%rooter      = .false._LK ! ATTN: this is to be set by the user at runtime, depending on the parallelism type.
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Synchronize all existing parallel images and return nothing.
    !>
    !>  \details
    !>  This is a static member of the [image_type](@ref pm_parallelism::image_type) class.<br>
    !>  This procedure synchronizes all images in the current communication world by,
    !>  <ol>
    !>      <li>    calling `sync all` in Coarray applications.
    !>      <li>    calling `mpi_barrier()` in MPI applications.
    !>      <li>    calling nothing in all other (including serial) applications.
    !>  </ol>
    !>
    !>  \warning
    !>  This routine contains global Coarray and MPI synchronization barriers and therefore, must be called by all processes in the current simulation.
    !>
    !>  \warning
    !>  The MPI library must be initialized and not finalized prior to calling this routine for MPI-parallel applications.
    !>
    !>  \interface{setImageSynced}
    !>  \code{.F90}
    !>
    !>      use pm_parallelism, only: setImageSynced()
    !>
    !>      call setImageSynced()
    !>
    !>  \endcode
    !>
    !>  \test
    !>  [test_pm_parallelism](@ref test_pm_parallelism)
    !>
    !>  \final{setImageSynced}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 4:13 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    subroutine setImageSynced()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setImageSynced
#endif
#if     CAF_ENABLED
        sync all
#elif   MPI_ENABLED
        block
            use mpi !, only: mpi_barrier, mpi_comm_world
            integer :: ierrMPI
            call mpi_barrier(mpi_comm_world, ierrMPI)
        end block
#elif   OMP_ENABLED
        !$omp barrier
#endif
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Finalize the current parallel simulation and return nothing.<br>
    !>
    !>  \details
    !>  This is a static member of the [image_type](@ref pm_parallelism::image_type) class.<br>
    !>  <ol>
    !>      <li>    This procedure facilitates a parallelism-agnostic way of finalizing serial and parallel applications by hiding the finalization step within this procedure.<br>
    !>      <li>    This procedure is primarily relevant to MPI applications where MPI needs to be finalized before the simulation stop.<br>
    !>      <li>    This procedure performs a `sync all` for Coarray applications and returns.<br>
    !>      <li>    This procedure performs nothing for other all applications (e.g., serial, OpenMP, ...).<br>
    !>  </ol>
    !>
    !>  \interface{setImageFinalized}
    !>  \code{.F90}
    !>
    !>      use pm_parallelism, only: setImageFinalized
    !>
    !>      call setImageFinalized()
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The MPI library must be initialized and not finalized prior to calling this routine for MPI-parallel applications.
    !>
    !>  \warning
    !>  The MPI communications will be shut down upon calling this routine and further inter-process communications will be impossible.
    !>
    !>  \test
    !>  [test_pm_parallelism](@ref test_pm_parallelism)
    !>
    !>  \final{setImageFinalized}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 4:13 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    subroutine setImageFinalized() ! LCOV_EXCL_LINE
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setImageFinalized
#endif
#if     CAF_ENABLED
        sync all
#elif   MPI_ENABLED
        use mpi !mpi_f08, only: mpi_comm_world, mpi_finalized, mpi_finalize, mpi_barrier
        implicit none
        integer :: ierrMPI
        logical :: isFinalized
        call mpi_finalized(isFinalized, ierrMPI)
        if (.not. isFinalized) then
            call mpi_barrier(mpi_comm_world, ierrMPI)
            call mpi_finalize(ierrMPI)
        end if
#elif   OMP_ENABLED
        !$omp barrier
#endif
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the ID of the current Coarray image / MPI process / OpenMP thread, all starting with `1`.<br>
    !>
    !>  \return
    !>  `imageID`   :   The output scalar `integer` of default kind \IK representing the ID of the current image/process/thread.<br>
    !>                  The output `imageID` is set to,
    !>                  <ol>
    !>                      <li>    The value returned by `this_image()` in Coarray applications,
    !>                      <li>    The value returned by `mpi_comm_rank() + 1` in MPI applications,
    !>                      <li>    The value returned by `omp_get_thread_num() + 1` in OpenMP applications,
    !>                      <li>    The value `1` in serial applications,
    !>                  </ol>
    !>
    !>  \interface{getImageID}
    !>  \code{.F90}
    !>
    !>      use pm_parallelism, only: getImageID
    !>      use pm_kind, only: IK
    !>      integer(IK) :: id
    !>
    !>      id = getImageID() ! The ID (starting at 1) of the current Coarray image, MPI process, or OpenMP thread.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The MPI library must be initialized and not finalized prior to calling this routine for MPI-parallel applications.<br>
    !>  Otherwise, the MPI library will automatically be initialized, which may fail if it has been already finalized.<br>
    !>
    !>  \test
    !>  [test_pm_parallelism](@ref test_pm_parallelism)
    !>
    !>  \final{getImageID}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 4:13 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    function getImageID() result(imageID)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getImageID
#endif
        integer(IK) :: imageID
#if     CAF_ENABLED
        imageID = this_image()
#elif   MPI_ENABLED
        block
            use mpi !mpi_f08, only : mpi_initialized, mpi_comm_world, mpi_comm_size, mpi_init
            integer :: rank, ierrMPI
            logical :: isinit, isfinit
            call mpi_initialized(isinit, ierrMPI)
            if (.not. isinit) then
                call mpi_finalized(isfinit, ierrMPI)
                if (isfinit) error stop MODULE_NAME//"@getImageID(): Error occurred. A finalized MPI library cannot be reinitialized."
                call mpi_init(ierrMPI)
            end if
            call mpi_comm_rank(mpi_comm_world, rank, ierrMPI)
            if (ierrMPI /= 0) error stop "Failed to fetch the MPI process counts."
            imageID = int(rank, IK) + 1_IK
        end block
#elif   OMP_ENABLED
        block
            use omp_lib, only: omp_get_thread_num
            imageID = int(omp_get_thread_num() + 1, IK)
        end block
#else
        imageID = 1_IK
#endif
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the number of available processes in the current parallel world communication.
    !>
    !>  \return
    !>  `imageCount`    :   The output scalar `integer` of default kind \IK representing the number of available processes in the current world communication.<br>
    !>                      The output `imageCount` is set to,<br>
    !>                      <ol>
    !>                          <li>    the value returned by the Fortran intrinsic `num_images()` in Coarray parallel applications,
    !>                          <li>    the value returned by `mpi_comm_size()` subroutine in MPI parallel applications,
    !>                          <li>    the value returned by `omp_get_num_threads()` in OpenMP parallel applications,
    !>                          <li>    the value `1` in serial applications,
    !>                      </ol>
    !>
    !>  \warning
    !>  For **OpenMP-enabled** ParaMonte library builds, this routine return the value of the OpenMP-intrinsic `omp_get_num_threads()` routine,
    !>  that is, the number of threads in the team that is executing the parallel region to which the routine region binds.<br>
    !>  If this generic interface is called from the sequential part of an OpenMP-enabled program, this routine returns `1`.<br>
    !>  To return a value larger than `1`, this routine must be called from within an OpenMP-enabled parallel region.<br>
    !>  The desired number of OpenMP threads can be set at runtime via [setImageCount()](@ref pm_parallelism::setImageCount),
    !>  or directly via the OpenMP intrinsic subroutine `omp_set_num_threads()`.<br>
    !>
    !>  \interface{getImageCount}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_parallelism, only: getImageCount
    !>      integer(IK) :: imageCount
    !>
    !>      imageCount = getImageCount()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  See the warnings associated with [getImageCountMPI](@ref pm_parallelism::getImageCountMPI).<br>
    !>
    !>  \test
    !>  [test_pm_parallelism](@ref test_pm_parallelism)
    !>
    !>  \final{getImageCount}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 4:13 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    function getImageCount() result(imageCount)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getImageCount
#endif
        integer(IK) :: imageCount
#if     CAF_ENABLED
        imageCount = num_images()
#elif   MPI_ENABLED
        imageCount = getImageCountMPI()
#elif   OMP_ENABLED
        imageCount = getImageCountOMP()
#else
        imageCount = 1_IK
#endif
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the number of available processes in the current **MPI-parallel** world communication.
    !>
    !>  \return
    !>  `imageCount`    :   The output scalar `integer` of default kind \IK representing the number of
    !>                      available processes in the current **MPI-parallel** world communication.<br>
    !>                      The output `imageCount` is set to,<br>
    !>                      <ol>
    !>                          <li>    the value returned by `mpi_comm_size()` subroutine in MPI parallel applications,
    !>                          <li>    the value `0` if the ParaMonte library build is not for MPI applications,
    !>                      </ol>
    !>
    !>  \interface{getImageCountMPI}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_parallelism, only: getImageCountMPI
    !>      integer(IK) :: imageCount
    !>
    !>      imageCount = getImageCountMPI()
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The MPI library must be initialized and not finalized prior to calling this routine for MPI-parallel applications.<br>
    !>  Otherwise, the MPI library will automatically be initialized, which may fail if it has been already finalized.<br>
    !>
    !>  \note
    !>  Th C-binding of this routine is required for automatic detection of the use of `mpiexec` launcher in dynamic programming languages.
    !>
    !>  \test
    !>  [test_pm_parallelism](@ref test_pm_parallelism)
    !>
    !>  \final{getImageCountMPI}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 4:13 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    function getImageCountMPI() result(imageCount) bind(C, name = "getImageCountMPI")
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getImageCountMPI
#endif
#if     MPI_ENABLED
        use mpi !mpi_f08, only : mpi_initialized, mpi_comm_world, mpi_comm_size, mpi_init
        integer :: nproc
        integer :: ierrMPI
        logical :: isinit
        integer(IK) :: imageCount
        imageCount = 0_IK
        call mpi_initialized(isinit, ierrMPI)
        if (ierrMPI /= 0) return ! LCOV_EXCL_LINE
        if (.not. isinit) then
            call mpi_init(ierrMPI) ! LCOV_EXCL_LINE
            if (ierrMPI /= 0) return ! LCOV_EXCL_LINE
        end if
        call mpi_comm_size(mpi_comm_world, nproc, ierrMPI)
        if (ierrMPI /= 0) return ! LCOV_EXCL_LINE
        imageCount = int(nproc, IK)
#else
        integer(IK) :: imageCount
        imageCount = 0_IK
#endif
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the number of available processes in the current **OpenMP-parallel** world communication.
    !>
    !>  \return
    !>  `imageCount`    :   The output scalar `integer` of default kind \IK representing the number of
    !>                      available processes in the current **OpenMP-parallel** world communication.<br>
    !>                      The output `imageCount` is set to,<br>
    !>                      <ol>
    !>                          <li>    the value returned by `mpi_comm_size()` subroutine in OpenMP parallel applications,
    !>                          <li>    the value `0` if the ParaMonte library build is not for OpenMP applications,
    !>                      </ol>
    !>
    !>  \interface{getImageCountOMP}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_parallelism, only: getImageCountOMP
    !>      integer(IK) :: imageCount
    !>
    !>      imageCount = getImageCountOMP()
    !>
    !>  \endcode
    !>
    !>  \test
    !>  [test_pm_parallelism](@ref test_pm_parallelism)
    !>
    !>  \final{getImageCountOMP}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 4:13 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    function getImageCountOMP() result(imageCount) bind(C, name = "getImageCountOMP")
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getImageCountOMP
#endif
#if     OMP_ENABLED
        use omp_lib, only: omp_get_num_threads
        integer(IK) :: imageCount
        imageCount = omp_get_num_threads()
#else
        integer(IK) :: imageCount
        imageCount = 0_IK
#endif
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Set the number the parallel **threads** for an OpenMP-enabled application.
    !>
    !>  \details
    !>  This generic interface is exclusively relevant to the OpenMP-enabled ParaMonte library builds and applications.<br>
    !>  This interface offers a convenient consistent cross-platform method of setting the number of runtime OpenMP threads for a parallel code section.<br>
    !>  It does so by calling the OpenMP intrinsic
    !>
    !>  \param[in]  count : The input scalar of type `integer` of default kind \IK containing the requested number of OpenMP threads.<br>
    !>                      This argument is relevant only if the OpenMP parallelism is enabled at the time of building the ParaMonte library.<br>
    !>                      Unlike Coarray or MPI parallelism paradigms, the number of images (threads) can be (re)set at runtime in OpenMP applications, hence the need for this argument.<br>
    !>                      If the input `count` is **non-positive**, this generic interface uses the output of OpenMP intrinsic `omp_get_num_procs()` as the number of threads.<br>
    !>                      This is unlike the default behavior of Intel and GNU OpenMP libraries which set the number of threads to `max(count, 1)`.<br>
    !>                      (**optional**, default = `omp_get_num_procs()`)
    !>
    !>  \devnote
    !>  Although this procedure exists for all parallelism builds of the ParaMonte library, it does nothing for non OpenMP-enabled library builds.<br>
    !>  Note that the number of images/processes for Coarray/MPI -enabled applications must be set by the user at the time of running the parallel application.<br>
    !>  This is unlike the OpenMP parallel applications where the number of threads can be (re)set at runtime during the program execution.<br>
    !>
    !>  \interface{setImageCount}
    !>  \code{.F90}
    !>
    !>      use pm_parallelism, only: setImageCount
    !>
    !>      call setImageCount(count = count)
    !>
    !>  \endcode
    !>
    !>  \final{setImageCount}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    subroutine setImageCount(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setImageCount
#endif
        integer(IK), intent(in), optional :: count
#if     OMP_ENABLED
        block
            use omp_lib, only: omp_get_num_procs, omp_set_num_threads
            integer :: count_def
            count_def = omp_get_num_procs()
            if (present(count)) then
                if (0_IK < count) count_def = int(count)
            end if
            call omp_set_num_threads(count_def)
        end block
#endif
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !>  \brief
        !>  Broadcast the error condition from all images/processes to all images/processes.<br>
        !>
        !>  \brief
        !>  This procedure exists primarily to aid the detection and graceful
        !>  handling of runtime errors during the testing of ParaMonte library samplers
        !>  or less frequently, when calling the samplers from dynamic programming languages or interactive
        !>  environments such as Jupyter which require graceful handling of any errors to avoid environment shutdown.<br>
        !>  The use of this procedure is essential to avoid runtime parallel deadlocks.<br>
        !>
        !>  \param[in]  failed  :   The input scalar `logical` of default kind \LK that should be `.true.` if and only if a fatal error has occurred.<br>
        !>                          The value of `failed` is typically the `occurred` component of an object of type [err_type](@ref pm_err::err_type).<br>
        !>                          On input, the specified value for `failed` is broadcast to all other active processes in the program.<br>
        !>
        !>  \return
        !>  `failedParallelism` :   The output scalar of type `logical` of default kind \LK that is `.true.` if and only if
        !>                          the input argument `failed` on **any** image/process in the current parallel simulation is `.true.`.<br>
        !>
        !>  \interface{isFailedImage}
        !>  \code{.F90}
        !>
        !>      use pm_kind, only: LK
        !>      use pm_parallelism, only: isFailedImage
        !>      logical(LK) :: failedParallelism
        !>      logical(LK) :: failed
        !>
        !>      failedParallelism = isFailedImage(failed)
        !>
        !>  \endcode
        !>
        !>  \warning
        !>  This subroutine must be called in parallel by **all** images or **none**.<br>
        !>
        !>  \warning
        !>  This function primarily exists for soft handling of fatal errors in parallel
        !>  testing mode and should therefore ideally not be used in production builds.<br>
        !>  If it is used repeatedly within an application, it can **significantly degrade** the parallel performance.<br>
        !>
        !>  \impure
        !>
        !>  \test
        !>  [test_pm_parallelism](@ref test_pm_parallelism)
        !>
        !>  \final{isFailedImage}
        !>
        !>  \author
        !>  \AmirShahmoradi, 9:49 PM Friday, March 1, 2013, Institute for Fusion Studies, The University of Texas Austin<br>
        function isFailedImage(failed) result(failedParallelism)
#if     __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
            !DEC$ ATTRIBUTES DLLEXPORT :: isFailedImage
#endif
            logical(LK), intent(in) :: failed
            logical(LK) :: failedParallelism
#if         CAF_ENABLED
            logical(LK), allocatable, save :: failure(:)[:]
            integer :: iid, imageCount
            imageCount = num_images()
            allocate(failure(imageCount)[*])
            failure(this_image()) = failed
            sync all
            do iid = 1, imageCount
                failedParallelism = failure(iid)[iid]
                if (failedParallelism) return
            end do
#elif       MPI_ENABLED
            logical, allocatable :: failure(:) ! This must be default kind.
            integer :: imageCount!, imageID ! This must be default kind.
            block
                use mpi !mpi_f08, only: mpi_comm_world, mpi_comm_rank, mpi_comm_size, mpi_allgather, mpi_logical
                use pm_arrayResize, only: setResized
                integer :: ierrMPI
               !call mpi_comm_rank(mpi_comm_world, imageID, ierrMPI)
                call mpi_comm_size(mpi_comm_world, imageCount, ierrMPI)
                call setResized(failure, int(imageCount, IK))
                call mpi_allgather  ( logical(failed) & ! LCOV_EXCL_LINE : send buffer
                                    , 1 & ! LCOV_EXCL_LINE : send count
                                    , mpi_logical & ! LCOV_EXCL_LINE : send datatype
                                    , failure(:) & ! LCOV_EXCL_LINE : receive buffer
                                    , 1 & ! LCOV_EXCL_LINE : receive count
                                    , mpi_logical & ! LCOV_EXCL_LINE : receive datatype
                                    , mpi_comm_world & ! LCOV_EXCL_LINE : comm
                                    , ierrMPI & ! LCOV_EXCL_LINE : error code
                                    )
                !call mpi_alltoall   ( err%occurred &    ! buffer_send   : The buffer containing the data that will be scattered to other processes.<br>
                !                    , 1 &               ! count_send    : The number of elements that will be sent to each process.<br>
                !                    , mpi_logical &     ! datatype_send : The type of one send buffer element.<br>
                !                    , failure &   ! buffer_recv   : The buffer in which store the gathered data.<br>
                !                    , imageCount &      ! count_recv    : The number of elements in the message to receive per process, not the total number of elements to receive from all processes altogether.<br>
                !                    , mpi_logical &     ! datatype_recv : The type of one receive buffer element.<br>
                !                    , mpi_comm_world &  ! communicator  : The communicator in which the all to all takes place.<br>
                !                    , ierrMPI & ! LCOV_EXCL_LINE : error code
                !                    )
                failedParallelism = logical(any(failure), LK)
            end block
#elif       OMP_ENABLED
            !$omp critical
            mv_failed = mv_failed .or. failed
            !$omp end critical
            !$omp barrier
            failedParallelism = mv_failed
            !$omp master
            mv_failed = .false._LK
            !$omp end master
#else
            failedParallelism = failed
#endif
        end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_parallelism ! LCOV_EXCL_LINE