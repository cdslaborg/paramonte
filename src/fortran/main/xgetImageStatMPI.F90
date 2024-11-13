!>  \brief
!>  This program contains instructions to determine the number of
!>  MPI-parallel processes and MPI image ranks used in MPI-parallelized programs.<br>
!>
!>  \details
!>  A primary goal for the design of this program is to allow automatic invocation of the MPI-parallel versions
!>  of the ParaMonte library from within dynamic programming language environments (e.g., MATLAB, Python, R, ...)
!>  without the use having to explicitly specify the MPI-parallelism or the specific MPI vendor used.<br>
!>  The goal is to allow to write scripts in dynamic languages that are parallelism-agnostic.<br>
!>
!>  Corresponding to each MPI-parallelized ParaMonte library, a binary for this program will also be generated
!>  that can be called from within dynamic languages using system-level applications to check if and which MPI
!>  library is being used on the command line.<br>
!>
!>  Multiple execution scenarios are possible:<br>
!>  <ol>
!>      <li>    If program executes successfully, each MPI process will output the following information:<br>
!>              \code{.sh}
!>                  imageCountRank<imageCount>imageCountRank<imageRank>imageCountRank
!>              \endcode
!>              where `<imageCount>` and `<imageRank>` are replaced with the actual inferred positive integer values at runtime.<br>
!>              This pattern is intentionally used to ensure easy parsing of the binary output.<br>
!>              A value of `1` for `imageCount` implies that either:<br>
!>              <ol>
!>                  <li>    the program has been called serially (without invoking the `mpiexec` launcher), or,<br>
!>                  <li>    the number of images (processes) specified with the MPI launcher `mpiexec` is `1`.<br>
!>              </ol>
!>              In either case, it is more sensible to simply call the serial version of the ParaMonte library,
!>              if available, which is always true for the final releases of the ParaMonte library for dynamic languages.<br>
!>      <li>    If program executes but an MPI-specific error occurs at runtime, the output message pattern will be the same as above.<br>
!>              However, the values are `0` for both the `<imageCount>` and `<imageRank>`, as in the following:<br>
!>              \code{.sh}
!>                  imageCountRank0imageCountRank0imageCountRank
!>              \endcode
!>      <li>    If the execution fails (e.g., when the corresponding MPI library is missing on the system),
!>              the process will most likely exist with a platform/process-dependent error message.<br>
!>              Parsing the output of the error message and searching for the keyword `imageCountRank`
!>              can determine whether a runtime error has occurred or not.<br>
!>  </ol>
!>
!>  \note
!>  The prefix `xget` in the program name stands for *execute* and *get*.<br>
!>
!>  \final{xgetImageStatMPI}
!>
!>  \author
!>  \AmirShahmoradi, 1:01 PM Sunday, November 10, 2024, Dallas, TX<br>
program xgetImageStatMPI
    use mpi !mpi_f08, only : mpi_initialized, mpi_comm_world, mpi_comm_size, mpi_init
    implicit none
    logical :: isinit
    logical :: isfinit
    integer :: ierrMPI
    integer :: imageRank = -1
    integer :: imageCount = 0
    imageCountMPI_block: block
        call mpi_finalized(isfinit, ierrMPI)
        if (isfinit) error stop "@xgetImageStatMPI(): Error occurred. A finalized MPI library cannot be reinitialized."
        call mpi_initialized(isinit, ierrMPI)
        if (ierrMPI /= 0) exit imageCountMPI_block
        if (.not. isinit) then
            call mpi_init(ierrMPI)
            if (ierrMPI /= 0) exit imageCountMPI_block
        end if
        call mpi_comm_size(mpi_comm_world, imageCount, ierrMPI)
        if (ierrMPI /= 0) imageCount = 0
        call mpi_comm_rank(mpi_comm_world, imageRank, ierrMPI)
        if (ierrMPI /= 0) imageRank = -1
        imageRank = imageRank + 1
    end block imageCountMPI_block
    ! All images report their rank and the MPI image count.
    ! The current bash process is responsible for capturing it.
    write(*, "(*(g0))") "imageCountRank", imageCount, "imageCountRank", imageRank, "imageCountRank"
    !call mpi_finalize(ierrMPI)
end program xgetImageStatMPI