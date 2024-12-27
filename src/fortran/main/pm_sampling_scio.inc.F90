    !%%%%%%%%%%%%%%%
#if namelist_ENABLED
    !%%%%%%%%%%%%%%%

    ! specbase
    namelist /NAMELIST/ description
    namelist /NAMELIST/ domain
    namelist /NAMELIST/ domainAxisName
    namelist /NAMELIST/ domainBallAvg
    namelist /NAMELIST/ domainBallCor
    namelist /NAMELIST/ domainBallCov
    namelist /NAMELIST/ domainBallStd
    namelist /NAMELIST/ domainCubeLimitLower
    namelist /NAMELIST/ domainCubeLimitUpper
    namelist /NAMELIST/ domainErrCount
    namelist /NAMELIST/ domainErrCountMax
    namelist /NAMELIST/ inputFileHasPriority
    namelist /NAMELIST/ outputChainFileFormat
    namelist /NAMELIST/ outputColumnWidth
    namelist /NAMELIST/ outputFileName
    namelist /NAMELIST/ outputStatus
    namelist /NAMELIST/ outputPrecision
    namelist /NAMELIST/ outputReportPeriod
    namelist /NAMELIST/ outputRestartFileFormat
    namelist /NAMELIST/ outputSampleSize
    namelist /NAMELIST/ outputSeparator
    namelist /NAMELIST/ outputSplashMode
    namelist /NAMELIST/ parallelism
    namelist /NAMELIST/ parallelismMpiFinalizeEnabled
    namelist /NAMELIST/ parallelismNumThread
    namelist /NAMELIST/ plang
    namelist /NAMELIST/ randomSeed
    namelist /NAMELIST/ sysInfoFilePath
    namelist /NAMELIST/ targetAcceptanceRate
    ! specmcmc
    namelist /NAMELIST/ outputChainSize
    namelist /NAMELIST/ outputSampleRefinementCount
    namelist /NAMELIST/ outputSampleRefinementMethod
    namelist /NAMELIST/ proposal
    namelist /NAMELIST/ proposalCor
    namelist /NAMELIST/ proposalCov
    namelist /NAMELIST/ proposalScale
    namelist /NAMELIST/ proposalStart
    namelist /NAMELIST/ proposalStartDomainCubeLimitLower
    namelist /NAMELIST/ proposalStartDomainCubeLimitUpper
    namelist /NAMELIST/ proposalStartRandomized
    namelist /NAMELIST/ proposalStd
    ! specdram
    namelist /NAMELIST/ proposalAdaptationBurnin
    namelist /NAMELIST/ proposalAdaptationCount
    namelist /NAMELIST/ proposalAdaptationCountGreedy
    namelist /NAMELIST/ proposalAdaptationPeriod
    namelist /NAMELIST/ proposalDelayedRejectionCount
    namelist /NAMELIST/ proposalDelayedRejectionScale
    ! specnest
    namelist /NAMELIST/ domainPartitionAdaptationCount
    namelist /NAMELIST/ domainPartitionAdaptationPeriod
    namelist /NAMELIST/ domainPartitionBiasCorrectionEnabled
    namelist /NAMELIST/ domainPartitionCountMax
    namelist /NAMELIST/ domainPartitionFactorExpansion
    namelist /NAMELIST/ domainPartitionFactorShrinkage
    namelist /NAMELIST/ domainPartitionKmeansClusterCountMax
    namelist /NAMELIST/ domainPartitionKmeansClusterSizeMin
    namelist /NAMELIST/ domainPartitionKmeansNormalizationEnabled
    namelist /NAMELIST/ domainPartitionKmeansNumFailMax
    namelist /NAMELIST/ domainPartitionKmeansNumRecursionMax
    namelist /NAMELIST/ domainPartitionKmeansNumTry
    namelist /NAMELIST/ domainPartitionKvolumeNumRecursionMax
    namelist /NAMELIST/ domainPartitionKvolumeWeightExponent
    namelist /NAMELIST/ domainPartitionMethod
    namelist /NAMELIST/ domainPartitionObject
    namelist /NAMELIST/ domainPartitionOptimizationScaleEnabled
    namelist /NAMELIST/ domainPartitionOptimizationShapeEnabled
    namelist /NAMELIST/ domainPartitionOptimizationShapeScaleEnabled
    namelist /NAMELIST/ domainPartitionScaleFactor
    namelist /NAMELIST/ domainSampler
    namelist /NAMELIST/ liveSampleSize
    namelist /NAMELIST/ tolerance

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   ParaDISE_ENABLED || ParaDRAM_ENABLED || ParaNest_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RETURN_IF_FAILED(LINE,FAILED,MSG) \
if (FAILED) then; \
err%occurred = .true._LK; \
err%msg = PROCEDURE_NAME//getFine(__FILE__, LINE)//SK_": "//trim(MSG); \
if (.not. opened) close(unit, iostat = err%stat); \
return; \
end if;
#if     ParaDISE_ENABLED
#define CFC_TYPE cfcdise_type
#define getErrChainRead getErrChainReadDISE
#define getErrChainWrite getErrChainWriteDISE
#define chainFileColName chainFileColNameDISE
#define isFailedChainResize isFailedChainResizeDISE
#elif   ParaDRAM_ENABLED
#define CFC_TYPE cfcdram_type
#define getErrChainRead getErrChainReadDRAM
#define getErrChainWrite getErrChainWriteDRAM
#define chainFileColName chainFileColNameDRAM
#define isFailedChainResize isFailedChainResizeDRAM
#elif   ParaNest_ENABLED
#define CFC_TYPE cfcnest_type
#define getErrChainRead getErrChainReadNest
#define getErrChainWrite getErrChainWriteNest
#define chainFileColName chainFileColNameNest
#define isFailedChainResize isFailedChainResizeNest
#else
#error  "Unrecognized interface."
#endif
    !>  \brief
    !>  Generate and return an object of type [err_type](@ref pm_err:err_type) whose components inform
    !>  the occurrence of an during reading and returning the contents of the specified input chain `file` of a [ParaMonte sampler](@ref pm_sampling).<br>
    !>
    !>  \details
    !>  This generic interface has two primary calling formats:<br>
    !>  <ol>
    !>      <li>    One in which all chain file contents are stored implicitly in the `chain` argument.
    !>      <li>    One in which all chain file contents are stored explicitly in the respectively arguments following `chain`.
    !>  </ol>
    !>  The latter exists solely because of the lack of PDT support in \gfortran{13.1} and earlier versions.<br>
    !>
    !>  \param[inout]   cfc                     :   The input/output scalar object that can be of,<br>
    !>                                              <ol>
    !>                                                  <li>    type [cfcdise_type](@ref pm_sampling_dise::cfcdise_type) indicating the chain to read is an output of the ParaDRAM sampler,<br>
    !>                                                  <li>    type [cfcdram_type](@ref pm_sampling_dram::cfcdram_type) indicating the chain to read is an output of the ParaDISE sampler,<br>
    !>                                                  <li>    type [cfcnest_type](@ref pm_sampling_nest::cfcnest_type) indicating the chain to read is an output of the ParaNest sampler,<br>
    !>                                              </ol>
    !>                                              If any optional arguments below are present, then `chain` has `intent(in)` and will remain untouched.<br>
    !>                                              If all optional arguments below are missing, then `chain` has `intent(inout)` and it components will be filled with chain file contents upon return.<br>
    !>  \param[in]      file                    :   The input scalar `character` of default kind \SK, containing the path to the chain file whose contents is to be read.<br>
    !>  \param[in]      form                    :   The input scalar `character` of default kind \SK, containing the `form` of the chain file whose contents is to be read.<br>
    !>                                              This argument is typically the value of `outputChainFileFormat` in the [base settings](@ref pm_sampling_base) of a ParaMonte sampler.<br>
    !>                                              But in general, the following values are accepted and recognized:<br>
    !>                                              <ol>
    !>                                                  <li>    The string `"unformatted"`, implying that the chain file form is `unformatted` (i.e., binary, non-humane-readable) whose contents must be read with `stream` access.<br>
    !>                                                  <li>    The string `"formatted"`, implying that the chain file form is `formatted` (i.e., ASCII, humane-readable).<br>
    !>                                                  <li>    The string `"compact"`, implying that the chain file form is `formatted` (i.e., ASCII, humane-readable).<br>
    !>                                                          Additionally, it implies that the individual states within the chain can be potentially unequally weighted.<br>
    !>                                                          This option is solely relevant to MCMC type of chain files (returned by [ParaDISE](pm_sampling::paradise_type) and [ParaDRAM](pm_sampling::paradram_type) samplers.<br>
    !>                                                          For all other samplers, it is functionally equivalent to specifying `"formatted"` for the input argument `form`.<br>
    !>                                                  <li>    The string `"verbose"`, implying that the chain file form is `formatted` (i.e., ASCII, humane-readable).<br>
    !>                                                          Additionally, any two consecutive states in the chain file that are identical will be merged and the corresponding `sampleWeight` will be incremented by one.<br>
    !>                                                          This option is solely relevant to MCMC type of chain files (returned by [ParaDISE](pm_sampling::paradise_type) and [ParaDRAM](pm_sampling::paradram_type) samplers.<br>
    !>                                                          For all other samplers, it is functionally equivalent to specifying `"formatted"` for the input argument `form`.<br>
    !>                                              </ol>
    !>                                              The specified value for `form` must **not** be `"unformatted"` (binary) if `chain` has `intent(inout)`.<br>
    !>                                              Furthermore, if the contents of the chain to read are in binary format, the kinds of the column storage space vectors,
    !>                                              particularly, for columns holding `real` value must be correct, otherwise, the IO process is bound to fail, as the procedure has no way of discerning the kinds of the values from the file.<br>
    !>  \param[in]      pre                     :   The input scalar `integer` of default kind \IK, containing the pre-specified chain size with which all output components used for storing the chain file contents must be allocated.<br>
    !>                                              This is relevant to the case where further samples are to be added to the components, for example, in a simulation restart.<br>
    !>                                              Specifying `pre` will also lead to faster chain read action, because it requires only one memory allocation per chain file column as opposed to repeated resizing.<br>
    !>                                              However, if specified, `pre` must be at least as large as the total number rows in the chain file (minus `1`, the header line).<br>
    !>                                              (**optional**. If missing, data will be written to the corresponding component of `chain`. It can be present **only if** the specified `chain` object contains component with the same name.)
    !>  \param[out]     pos                     :   The output `allocatable` vector of the same size as the output (weighted) chain (that is `nsam`) of type `integer` of default **processor** kind (which may not necessarily be the same as \IK),
    !>                                              containing the positions of the beginning of the row records of the binary chain file as returned by the `pos` specifier of the Fortran intrinsic `read()` statement.<br>
    !>                                              Note that `pos` does not contain the starting position of the file header as it is `1` by definition.<br>
    !>                                              This vector can be later used for repositioning the chain file for subsequent writes to the file.<br>
    !>                                              (**optional**. Its presence is relevant only if the input `form` is `"binary"` or `"unformatted"`.)
    !>
    !>  \return
    !>  `err`                                   :   The output scalar of type [err_type](@ref pm_err::err_type) whose component `occur` is set to the `.true.` **if and only if** the
    !>                                              procedure fails to accomplish the task, in which case, the `msg` component of `err` is set to a description of the error origin.<br>
    !>                                              On output, the `msg` component of `err` is an empty string if the task is accomplished successfully.<br>
    !>                                              Possible causes of read action failure include (but are not limited to):<br>
    !>                                              <ol>
    !>                                                  <li>    The file header contents have been changed.
    !>                                                  <li>    The file is read with a wrong specified `form`.
    !>                                                  <li>    The file contents have been compromised by the user.
    !>                                              </ol>
    !>
    !>  \interface{getErrChainRead}
    !>  \code{.F90}
    !>
    !>      use pm_sampling_scio, only: getErrChainRead
    !>      use pm_err, only: err_type
    !>      type(err_type) :: err
    !>
    !>      err = getErrChainRead(chain, file, form, pre = pre, pos = pos)
    !>      !
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \devnote
    !>  Upon entry, it file is closed already, it will be opened and rewound and will be closed upon return.<br>
    !>  If file is open, it is closed and reopened and rewound, however, it will remain open and repositioned correctly such that further write to the file can resume upon return.<br>
    !>
    !>  \example{getErrChainRead}
    !>  \include{lineno} example/pm_sampling_scio/getErrChainRead/main.F90
    !>  \compilef{getErrChainRead}
    !>  \output{getErrChainRead}
    !>  \include{lineno} example/pm_sampling_scio/getErrChainRead/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampling_scio](@ref test_pm_sampling_scio)
    !>
    !>  \todo
    !>  \pvhigh
    !>  The madness seen in the current interface of the procedures of this generic interface is due to the lack of PDT support in \gfortran{13.1}.<br>
    !>  Ideally, the derived type [chainParaDRAM_type](@ref pm_sampling_scio::chainParaDRAM_type) must be converted to a parameterized derived type (PDT) once supported by \gfortran.<br>
    !>  Consequently, all chain contents components that are explicitly passed to the procedures can be implicitly returned as components of the output `chain` argument.<br>
    !>
    !>  \final{getErrChainRead}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas Austin<br>
    function getErrChainRead(cfc, file, form, pre, pos) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErrChainRead
#endif
        use pm_io, only: setContentsFrom
        use pm_strASCII, only: getStrLower
        !use pm_arrayRebind, only: setRebound
        use pm_arrayResize, only: setResized
        use pm_arrayFind, only: getCountLoc
        use pm_arraySplit, only: setSplit
        use pm_arrayFind, only: discrete
        use pm_io, only: getCountRecord
        use pm_io, only: setRecordFrom
        use pm_sysPath, only: isExtant
        use pm_val2str, only: getStr
        use pm_io, only: isOpen

        type(CFC_TYPE), intent(inout) :: cfc
        character(*, SK), intent(in) :: file, form
        integer, intent(out), allocatable, optional :: pos(:)
        integer(IK), intent(in), optional :: pre
        type(CFC_TYPE) :: row
        type(err_type) :: err

        integer :: unit, posCurrent
        character(:, SK), allocatable :: record
        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//SK_"@getErrChainRead()"
        logical(LK) :: opened, isBinary, isASCII, isVerbose
        integer(IK) :: lenrec, ibeg, iend, ndim, pre_def, nsep, seplen
#if     ParaDISE_ENABLED || ParaDRAM_ENABLED
        logical(LK) :: isPastFirstRead
        isPastFirstRead = .false._LK
#elif   !ParaNest_ENABLED
#error  "Unrecognized interface."
#endif
        record = getStrLower(form)
        isVerbose = record == SK_"verbose"
        isBinary = record == SK_"unformatted" .or. record == SK_"binary"
        isASCII = record == SK_"formatted" .or. record == SK_"compact" .or. record == SK_"ascii" .or. isVerbose
        err = err_type(msg = repeat(SK_" ", 255))

        ! First check if file is open and close it.

        inquire(file = file, exist = err%occurred, opened = opened, number = unit, iostat = err%stat, iomsg = err%msg)
        RETURN_IF_FAILED(__LINE__,err%stat /= 0 .or. .not. err%occurred,err%msg)
        !if (opened) close(unit, iostat = err%stat, iomsg = err%msg)
        !RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
        err%occurred = .false._LK ! THIS MUST BE HERE.

        ! open file with the right access and read the header.

        if (.not. opened) then
            if (isBinary) then
                open(newunit = unit, file = file, form = "unformatted", access = "stream", status = "old", iostat = err%stat, iomsg = err%msg SHARED)
            elseif (isASCII) then
                open(newunit = unit, file = file, status = "old", iostat = err%stat, iomsg = err%msg SHARED)
            end if
            RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
        end if

        if (isBinary) then
            ! header is assumed to terminate with a null character in unformatted files.
            call setResized(record, 132_IK)
            lenrec = 0_IK
            do
                lenrec = lenrec + 1_IK
                if (len(record, IK) < lenrec) call setResized(record)
                read(unit, pos = lenrec, iostat = err%stat, iomsg = err%msg) record(lenrec : lenrec)
                RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
                if (record(lenrec : lenrec) /= achar(0)) cycle
                lenrec = lenrec - 1_IK
                exit
            end do
        elseif (isASCII) then
            call setRecordFrom(unit, record, iostat = err%stat, iomsg = err%msg, ub = lenrec)
            RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
        else
            err%occurred = .true._LK
            err%msg = PROCEDURE_NAME//SK_"Unrecognized file `form`: """//form//SK_""""
            return
        end if

        ! Parse the header to infer the separator and the column names.

        ibeg = len_trim(chainFileColName(1)) + 1 ! beginning of the separator.
        seplen = index(record(ibeg : lenrec), trim(chainFileColName(2)), kind = IK) - 1_IK
        cfc%sep = record(ibeg : ibeg + seplen - 1)
        cfc%header = record(1 : lenrec)
        RETURN_IF_FAILED(__LINE__,seplen < 1_IK,SK_"Failed to infer the chain file column separator: """//cfc%header//SK_"""")
        nsep = getCountLoc(cfc%header, cfc%sep, discrete, blindness = len(cfc%sep, IK))
        ndim = nsep - size(chainFileColName, 1, IK) + 1_IK
        RETURN_IF_FAILED(__LINE__,ndim < 1_IK,SK_"Failed to infer the chain file column separator: """//cfc%header//SK_"""")

        ! \todo
        ! Additional check could be performed here to ensure the first default column names are indeed correctly inferred from the chain file.

        ! Allocate components.

        pre_def = 100000_IK
        if (present(pre)) pre_def = pre
        RETURN_IF_FAILED(__LINE__,isFailedChainResize(row, ndim, 1_IK, err%msg),SK_"Failed to allocate storage for the chain row.")
        RETURN_IF_FAILED(__LINE__,isFailedChainResize(cfc, ndim, pre_def, err%msg),SK_"Failed to allocate storage for the chain columns.")
        if (present(pos)) then
            call setResized(pos, pre_def, failed = err%occurred, errmsg = err%msg)
            RETURN_IF_FAILED(__LINE__,err%occurred,err%msg)
        end if

        cfc%nsam = 0_IK
        loopOverRows: do

            blockFileForm: if (isBinary) then
                ! \todo
                ! The following reading approach cannot capture end of record errors in binary mode.
                ! If the record is incomplete, an end of file error will be returned instead.
                ! As such, restart with compromised binary chain file may fail.
                ! This could be fixed by separating the individual field reads from each other and
                ! setting end of file error for all fields (except the last) to `iostat_eor` for later correction outside the read block.
                ! This was currently deemed too much work with little obvious gain. The future developer can decide to proceed with this or not if needed.
                ! There may be a remote possibility of natural incomplete record to occur under natural premature simulation stop (without manual chain file compromise).
                ! But the odds are extremely low and the effort to check for this not worth currently.
                ! An alternative solution is to always ignore the last chain row read which is followed by an end of file record.
                ! But this has to be done within the calling routine and not here to guarantee the output of this procedure is consistent for all chain form types.
                ! Update Nov 2023: After extensive tests, the current approach appears to work even for compromised binary chain files.
#if             ParaDRAM_ENABLED || ParaDISE_ENABLED
                inquire(unit = unit, pos = posCurrent, iostat = err%stat, iomsg = err%msg)
                RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
                read(unit, iostat = err%stat, iomsg = err%msg   ) row%processID             (1) & ! LCOV_EXCL_LINE
                                                                , row%delayedRejectionStage (1) & ! LCOV_EXCL_LINE
                                                                , row%meanAcceptanceRate    (1) & ! LCOV_EXCL_LINE
                                                                , row%proposalAdaptation    (1) & ! LCOV_EXCL_LINE
                                                                , row%burninLocation        (1) & ! LCOV_EXCL_LINE
                                                                , row%sampleWeight          (1) & ! LCOV_EXCL_LINE
                                                                , row%sampleLogFunc         (1) & ! LCOV_EXCL_LINE
                                                                , row%sampleState           (1:ndim,1)
#elif           ParaNest_ENABLED
                read(unit, iostat = err%stat, iomsg = err%msg   ) row%processID             (1) & ! LCOV_EXCL_LINE
                                                                , row%meanAcceptanceRate    (1) & ! LCOV_EXCL_LINE
                                                                , row%domainLogVol          (1) & ! LCOV_EXCL_LINE
                                                                , row%funcLogIntegral       (1) & ! LCOV_EXCL_LINE
                                                                , row%logMaxRelErr          (1) & ! LCOV_EXCL_LINE
                                                                , row%sampleLogWeight       (1) & ! LCOV_EXCL_LINE
                                                                , row%sampleLogFunc         (1) & ! LCOV_EXCL_LINE
                                                                , row%sampleState           (1:ndim,1)
#else
#error          "Unrecognized interface."
#endif
            else blockFileForm ! only if compact or verbose.
                ! Read the row as a string record.
                ! We do so, because we do not know a priori what sort of `sep` we are dealing with.
                call setRecordFrom(unit, record, iostat = err%stat, iomsg = err%msg, ub = lenrec)
                blockParseRecord: if (err%stat == 0) then
                    ! ibeg always refers to the position of the beginning of a column.
                    iend = 1_IK - seplen
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%processID            ; if (iserf(err%stat)) exit blockParseRecord
#if                 ParaDRAM_ENABLED || ParaDISE_ENABLED
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%delayedRejectionStage; if (iserf(err%stat)) exit blockParseRecord
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%meanAcceptanceRate   ; if (iserf(err%stat)) exit blockParseRecord
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%proposalAdaptation   ; if (iserf(err%stat)) exit blockParseRecord
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%burninLocation       ; if (iserf(err%stat)) exit blockParseRecord
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%sampleWeight         ; if (iserf(err%stat)) exit blockParseRecord
#elif               ParaNest_ENABLED
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%meanAcceptanceRate   ; if (iserf(err%stat)) exit blockParseRecord
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%domainLogVol         ; if (iserf(err%stat)) exit blockParseRecord
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%funcLogIntegral      ; if (iserf(err%stat)) exit blockParseRecord
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%logMaxRelErr         ; if (iserf(err%stat)) exit blockParseRecord
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%sampleLogWeight      ; if (iserf(err%stat)) exit blockParseRecord
#else
#error              "Unrecognized interface."
#endif
                    ibeg = iend + seplen; iend = ibeg + index(record(ibeg:lenrec), cfc%sep, kind = IK) - 1; read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%sampleLogFunc        ; if (iserf(err%stat)) exit blockParseRecord
                    ibeg = iend + seplen; iend = lenrec;                                                    read(record(ibeg:iend), *, iostat = err%stat, iomsg = err%msg) row%sampleState          ; if (iserf(err%stat)) exit blockParseRecord
                end if blockParseRecord
            end if blockFileForm

            blockParseSuccess: if (err%stat == 0) then
                cfc%nsam = cfc%nsam + 1_IK
                if (pre_def < cfc%nsam) then
                    pre_def = pre_def * 2_IK
                    RETURN_IF_FAILED(__LINE__,present(pre),SK_"The chain file contains more rows ("//getStr(cfc%nsam)//SK_") than the specified input `pre` ("//getStr(pre_def/2)//SK_"). If you are unsure of the right size, do not specify the optional input argument `pre`.")
                    RETURN_IF_FAILED(__LINE__,isFailedChainResize(cfc, ndim, pre_def, err%msg),SK_"Failed to allocate storage for the chain columns.")
                    if (present(pos)) then
                        call setResized(pos, pre_def, failed = err%occurred, errmsg = err%msg)
                        RETURN_IF_FAILED(__LINE__,err%occurred,err%msg)
                    end if
                end if
                if (present(pos)) pos(cfc%nsam) = posCurrent
#if             ParaDRAM_ENABLED || ParaDISE_ENABLED
                if (isVerbose) then
                    if (isPastFirstRead) then
                        if (1_IK < cfc%nsam) then
                            if (all(row%sampleState(1 : ndim, 1) == cfc%sampleState(1 : ndim, cfc%nsam - 1))) then
                                row%proposalAdaptation(1) = max(row%proposalAdaptation(1), cfc%proposalAdaptation(cfc%nsam - 1))
                                row%sampleWeight(1) = row%sampleWeight(1) + cfc%sampleWeight(cfc%nsam - 1)
                                cfc%nsam = cfc%nsam - 1_IK
                            end if
                        end if
                    else ! happens once only after reading the first line (after header) in the chain file.
                        isPastFirstRead = .true._LK
                    end if
                end if
                cfc%processID             (cfc%nsam)            = row%processID             (1)
                cfc%delayedRejectionStage (cfc%nsam)            = row%delayedRejectionStage (1)
                cfc%meanAcceptanceRate    (cfc%nsam)            = row%meanAcceptanceRate    (1)
                cfc%proposalAdaptation    (cfc%nsam)            = row%proposalAdaptation    (1)
                cfc%burninLocation        (cfc%nsam)            = row%burninLocation        (1)
                cfc%sampleWeight          (cfc%nsam)            = row%sampleWeight          (1)
                cfc%sampleLogFunc         (cfc%nsam)            = row%sampleLogFunc         (1)
                cfc%sampleState           (1 : ndim, cfc%nsam)  = row%sampleState           (1:ndim, 1)
#elif           ParaNest_ENABLED
                cfc%processID             (cfc%nsam)            = row%processID             (1)
                cfc%meanAcceptanceRate    (cfc%nsam)            = row%meanAcceptanceRate    (1)
                cfc%domainLogVol          (cfc%nsam)            = row%domainLogVol          (1)
                cfc%funcLogIntegral       (cfc%nsam)            = row%funcLogIntegral       (1)
                cfc%logMaxRelErr          (cfc%nsam)            = row%logMaxRelErr          (1)
                cfc%sampleLogWeight       (cfc%nsam)            = row%sampleLogWeight       (1)
                cfc%sampleLogFunc         (cfc%nsam)            = row%sampleLogFunc         (1)
                cfc%sampleState           (1 : ndim, cfc%nsam)  = row%sampleState           (1:ndim, 1)
#else
#error          "Unrecognized interface."
#endif
                cycle loopOverRows
            end if blockParseSuccess

            exit loopOverRows

        end do loopOverRows

        ! last read is always end-of-file (normal) or end-of-record error (corrupt file).

        RETURN_IF_FAILED(__LINE__,.not. iserf(err%stat),err%msg)

        ! End of record condition requires ignoring the last read sample. This gets tricky for verbose chain where the IO could fail amid writing the same unweighted state repeatedly.

        if (is_iostat_eor(err%stat)) then
            err%iswarned = .true._LK
            err%msg = SKG_"An end-of-record condition occurred while parsing the contents of the chain file at chain row "//getStr(cfc%nsam)//SKG_": "//trim(err%msg)//& ! LCOV_EXCL_LINE
            SKG_". Assuming the previous chain sampleState as the last in the chain file with the corresponding values: "
            if (0 == cfc%nsam) then
                err%msg = err%msg//SKG_"Nothing. The first row of the chain file is corrupt (incomplete)."
            else
                err%msg = err%msg//getStr(cfc%sampleState(1 : ndim, cfc%nsam - 1))
                ! return the position to the end of the last successful row read.
                if (opened) then
                    if (isBinary) then
                        read(unit, pos = posCurrent, iostat = err%stat, iomsg = err%msg)
                        RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
                    else
                        backspace(unit, iostat = err%stat, iomsg = err%msg)
                        RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
#if                     ParaDRAM_ENABLED || ParaDISE_ENABLED
                        if (isVerbose) then
                            do ibeg = 1, cfc%sampleWeight(cfc%nsam)
                                backspace(unit, iostat = err%stat, iomsg = err%msg)
                                RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
                            end do
                        end if
#elif                   !ParaNest_ENABLED
#error                  "Unrecognized interface."
#endif
                    end if
                end if
            end if
        elseif (is_iostat_end(err%stat)) then
            if (isASCII) then
                ! End of file condition for formatted files requires a backspace to go before the end of file record.
                backspace(unit, iostat = err%stat, iomsg = err%msg)
            elseif (isBinary) then
                ! Ignore the last read row before hitting end of file, as the record may have been incomplete.
                read(unit, pos = posCurrent, iostat = err%stat, iomsg = err%msg)
            end if
            RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
        end if

        ! Clean up the output.

        if (.not. present(pre)) then
            RETURN_IF_FAILED(__LINE__,isFailedChainResize(cfc, ndim, cfc%nsam, err%msg),SK_"Failed to resize the storage for the chain columns.")
            if (present(pos)) then
                call setResized(pos, cfc%nsam, failed = err%occurred, errmsg = err%msg)
                RETURN_IF_FAILED(__LINE__,err%occurred,err%msg)
            end if
        end if

        if (.not. opened) then
            close(unit, iostat = err%stat, iomsg = err%msg)
            RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
        end if

    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !   \todo
    !   The performance of the generic interface might slightly (perhaps negligibly) improve
    !   if the array initializations are only done for the left-over elements when the input argument `pre` is specified.
    !   For now, code clarity by avoiding code duplication wins over the slight performance improvement.
    function isFailedChainResize(cfc, ndim, nrow, errmsg) result(failed)
        use pm_arrayRefill, only: setRefilled
        real(RKG), parameter :: NULL_RK = -huge(0._RKG)
        integer(IK), parameter :: NULL_IK = -huge(0_IK)
        type(CFC_TYPE), intent(inout) :: cfc
        integer(IK), intent(in) :: ndim, nrow
        character(*, SK), intent(out) :: errmsg
        logical(LK) :: failed
        call setRefilled(cfc%processID              , NULL_IK  , nrow           , failed = failed, errmsg = errmsg)
        call setRefilled(cfc%meanAcceptanceRate     , NULL_RK  , nrow           , failed = failed, errmsg = errmsg)
#if     ParaDISE_ENABLED || ParaDRAM_ENABLED
        call setRefilled(cfc%delayedRejectionStage  , NULL_IK  , nrow           , failed = failed, errmsg = errmsg)
        call setRefilled(cfc%proposalAdaptation     , NULL_RK  , nrow           , failed = failed, errmsg = errmsg)
        call setRefilled(cfc%burninLocation         , NULL_IK  , nrow           , failed = failed, errmsg = errmsg)
        call setRefilled(cfc%sampleWeight           , NULL_IK  , nrow           , failed = failed, errmsg = errmsg)
#elif   ParaNest_ENABLED
        call setRefilled(cfc%domainLogVol           , NULL_RK  , nrow           , failed = failed, errmsg = errmsg)
        call setRefilled(cfc%funcLogIntegral        , NULL_RK  , nrow           , failed = failed, errmsg = errmsg)
        call setRefilled(cfc%logMaxRelErr           , NULL_RK  , nrow           , failed = failed, errmsg = errmsg)
        call setRefilled(cfc%sampleLogWeight        , NULL_RK  , nrow           , failed = failed, errmsg = errmsg)
#else
#error  "Unrecognized interface."
#endif
        call setRefilled(cfc%sampleLogFunc          , NULL_RK  , nrow           , failed = failed, errmsg = errmsg)
        call setRefilled(cfc%sampleState            , NULL_RK  , [ndim, nrow]   , failed = failed, errmsg = errmsg)
    end function

    !>  \brief
    !>  Write the chain properties to the chain file.
    !>
    !>  \param[in]  cfc                         :   The object of class [ChainFileContents_type](@ref ChainFileContents_type).
    !>  \param[in]  ibeg                        :   The beginning index of the compact chain beyond which the elements of the chain will be written to the output file.
    !>  \param[in]  iend                        :   The ending index of the compact chain below which the elements of the chain will be written to the output file.
    !>  \param[in]  unit                        :   The unit ID of the chain file to which the header should be written.
    !>  \param[in]  form                        :   The file format of the chain file (`"binary"` vs. `"compact"` vs. `"verbose"`).
    !>  \param[in]  format                      :   The Fortran IO formatting string to be used to read the contents of the chain file (**optional**).
    !>                                              This argument is only required with a non-binary chain file, i.e., when `isBinary = .false.`.
    !> \param[in]   proposalAdaptationPeriod    :   The adaptive update period (**optional**). It must be provided if `form = "verbose"`.
    !>
    !>  \warning
    !>  The file header is written only if the input is not already open.
    !>
    !>  \devnote
    !>  Upon entry, it file is closed already, it will be opened and rewound and will be closed upon return.
    !>  If file is open, it is closed and reopened and rewound, however, it will remain open at the last position upon return.
    function getErrChainWrite(cfc, file, form, format, ibeg, iend, proposalAdaptationPeriod) result(err)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErrChainWrite
#endif
        use pm_strASCII, only: setStrLower
        type(CFC_TYPE)  , intent(in)            :: cfc
        character(*, SK), intent(in)            :: file, form
        integer(IK)     , intent(in), optional  :: ibeg, iend
        character(*, SK), intent(in), optional  :: format
        integer(IK)     , intent(in), optional  :: proposalAdaptationPeriod
        type(err_type)                          :: err

        character(*, SK), parameter             :: PROCEDURE_NAME = MODULE_NAME//SK_"@getErrChainWrite()"
        logical(LK)                             :: opened, isBinary, isCompact, isVerbose, isASCII!, exist
        integer(IK)                             :: i, ibeg_def, iend_def
        character((len(form, IK)), SK)          :: lcform
        integer                                 :: unit

        ibeg_def = 1_IK
        iend_def = cfc%nsam
        if (present(ibeg)) ibeg_def = ibeg
        if (present(iend)) iend_def = iend
        RETURN_IF_FAILED(__LINE__,ibeg_def < 1_IK .or. iend_def < ibeg_def,SK_": The condition `0 < ibeg <= iend <= cfc%nsam` must hold. ibeg, iend, cfc%nsam = "//getStr([ibeg_def, iend_def, cfc%nsam]))

        isBinary = .false._LK
        isCompact = .false._LK
        isVerbose = .false._LK
        lcform = form
        call setStrLower(lcform)
        isVerbose = lcform == SK_"verbose"
        isBinary = lcform == SK_"unformatted" .or. lcform == SK_"binary"
        isASCII = lcform == SK_"formatted" .or. lcform == SK_"compact" .or. isVerbose
        RETURN_IF_FAILED(__LINE__,isASCII .and. .not. present(format),SK_": The optional argument `format` must be present when the file is formatted. form = """//form//SK_"""")
        RETURN_IF_FAILED(__LINE__,isVerbose .and. .not. present(proposalAdaptationPeriod),SK_": The optional argument `proposalAdaptationPeriod` must be present when the file is verbose. form = """//form//SK_"""")
#if     ParaNest_ENABLED
        isCompact = isCompact .or. isVerbose
        isVerbose = .false._LK
#endif
        err%msg = repeat(SKG_" ", 255)
        inquire(file = file, opened = opened, iostat = err%stat, iomsg = err%msg)!, exist = exist
        RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)!.not. exist .or.

        if (.not. opened) then
            if (isBinary) then
                open(newunit = unit, file = file, form = "unformatted", access = "stream", iostat = err%stat, iomsg = err%msg SHARED)
            else
                open(newunit = unit, file = file, form = "formatted", access = "sequential", iostat = err%stat, iomsg = err%msg SHARED)
            end if
            RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
            if (isBinary) then
                write(unit, iostat = err%stat, iomsg = err%msg) cfc%header//achar(0)
            else
                write(unit, "(A)", iostat = err%stat, iomsg = err%msg) cfc%header
            end if
            RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
        end if

        blockWrite: if (isBinary) then
            do i = ibeg_def, iend_def
                write(unit, iostat = err%stat, iomsg = err%msg  ) cfc%processID(i) &
#if                                                             ParaDISE_ENABLED || ParaDRAM_ENABLED
                                                                , cfc%delayedRejectionStage(i) &
                                                                , cfc%meanAcceptanceRate(i) &
                                                                , cfc%proposalAdaptation(i) &
                                                                , cfc%burninLocation(i) &
                                                                , cfc%sampleWeight(i) &
#elif                                                           ParaNest_ENABLED
                                                                , cfc%meanAcceptanceRate(i) &
                                                                , cfc%domainLogVol(i) &
                                                                , cfc%funcLogIntegral(i) &
                                                                , cfc%sampleLogWeight(i) &
#else
#error                                                          "Unrecognized interface."
#endif
                                                                , cfc%sampleLogFunc(i) &
                                                                , cfc%sampleState(:, i)
                RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
            end do
        elseif (isASCII) then blockWrite
            do i = ibeg_def, iend_def
                write(unit, format, iostat = err%stat, iomsg = err%msg  ) cfc%processID(i) &
#if                                                                     ParaDISE_ENABLED || ParaDRAM_ENABLED
                                                                        , cfc%delayedRejectionStage(i) &
                                                                        , cfc%meanAcceptanceRate(i) &
                                                                        , cfc%proposalAdaptation(i) &
                                                                        , cfc%burninLocation(i) &
                                                                        , cfc%sampleWeight(i) &
#elif                                                                   ParaNest_ENABLED
                                                                        , cfc%meanAcceptanceRate(i) &
                                                                        , cfc%domainLogVol(i) &
                                                                        , cfc%funcLogIntegral(i) &
                                                                        , cfc%sampleLogWeight(i) &
#else
#error                                                                  "Unrecognized interface."
#endif
                                                                        , cfc%sampleLogFunc(i) &
                                                                        , cfc%sampleState(:, i)
                RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
            end do
        else
#if         ParaDISE_ENABLED || ParaDRAM_ENABLED
            if (isVerbose) then
                block
                    integer(IK) :: j, counter
                    real(RKG) :: adaptation
                    counter = ibeg_def
                    do i = ibeg_def, iend_def
                        do j = 1, cfc%sampleWeight(i)
                            if (mod(counter, proposalAdaptationPeriod) == 0_IK) then
                                adaptation = cfc%proposalAdaptation(i)
                            else
                                adaptation = 0._RKG
                            end if
                            write(unit, format, iostat = err%stat, iomsg = err%msg  ) cfc%processID(i) &
                                                                                    , cfc%delayedRejectionStage(i) &
                                                                                    , cfc%meanAcceptanceRate(i) &
                                                                                    , cfc%proposalAdaptation(i) &
                                                                                    , cfc%burninLocation(i) &
                                                                                    , 1_IK &
                                                                                    , cfc%sampleLogFunc(i) &
                                                                                    , cfc%sampleState(:, i)
                            RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)
                            counter = counter + 1
                        end do
                    end do
                    exit blockWrite
                end block
            end if
#endif
            RETURN_IF_FAILED(__LINE__,.true.,SK_": Unrecognized input value for the input argument `form`. form = """//form//SK_"""")
        end if blockWrite

        flush(unit)
        if (.not. opened) close(unit, iostat = err%stat, iomsg = err%msg)
        RETURN_IF_FAILED(__LINE__,err%stat /= 0,err%msg)

    end function getErrChainWrite

#undef  isFailedChainResize
#undef  RETURN_IF_FAILED
#undef  chainFileColName
#undef  getErrChainWrite
#undef  getErrChainRead
#undef  CFC_TYPE

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif