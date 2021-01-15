! The MACRO variable "NAMELIST" will have to be replaced by the sampler name: ParaDRAM, ParaDISE, ParaXXXX

        ! ParaMonte variables
        namelist /NAMELIST/ sampleSize
        namelist /NAMELIST/ randomSeed
        namelist /NAMELIST/ description
        namelist /NAMELIST/ outputFileName
        namelist /NAMELIST/ outputDelimiter
        namelist /NAMELIST/ chainFileFormat
        namelist /NAMELIST/ variableNameList
        namelist /NAMELIST/ restartFileFormat
        namelist /NAMELIST/ outputColumnWidth
        namelist /NAMELIST/ overwriteRequested
        namelist /NAMELIST/ outputRealPrecision
        namelist /NAMELIST/ silentModeRequested
        namelist /NAMELIST/ domainLowerLimitVec
        namelist /NAMELIST/ domainUpperLimitVec
        namelist /NAMELIST/ parallelizationModel
        namelist /NAMELIST/ progressReportPeriod
        namelist /NAMELIST/ inputFileHasPriority
        namelist /NAMELIST/ targetAcceptanceRate
        namelist /NAMELIST/ mpiFinalizeRequested
        namelist /NAMELIST/ maxNumDomainCheckToWarn
        namelist /NAMELIST/ maxNumDomainCheckToStop
        namelist /NAMELIST/ systemInfoFilePath
        namelist /NAMELIST/ interfaceType

        ! ParaMCMC variables

        namelist /NAMELIST/ chainSize
        namelist /NAMELIST/ startPointVec
        namelist /NAMELIST/ proposalStartStdVec
        namelist /NAMELIST/ proposalStartCorMat
        namelist /NAMELIST/ proposalStartCovMat
        namelist /NAMELIST/ sampleRefinementCount
        namelist /NAMELIST/ sampleRefinementMethod
        namelist /NAMELIST/ randomStartPointRequested
        namelist /NAMELIST/ randomStartPointDomainLowerLimitVec
        namelist /NAMELIST/ randomStartPointDomainUpperLimitVec

#if defined PARADRAM || defined PARADISE
        namelist /NAMELIST/ scaleFactor
        namelist /NAMELIST/ proposalModel
#endif

        ! ParaDRAM variables

        namelist /NAMELIST/ greedyAdaptationCount
        namelist /NAMELIST/ delayedRejectionCount
        namelist /NAMELIST/ burninAdaptationMeasure
        namelist /NAMELIST/ delayedRejectionScaleFactorVec

#if defined PARADRAM || defined PARADISE
        namelist /NAMELIST/ adaptiveUpdateCount
        namelist /NAMELIST/ adaptiveUpdatePeriod
#endif

        ! ParaNest variables

        namelist /NAMELIST/ tightness
        namelist /NAMELIST/ tolerance
        namelist /NAMELIST/ liveSampleSize
        namelist /NAMELIST/ inclusionFraction
        namelist /NAMELIST/ mahalSqWeightExponent
        namelist /NAMELIST/ stabilizationRequested
        namelist /NAMELIST/ MaxAllowedKmeansFailure
        namelist /NAMELIST/ maxAllowedMinVolFailure
        namelist /NAMELIST/ maxKvolumeLoopRecursion

#if defined PARANEST
        namelist /NAMELIST/ scaleFactor
        namelist /NAMELIST/ proposalModel
        namelist /NAMELIST/ adaptiveUpdateCount
        namelist /NAMELIST/ adaptiveUpdatePeriod
#endif
