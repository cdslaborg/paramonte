function namelist = getInputFile(self)

    isParaDRAM = strcmp(self.methodName,"ParaDRAM");

    if isempty(self.inputFile)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% begin namelist generation from arguments
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SpecBase = SpecBaseVerify(self.objectName);
        SpecMCMC = SpecMCMCVerify(self.objectName);
        SpecDRAM = SpecDRAMVerify(self.objectName);
        namelist = "";

        % set up outputFileName

        if isempty(self.spec.outputFileName)
            self.spec.outputFileName = fullfile( string(pwd) , self.genOutputFileName() );
        else
            self.spec.outputFileName = convertStringsToChars(SpecBase.outputFileName(self.spec.outputFileName));
            if endsWith(self.spec.outputFileName,'\') || endsWith(self.spec.outputFileName,"/")
                self.spec.outputFileName = fullfile( string(GetFullPath(self.spec.outputFileName,'lean')) , self.genOutputFileName() );
            end
            self.spec.outputFileName = string(outputFileName);
        end

        % ParaMonte variables

        if  ~isempty(self.spec.sampleSize                          );  namelist = namelist + SpecBase.sampleSize                           (self.spec.sampleSize                         ); end
        if  ~isempty(self.spec.randomSeed                          );  namelist = namelist + SpecBase.randomSeed                           (self.spec.randomSeed                         ); end
        if  ~isempty(self.spec.description                         );  namelist = namelist + SpecBase.description                          (self.spec.description                        ); end
        if  ~isempty(self.spec.outputFileName                      );  namelist = namelist + SpecBase.outputFileName                       (self.spec.outputFileName                     ); end
        if  ~isempty(self.spec.outputDelimiter                     );  namelist = namelist + SpecBase.outputDelimiter                      (self.spec.outputDelimiter                    ); end
        if  ~isempty(self.spec.chainFileFormat                     );  namelist = namelist + SpecBase.chainFileFormat                      (self.spec.chainFileFormat                    ); end
        if  ~isempty(self.spec.variableNameList                    );  namelist = namelist + SpecBase.variableNameList                     (self.spec.variableNameList                   ); end
        if  ~isempty(self.spec.restartFileFormat                   );  namelist = namelist + SpecBase.restartFileFormat                    (self.spec.restartFileFormat                  ); end
        if  ~isempty(self.spec.outputColumnWidth                   );  namelist = namelist + SpecBase.outputColumnWidth                    (self.spec.outputColumnWidth                  ); end
        if  ~isempty(self.spec.outputRealPrecision                 );  namelist = namelist + SpecBase.outputRealPrecision                  (self.spec.outputRealPrecision                ); end
        if  ~isempty(self.spec.silentModeRequested                 );  namelist = namelist + SpecBase.silentModeRequested                  (self.spec.silentModeRequested                ); end
        if  ~isempty(self.spec.domainLowerLimitVec                 );  namelist = namelist + SpecBase.domainLowerLimitVec                  (self.spec.domainLowerLimitVec                ); end
        if  ~isempty(self.spec.domainUpperLimitVec                 );  namelist = namelist + SpecBase.domainUpperLimitVec                  (self.spec.domainUpperLimitVec                ); end
        if  ~isempty(self.spec.parallelizationModel                );  namelist = namelist + SpecBase.parallelizationModel                 (self.spec.parallelizationModel               ); end
        if  ~isempty(self.spec.progressReportPeriod                );  namelist = namelist + SpecBase.progressReportPeriod                 (self.spec.progressReportPeriod               ); end
        if  ~isempty(self.spec.targetAcceptanceRate                );  namelist = namelist + SpecBase.targetAcceptanceRate                 (self.spec.targetAcceptanceRate               ); end
        if  ~isempty(self.spec.mpiFinalizeRequested                );  namelist = namelist + SpecBase.mpiFinalizeRequested                 (self.spec.mpiFinalizeRequested               ); end
        if  ~isempty(self.spec.maxNumDomainCheckToWarn             );  namelist = namelist + SpecBase.maxNumDomainCheckToWarn              (self.spec.maxNumDomainCheckToWarn            ); end
        if  ~isempty(self.spec.maxNumDomainCheckToStop             );  namelist = namelist + SpecBase.maxNumDomainCheckToStop              (self.spec.maxNumDomainCheckToStop            ); end
        if  isParaDRAM
        if  ~isempty(self.spec.chainSize                           );  namelist = namelist + SpecMCMC.chainSize                            (self.spec.chainSize                          ); end
        if  ~isempty(self.spec.scaleFactor                         );  namelist = namelist + SpecMCMC.scaleFactor                          (self.spec.scaleFactor                        ); end
        if  ~isempty(self.spec.startPointVec                       );  namelist = namelist + SpecMCMC.startPointVec                        (self.spec.startPointVec                      ); end
        if  ~isempty(self.spec.proposalModel                       );  namelist = namelist + SpecMCMC.proposalModel                        (self.spec.proposalModel                      ); end
        if  ~isempty(self.spec.proposalStartCovMat                 );  namelist = namelist + SpecMCMC.proposalStartCovMat                  (self.spec.proposalStartCovMat                ); end
        if  ~isempty(self.spec.proposalStartCorMat                 );  namelist = namelist + SpecMCMC.proposalStartCorMat                  (self.spec.proposalStartCorMat                ); end
        if  ~isempty(self.spec.proposalStartStdVec                 );  namelist = namelist + SpecMCMC.proposalStartStdVec                  (self.spec.proposalStartStdVec                ); end
        if  ~isempty(self.spec.sampleRefinementCount               );  namelist = namelist + SpecMCMC.sampleRefinementCount                (self.spec.sampleRefinementCount              ); end
        if  ~isempty(self.spec.sampleRefinementMethod              );  namelist = namelist + SpecMCMC.sampleRefinementMethod               (self.spec.sampleRefinementMethod             ); end
        if  ~isempty(self.spec.randomStartPointRequested           );  namelist = namelist + SpecMCMC.randomStartPointRequested            (self.spec.randomStartPointRequested          ); end
        if  ~isempty(self.spec.randomStartPointDomainLowerLimitVec );  namelist = namelist + SpecMCMC.randomStartPointDomainLowerLimitVec  (self.spec.randomStartPointDomainLowerLimitVec); end
        if  ~isempty(self.spec.randomStartPointDomainUpperLimitVec );  namelist = namelist + SpecMCMC.randomStartPointDomainUpperLimitVec  (self.spec.randomStartPointDomainUpperLimitVec); end
        end
        if  isParaDRAM
        if  ~isempty(self.spec.adaptiveUpdateCount                 );  namelist = namelist + SpecDRAM.adaptiveUpdateCount                  (self.spec.adaptiveUpdateCount                ); end
        if  ~isempty(self.spec.adaptiveUpdatePeriod                );  namelist = namelist + SpecDRAM.adaptiveUpdatePeriod                 (self.spec.adaptiveUpdatePeriod               ); end
        if  ~isempty(self.spec.greedyAdaptationCount               );  namelist = namelist + SpecDRAM.greedyAdaptationCount                (self.spec.greedyAdaptationCount              ); end
        if  ~isempty(self.spec.delayedRejectionCount               );  namelist = namelist + SpecDRAM.delayedRejectionCount                (self.spec.delayedRejectionCount              ); end
        if  ~isempty(self.spec.burninAdaptationMeasure             );  namelist = namelist + SpecDRAM.burninAdaptationMeasure              (self.spec.burninAdaptationMeasure            ); end
        if  ~isempty(self.spec.delayedRejectionScaleFactorVec      );  namelist = namelist + SpecDRAM.delayedRejectionScaleFactorVec       (self.spec.delayedRejectionScaleFactorVec     ); end
        end

        namelist = "&" + self.methodName + " " + namelist + SpecBase.interfaceType() + "/";

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% end namelist generation from arguments
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else

        if exist(self.inputFile,"file")
            namelist = string(GetFullPath(self.inputFile,'lean'));
        else
            namelist = string(self.inputFile);
        end
        if ~self.mpiEnabled
            self.Err.msg = "Input namelist file is given by the user." + newline + ...
                           "All simulation specifications will be read from the input file.";
            self.Err.warn();
        end
        
    end

end