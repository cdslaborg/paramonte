function Contents = readCluster(file,delim)

    if nargin==1
        delim = ",";
    elseif nargin~=2
        error("The numbber of input arguments must be either one (file) or two (file,type).");
    end

    fileContentsString = strrep(fileread(file),char(13),'');
    LineList = strsplit(fileContentsString,char(10));
    lenLineList = length(LineList);

    %%%% First, determine the data labels

    LabList = strings(lenLineList, 1);
    uniqueLabelCount = 0;
    LabelCount = struct();
    for iline = 1:lenLineList
        if length(LineList{iline}) > 0 && isletter(LineList{iline}(1))
            if sum(contains(LabList,LineList{iline}))
                LabelCount.(LabList(uniqueLabelCount)) = LabelCount.(LabList(uniqueLabelCount)) + 1;
            else
                uniqueLabelCount = uniqueLabelCount + 1;
                LabList(uniqueLabelCount) = string(LineList{iline});
                LabelCount.(LabList(uniqueLabelCount)) = 0;
            end
        end
    end
    LabList = unique(LabList(1:uniqueLabelCount));
    lenLabList = length(LabList);

    %%%% Ensure unique labels appear with equal frequency.

    LabelCountVec = zeros(lenLabList,1);
    for ilab = 1:lenLabList
        LabelCountVec(ilab) = LabelCount.(LabList(ilab));
    end

    if any(LabelCountVec ~= LabelCountVec(1))
        disp(LabelCount);
        disp(LabelCountVec);
        error("Inconciststent number of labels found in the input file.");
    end

    Contents.count = count(fileContentsString, LabList(1));
    Contents.case = cell(Contents.count,1);

    iline = 0;
    for icase = 1:Contents.count

        for ilab = 1:lenLabList

            iline = iline + 1;

            for jlab = 1:lenLabList

                label = LabList(jlab);
                if strcmp(LineList{iline}, label)

                    iline = iline + 1;

                    if sum(label == ["nd", "np", "nc", "nt", "niter", "nfail"])
                        format = "%d";
                    elseif sum(label == ["dist"])
                        format = "%s";
                    else
                        format = "%f";
                    end

                    dumCell = textscan(LineList{iline}, format, "Delimiter", delim);

                    if strcmp(label, "CholeskyLower")

                        CholeskyLower = zeros(Contents.case{icase}.nd, Contents.case{icase}.nd);
                        CholeskyUpper = zeros(Contents.case{icase}.nd, Contents.case{icase}.nd);
                        Contents.case{icase}.CovMatUpper = zeros(Contents.case{icase}.nd, Contents.case{icase}.nd, Contents.case{icase}.nc);
                        icount = 0;
                        for ic = 1:Contents.case{icase}.nc
                            CholeskyLower(:,:) = 0.;
                            CholeskyUpper(:,:) = 0.;
                            for j = 1:Contents.case{icase}.nd
                                for i = j:Contents.case{icase}.nd
                                    icount = icount + 1;
                                    CholeskyLower(i,j) = dumCell{1}(icount);
                                    CholeskyUpper(j,i) = dumCell{1}(icount);
                                end
                            end
                            Contents.case{icase}.CovMatUpper(:,:,ic) = CholeskyLower * CholeskyUpper;
                        end

                    elseif strcmp(label, "Center") % assumes (nd,nc) shape

                        Contents.case{icase}.(label) = reshape( dumCell{1} , [Contents.case{icase}.nd, Contents.case{icase}.nc] );

                    elseif strcmp(label, "Point") || strcmp(label, "HubPoint") % assumes (nd,np) shape

                        npdum = length(dumCell{1}) / Contents.case{icase}.nd;
                        Contents.case{icase}.(label) = reshape( dumCell{1}, [Contents.case{icase}.nd, npdum] )';

                    elseif strcmp(label, "MahalSq") % assumes (np,nc) shape

                        Contents.case{icase}.(label) = reshape( dumCell{1}, [Contents.case{icase}.np, Contents.case{icase}.nc] );

                    else

                        Contents.case{icase}.(label) = dumCell{1};

                    end

                    break

                end

            end % jlab

        end % ilab

    end % icase

end % function