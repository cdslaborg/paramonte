clear all;
close all;

file = "../output/Test_MinVolPartition_mod@test_runMinVolPartition_2.1.txt";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contents = strrep(fileread(file),char(13),'');
lineList = strsplit(contents,char(10));
lineListLen = length(lineList);

proposalUpdates.count = count(contents, "ClusterSize");
proposalUpdates.update = cell(proposalUpdates.count,1);

iline = 0;
for iupdate = 1:proposalUpdates.count

    %%%% nd, np, nc

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%d", "Delimiter", ",");
    proposalUpdates.update{iupdate}.nd = dumCell{1}(1);
    proposalUpdates.update{iupdate}.np = dumCell{1}(2);
    proposalUpdates.update{iupdate}.nc = dumCell{1}(3);

    %%%% ClusterSize

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%d", "Delimiter", ",");
    proposalUpdates.update{iupdate}.Size = dumCell{1};

    %%%% ClusterCenter

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%f", "Delimiter", ",");
    proposalUpdates.update{iupdate}.Center = reshape( dumCell{1} , [proposalUpdates.update{iupdate}.nd, proposalUpdates.update{iupdate}.nc] );

    %%%% ClusterLogVolume

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%f", "Delimiter", ",");
    proposalUpdates.update{iupdate}.LogVolume = dumCell{1};

    %%%% ClusterCholeskyLower

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%f", "Delimiter", ",");
    CholeskyLower = zeros(proposalUpdates.update{iupdate}.nd, proposalUpdates.update{iupdate}.nd);
    CholeskyUpper = zeros(proposalUpdates.update{iupdate}.nd, proposalUpdates.update{iupdate}.nd);
    proposalUpdates.update{iupdate}.CovMatUpper = zeros ( proposalUpdates.update{iupdate}.nd ...
                                                        , proposalUpdates.update{iupdate}.nd ...
                                                        , proposalUpdates.update{iupdate}.nc ...
                                                        );
    icount = 0;
    for ic = 1:proposalUpdates.update{iupdate}.nc
        CholeskyLower(:,:) = 0.;
        CholeskyUpper(:,:) = 0.;
        for j = 1:proposalUpdates.update{iupdate}.nd
            for i = j:proposalUpdates.update{iupdate}.nd
                icount = icount + 1;
                CholeskyLower(i,j) = dumCell{1}(icount);
                CholeskyUpper(j,i) = dumCell{1}(icount);
            end
        end
        proposalUpdates.update{iupdate}.CovMatUpper(:,:,ic) = CholeskyLower * CholeskyUpper;
    end

    %%%% Point

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%f", "Delimiter", ",");
    proposalUpdates.update{iupdate}.Point = reshape( dumCell{1}, [proposalUpdates.update{iupdate}.nd, proposalUpdates.update{iupdate}.np] )';

    %%%% Membership

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%d", "Delimiter", ",");
    proposalUpdates.update{iupdate}.Membership = dumCell{1};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% plot ellipsoids
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure; hold on; box on;

    plot( proposalUpdates.update{iupdate}.Point(:,1) ...
        , proposalUpdates.update{iupdate}.Point(:,2) ...
        , "." ...
        );

    % get ellipsoid boundary

    for ic = 1:proposalUpdates.update{iupdate}.nc
        bcrd = getEllipsoidBoundary ( proposalUpdates.update{iupdate}.CovMatUpper(:,:,ic) ... covMat
                                    , proposalUpdates.update{iupdate}.Center(:,ic) ... meanVec
                                    , 50 ... npoint
                                    );
        plot( bcrd(:,1) ...
            , bcrd(:,2) ...
            , "color", "red" ...
            );
    end

    hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

