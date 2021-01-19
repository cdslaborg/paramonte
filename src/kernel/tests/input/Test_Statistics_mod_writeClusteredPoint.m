%clear all;
%close all;

file = "../output/Test_MinVolPartition_mod@test_runMinVolPartition_2.ClusteredPoint.1.txt";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contents = strrep(fileread(file),char(13),'');
lineList = strsplit(contents,char(10));
lineListLen = length(lineList);

ClusteredPoint.count = count(contents, "CholeskyLower");
ClusteredPoint.update = cell(ClusteredPoint.count,1);

iline = 0;
for iupdate = 1:ClusteredPoint.count

    %%%% dist

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%s", "Delimiter", ",");
    ClusteredPoint.update{iupdate}.dist = dumCell{1};

    %%%% nd, np, nc

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%d", "Delimiter", ",");
    ClusteredPoint.update{iupdate}.nd = dumCell{1}(1);
    ClusteredPoint.update{iupdate}.np = dumCell{1}(2);
    ClusteredPoint.update{iupdate}.nc = dumCell{1}(3);

    %%%% Size

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%d", "Delimiter", ",");
    ClusteredPoint.update{iupdate}.Size = dumCell{1};

    %%%% Center

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%f", "Delimiter", ",");
    ClusteredPoint.update{iupdate}.Center = reshape( dumCell{1} , [ClusteredPoint.update{iupdate}.nd, ClusteredPoint.update{iupdate}.nc] );

    %%%% LogVolume

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%f", "Delimiter", ",");
    ClusteredPoint.update{iupdate}.LogVolume = dumCell{1};

    %%%% CholeskyLower

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%f", "Delimiter", ",");
    CholeskyLower = zeros(ClusteredPoint.update{iupdate}.nd, ClusteredPoint.update{iupdate}.nd);
    CholeskyUpper = zeros(ClusteredPoint.update{iupdate}.nd, ClusteredPoint.update{iupdate}.nd);
    ClusteredPoint.update{iupdate}.CovMatUpper = zeros ( ClusteredPoint.update{iupdate}.nd ...
                                                        , ClusteredPoint.update{iupdate}.nd ...
                                                        , ClusteredPoint.update{iupdate}.nc ...
                                                        );
    icount = 0;
    for ic = 1:ClusteredPoint.update{iupdate}.nc
        CholeskyLower(:,:) = 0.;
        CholeskyUpper(:,:) = 0.;
        for j = 1:ClusteredPoint.update{iupdate}.nd
            for i = j:ClusteredPoint.update{iupdate}.nd
                icount = icount + 1;
                CholeskyLower(i,j) = dumCell{1}(icount);
                CholeskyUpper(j,i) = dumCell{1}(icount);
            end
        end
        ClusteredPoint.update{iupdate}.CovMatUpper(:,:,ic) = CholeskyLower * CholeskyUpper;
    end

    %%%% Point

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%f", "Delimiter", ",");
    ClusteredPoint.update{iupdate}.Point = reshape( dumCell{1}, [ClusteredPoint.update{iupdate}.nd, ClusteredPoint.update{iupdate}.np] )';

    %%%% Membership

    iline = iline + 2;
    dumCell = textscan(lineList{iline}, "%d", "Delimiter", ",");
    ClusteredPoint.update{iupdate}.Membership = dumCell{1};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% plot ellipsoids
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure; hold on; box on;

    plot( ClusteredPoint.update{iupdate}.Point(:,1) ...
        , ClusteredPoint.update{iupdate}.Point(:,2) ...
        , "." ...
        );

    % get ellipsoid boundary

    for ic = 1:ClusteredPoint.update{iupdate}.nc
        bcrd = getEllipsoidBoundary ( ClusteredPoint.update{iupdate}.CovMatUpper(:,:,ic) ... covMat
                                    , ClusteredPoint.update{iupdate}.Center(:,ic) ... meanVec
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

