clear all;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the actual clusters

%Test_Statistics_mod_writeClusteredPoint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methodList = ["", "Kmeans", "MaxDen", "MinVol"];
fileTemplate = "../output/Test_Partition__SAMPLING_METHOD___mod@test_runPartition_2@__FILE_TYPE__@1.txt";

iline = 0;

for method = methodList

    file = strrep(fileTemplate,"__SAMPLING_METHOD__",method);
    file = strrep(file,"__FILE_TYPE__","Partition");
    partition = readPartition(file,"p");

    file = strrep(fileTemplate,"__SAMPLING_METHOD__",method);
    file = strrep(file,"__FILE_TYPE__","ClusteredPoint");
    cluspointExists = isfile(file);
    if cluspointExists
        cluspoint = readPartition(file,"cp");
    end

    for icase = 1:partition.count

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% plot ellipsoids
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        id1 = 1;
        id2 = 2;
        dum = 1:partition.case{icase}.nd;
        mask = dum(dum==id1 | dum==id2);

        figure; hold on; box on;

        istart = 1;
        iend = partition.case{icase}.Size(1);
        for ic = 1:partition.case{icase}.nc
            if ic > 1
                istart = istart + partition.case{icase}.Size(ic-1);
                iend = iend + partition.case{icase}.Size(ic);
            end
            plot( partition.case{icase}.Point(istart:iend,id1) ...
                , partition.case{icase}.Point(istart:iend,id2) ...
                , "." ...
                , 'markersize', 20 ...
                );
        end

        % get the actual ellipsoid boundary

        if cluspointExists
            for ic = 1:cluspoint.case{icase}.nc
                bcrd = getEllipsoidBoundary ( cluspoint.case{icase}.CovMatUpper(mask,mask,ic) ... covMat
                                            , cluspoint.case{icase}.Center(mask,ic) ... meanVec
                                            , 100 ... npoint
                                            );
                plot( bcrd(:,1) ...
                    , bcrd(:,2) ...
                    , "color", "green" ...
                    );
            end
        end

        % get the predicted ellipsoid boundary

        for ic = 1:partition.case{icase}.nc
            bcrd = getEllipsoidBoundary ( partition.case{icase}.CovMatUpper(mask,mask,ic) ... covMat
                                        , partition.case{icase}.Center(mask,ic) ... meanVec
                                        , 100 ... npoint
                                        );
            plot( bcrd(:,1) ...
                , bcrd(:,2) ...
                , "color", "red" ...
                );
        end

        title(method);

        hold off;

    end

end
