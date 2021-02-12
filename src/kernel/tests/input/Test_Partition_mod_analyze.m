clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the actual clusters

%Test_Statistics_mod_writeClusteredPoint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% things to plot

dbscanEnabled = true;
pointsEnabled = 1;
innestPointEnabled = 1;
originalBoundsEnabled = 1;
predictedBoundsEnabled = 1;
distMatEnabled = 1;

methodList = ["", "Kmeans", "MaxDen", "MinVol"];
fileTemplate = "../output/Test_Partition__SAMPLING_METHOD___mod@test_runPartition_2@__FILE_TYPE__@1.txt";

iline = 0;

global DistSqSorted np LogFac

for method = methodList

    file = strrep(fileTemplate,"__SAMPLING_METHOD__",method);
    file = strrep(file,"__FILE_TYPE__","Partition");
    %partition = readPartition(file,"p");
    partition = readCluster(file);

    file = strrep(fileTemplate,"__SAMPLING_METHOD__",method);
    file = strrep(file,"__FILE_TYPE__","ClusteredPoint");
    cluspointExists = isfile(file);
    if cluspointExists
        %cluspoint = readPartition(file,"cp");
        cluspoint = readCluster(file);
    end

    for icase = 1:partition.count

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% plot ellipsoids
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        id1 = 1;
        id2 = 2;
        dum = 1:partition.case{icase}.nd;
        mask = dum(dum==id1 | dum==id2);

        ProxyCenterIndex = dsearchn ( partition.case{icase}.Point ...
                                    , cluspoint.case{1}.Center' ...
                                    );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % knn analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if distMatEnabled

            nd = double(partition.case{icase}.nd);
            np = double(partition.case{icase}.np);
            DistSq = squareform(pdist(partition.case{icase}.Point));
            DistSqSorted = sort(DistSq,"ascend");
            DistSqCDF = zeros(np,np);
            npIndexMinSumDistSqCDF = 1; %floor(sqrt(np));

            % optimize the density

            LogFac = zeros(np,1);
            LogFac(1) = 0.;
            for ip = 2:np
                LogFac(ip) = LogFac(ip-1) + log(ip);
            end

            % Identify hubs

            [~,MinDistIndex] = min(DistSq);
            %MinDistIndex = [MinDistIndex, 1:np];
            [UniqueValue,ia,ic] = unique(MinDistIndex);
            UniqueCount = accumarray(ic,1);
            %UniqueValueCount = [UniqueValue, UniqueCount];
            max(UniqueCount)

            figure; hold on; box on;

            % plot CDFs

            skip = 1;
            for ip = 1:skip:np
                DistSqCDF(1:np,ip) = cumsum(DistSqSorted(1:np,ip));
                DistSqCDF(1:np,ip) = DistSqCDF(1:np,ip) / DistSqCDF(np,ip);
                p = plot( DistSqCDF(1:np,ip) ...
                        , "linewidth", 0.1 ...
                        , "color", [255 0 0 25] / 255 ...
                        );
                %p.Color(4) = 0.1;
            end

            % plot the reference x^2 curve

            x = 0:0.5:np;
            plot( x ...
                , x.^2 / np^2 ...
                , "--" ...
                , "linewidth", 2 ...
                , "color", "black" ...
                );

            % plot CDFs of closest points to centers

            LineWidth = max(1, 3 * cluspoint.case{1}.Size / max(cluspoint.case{1}.Size));
            count = 0;
            for ip = ProxyCenterIndex
                count = count + 1;
                plot( DistSqCDF(1:np,ip) ...
                    , "-" ...
                    , "linewidth", LineWidth(count) ...
                    , "color", "black" ...
                    );
            end

            % Find and plot the point with the least sum(DistSqCDF)

            SumDistSqCDF = zeros(np);
            for ip = 1:np
                SumDistSqCDF(ip) = sum(DistSqCDF(1:np,ip));
            end
            [~, IndexMinSumDistSqCDF] = sort(SumDistSqCDF);
            plot( DistSqCDF(1:np,IndexMinSumDistSqCDF(1)) ...
                , "-" ...
                , "linewidth", 2.5 ...
                , "color", "blue" ...
                );
            
            hold off;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot points and clusters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if pointsEnabled || (innestPointEnabled && cluspointExists) || distMatEnabled || (originalBoundsEnabled && cluspointExists) || dbscanEnabled || (predictedBoundsEnabled && cluspointExists)

            figure; hold on; box on;

            if pointsEnabled

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
                        , 'markersize', 15 ...
                        );
                end

                for ip = ProxyCenterIndex
                    plot( partition.case{icase}.Point(ip,id1) ...
                        , partition.case{icase}.Point(ip,id2) ...
                        , "." ...
                        , 'markersize', 40 ...
                        , "color", "black" ...
                        );
                end

                %hsizeMax = cluspoint.case{icase}.HubEdgeCount(1);
                for ih = 1:1 %length(cluspoint.case{icase}.HubPoint(:,id1))
                    %if cluspoint.case{icase}.HubEdgeCount(ih) == hsizeMax
                    if cluspoint.case{icase}.HubEdgeCount(ih) > 1 %cluspoint.case{icase}.nd
                        plot( cluspoint.case{icase}.HubPoint(ih,id1) ...
                            , cluspoint.case{icase}.HubPoint(ih,id2) ...
                            , "*" ...
                            , 'markersize', 30 ...
                            , "color", "black" ...
                            );
                    else
                        break
                    end
                end

            end


            % plot the innestPoint

            if innestPointEnabled && cluspointExists
                plot( cluspoint.case{1}.InnestPoint(1) ...
                    , cluspoint.case{1}.InnestPoint(2) ...
                    , ".", "markerSize", 40 ...
                    , "color", "red" ...
                    )
            end

            % plot the point witth the least sum(DistSqCDF)

            if distMatEnabled
                plot( partition.case{icase}.Point(IndexMinSumDistSqCDF(1:npIndexMinSumDistSqCDF),id1) ...
                    , partition.case{icase}.Point(IndexMinSumDistSqCDF(1:npIndexMinSumDistSqCDF),id2) ...
                    , "." ...
                    , 'markersize', 50 ...
                    , "color", "blue" ...
                    );
            end

            % get the actual ellipsoid boundary

            if originalBoundsEnabled && cluspointExists
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

            % add dbscan IDs

            if dbscanEnabled
                idx = dbscan( cluspoint.case{1}.Point ...
                            , 0.5 ...
                            , 3 ...
                            ..., "distance", 'mahalanobis' ...
                            ..., "distance", 'correlation' ...
                            );
                gscatter(cluspoint.case{1}.Point(:,1),cluspoint.case{1}.Point(:,2),idx);
                disp([sum(idx<0), length(cluspoint.case{1}.Point)]);
            end

            % get the predicted ellipsoid boundary

            if predictedBoundsEnabled && cluspointExists
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
            end

            title(method);

            hold off;

        end

    end

end
