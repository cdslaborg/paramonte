%clear all;
%close all;

outDir = "../output/";

KmeansFileList = dir("../output/Test_Kmeans_mod*");

lenKmeansFileList = length(KmeansFileList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifile = 1:lenKmeansFileList

    if ~KmeansFileList(ifile).isdir

        Kmeans = readCluster( fullfile(KmeansFileList(ifile).folder,KmeansFileList(ifile).name) );

        for icase = 1:Kmeans.count

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% plot ellipsoids
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            id1 = 1;
            id2 = 2;
            dum = 1:Kmeans.case{icase}.nd;
            mask = dum(dum==id1 | dum==id2);

            figure; hold on; box on;

            for ic = 1:Kmeans.case{icase}.nc
                Mask = Kmeans.case{icase}.Membership == ic;
                plot( Kmeans.case{icase}.Point(id1, Mask) ...
                    , Kmeans.case{icase}.Point(id2, Mask) ...
                    , "." ...
                    , "markersize", 20 ...
                    );
            end

            % Plot the cluster centers

            plot( Kmeans.case{icase}.Center(id1,:) ...
                , Kmeans.case{icase}.Center(id2,:) ...
                , "." ...
                , "color", "black" ...
                , "markersize", 20 ...
                );

            % get the actual ellipsoid boundary

            if isfield(Kmeans.case{icase},"CovMatUpper")
                for ic = 1:Kmeans.case{icase}.nc
                    bcrd = getEllipsoidBoundary ( Kmeans.case{icase}.CovMatUpper(mask,mask,ic) * Kmeans.case{icase}.ScaleFactorSq(ic) ... covMat
                                                , Kmeans.case{icase}.Center(mask,ic) ... meanVec
                                                , 100 ... npoint
                                                );
                    plot( bcrd(:,1) ...
                        , bcrd(:,2) ...
                        , "color", "black" ...
                        );
                end
            end

            hold off;

        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

