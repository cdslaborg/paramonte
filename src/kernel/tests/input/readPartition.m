function partition = readPartition(file,type,delim)

    if nargin==2
        delim = ",";
    elseif nargin~=3
        error("The numbber of input arguments must be either two (file,type) or three (file,type,delim).");
    end

    contents = strrep(fileread(file),char(13),'');
    lineList = strsplit(contents,char(10));
    lineListLen = length(lineList);

    partition.count = count(contents, "CholeskyLower");
    partition.case = cell(partition.count,1);

    iline = 0;
    for icase = 1:partition.count

        if strcmpi(type,"cp")

            %%%% dist

            iline = iline + 2;
            dumCell = textscan(lineList{iline}, "%s", "Delimiter", ",");
            partition.case{icase}.dist = dumCell{1};

        end

        %%%% nd, np, nc

        iline = iline + 2;
        dumCell = textscan(lineList{iline}, "%d", "Delimiter", delim);
        partition.case{icase}.nd = dumCell{1}(1);
        partition.case{icase}.np = dumCell{1}(2);
        partition.case{icase}.nc = dumCell{1}(3);

        %%%% Size

        iline = iline + 2;
        dumCell = textscan(lineList{iline}, "%d", "Delimiter", delim);
        partition.case{icase}.Size = dumCell{1};

        %%%% Center

        iline = iline + 2;
        dumCell = textscan(lineList{iline}, "%f", "Delimiter", delim);
        partition.case{icase}.Center = reshape( dumCell{1} , [partition.case{icase}.nd, partition.case{icase}.nc] );

        %%%% LogVolume

        iline = iline + 2;
        dumCell = textscan(lineList{iline}, "%f", "Delimiter", delim);
        partition.case{icase}.LogVolume = dumCell{1};

        %%%% CholeskyLower

        iline = iline + 2;
        dumCell = textscan(lineList{iline}, "%f", "Delimiter", delim);
        CholeskyLower = zeros(partition.case{icase}.nd, partition.case{icase}.nd);
        CholeskyUpper = zeros(partition.case{icase}.nd, partition.case{icase}.nd);
        partition.case{icase}.CovMatUpper = zeros   ( partition.case{icase}.nd ...
                                                    , partition.case{icase}.nd ...
                                                    , partition.case{icase}.nc ...
                                                    );
        icount = 0;
        for ic = 1:partition.case{icase}.nc
            CholeskyLower(:,:) = 0.;
            CholeskyUpper(:,:) = 0.;
            for j = 1:partition.case{icase}.nd
                for i = j:partition.case{icase}.nd
                    icount = icount + 1;
                    CholeskyLower(i,j) = dumCell{1}(icount);
                    CholeskyUpper(j,i) = dumCell{1}(icount);
                end
            end
            partition.case{icase}.CovMatUpper(:,:,ic) = CholeskyLower * CholeskyUpper;
        end

        %%%% Point

        iline = iline + 2;
        dumCell = textscan(lineList{iline}, "%f", "Delimiter", delim);
        partition.case{icase}.Point = reshape( dumCell{1}, [partition.case{icase}.nd, partition.case{icase}.np] )';

        %%%% Membership

        iline = iline + 2;
        dumCell = textscan(lineList{iline}, "%d", "Delimiter", delim);
        partition.case{icase}.Membership = dumCell{1};

    end

end