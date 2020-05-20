%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonte: plain powerful parallel Monte Carlo library.
%
%   Copyright (C) 2012-present, The Computational Data Science Lab
%
%   This file is part of the ParaMonte library.
%
%   ParaMonte is free software: you can redistribute it and/or modify it 
%   under the terms of the GNU Lesser General Public License as published 
%   by the Free Software Foundation, version 3 of the License.
%
%   ParaMonte is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with the ParaMonte library. If not, see, 
%
%       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
%
%   ACKNOWLEDGMENT
%
%   As per the ParaMonte library license agreement terms, 
%   if you use any parts of this library for any purposes, 
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [colnames, colindex] = getColNamesIndex( dfColumns, columns )

    allColumnsNeeded = false;
    columnsLen = getVecLen(columns);
    if columnsLen==0
        allColumnsNeeded = true;
    elseif isa(columns,"numeric")
        columns = [columns];
    elseif columnsLen>0
        columns = string(columns);
    else
        error   ( "The input argument 'columns' must be a list whose elements are all" + newline ...
                + "    1.   string-valued, each representing the name of the column from the dataframe, or," + newline ...
                + "    2.   integer-valued, each representing the index of the column from the dataframe." + newline ...
                + "You have entered:" + newline ...
                + string(columns) + newline ...
                );
    end

    dfColumnsLen = length(dfColumns);
    for j = dfColumnsLen:-1:1
        dfColumnsString(j) = string(dfColumns{j});
    end
    if allColumnsNeeded
        colnames = dfColumnsString;
        colindex = 1:1:length(colnames);
    else
        for i = columnsLen:-1:1
            errorOccurred = true;
            if isa(columns(i),"string") || isa(columns(i),"char")
                if isNumericString(columns(i))
                    colindex(i) = str2double(columns(i));
                    colnames(i) = dfColumnsString(colindex(i));
                    errorOccurred = false;
                else % is alphabetical string
                    %colnames(i) = columns(i);
                    colnameLower = columns(i);
                    for j = 1:dfColumnsLen
                        if strcmpi(colnameLower,dfColumnsString(j))
                            errorOccurred = false;
                            colnames(i) = dfColumnsString(j);
                            colindex(i) = j;
                            break;
                        end
                    end
                    if errorOccurred; break; end
                end
            elseif isa(columns(i),"numeric")
                colindex(i) = columns(i);
                colnames(i) = string(dfColumns{colindex(i)});
                errorOccurred = false;
            end
            if errorOccurred; break; end
        end

        if errorOccurred
            dummy = "";
            for j = 1:dfColumnsLen
                dummy = dummy + " """ + string(dfColumns{j}) + """ ";
            end
            error   ( newline ...
                    + "The input element #" + string(i) + " of the ""columns"" variable does not match any elements of the input column names:" + newline ...
                    + "    columns = " + string(columns(i)) + newline ...
                    + "The input column names of the data table: " + newline ...
                    + dummy ...
                    + newline ...
                    );
        end
    
    end

end
