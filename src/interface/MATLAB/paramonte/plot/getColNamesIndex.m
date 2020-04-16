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
            for i = 1:dfColumnsLen
                dummy = dummy + " """ + string(dfColumns{i}) + """ ";
            end
            error   ( "The " + string(i) + "th input element of the ""columns"" variable is neither numeric, nor char or string:" + newline ...
                        + "    columns = " + string(columns(i)) + newline ...
                        + "The input column names of the data table: " + newline ...
                        + dummy + newline ...
                        );
        end
    
    end

end
