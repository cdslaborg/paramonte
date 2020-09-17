function maxLogFunc = getMaxLogFunc(dataFrame, column)
    %
    %   Returns a structure with components containing the properties of the 
    %   input ``dataFrame``, corresponding to the mode (maximum) of the data 
    %   in the column of the dataFrame that is identified by ``column``. 
    %
    %       **Parameters**
    %
    %           dataFrame
    %
    %               A MATLAB Table containing the output sample data 
    %               from a ParaMonte simulation.
    %
    %           column
    %
    %               A string containing the name of the column of the input 
    %               ``dataFrame`` that contains values of the objective 
    %               function (or its logarithm). The default value is
    %               "SampleLogFunc".
    %
    %       **Returns**
    %
    %           maxLogFunc
    %
    %               A structure with the following components:
    %
    %                   idrow
    %
    %                       The ID of the row of ``dataFrame`` corresponding 
    %                       to the mode of ``column``, 
    %
    %                   value
    %
    %                       The value of ``column`` at maximum, 
    %
    %                   dfrow
    %
    %                       The entire row of ``dataFrame`` corresponding to 
    %                       the mode of ``column``, 
    %
    %                   state
    %
    %                       The location (state) within the domain of the objective 
    %                       function where the maximum of ``column`` occurs.
    %
    if nargin<2; column = "SampleLogFunc"; end
    for i = 1:length(dataFrame.Properties.VariableNames)
        if strcmpi(dataFrame.Properties.VariableNames{i},column)
            column = dataFrame.Properties.VariableNames{i};
            break;
        end
    end
    offset = i + 1;
    maxLogFunc = struct();
    [max_val,max_idx] = max(dataFrame.(column));
    maxLogFunc.idrow = max_idx;
    maxLogFunc.value = max_val;
    maxLogFunc.dfrow = dataFrame{maxLogFunc.idrow,:};
    maxLogFunc.state = dataFrame{maxLogFunc.idrow,offset:end};
end