function result = genOutputFileName(methodName)
    % generate output filename
    result = methodName + '_run_' + datestr(now, 'ddmmyy') + '_' + datestr(now, 'HHMMSS_FFF');
end
