function result = genOutputFileName(self)
    % generate output filename
    result = self.methodName + "_run_" + datestr(now, "yyyymmdd") + "_" + datestr(now, "HHMMSS_FFF");
end
