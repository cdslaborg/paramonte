function versionTriplet = getVersionTriplet(versionDumpString)
    %
    % Take an input version string like, "1.1.1" and return an integer triplet list.
    %
    versionTriplet = str2double(string(strsplit(versionDumpString,".")));
end
