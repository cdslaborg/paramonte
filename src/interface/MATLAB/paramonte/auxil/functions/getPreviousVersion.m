function previousVerionString = getPreviousVersion(currentVerionString)
    %
    % Take an input version string like, "1.1.1" and return another string representing the version before the input version, like, 1.1.0.
    %
    currentVerionTriplet = getVersionTriplet(currentVerionString);
    previousVerionString = [];
    index = 3;
    while true
        index = index - 1;
        if index <= 0
            return;
        else
            if currentVerionTriplet(index) > 0
                previousVerionTriplet = currentVerionTriplet;
                previousVerionTriplet(index) = previousVerionTriplet(index) - 1;
                previousVerionString = join(string(previousVerionTriplet),".");
                return;
            end
        end
    end
end