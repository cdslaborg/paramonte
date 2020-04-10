function Quantile = getQuantile(np,nq,SortedQuantileProbability,Point,Weight,sumWeight)

    Quantile    = zeros(nq,1);

    Indx        = zeros(np,1);
    SortedQuantileDensity = zeros(nq,1);

    iq              = 1;
    Quantile        = 0;
    probability     = 0;
    weightCounter   = 0;

    Indx = indexArray(np, Point);

    if ~isempty(sumWeight)
        SortedQuantileDensity = round(SortedQuantileProbability * sumWeight);
        for ip = 1 : np % loopWeighted
            for iw = 1 : Weight(Indx(ip))
                weightCounter = weightCounter + 1;
                if weightCounter >= SortedQuantileDensity(iq)
                    Quantile(iq) = Point(Indx(ip));
                    iq = iq + 1;
                    if iq > nq, break; end
                end
            end
        end % loopWeighted
    else
        SortedQuantileDensity = round(SortedQuantileProbability * np);
        for ip = 1 : np % loopNonWeighted
            if ip >= SortedQuantileDensity(iq)
                Quantile(iq) = Point(Indx(ip));
                iq = iq + 1;
                if iq>nq, break; end
            end
        end % loopNonWeighted
    end

end