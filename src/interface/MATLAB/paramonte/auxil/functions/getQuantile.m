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