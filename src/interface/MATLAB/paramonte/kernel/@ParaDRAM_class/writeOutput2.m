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
function writeOutput2(self)
    
    % if new point has been sampled, write the previous sampled point to output file

    data = [""];

    if self.SpecBase.chainFileFormat.isCompact

        data =  [ self.Chain.ProcessID      ...
                ; self.Chain.DelRejStage    ...
                ; self.Chain.MeanAccRate    ...
                ; self.Chain.Adaptation     ...
                ; self.Chain.BurninLoc      ...
                ; self.Chain.Weight         ...
                ; self.Chain.LogFunc        ...
                ; self.Chain.State          ...
                ] ;
        fprintf(self.ChainFile.unit, self.ChainFile.format, data);

    elseif self.SpecBase.chainFileFormat.isVerbose

        count = ones(1, length(self.Chain.Weight));
        data = repelem( [ self.Chain.ProcessID      ...
                        ; self.Chain.DelRejStage    ...
                        ; self.Chain.MeanAccRate    ...
                        ; self.Chain.Adaptation     ...
                        ; self.Chain.BurninLoc      ...
                        ; count                     ...
                        ; self.Chain.LogFunc        ...
                        ; self.Chain.State        ] ...
                        , 1, self.Chain.Weight      ...
                        ) ;
        fprintf(self.ChainFile.unit, self.ChainFile.format, data);

    end

    fprintf(self.ChainFile.unit, "\n\t\twriteOutput2");

end