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
classdef logfunc

    properties (Constant)
        NDIM = 2; % The number of dimensions of the domain of the objective function.
    end

    % define the natural logarithm of the objective function: getLogFunc

    methods (Static)
        function logFuncVal = get(point)
            % Return the negative natural logarithm of the 2-dimensional Himmelblau's function.
            % Reference: https://en.wikipedia.org/wiki/Himmelblau%27s_function
            %
            % note that the input point is an array of NDIM rows and 1 column.
            % Therefore, the MEAN array must have the same dimensions.

            % the logarithm of objective function: log(ObjectiveFunction)

            logFuncVal = -log( (point(1)^2 + point(2) - 11)^2 + (point(1) + point(2)^2 - 7)^2 );

        end
    end

end
