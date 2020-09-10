%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a 
%%%%   copy of this software and associated documentation files (the "Software"), 
%%%%   to deal in the Software without restriction, including without limitation 
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%%%%   and/or sell copies of the Software, and to permit persons to whom the 
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be 
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of 
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
%%%%   work (education/research/industry/development/...) by citing the ParaMonte 
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2009, Joseph Kirk
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmap = cold(m)
%COLD    Black-Blue-Cyan-White Color Map
%   COLD(M) returns an M-by-3 matrix containing a "cold" colormap
%   COLD, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   Example:
%       imagesc(peaks(500))
%       colormap(cold); colorbar
%
%   Example:
%       load topo
%       imagesc(0:360,-90:90,topo), axis xy
%       colormap(cold); colorbar
%
% See also: hot, cool, jet, hsv, gray, copper, bone, vivid
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.0
% Date: 04/21/09

if nargin < 1
    m = size(get(gcf,'colormap'),1);
end
n = fix(3/8*m);

r = [zeros(2*n,1); (1:m-2*n)'/(m-2*n)];
g = [zeros(n,1); (1:n)'/n; ones(m-2*n,1)];
b = [(1:n)'/n; ones(m-n,1)];

cmap = [r g b];



