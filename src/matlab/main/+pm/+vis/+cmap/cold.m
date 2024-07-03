%>  \brief
%>  Return a black-blue-cyan-white colormap matrix
%>  of shape ``(nell, 3)`` containing a **cold** colormap.<br>
%>
%>  \param[in]  nell    :   The input scalar MATLAB integer representing
%>                          the number of rows of the output colormap matrix.<br>
%>                          (**optional**, default = ``size(get(gcf, 'colormap'), 1)``)
%>
%>  \return
%>  ``cmap``            :   The output black-blue-cyan-white colormap matrix
%>                          of shape ``[nell, 3]`` containing a **cold** colormap.<br>
%>
%>  \interface{cold}
%>  \code{.m}
%>
%>      cmap = pm.vis.cmap.cold()
%>      cmap = pm.vis.cmap.cold([])
%>      cmap = pm.vis.cmap.cold(nell)
%>
%>  \endcode
%>
%>  \example{cold}
%>
%>      imagesc(peaks(500));
%>      colormap(pm.vis.cmap.cold());
%>      colorbar;
%>
%>      load topo
%>      imagesc(0:360, -90:90, topo); axis("xy");
%>      colormap(pm.vis.cmap.cold());
%>      colorbar;
%>
%>  See also: hot, cool, jet, hsv, gray, copper, bone, vivid
%>
%>  \final{cold}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:42 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
%>      Originally by:
%>
%>          Author: Joseph Kirk
%>          Email: jdkirk630@gmail.com
%>          Release: 1.0
%>          Date: 04/21/09
%>  
%>          Copyright (c) 2009, Joseph Kirk
%>          All rights reserved.
%>  
%>          Redistribution and use in source and binary forms, with or without
%>          modification, are permitted provided that the following conditions are
%>          met:
%>  
%>              * Redistributions of source code must retain the above copyright
%>                notice, this list of conditions and the following disclaimer.
%>              * Redistributions in binary form must reproduce the above copyright
%>                notice, this list of conditions and the following disclaimer in
%>                the documentation and/or other materials provided with the distribution
%>  
%>          THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%>          AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%>          IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%>          ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%>          LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%>          CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%>          SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%>          INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%>          CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%>          ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%>          POSSIBILITY OF SUCH DAMAGE.
%>
function cmap = cold(nell)
    if nargin < 1
        nell = [];
    end
    if  isempty(nell)
        nell = size(get(gcf, 'colormap'), 1);
    end
    n = fix(3 / 8 * nell);
    r = [zeros(2 * n, 1); (1 : nell - 2 * n)' / (nell - 2 * n)];
    g = [zeros(n, 1); (1 : n)' / n; ones(nell - 2 * n, 1)];
    b = [(1 : n)' / n; ones(nell - n, 1)];
    cmap = [r, g, b];
end