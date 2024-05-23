%>  \brief
%>  Return a blue-red colormap matrix
%>  of shape ``(nell, 3)`` containing a **redblue** colormap.
%>
%>  \details
%>  The colors begin with bright blue, range through shades of
%>  blue to white, and then through shades of red to bright red.
%>  The output is the same length as the current figure colormap.
%>  If no figure exists, MATLAB will create one.
%>
%>  \param[in]  nell    :   The input scalar MATLAB integer representing
%>                          the number of rows of the output colormap matrix.
%>                          (**optional**, default = ``size(get(gcf, 'colormap'), 1)``)
%>
%>  \return
%>  `cmap`              :   The output blue-red colormap matrix
%>                          of shape ``[nell, 3]`` containing a **redblue** colormap.
%>
%>  \interface{redblue}
%>  \code{.m}
%>
%>      cmap = pm.vis.cmap.redblue()
%>      cmap = pm.vis.cmap.redblue([])
%>      cmap = pm.vis.cmap.redblue(nell)
%>
%>  \endcode
%>
%>  \example{redblue}
%>
%>      imagesc(peaks(500));
%>      colormap(pm.vis.cmap.redblue());
%>      colorbar;
%>
%>      load topo
%>      imagesc(0:360, -90:90, topo); axis("xy");
%>      colormap(pm.vis.cmap.redblue());
%>      colorbar;
%>
%>  See also: hsv, gray, hot, bone, copper, pink, flag, colormap, rgbplot.
%>
%>  \final{redblue}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:44 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
%>      Originally by:
%>
%>          Adam Auton, 9th October 2009
%>          Copyright (c) 2009, Adam Auton
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
function cmap = redblue(nell)
    if nargin < 1
        nell = size(get(gcf, 'colormap'), 1);
    end
    if (mod(nell, 2) == 0)
        % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
        nellHalf = nell * 0.5;
        r = (0 : nellHalf - 1)' / max(nellHalf - 1, 1);
        g = r;
        r = [r; ones(nellHalf, 1)];
        g = [g; flipud(g)];
        b = flipud(r);
    else
        % From [0 0 1] to [1 1 1] to [1 0 0];
        nellHalf = floor(nell * 0.5);
        r = (0 : nellHalf - 1)' / max(nellHalf, 1);
        g = r;
        r = [r; ones(nellHalf + 1, 1)];
        g = [g; 1; flipud(g)];
        b = flipud(r);
    end
    cmap = [r, g, b];
end