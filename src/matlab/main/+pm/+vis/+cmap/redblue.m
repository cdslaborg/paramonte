%>  \brief
%>  Return a blue-red colormap matrix
%>  of shape ``(nell, 3)`` containing a **redblue** colormap.<br>
%>
%>  \details
%>  The colors begin with bright blue, range through shades of
%>  blue to white, and then through shades of red to bright red.<br>
%>  The output is the same length as the current figure colormap.<br>
%>  If no figure exists, MATLAB will create one.<br>
%>
%>  \param[in]  nell    :   The input scalar MATLAB integer representing
%>                          the number of elements (rows) of the output colormap matrix.<br>
%>                          (**optional**, default = ``size(get(gcf, 'colormap'), 1)``)
%>
%>  \return
%>  ``cmap``            :   The output blue-red colormap matrix
%>                          of shape ``[nell, 3]`` containing a **redblue** colormap.<br>
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
%>  \include{lineno} example/vis/cmap/redblue/main.m
%>  \vis{redblue}
%>  \image html example/vis/cmap/redblue/redblue.1.png width=700
%>  \image html example/vis/cmap/redblue/redblue.2.png width=700
%>
%>  \see
%>  MATLAB intrinsic functions ``hsv()``, ``gray()``, ``hot()``, ``bone()``, ``copper()``, ``pink()``, ``flag()``, ``colormap()``, ``rgbplot()``.<br>
%>
%>  \final{redblue}
%>
%>  Adam Auton, 9th October 2009
%>  Copyright (c) 2009, Adam Auton
%>  All rights reserved.
%>
%>  Redistribution and use in source and binary forms, with or without
%>  modification, are permitted provided that the following conditions are
%>  met:
%>
%>  * Redistributions of source code must retain the above copyright
%>    notice, this list of conditions and the following disclaimer.
%>  * Redistributions in binary form must reproduce the above copyright
%>    notice, this list of conditions and the following disclaimer in
%>    the documentation and/or other materials provided with the distribution
%>
%>  This software is provided by the copyright holders and contributors "as is"
%>  and any express or implied warranties, including, but not limited to, the
%>  implied warranties of merchantability and fitness for a particular purpose
%>  are disclaimed. in no event shall the copyright owner or contributors be
%>  liable for any direct, indirect, incidental, special, exemplary, or
%>  consequential damages (including, but not limited to, procurement of
%>  substitute goods or services; loss of use, data, or profits; or business
%>  interruption) however caused and on any theory of liability, whether in
%>  contract, strict liability, or tort (including negligence or otherwise)
%>  arising in any way out of the use of this software, even if advised of the
%>  possibility of such damage.
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:44 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
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