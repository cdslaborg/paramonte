%>  \brief
%>  Return a black-blue-cyan-white colormap matrix
%>  of shape ``(nell, 3)`` containing a **cold** colormap.<br>
%>
%>  \param[in]  nell    :   The input scalar MATLAB integer representing
%>                          the number of elements (rows) of the output colormap matrix.<br>
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
%>  \see
%>  MATLAB intrinsic functions ``hot()``, ``cool()``, ``jet()``, ``hsv()``, ``gray()``, ``copper()``, ``bone()``, ``vivid()``.<br>
%>
%>  \example{cold}
%>  \include{lineno} example/vis/cmap/cold/main.m
%>  \vis{cold}
%>  \image html example/vis/cmap/cold/cold.1.png width=700
%>  \image html example/vis/cmap/cold/cold.2.png width=700
%>
%>  \final{cold}
%>
%>  Author: Joseph Kirk
%>  Email: jdkirk630@gmail.com
%>  Release: 1.0
%>  Date: 04/21/09
%>
%>  Copyright (c) 2009, Joseph Kirk
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
%>  \JoshuaOsborne, May 21 2024, 7:42 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
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