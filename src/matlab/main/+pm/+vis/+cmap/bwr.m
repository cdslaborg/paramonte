%>  \brief
%>  Return a blue-white-red colormap matrix
%>  of shape ``(nell, 3)`` containing the RGB values with
%>  white corresponding to the CAXIS value closest to zero.<br>
%>
%>  \note
%>  This colormap is most useful for images and surface plots
%>  with positive and negative values.<br>
%>
%>  \param[in]  nell    :   The input scalar MATLAB integer representing
%>                          the number of elements (rows) of the output colormap matrix.<br>
%>                          (**optional**, default = ``size(get(gcf, 'colormap'), 1)``)
%>  \param[in]  clim    :   The input vector length ``2`` of MATLAB doubles
%>                          representing the color limits of the output colormap.<br>
%>                          (**optional**, default = ``get(gca, 'CLim')``)
%>
%>  \return
%>  ``cmap``            :   The output blue-white-red colormap matrix
%>                          of shape ``[nell, 3]`` containing a **cold** colormap.<br>
%>
%>  \interface{bwr}
%>  \code{.m}
%>
%>      cmap = pm.vis.cmap.cold()
%>      cmap = pm.vis.cmap.cold([])
%>      cmap = pm.vis.cmap.cold(nell)
%>
%>  \endcode
%>
%>  \see
%>  MATLAB intrinsic functions ``hsv()``, ``hot()``, ``cool()``, ``bone()``, ``copper()``, ``pink()``, ``flag()``, ``colormap()``, ``rgbplot()``.<br>
%>
%>  \example{bwr}
%>  \include{lineno} example/vis/cmap/bwr/main.m
%>  \vis{bwr}
%>  \image html example/vis/cmap/bwr/bwr.1.png width=700
%>  \image html example/vis/cmap/bwr/bwr.2.png width=700
%>  \image html example/vis/cmap/bwr/bwr.3.png width=700
%>  \image html example/vis/cmap/bwr/bwr.4.png width=700
%>
%>  \final{bwr}
%>
%>  Copyright (c) 2009, Nathan Childress
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
%>  \JoshuaOsborne, May 21 2024, 7:40 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function newmap = bwr(nell, clim)

    if nargin < 1 || isempty(nell)
       nell = size(get(gcf,'colormap'), 1);
    end
    bottom = [0, 0, 0.5];
    botmiddle = [0, 0.5, 1];
    topmiddle = [1, 0, 0];
    middle = [1, 1, 1];
    top = [0.5, 0, 0];

    % Find middle.

    try
        lims = get(gca, 'CLim');
    catch
        lims = clim;
    end

    % Find ratio of negative to positive.

    if (lims(1) < 0) & (lims(2) > 0)

        % It has both negative and positive
        % Find ratio of negative to positive
        ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
        neglen = round(nell*ratio);
        poslen = nell - neglen;

        % Just negative

        new = [bottom; botmiddle; middle];
        len = length(new);
        oldsteps = linspace(0, 1, len);
        newsteps = linspace(0, 1, neglen);
        newmap1 = zeros(neglen, 3);

        for i = 1 : 3
            % Interpolate over RGB spaces of colormap
            newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
        end

        % Just positive.

        new = [middle; topmiddle; top];
        len = length(new);
        oldsteps = linspace(0, 1, len);
        newsteps = linspace(0, 1, poslen);
        newmap = zeros(poslen, 3);

        for i = 1 : 3
            % Interpolate over RGB spaces of colormap
            newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
        end

        % And put them together.

        newmap = [newmap1; newmap];

    elseif lims(1) >= 0

        % Just positive.

        new = [middle; topmiddle; top];
        len = length(new);
        oldsteps = linspace(0, 1, len);
        newsteps = linspace(0, 1, nell);
        newmap = zeros(nell, 3);

        for i = 1 : 3
            % Interpolate over RGB spaces of colormap
            newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
        end

    else

        % Just negative.

        new = [bottom; botmiddle; middle];
        len = length(new);
        oldsteps = linspace(0, 1, len);
        newsteps = linspace(0, 1, nell);
        newmap = zeros(nell, 3);

        for i=1:3
            % Interpolate over RGB spaces of colormap
            newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
        end

    end
    %
    % nell = 64;
    % new = [bottom; botmiddle; middle; topmiddle; top];
    % % x = 1:nell;
    %
    % oldsteps = linspace(0, 1, 5);
    % newsteps = linspace(0, 1, nell);
    % newmap = zeros(nell, 3);
    %
    % for i=1:3
    %     % Interpolate over RGB spaces of colormap
    %     newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    % end
    %
    % % set(gcf, 'colormap', newmap), colorbar

end