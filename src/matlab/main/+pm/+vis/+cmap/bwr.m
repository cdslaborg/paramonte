function newmap = bwr(nell, clim)
    %
    %   Return a blue-white-red colormap matrix
    %   of shape ``(nell, 3)`` containing the RGB values with
    %   white corresponding to the CAXIS value closest to zero.
    %
    %   \note
    %
    %       This colormap is most useful for images and surface plots
    %       with positive and negative values.
    %
    %   Parameters
    %   ----------
    %
    %       nell
    %
    %           The input scalar MATLAB integer representing
    %           the number of rows of the output colormap matrix.
    %           (**optional**, default = ``size(get(gcf, 'colormap'), 1)``)
    %
    %       clim
    %
    %           The input vector length ``2`` of MATLAB doubles
    %           representing the color limits of the output colormap.
    %           (**optional**, default = ``get(gca, 'CLim')``)
    %
    %   Returns
    %   -------
    %
    %       cmap
    %
    %           The output blue-white-red colormap matrix
    %           of shape ``[nell, 3]`` containing a **cold** colormap.
    %
    %   Interface
    %   ---------
    %
    %       cmap = pm.vis.cmap.cold()
    %       cmap = pm.vis.cmap.cold([])
    %       cmap = pm.vis.cmap.cold(nell)
    %
    %   Example
    %   -------
    %
    %       figure;
    %       imagesc(peaks(500));
    %       colormap(pm.vis.cmap.bwr(256));
    %       colorbar;
    %
    %       figure;
    %       imagesc(peaks(500), [0, 8]);
    %       colormap(pm.vis.cmap.bwr());
    %       colorbar;
    %
    %       figure;
    %       imagesc(peaks(250), [-6, 0]);
    %       colormap(pm.vis.cmap.bwr());
    %       colorbar;
    %
    %       figure;
    %       surf(peaks);
    %       colormap(pm.vis.cmap.bwr());
    %       axis("tight");
    %
    %   See also hsv, hot, cool, bone, copper, pink, flag, colormap, rgbplot.
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    %       Originally by:
    %
    %           Copyright (c) 2009, Nathan Childress
    %           All rights reserved.
    %
    %           Redistribution and use in source and binary forms, with or without
    %           modification, are permitted provided that the following conditions are
    %           met:
    %
    %               * Redistributions of source code must retain the above copyright
    %                 notice, this list of conditions and the following disclaimer.
    %               * Redistributions in binary form must reproduce the above copyright
    %                 notice, this list of conditions and the following disclaimer in
    %                 the documentation and/or other materials provided with the distribution
    %
    %           THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    %           AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    %           IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    %           ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    %           LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    %           CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    %           SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    %           INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    %           CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    %           ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    %           POSSIBILITY OF SUCH DAMAGE.
    %
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