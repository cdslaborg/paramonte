close all
clear all
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

nrow = 5;
ncol = 5;
tl = pm.vis.TilingLayout(nrow, ncol, "cbarh", {"enabled", true}, "cbarv", {"enabled", true});
fig = pm.vis.figure.Figure("figure", {"position", tl.position});
fig.make();
ax = struct();
ax.main = axes("position", tl.tiling.position, "visible", "off");
ax.tile = cell(nrow, ncol);
for icol = 1 : ncol
    for irow = 1 : nrow
        ax.tile{irow, icol} = axes("position", tl.tiling.tile.position(:, irow, icol), "visible", "on", "box", "on", "FontSize", 10);
        if  1 < icol
            ax.tile{irow, icol}.YTickLabel = [];
        end
        if  irow < nrow
            ax.tile{irow, icol}.XTickLabel = [];
        end
    end
end
if  tl.cbarv.enabled
    %cv = colorbar(ax.main, "location", "north");%, "position", tl.cbarv.position);
    cv = colorbar(ax.main, "position", tl.cbarv.position);
    xlabel(cv, "Vertical Label", "fontsize", 13);
end
if  tl.cbarv.enabled
    ch = colorbar(ax.main, "position", tl.cbarh.position, "location", "north");
    xlabel(ch, "Horizontal Label", "fontsize", 13);
end
xlabel(ax.main, "X-Axis Label", "fontsize", 13, "visible", "on");
ylabel(ax.main, "Y-Axis Label", "fontsize", 13, "visible", "on");
pm.vis.figure.savefig("TilingLayout.1.png", "-m3");