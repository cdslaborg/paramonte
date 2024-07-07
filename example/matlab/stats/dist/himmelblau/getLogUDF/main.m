cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

range = -6 : 0.01 : 6;
[x, y] = meshgrid(range, range);
z = exp(pm.stats.dist.himmelblau.getLogUDF(x, y, 1));
for dim = 2 : 3
    figure("color", "white");
    if  dim == 2
        contourf(x, y, log(z), 50, "EdgeColor", "none");
    else
        surf(x, y, z, log(z), "EdgeColor", "none");
        set(gca, "zscale", "log");
        zlabel("Z");
    end
    xlabel("X");
    ylabel("Y");
    cbar = colorbar();
    ylabel(cbar, "Function Value")
    set(gca, "fontsize", 13);
    pm.vis.figure.savefig("himmelblau.getLogUDF." + string(dim) + "d.png", "-m3");
end