cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

range = -6 : 0.01 : 6;
[x, y] = meshgrid(range, range);
z = pm.math.func.himmelblau.getFunc(x, y);
for dim = 2 : 3
    figure("color", "white");
    if  dim == 2
        contourf(x, y, z, 50);
    else
        surf(x, y, z, "EdgeColor", "none");
        zlabel("Z");
    end
    xlabel("X");
    ylabel("Y");
    cbar = colorbar();
    ylabel(cbar, "Function Value")
    pm.vis.figure.savefig("himmelblau.getFunc." + string(dim) + "d.png", "-m3");
end