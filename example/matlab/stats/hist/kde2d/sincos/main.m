cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

%%%%
%%%% Generate a Sinusoidal density.
%%%%

resol = 1000;
meshx = 2 * cos(rand(resol, 1) * 10 * pi) + randn(resol, 1) / 3;%rand(1000, 1); 
meshy = 2 * sin(rand(resol, 1) * 10 * pi) + randn(resol, 1) / 3;
data = [meshx, meshy];

[bandwidth, density, meshx, meshy] = pm.stats.hist.kde2d(data);

% Plot the data and the density estimate.
figure("color", "white");
surf(meshx, meshy, density, 'LineStyle', 'none'); 
%colorbar();
%view([0, 70]); hold on; alpha(.8);
%plot(data(:, 1), data(:, 2), 'w.', 'MarkerSize', 5);
pm.vis.figure.savefig("pm.stats.hist.kde2d.sincos.png", "-m4");