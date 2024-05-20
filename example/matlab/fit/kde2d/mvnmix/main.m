cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

%%%%
%%%% Generate a Gaussian mixture with distant modes.
%%%%

data = [ randn(500, 2) ...
       ; randn(500, 1) + 3.5, randn(500, 1) ...
       ];

[bandwidth, density, meshx, meshy] = pm.fit.kde2d(data);

% Plot the data and the density estimate.
contour3(meshx, meshy, density, 50); hold on;
plot(data(:, 1), data(:, 2), 'r.', 'MarkerSize', 5);
pm.vis.figure.savefig("kde2d.png", "-m4");