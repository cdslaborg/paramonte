cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

%%%%
%%%% Generate a Gaussian mixture with distant modes.
%%%%

data =  [ randn(100, 1), randn(100, 1) / 4 ...
        ; randn(100, 1) + 18, randn(100,1) ...
        ; randn(100, 1) + 15, randn(100, 1)/ 2 - 18 ...
        ];

[bandwidth, density, meshx, meshy] = pm.fit.kde2d(data);

% Plot the data and the density estimate.
figure("color", "white");
contour3(meshx, meshy, density, 50);
%colormap("hot"); hold on; alpha(.8); view([0, 60]);
%plot(data(:, 1), data(:, 2), 'w.', 'MarkerSize', 5);
pm.vis.figure.savefig("pm.fit.kde2d.mvnmixfar.png", "-m4");