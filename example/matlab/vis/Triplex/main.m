close all
clear all
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

sampler = pm.sampling.Paradram();
sampler.spec.outputRestartFileFormat = "ascii";
sampler.spec.outputChainSize = 10000;
sampler.spec.outputFileName = "./test";
sampler.spec.outputStatus = "retry";
sampler.run(@(x) -sum(x.^2), 4);
restart = sampler.readRestart();
sample = sampler.readSample();

%%%%

cols = sample{1}.slfc + [1 : sample{1}.ndim];
kws = {"colorbar", {"enabled", false}};
tx = pm.vis.Triplex ( pm.vis.SubplotScatter(sample{1}.df, "colx", cols, "coly", cols, kws{:}) ...
                    , pm.vis.SubplotHistogram(sample{1}.df, "colx", cols) ...
                    , pm.vis.SubplotContour(sample{1}.df, "colx", cols, "coly", cols, kws{:}) ...
                    );
%tx.layout.reset();
%tx.layout.cbarh.enabled = true;
%tx.layout.cbarv.enabled = true;
%tx.layout.cbarh.position = [];
%tx.layout.cbarv.position = [];
%tx.layout.tiling.tile.width = [];
tx.make();
pm.vis.figure.savefig("Triplex.1.png", "-m3");

%%%%
%%%%

url = "https://raw.githubusercontent.com/cdslaborg/paramontex/fbeca6745684c798ff28c1bf57cfae0c190db478/MATLAB/mlx/sampling_multivariate_normal_distribution_via_paradram/out/mvn_serial_process_1_chain.txt";
chain = sampler.readChain(url);
cols = chain{1}.slfc + [1 : chain{1}.ndim];

%%%%

tx = pm.vis.Triplex ( pm.vis.SubplotLineScatter(chain{1}.df, "colx", cols, "coly", cols, "colc", chain{1}.slfc, "colorbar", {"enabled", false}, "colormap", {"map", "autumn"}) ...
                    , pm.vis.SubplotHistogram(chain{1}.df, "colx", cols, "histogram", {"faceColor", "red"}) ...
                    , pm.vis.SubplotContour3(chain{1}.df, "colx", cols, "coly", cols, "colorbar", {"enabled", false}) ...
                    );
tx.make();
pm.vis.figure.savefig("Triplex.2.png", "-m3");

tx.fout.cbarv.Label.Interpreter = "tex";
tx.fout.cbarv.Label.FontSize = 13;
tx.fout.cbarv.Label.String = "Log_e ( Probbability Density Function of the MVN distribution )";
pm.vis.figure.savefig("Triplex.3.png", "-m3");

%%%%

tx = pm.vis.Triplex ( pm.vis.SubplotLineScatter(chain{1}.df, "colx", cols, "coly", cols, "colc", chain{1}.slfc, "colorbar", {"enabled", false}) ...
                    , pm.vis.SubplotHistfit(chain{1}.df, "colx", cols) ...
                    , [] ...pm.vis.SubplotContour3(chain{1}.df, "colx", cols, "coly", cols, "colorbar", {"enabled", false}) ...
                    );
tx.layout.tiling.position = [.1, nan, nan, nan]; % allow more room for axes labels.
tx.make();
pm.vis.figure.savefig("Triplex.4.png", "-m3");

%%%%

tx = pm.vis.Triplex ( pm.vis.SubplotLineScatter(chain{1}.df, "colx", cols, "coly", cols, "colormap", {"enabled", false}) ...
                    , pm.vis.SubplotHistfit(chain{1}.df, "colx", cols) ...
                    , pm.vis.SubplotContour3(chain{1}.df, "colx", cols, "coly", cols, "colormap", {"enabled", false}) ...
                    );
tx.layout.tiling.position = [.1, nan, nan, nan]; % allow more room for axes labels.
tx.make();
pm.vis.figure.savefig("Triplex.5.png", "-m3");

%%%%

tx = pm.vis.Triplex ( [] ...
                    , pm.vis.SubplotHistfit(chain{1}.df, "colx", cols) ...
                    , pm.vis.SubplotScatter(chain{1}.df, "colx", cols, "coly", cols, "colorbar", {"enabled", false}) ...
                    );
tx.make();
pm.vis.figure.savefig("Triplex.6.png", "-m3");

%%%%

tx = pm.vis.Triplex ( pm.vis.SubplotLineScatter(chain{1}.df, "colx", cols, "coly", cols, "colormap", {"enabled", false}) ...
                    , [] ...
                    , pm.vis.SubplotContourf(chain{1}.df, "colx", cols, "coly", cols, "colorbar", {"enabled", false}, "colormap", {"enabled", true}) ...
                    );
tx.make();
pm.vis.figure.savefig("Triplex.7.png", "-m3");