d = importdata("../output/Test_Knn_mod@test_getLogDensity_1@1.txt");

Diff = d(:,2) - d(:,1);
meanDiff = mean(Diff);
fitLine = @(x) x; % + meanDiff;
x = min(d(:,1)) : 0.1 : max(d(:,1));

figure; hold on; box on;
plot(d(:,1), d(:,2), ".");
plot(x, fitLine(x), "-", "linewidth", 2);
hold off;

figure; hold on; box on;
histfit(Diff)
hold off;
