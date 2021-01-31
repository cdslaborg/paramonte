clear all;
close all;

d = importdata("../output/Test_MultiSkewNorm_mod@test_getRandMSN_1@1.txt");

figure;

plot( d(:,1) ...
    , d(:,2) ...
    , '.' ...
    )
