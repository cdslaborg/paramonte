d = importdata("WeightedData.txt");
sumWeight = sum(d(:,1));
Verbose = zeros(sumWeight,1);
counter = 0;
for i = 1:length(d(:,1))
    for j = 1:d(i,1)
        counter = counter + 1;
        Verbose(counter) = d(i,2);
    end
end
acfCompact = autocorr( d(:,2), length(d(:,2))-1);
acfVerbose = autocorr( Verbose, length(Verbose)-1);

Mean = dot(d(:,1), d(:,2)) / sumWeight;
NormedData = d(:,2) - Mean;
