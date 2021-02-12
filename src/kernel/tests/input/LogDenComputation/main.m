LogRad.min = log(0.02);
LogRad.max = log(0.2);
LogRad.del = (LogRad.max - LogRad.min) / 199;
LogRad.Value = LogRad.min:LogRad.del:LogRad.max;
LogRad.count = length(LogRad.Value);

LogDen.min = log(0.000001);
LogDen.max = log(10000);
LogDen.del = (LogDen.max - LogDen.min) / 99;
LogDen.Value = LogDen.min:LogDen.del:LogDen.max;
LogDen.count = length(LogDen.Value);

figure; hold on;
LogLike = zeros(LogRad.count,LogDen.count);
for id = 1:LogDen.count
    id
    for ir = 1:LogRad.count
        param = [ LogDen.Value(id), LogRad.Value(ir) ];
        LogLike(ir,id) = getNegLogLike(param);
    end
    plot( exp(LogRad.Value) ...
        , LogLike(:,id) ...
        );
end
set(gca, "xscale", "log", "yscale", "linear")
hold off;

figure; hold on;
[X,Y] = meshgrid(LogRad.Value, LogDen.Value);
surf(log(LogLike)');
hold off;
