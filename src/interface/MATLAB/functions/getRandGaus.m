function getRandGaus = getRandGaus()
    % This code returns a normally distributed deviate with zero mean and unit variance.

    persistent iset
    persistent gset

    if isempty(iset), iset = 0; end

    if iset == 0
        while true
            vec = rand(1, 2);
            vec = 2*vec - 1;
            rsq = vec(1,1).^2 + vec(1,2).^2;
            if rsq > 0 && rsq < 1, break; end
        end
        fac = sqrt(-2 * log(rsq)/rsq);
        gset = vec(1)*fac;
        getRandGaus = vec(2)*fac;
        iset = 1;
    else
        getRandGaus = gset;
        iset = 0;
    end

end