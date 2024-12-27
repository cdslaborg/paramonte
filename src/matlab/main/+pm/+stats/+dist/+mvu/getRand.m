%>  \brief
%>  Return a (set of) multivariate Uniform random vector(s),
%>  from within a hyper-ellipsoidal domain.<br>
%>
%>  \details
%>  The returned random vectors are uniformly distributed
%>  within the hyper-ellipsoidal domain of the Uniform distribution.<br>
%>
%>  \param[in]  mean    :   The input vector of MATLAB ``real``,
%>                          representing the mean of a Multivariate Uniform
%>                          distribution in ``size(mean)`` dimensional space.<br>
%>                          (**optional**. default = ``[]``. It must be present if ``cholow`` is missing.)
%>  \param[in]  cholow  :   The input square matrix of MATLAB ``real``,
%>                          representing the lower-triangle of the Cholesky
%>                          factorization of the Gramian matrix of the target
%>                          Multivariate Uniform distribution in ``numel(mean)`` dimensional space.<br>
%>                          This argument can be obtained by passing the Gramian matrix ``gramian``
%>                          of the distribution to the MATLAB intrinsic function ``chol(gramian, "lower")``.<br>
%>                          (**optional**. default = ``[]``. It must be present if ``mean`` is missing.)
%>  \param[in]  s1      :   The input scalar MATLAB whole-number,
%>                          representing the number of random vectors to generate.<br>
%>                          (**optional**. default = ``1``)
%>
%>  \return
%>  ``rand``            :   The output vector of MATLAB ``real`` of
%>                          shape ``(numel(mean), 1)`` containing a random vector
%>                          from the specified Multivariate Uniform distribution.<br>
%>
%>  \interface{getRand}
%>  \code{.m}
%>
%>      rand = pm.stats.dist.mvu.getRand(mean)
%>      rand = pm.stats.dist.mvu.getRand([], cholow)
%>      rand = pm.stats.dist.mvu.getRand(mean, cholow)
%>      rand = pm.stats.dist.mvu.getRand([], cholow, s1)
%>      rand = pm.stats.dist.mvu.getRand(mean, cholow, s1)
%>
%>  \endcode
%>
%>  \example{getRand}
%>
%>  \example{getRand}
%>  \include{lineno} example/stats/dist/mvu/getRand/main.m
%>  \vis{getRand}
%>  \image html example/stats/dist/mvu/getRand/mvu.getRand.hist1.png width=700
%>  \image html example/stats/dist/mvu/getRand/mvu.getRand.hist2.png width=700
%>  \image html example/stats/dist/mvu/getRand/mvu.getRand.scatter1.png width=700
%>
%>  \final{getRand}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 4:03 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function rand = getRand(mean, cholow, s1)
    if  nargin < 3
        s1 = 1;
    end
    if  1 < nargin
        if ~isempty(cholow)
            ndiminv = 1 / size(cholow, 1);
            rand = randn(size(cholow, 1), s1);
            for irand = 1 : s1
                rand(:, irand) = cholow * rand(:, irand) * (unifrnd(0, 1)^ndiminv / norm(rand(:, irand)));
            end
        else
            ndiminv = 1 / numel(mean);
            rand = randn(numel(mean), s1);
            for irand = 1 : s1
                rand(:, irand) = rand(:, irand) * (unifrnd(0, 1)^ndiminv / norm(rand(:, irand)));
            end
        end
    else
        ndiminv = 1 / numel(mean);
        rand = randn(numel(mean), s1);
        for irand = 1 : s1
            rand(:, irand) = rand(:, irand) * (unifrnd(0, 1)^ndiminv / norm(rand(:, irand)));
        end
    end
    if ~isempty(mean)
        rand = bsxfun(@plus, mean(:), rand);
    end
end