%>  \brief
%>  Return a (set of) multivariate Normal random vector(s).
%>
%>  \param[in]  mean    :   The input vector of MATLAB ``real``,
%>                          representing the mean of a Multivariate Normal
%>                          distribution in ``size(mean)`` dimensional space.
%>                          (**optional**. default = ``[]``. It must be present if ``cholow`` is missing.)
%>  
%>  \param[in]  cholow  :   The input square matrix of MATLAB ``real``,
%>                          representing the lower-triangle of the Cholesky
%>                          factorization of the covariance matrix of the target 
%>                          Multivariate Normal distribution in ``numel(mean)`` dimensional space.
%>                          This argument can be obtained by passing the covariance matrix ``covmat``
%>                          of the distribution to the MATLAB intrinsic function ``chol(covmat, "lower")``.
%>                          (**optional**. default = ``[]``. It must be present if ``mean`` is missing.)
%>  
%>  \param[in]  s1      :   The input vector of MATLAB ``real``,
%>                          representing the mean of a Multivariate Normal
%>                          distribution in ``size(mean)`` dimensional space.
%>                          (**optional**. default = ``1``)
%>  
%>  \return
%>  `rand`              :   The output vector of MATLAB ``real`` of
%>                          size ``(size(mean), 1)`` containing a random vector
%>                          from the specified Multivariate Normal distribution.
%>
%>  \interface{rand}
%>  \code{.m}
%>
%>      rand = pm.stats.dist.mvn.rand(mean)
%>      rand = pm.stats.dist.mvn.rand([], cholow)
%>      rand = pm.stats.dist.mvn.rand(mean, cholow)
%>      rand = pm.stats.dist.mvn.rand([], cholow, s1)
%>      rand = pm.stats.dist.mvn.rand(mean, cholow, s1)
%>
%>  \endcode
%>
%>  \example{getBorder}
%>
%>      rand = pm.stats.dist.mvn.rand(zeros(2, 1))
%>      histogram(pm.stats.dist.mvn.rand(zeros(2, 1), [], 10000))
%>      histogram(pm.stats.dist.mvn.rand([-3, 3], chol([1 .5; .5, 1], "lower"), 10000))
%>      rand = pm.stats.dist.mvn.rand([-3, 3], chol([1 .9; .9, 1], "lower"), 10000); scatter(rand(1, :), rand(2, :));
%>
%>  \final{rand}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:06 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function val = rand(mean, cholow, s1)
    if  nargin < 3
        s1 = 1;
    end
    if  1 < nargin
        if ~isempty(cholow)
            if 1 < s1
                val = randn(size(cholow, 1), s1);
                for irand = 1 : s1
                    val(:, irand) = cholow * val(:, irand);
                end
            else
                val = cholow * randn(size(cholow, 1), s1);
            end
        else
            val = randn(numel(mean), s1);
        end
    else
        val = randn(numel(mean), s1);
    end
    if ~isempty(mean)
        val = bsxfun(@plus, mean(:), val);
    end
end