%>  \brief
%>  Return the corresponding natural logarithm(s) of Probability Density Function (PDF)
%>  of the input ``ndim``-dimensional multivariate Normal random vector(s).<br>
%>
%>  \details
%>  This log-PDF function can be potentially faster and more flexible than the intrinsic MATLAB
%>  equivalent because it optionally accepts vectors of MVN random variates and returns the
%>  natural logarithm of the density function which avoids a costly exponentiation.<br>
%>
%>  \param[in]  x       :   The input vector of shape ``[ndim, nvec]`` of MATLAB ``real``,
%>                          representing the set of ``nvec`` MVN vectors from the ``ndim``-dimensional MVN space.<br>
%>  \param[in]  mean    :   The input vector of shape ``[ndim, 1]`` of MATLAB ``real``,
%>                          representing the mean of a Multivariate Normal
%>                          distribution in ``size(mean)`` dimensional space.<br>
%>                          If the input value is empty, it is set to the default value.<br>
%>                          (**optional**. default = ``zeros(ndim, 1)``.)
%>  \param[in]  invcov  :   The input positive-definite square matrix of shape ``[ndim, ndim]`` of MATLAB ``real``,
%>                          representing the inverse of the covariance matrix of the target Multivariate
%>                          Normal distribution in ``numel(mean)`` dimensional space.<br>
%>                          If the input value is empty, it is set to the default value.<br>
%>                          Note that the inverse covariance matrix can be readily obtained from
%>                          its lower or upper Cholesky factorization (``cholow`` or ``choupp``) via
%>                          [pm.matlab.inv(cholow)](@ref inv) or [pm.matlab.inv(choupp')](@ref inv).<br>
%>                          (**optional**. default = ``eye(ndim, ndim)``.)
%>
%>  \return
%>  ``logPDF``          :   The output column vector of MATLAB ``real`` of shape ``[nvec, 1]``
%>                          containing the natural logarithm(s) of the input MVN vector(s).<br>
%>
%>  \interface{getLogPDF}
%>  \code{.m}
%>
%>      logPDF = pm.stats.dist.mvn.getLogPDF(x)
%>      logPDF = pm.stats.dist.mvn.getLogPDF(x, mean)
%>      logPDF = pm.stats.dist.mvn.getLogPDF(x, [], invcov)
%>      logPDF = pm.stats.dist.mvn.getLogPDF(x, mean, invcov)
%>
%>  \endcode
%>
%>  \example{getLogPDF}
%>  \include{lineno} example/stats/dist/mvn/getLogPDF/main.m
%>  \output{getLogPDF}
%>  \include{lineno} example/stats/dist/mvn/getLogPDF/main.m
%>  \vis{getLogPDF}
%>  \image html example/stats/dist/mvn/getLogPDF/mvn.getLogPDF.scatter.3d.png width=700
%>
%>  \final{getLogPDF}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:06 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function logPDF = getLogPDF(x, mean, invcov)
    ndim = size(x, 1);
    nvec = size(x, 2);
    if  nargin < 3
        invcov = [];
    end
    if  nargin < 2
        mean = [];
    end
    if  isempty(invcov)
        logDetInvCov = 0.;
        invcov = eye(ndim, ndim);
    else
        logDetInvCov = log(det(invcov));
    end
    if  isempty(mean)
        mean = zeros(ndim, 1);
    else
        mean = mean(:);
    end
    logPDF = zeros(nvec, 1);
    log2pi = 1.837877066409345; % log(2 * pi)
    logNormFac = -0.5 * (ndim * log2pi - logDetInvCov);
    for ivec = 1 : nvec
        xnormed = x(:, ivec) - mean;
        logPDF(ivec) = logNormFac - .5 * (xnormed' * invcov * xnormed);
    end
end