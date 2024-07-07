%>  \brief
%>  Generate and return a random positive-definite (correlation or covariance) matrix using the Gram method.<br>
%>
%>  \details
%>  The Gram method of generating random positive-definite square matrices is based on the
%>  observation that every real positive definite matrix \f$M\f$ has a Cholesky factorization
%>  \f{equation}{
%>      M = LL*
%>  \f}
%>  where \f$L\f$ is a uniquely defined lower triangular matrix with positive diagonal entries.<br>
%>  Therefore, \f$M\f$ can be constructed from a given random \f$L\f$.<br>
%>  The **Gram method** is fast, however, the resulting matrix \f$M\f$ does not possess any particular structure.<br>
%>  because it uses the Cholesky factorization of the distribution covariance matrix.<br>
%>
%>  \param[in]  ndim        :   The input positive scalar MATLAB whole number(``integer``),
%>                              representing the rank of the matrix (the number of dimensions) of shape ``(ndim, ndim)``.<br>
%>                              (**optional**. It must be present **if and only if** the input ``scale`` argument is missing or is a scalar.)
%>  \param[in]  scale       :   The input scalar or `contiguous` vector of size `ndim` of type ``real``, representing the scale
%>                              of the matrix (e.g., the standard deviation of a covariance matrix) along each dimension.<br>
%>                              (**optional**. default = ``1``. It can be present **if and only if** it is a scalar or is a vector of size ``ndim``.)
%>
%>  \return
%>  `rand`                  :   The output matrix of shape ``(1:ndim, 1:ndim)`` of type MATLAB ``double``,
%>                              containing a random positive-definite matrix.<br>
%>                              If the optional input scale is missing, then the output ``rand`` is a correlation matrix.<br>
%>
%>  \interface{getRand}
%>  \code{.F90}
%>
%>      rand(1:ndim, 1:ndim) = pm.stats.dist.cov.getRand(ndim)
%>      rand(1:ndim, 1:ndim) = pm.stats.dist.cov.getRand(ndim, scale(1:ndim))
%>
%>  \endcode
%>
%>  \warning
%>  The condition `all([0 < scale])` must hold for the corresponding input arguments.<br>
%>
%>  \example{getRand}
%>  \include{lineno} example/stats/dist/cov/getRand/main.m
%>  \output{getRand}
%>  \include{lineno} example/stats/dist/cov/getRand/main.out.m
%>
%>  \final{getRand}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, July 6 2024, 7:07 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
function rand = getRand(ndim, scale)
    rand = zeros(ndim, ndim);
    for idim = 1 : ndim
        while true
            rand(1 : idim, idim) = unifrnd(-1, 1, idim, 1);
            %%%% \todo
            %%%% The following is redundant for the element and can be optimized.
            normfac = sqrt(dot(rand(1 : idim, idim), rand(1 : idim, idim)));
            if  normfac ~= 0
                break;
            end
        end
        normfac = 1 / normfac;
        rand(1 : idim, idim) = rand(1 : idim, idim) * normfac;
    end
    rand = transpose(rand) * rand;
    %%%%
    %%%%    Set the diagonals to unity to create correlation matrix.
    %%%%    Perhaps not really essential as it is guaranteed to be 1, at least theoretically.
    %%%%
    %rand = rand - diag(diag(rand)) + eye(ndim, ndim);
    %%%%
    %%%% Rescale if requested.
    %%%%
    if  1 < nargin
        if  numel(scale) == 1
            scalesq = scale^2;
            rand = rand * scalesq;
        else
            for jdim = 1 : ndim
                for idim = 1 : jdim - 1
                    rand(idim, jdim) = rand(idim, jdim) * scale(idim) * scale(jdim);
                end
                rand(jdim, jdim) = rand(jdim, jdim) * scale(jdim)^2;
                rand(jdim, idim) = rand(idim, jdim);
            end
        end
    end
end