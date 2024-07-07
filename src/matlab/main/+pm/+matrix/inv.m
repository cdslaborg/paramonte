%>  \brief
%>  Generate and return the inverse of a positive-definite matrix
%>  whose Lower-triangular Cholesky factorization is given as input.<br>
%>
%>  \param[in]  cholow  :   The input square matrix representing the lower-triangular
%>                          Cholesky factorization whose inverse must be computed.<br>
%>
%>  \return
%>  ``invmat``          :   The output square positive-definite matrix, representing the inverse of
%>                          the matrix whose lower-triangular Cholesky factorization is given as input.<br>
%>
%>  \interface{inv}
%>  \code{.m}
%>
%>      invmat = pm.matrix.inv(cholow);
%>
%>  \endcode
%>
%>  \see
%>  See *Matrix Inversion Using Cholesky Decomposition*, Aravindh Krishnamoorthy, Deepak Menon, (arXiv:1111.4144).<br>
%>
%>  \example{inv}
%>  \include{lineno} example/matrix/inv/main.m
%>  \output{inv}
%>  \include{lineno} example/matrix/inv/main.out.m
%>
%>  \final{inv}
%>
%>  This function builds upon the work of Aravindh Krishnamoorthy (2024).
%>  [Matrix Inversion using Cholesky Decomposition](https://www.mathworks.com/matlabcentral/fileexchange/41957-matrix-inversion-using-cholow-decomposition),
%>  MATLAB Central File Exchange.<br>
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:06 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function invmat = inv(cholow)
    %%%%
    %%%% Initializations.
    %%%%
    ndim = size(cholow, 1);
    invmat = zeros(ndim, ndim);
    %%%%
    %%%% Construct the auxillary diagonal matrix ``auxdia = 1 / rii``.
    %%%%
    auxdia = inv(diag(diag(cholow)));
    %%%%
    %%%% Compute the inverse.
    %%%%
    for j = ndim : -1 : 1
        for i = j : -1 : 1
            invmat(i, j) = auxdia(i, j) - cholow(i + 1 : end, i)' * invmat(i + 1 : end, j);
            invmat(i, j) = invmat(i, j) / cholow(i, i);
            %%%%
            %%%% Write out the symmetric element.
            %%%%
            invmat(j, i) = conj(invmat(i, j));
        end
    end
end