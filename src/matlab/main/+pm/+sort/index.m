%>  \brief
%>  Given ``array(1 : n)``, return the array index ``indx(1 : n)``
%>  such that ``array(indx)`` is in ascending order.<br>
%>
%>  \details
%>  This is a simple convenience wrapper around
%>  the MATLAB intrinsic function ``sort()``.<br>
%>
%>  \note
%>  Beware the output ``indx`` is **not** the
%>  rank of the elements of the input array.<br>
%>
%>  \param[in]  array   :   The input MATLAB vector of sortable values
%>                          (that can be passed to the MATLAB intrinsic function ``sort()``).<br>
%>
%>  \return
%>  ``indx``            :   The output MATLAB integer vector of the same size as the
%>                          input ``array`` such that ``array(indx)`` is in ascending order.<br>
%>
%>  \interface{index}
%>  \code{.m}
%>
%>      indx = pm.sort.index(array)
%>
%>  \endcode
%>
%>  \example{index}
%>  \include{lineno} example/sort/index/main.m
%>  \output{index}
%>  \include{lineno} example/sort/index/main.out.m
%>
%>  \final{index}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 3:55 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function indx = index(array)
    [~, indx] = sort(array);
end