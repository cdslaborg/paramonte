%>  \brief
%>  Return a scalar MATLAB string containing the
%>  MATLAB release version, year, or season as requested.
%>
%>  \param[in]  type    :   The input scalar MATLAB string that can be either:
%>                          <ol>
%>                              <li>    the string ``"year"``, implying that the MATLAB release year must be returned.<br>
%>                              <li>    the string ``"season"`` or ``"minor"``, implying that the MATLAB release season must be returned.<br>
%>                          </ol>
%>                          (**optional**. default = ``""`` implying the full version.)
%>
%>  \return
%>  ``matlabRelease``   :   The output scalar MATLAB string containing the
%>                          MATLAB release version, year, or season as requested.
%>
%>  \interface{release}
%>  \code{.m}
%>
%>      matlabRelease = pm.matlab.release()
%>      matlabRelease = pm.matlab.release(type)
%>
%>  \endcode
%>
%>  \example{release}
%>  \include{lineno} example/matlab/release/main.m
%>  \output{release}
%>  \include{lineno} example/matlab/release/main.out.m
%>
%>  \final{release}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:19 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function matlabRelease = release(type)
    if nargin < 1
        type = "";
    end
    % WARNING: Do NOT convert the input char '-release' to string or the output will be wrong.
    matlabRelease = version('-release');
    if string(type) == "year"
        matlabRelease = matlabRelease(1:4);
    elseif string(type) == "season" || string(type) == "minor"
        matlabRelease = matlabRelease(5:5);
    elseif type ~= ""
        disp("type = ");
        disp(type);
        error("Invalid input value for the argument ``type`` of getMatlabRelease().");
    end
    matlabRelease = string(matlabRelease);
end