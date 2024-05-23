%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if the
%>  user response to a question on the MATLAB command prompt
%>  is ``Y`` or ``y``, representing YES, otherwise ``false``
%>  if the user response is ``N`` or ``n``.
%>
%>  \details
%>  Continues asking for as long as the user response
%>  is not among the accepted answers mentioned above.
%>
%>  \param[in]  question    :   The input scalar MATLAB string containing the
%>                              question to be displayed on the MATLAB prompt.
%>
%>  \return
%>  `itis`                  :   The output scalar MATLAB logical that is ``true`` if the
%>                              user response to a question on the MATLAB command prompt
%>                              is ``Y`` or ``y``, representing YES, otherwise `false`
%>                              if the user response is ``N`` or ``n``.
%>
%>  \interface{isyes}
%>  \code{.m}
%>
%>      pm.io.isyes()
%>      pm.io.isyes(question)
%>
%>  \endcode
%>  \final{isyes}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 6:00 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = isyes(question)
    if nargin == 0
        msg = "";
    elseif nargin == 1
        msg = string(question);
    end
    while true
        % WARNING: MATLAB 2018a compatibility: do NOT convert the arguments to `input()` to string. Arguments must be chars.
        answer = lower(string(input(msg{:}, 's')));
        if strcmp(answer, "y")
            itis = true;
            break;
        elseif strcmp(answer, "n")
            itis = false;
            break;
        else
            disp("Invalid answer. Please enter either y or n, then press enter.");
            continue;
        end
    end
end