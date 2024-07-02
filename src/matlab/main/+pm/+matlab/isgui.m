%>  \brief
%>  Return a scalar MATLAB logical ``true`` if and only if the
%>  current environment is the MATLAB Graphical User Interface (GUI),
%>  otherwise, return ``false`` implying that MATLAB engine is active
%>  without the GUI (e.g., called from the system shell command line).
%>
%>  \details
%>  This functionality is important for certain
%>  library features that require MATLAB GUI.
%>
%>  \warning
%>  This function relies on functionalities
%>  that are supported only in MATLAB ``> 2019a``.
%>
%>  \return
%>  ``itis``    :   The output scalar MATLAB logical that is ``true``
%>                  if and only if ParaMonte MATLAB library is being called
%>                  from MATLAB Graphical User Interface (GUI), otherwise, ``false``.
%>
%>  \interface{isgui}
%>  \code{.m}
%>
%>      itis = pm.matlab.isgui()
%>
%>  \endcode
%>
%>  \example{isgui}
%>  \include{lineno} example/matlab/isgui/main.m
%>  \output{isgui}
%>  \include{lineno} example/matlab/isgui/main.out.m
%>
%>  \final{isgui}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:39 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = isgui()
    persistent itis_persistent
    if isempty(itis_persistent)
        try
            % usejava() is sensitive to char vs. string input arguments. always input char.
            itis_persistent = usejava('desktop') && feature('ShowFigureWindows') && ~batchStartupOptionUsed; % batchStartupOptionUsed is introduced in R2019a and not supported in older versions of MATLAB
        catch
            itis_persistent = usejava('desktop') && feature('ShowFigureWindows');
        end
    end
    itis = itis_persistent;
end