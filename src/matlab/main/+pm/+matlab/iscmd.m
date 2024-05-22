%>  \brief
%>  Return a scalar MATLAB logical ``true`` if and only if the
%>  the MATLAB binary is being called from the shell command-line,
%>  in ``batch`` mode (without Graphical User Interface (GUI) interface).
%>  Otherwise, return ``false`` implying that MATLAB GUI engine is active.
%>
%>  \details
%>  This functionality is important for certain library
%>  features that require MATLAB to be called from the
%>  command-line, e.g., MPI-parallel functionalities.
%>
%>  \warning
%>  This function relies on functionalities
%>  that are supported only in MATLAB > 2019a.
%>
%>  \param[in]  `None`
%>
%>  \return
%>  `itis`  :   The output scalar MATLAB logical that is ``true``
%>              if and only if ParaMonte MATLAB library is being called
%>              from MATLAB Graphical User Interface (GUI), otherwise, ``false``.
%>
%>  \interface{iscmd}
%>  \code{.m}
%>
%>      itis = pm.matlab.iscmd()
%>
%>  \endcode
%>  \final{iscmd}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 11:38 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function itis = iscmd()
    try
        %   This error catching can be removed once
        %   no MATLAB < R2019a is expected to call it.
        itis = batchStartupOptionUsed; % batchStartupOptionUsed is introduced in R2019a and not supported in older versions of MATLAB
    catch
        try
            itis = ~pm.matlab.isgui();
        catch
            itis = isdeployed();
        end
    end
end