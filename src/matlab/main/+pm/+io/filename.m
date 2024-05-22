%>  \brief
%>  Return a scalar MATLAB string containing
%>  an ``outputFileName`` base simulation specification
%>  of the ParaMonte samplers based on the current date and time.
%>
%>  \note
%>  This functionality is primarily used by the ParaMonte
%>  samplers to generate random output file base name for
%>  situations where the use does not specify a name.
%>
%>  \param[in]  prefix  :   The input scalar MATLAB string,
%>                          containing the desired prefix, e.g.,
%>                          ParaDRAM, ParaDISE, ParaNest, ...
%>                          (**optional**, default = ``""``)
%>
%>  \return
%>  `str`               :   The output scalar MATLAB string containing
%>                          an output file name based on the current date and time.<br>
%>                          The following represents the template used
%>                          across all programming language environments:<br>
%>                          <pre>
%>                              YYYYMMDD_HHMMSS_MMM
%>                          or
%>                              prefix_YYYYMMDD_HHMMSS_MMM
%>                          where,
%>                              1.  ``prefix``  is replaced with the input ``prefix`` name.
%>                              2.  ``YYYY``    is replaced with the current Gregorian year.
%>                              2.  ``MM``      is replaced with the current month of year.
%>                              2.  ``DD``      is replaced with the current day of month.
%>                              2.  ``HH``      is replaced with the current hour of day.
%>                              2.  ``MM``      is replaced with the current minute of hour.
%>                              2.  ``SS``      is replaced with the current second of minute.
%>                              2.  ``MMM``     is replaced with the current millisecond of second.
%>                          </pre>
%>
%>  \interface{filename}
%>  \code{.m}
%>
%>      str = pm.io.filename(prefix)
%>
%>  \endcode
%>  \final{filename}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 6:13 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = filename(prefix)
    if nargin < 1
        prefix = [];
    end
    str = datestr(now, 'yyyymmdd_HHMMSS_FFF');
    % MATLAB recommends replacing the above with the following.
    % The plan was to keep the old syntax as it is not clear 
    % when the new syntax was instroduces.
    %str = string(datetime('now', 'format', 'yyyyMMdd_HHmmss_SSS'))
    if ~isempty(prefix)
        str = string(prefix) + "_" + str;
    end
end