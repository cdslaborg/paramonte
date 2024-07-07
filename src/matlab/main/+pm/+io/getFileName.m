%>  \brief
%>  Return a scalar MATLAB string containing
%>  a file name based on the current date and time.
%>
%>  \details
%>  This functionality is primarily used by the ParaMonte
%>  samplers to generate random output file base name for
%>  situations where the user does not specify a name.
%>
%>  \param[in]  prefix  :   The input scalar MATLAB string,
%>                          containing the desired prefix, e.g.,
%>                          ParaDRAM, ParaDISE, ParaNest, ...<br>
%>                          (**optional**, default = ``""``)
%>
%>  \return
%>  ``str``             :   The output scalar MATLAB string containing
%>                          an output file name based on the current date and time.<br>
%>                          The following represents the template used across
%>                          all programming language environments:<br>
%>                          \code{.m}
%>                              YYYYMMDD_HHMMSS_MMM
%>                          \endcode
%>                          or
%>                          \code{.m}
%>                              prefix_YYYYMMDD_HHMMSS_MMM
%>                          \endcode
%>                          where,
%>                          <ol>
%>                              <li>    ``prefix``  is replaced with the input ``prefix`` name.
%>                              <li>    ``YYYY``    is replaced with the current Gregorian year.
%>                              <li>    ``MM``      is replaced with the current month of year.
%>                              <li>    ``DD``      is replaced with the current day of month.
%>                              <li>    ``HH``      is replaced with the current hour of day.
%>                              <li>    ``MM``      is replaced with the current minute of hour.
%>                              <li>    ``SS``      is replaced with the current second of minute.
%>                              <li>    ``MMM``     is replaced with the current millisecond of second.
%>                          </ol>
%>
%>  \interface{getFileName}
%>  \code{.m}
%>
%>      str = pm.io.getFileName()
%>      str = pm.io.getFileName([])
%>      str = pm.io.getFileName(prefix)
%>
%>  \endcode
%>
%>  \example{getFileName}
%>  \include{lineno} example/io/getFileName/main.m
%>  \output{getFileName}
%>  \include{lineno} example/io/getFileName/main.out.m
%>
%>  \final{getFileName}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 6:13 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function str = getFileName(prefix)
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