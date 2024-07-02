%>  \brief
%>  Return a scalar MATLAB string containing
%>  the current ParaMonte MATLAB library version number,
%>  or the previous version from the specified current version.
%>
%>  \param[in]  current :   The input scalar MATLAB string containing the
%>                          current Semantic version of the library in
%>                          triplet format (e.g., ``"1.2.3"``).<br>
%>
%>  \return
%>  `vernum`            :   The output scalar MATLAB string containing either:
%>                          <ol>
%>                              <li>    the ParaMonte MATLAB version if
%>                                      the input argument ``current`` is missing.<br>
%>                              <li>    the semantically-previous version of the library
%>                                      if the input argument ``current`` is present.<br>
%>                                      The output version **may not necessarily exist**
%>                                      as a version of the ParaMonte library.<br>
%>                          </ol>
%>
%>  \interface{version}
%>  \code{.m}
%>
%>      vernum = pm.lib.version()
%>      vernum = pm.lib.version(current)
%>
%>  \endcode
%>
%>  \example{version}
%>  \include{lineno} example/lib/version/main.m
%>  \output{version}
%>  \include{lineno} example/lib/version/main.out.m
%>
%>  \final{version}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 8:13 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function vernum = version(current)
    persistent vernum_persistent;
    if  nargin < 1
        if isempty(vernum_persistent)
            try
                fid = fopen(fullfile(pm.lib.path.auxil(), "VERSION.md"));
                vernum_persistent = string(strip(fgetl(fid)));
                fclose(fid);
            catch
                vernum_persistent = "UNKNOWN";
            end
        end
        vernum = vernum_persistent;
    else
        currentVerionTriplet = pm.str.split2real(current);
        vernum = [];
        index = 4;
        while true
            index = index - 1;
            if index <= 0
                return;
            else
                if currentVerionTriplet(index) > 0
                    previousVerionTriplet = currentVerionTriplet;
                    previousVerionTriplet(index) = previousVerionTriplet(index) - 1;
                    vernum = join(string(previousVerionTriplet),".");
                    return;
                end
            end
        end
    end
end