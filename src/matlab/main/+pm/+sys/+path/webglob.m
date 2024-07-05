%>  \brief
%>  Return a scalar MATLAB string or vector of MATLAB strings
%>  containing the fully-resolved paths matching the input pattern.<br>
%>
%>  \details
%>  This function is very similar to [pm.sys.path.glob](@ref glob).<br>
%>  However, unlike [pm.sys.path.glob](@ref glob), if the input
%>  pattern matches a World Wide Web link, it will also
%>  download the file to a temporary path on the system
%>  and the temporary download path as the output.<br>
%>
%>  \param[in]  pattern :   The input scalar MATLAB string containing either:<br>
%>                          <ol>
%>                              <li>    the pattern to search for paths on the current system.<br>
%>                                      Wildcards may be used for basenames and for the directory parts.<br>
%>                                      If pattern contains directory parts, then these will
%>                                      be included in the output ``pathList``.<br>
%>                                      Following wildcards can be used:<br>
%>                                      <ol>
%>                                          <li>    ``*`` match zero or more characters.<br>
%>                                          <li>    ``?`` match any single character.<br>
%>                                          <li>    ``[ab12]`` match one of the specified characters.<br>
%>                                          <li>    ``[^ab12]`` match none of the specified characters.<br>
%>                                          <li>    ``[a-z]`` match one character in range of characters.<br>
%>                                          <li>    ``{a,b,c}`` matches any one of strings ``a``, ``b``, or ``c``.<br>
%>                                          <li>    All above wildcards do not match a file separator.<br>
%>                                          <li>    ``**`` match zero or more characters including file separators.<br>
%>                                                  This can be used to match zero or more directory parts
%>                                                  and will recursively list matching names.<br>
%>                                      </ol>
%>                              <li>    the weblink to download and save locally on the system temporary folder.<br>
%>                          </ol>
%>
%>  \param[in]  anycase :   The input scalar MATLAB logical.<br>
%>                          If ``true``, the search will be case-sensitive.<br>
%>                          If ``false``, the search will be case-insensitive.<br>
%>                          On Windows, ``anycase`` is always reset to ``true`` even if user-specified.<br>
%>                          (**optional**. default = ``false`` on Unix and ``true`` on Windows.)
%>
%>  \return
%>  ``pathList``        :   The output MATLAB cell array of strings containing the files
%>                          or directories that match the path specified by string ``pattern``.<br>
%>  ``isdirList``       :   The output MATLAB cell array of the same size as ``pathList``,
%>                          each element of which is a MATLAB logical value that is ``true`` if
%>                          and only if the corresponding element of ``pathList`` is a directory.<br>
%>
%>  \warning
%>  Symbolic linked directories or junctions may cause an infinite loop when using the ``**``.<br>
%>
%>  \interface{webglob}
%>  \code{.m}
%>
%>      [pathList, isdirList] = pm.sys.path.webglob(pattern)
%>      [pathList, isdirList] = pm.sys.path.webglob(pattern, anycase)
%>
%>  \endcode
%>
%>  \example{webglob}
%>  \include{lineno} example/sys/path/webglob/main.m
%>  \output{webglob}
%>  \include{lineno} example/sys/path/webglob/main.out.m
%>
%>  \final{webglob}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:10 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function [pathList, isdirList] = webglob(pattern, anycase)
    if nargin < 2
        anycase = [];
    end
    [pathList, isdirList] = pm.sys.path.glob(pattern, anycase);
    if  isempty(pathList)
        %%%%
        %%%% check if the input path is a weblink.
        %%%%
        if pm.web.isurl(pattern)
            %%%%
            %%%% Extract the weblink file name.
            %%%%
            parts = split(pattern, "/");
            if parts(end) ~= ""
                pathtmp = fullfile(tempdir(), parts(end));
            else
                pathtmp = tempname();
            end
            pathList = [string(websave(pathtmp, pattern))];
        end
    end
end