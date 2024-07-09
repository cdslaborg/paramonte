%>  \brief
%>  Return a scalar MATLAB string containing
%>  the current ParaMonte generic or MATLAB library version numbers,
%>  or optionally, the major, minor, or the patch version number.
%>
%>  \param[in]  lang    :   The input scalar MATLAB string that can be either:<br>
%>                          <ol>
%>                              <li>    the value ``"generic"``, indicating
%>                                      the generic library version found in the root directory of
%>                                      the [ParaMonte library repository](https://github.com/cdslaborg/paramonte).<br>
%>                              <li>    the value ``"matlab"``, indicating
%>                                      the MATLAB library version found in the ``src/matlab`` subdirectory in
%>                                      the [ParaMonte library repository](https://github.com/cdslaborg/paramonte).<br>
%>                          </ol>
%>                          The input version must be in triplet format (e.g., ``"1.2.3"``).<br>
%>                          If the input value is empty ``[]``, the default value is used.<br>
%>                          (**optional**, default = ``"matlab"``)
%>  \param[in]  type    :   The input scalar MATLAB string that can be either:<br>
%>                          <ol>
%>                              <li>    the value ``"all"``, indicating
%>                                      the full version number output in the format ``"major.minor.patch"``.<br>
%>                              <li>    the value ``"major"``, indicating
%>                                      the full version number output in the format ``"major"``.<br>
%>                              <li>    the value ``"minor"``, indicating
%>                                      the full version number output in the format ``"minor"``.<br>
%>                              <li>    the value ``"patch"``, indicating
%>                                      the full version number output in the format ``"patch"``.<br>
%>                          </ol>
%>                          If the input value is empty ``[]``, the default value is used.<br>
%>                          (**optional**, default = ``"all"``)
%>
%>  \return
%>  ``vernum``          :   The output scalar MATLAB string containing the requested
%>                          current generic or language-specific ParaMonte library version.<br>
%>
%>  \interface{version}
%>  \code{.m}
%>
%>      vernum = pm.lib.version()
%>      vernum = pm.lib.version([])
%>      vernum = pm.lib.version(lang)
%>
%>      vernum = pm.lib.version([], [])
%>      vernum = pm.lib.version([], type)
%>      vernum = pm.lib.version(lang, [])
%>      vernum = pm.lib.version(lang, type)
%>
%>  \endcode
%>
%>  \warning
%>  If the routine fails to find the respective VERSION files within the local ParaMonte library,
%>  it will return the value ``"UNKNOWN"`` for the output ``vernum``.
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
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function vernum = version(lang, type)

    if  nargin < 2
        type = [];
    end
    if  nargin < 1
        lang = [];
    end
    if ~isempty(type)
        type = string(lower(type));
    else
        type = "all";
    end
    if ~isempty(lang)
        lang = string(lower(lang));
    else
        lang = "matlab";
    end

    % persistent vernum_persistent;
    % if  isempty(vernum_persistent)
    %     vernum_persistent = struct();
    % end
    %
    % if ~isfield(vernum_persistent, lang)
    %     try
    %         fid = fopen(fullfile(pm.lib.path.auxil(), "VERSION." + lang + ".md"));
    %         vernum_persistent.(lang) = string(strip(fgetl(fid)));
    %         fclose(fid);
    %     catch
    %         vernum_persistent.(lang) = "UNKNOWN";
    %     end
    % end
    %
    % vernum = vernum_persistent.(lang);

    try
        version_file = fullfile(pm.lib.path.auxil(), "VERSION." + lang + ".md");
        fid = fopen(version_file);
        vernum = string(strip(fgetl(fid)));
        fclose(fid);
    catch
        vernum = "UNKNOWN";
    end

    if ~strcmp(type, "all")
        triplet = strsplit(vernum, ".");
        if  strcmp(type, "major")
            vernum = string(triplet(1));
        elseif strcmp(type, "minor") && numel(triplet) > 1
            vernum = string(triplet(2));
        elseif strcmp(type, "patch") && numel(triplet) > 2
            vernum = string(triplet(3));
        end
    end

end