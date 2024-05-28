%>  \brief
%>  Return a scalar MATLAB logical that is ``true`` if and
%>  only if the deletion of the input file fails,
%>  otherwise, return ``false``.
%>
%>  \details
%>  This function is a simple wrapper around the MATLAB function ``delete()``
%>  It is primarily meant to bring consistency to messaging file deletions
%>  if they fail. Such consistency is particularly required on Windows
%>  systems as the OS tends to lock file ownership to one application
%>  causing deletions to fail frequently.
%>
%>  \param[in]  file    :   The input scalar MATLAB string,
%>                          containing the file path to be deleted.
%>  
%>  \param[in]  desc    :   The input scalar MATLAB string, containing a descriptive message
%>                          to be printed on the MATLAB command line if the deletion task fails.
%>                          The input ``desc``, if not empty, will be added to the following template
%>                          before being displayed:<br>
%>                              "Failed to delete <desc> from the local disk. File may be protected."<br>
%>                          If the input ``desc`` is empty ``[]`` or empty string ``""``, a default
%>                          value ``"the requested file"`` for ``desc`` will be added to the template.<br>
%>                          (**optional**. If missing, no warning will be displayed upon failure.)
%>
%>  \return
%>  `failed`            :   The output scalar MATLAB logical that is ``true`` if and
%>                          only if the deletion of the input file fails,
%>                          otherwise, return ``false``.
%>
%>  \interface{rmfile}
%>  \code{.m}
%>
%>      failed = pm.sys.path.rmfile(file)
%>      failed = pm.sys.path.rmfile(file, desc)
%>
%>  \endcode
%>  \final{rmfile}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:28 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function failed = rmfile(file, desc)
    try
        delete(file);
    catch me
        if nargin == 2
            if isempty(desc) || desc == ""
                desc = "the requested file";
            end
            warning ( newline ...
                    + string(me.identifier) + " : " + string(me.message) + newline ...
                    + "Failed to delete " + string(desc) + " from the local disk." + newline ...
                    + "An application may possess the file ownership, particularly on Windows." + newline ...
                    + "You can try manually deleting it at:" + newline ...
                    + newline ...
                    + pm.io.tab + string(file) + newline ...
                    + newline ...
                    );
        end
    end
end