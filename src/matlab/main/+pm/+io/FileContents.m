%>  \brief
%>  This is the base class for generating objects
%>  that contain the contents of a given file.<br>
%>
%>  \details
%>  This class is meant to be primarily internally used
%>  by the ParaMonte library routines (e.g., samplers).<br>
%>
%>  \devnote
%>  The ``handle`` superclass of this class
%>  is critical for the class functionality.<br>
%>  See the documentation of the class constructor.<br>
%>
%>  \note
%>  See below for information on class attributes (properties).<br>
%>
%>  \note
%>  See below for information on the methods.<br>
%>
%>  \final{FileContents}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 6:03 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef FileContents < pm.matlab.Handle

    properties(Access = public)
        %>
        %>  silent  :   The scalar MATLAB logical (Boolean) indicator which is ``false`` by default.<br>
        %>              If it is set to ``true``, it will silence all output postprocessing messages.<br>
        %>
        silent = false;
        %>
        %>  file    :   The scalar MATLAB string containing the path to the file whose contents are read.
        %>
        file = "";
    end

    properties(Hidden)
        weblinks;
        spinner;
        timer;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Return a scalar object of class [pm.io.FileContents](@ref FileContents).
        %>
        %>  \details
        %>  This is the constructor of the class [pm.io.FileContents](@ref FileContents).<br>
        %>  It merely serves as the blueprint for the IO subclasses
        %>  accessible to the end users.<br>
        %>
        %>  \param[in]  file    :   The input scalar MATLAB string
        %>                          containing the path to an external file.
        %>
        %>  \param[in]  silent  :   The input scalar MATLAB logical.<br>
        %>                          If ``true``, all descriptive messages will be suppressed.<br>
        %>                          Setting this option to ``false`` is particularly useful
        %>                          in MPI-parallel simulations.<br>
        %>                          (**optional**, default = ``false``)
        %>
        %>  \return
        %>  `self`              :   The output scalar object of class [pm.io.FileContents](@ref FileContents).
        %>
        %>  \interface{FileContents}
        %>  \code{.m}
        %>
        %>      contents = pm.io.FileContents(file)
        %>      contents = pm.io.FileContents(file, [])
        %>      contents = pm.io.FileContents(file, silent)
        %>
        %>  \endcode
        %>
        %>  \final{FileContents}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:05 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = FileContents(file, silent)
            if  nargin < 2
                silent = [];
            end
            if  nargin < 1
                file = [];
            end
            if ~isempty(silent)
                self.silent = silent;
            end
            if ~isempty(file)
                self.file = string(file);
            end
            if ~self.silent
                disp("reading file: """ + self.file + """");
            end
            self.timer = pm.timing.Timer();
            self.spinner = pm.timing.Spinner();
            self.weblinks = pm.lib.weblinks();
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public, Hidden)

        %>  \brief
        %>  Display the input warning message about the line
        %>  number ``line`` of the file whose contents are read and return nothing.
        %>
        %>  \details
        %>  This is a ``Hidden`` method of the class [pm.io.FileContents](@ref FileContents).<br>
        %>
        %>  \param[in]  line    :   The input scalar MATLAB string or whole number,
        %>                          representing the line number within the file about which the warning message should be printed.<br>
        %>                          (**optional**, default = ``"UNKNOWN"``)
        %>
        %>  \param[in]  msg     :   The input scalar MATLAB string containing a message to display on the MATLAB console.<br>
        %>                          (**optional**, default = ``"done in " + sprintf("%.6f", string(self.timer.del())) + " seconds."``)
        %>
        %>  \interface{warn}
        %>  \code{.m}
        %>
        %>      fc = pm.io.FileContents(file)
        %>      fc.warn(line, msg)
        %>      fc.warn([], msg)
        %>      fc.warn([], [])
        %>      fc.warn()
        %>
        %>  \endcode
        %>
        %>  \final{warn}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:07 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function warn(self, line, msg)
            if  nargin < 3
                msg = [];
            end
            if  nargin < 2
                line = [];
            end
            if ~isempty(msg)
                msg = string(msg) + newline;
            else
                msg = "";
            end
            if ~isempty(line)
                line = string(line);
            else
                line = "UNKNOWN";
            end
            warning ( newline ...
                    + "The structure of the input file:" + newline ...
                    + newline ...
                    + pm.io.tab + self.file + newline ...
                    + newline ...
                    + "appears compromised around line: " + line + newline ...
                    + msg ...
                    + "The file parsing will proceed with no guarantee of success." + newline ...
                    + newline ...
                    );
        end

        %>  \brief
        %>  Display the input final message and return nothing.
        %>
        %>  \details
        %>  This is a ``Hidden`` method of the class [pm.io.FileContents](@ref FileContents).<br>
        %>
        %>  \param[in]  msg     :   The input scalar MATLAB string containing a
        %>                          message to display on the MATLAB console.<br>
        %>                          (**optional**, default = ``"done in " + sprintf("%.6f", string(self.timer.del())) + " seconds."``)
        %>
        %>  \param[in]  advance :   The input scalar MATLAB ``logical``.<br>
        %>                          If ``true``, an end of line character will be added at the end of the printed message.<br>
        %>                          (**optional**, default = ``true``)
        %>
        %>  \interface{checkpoint}
        %>  \code{.m}
        %>
        %>      fc = pm.io.FileContents(file)
        %>      fc.checkpoint(msg, advance)
        %>      fc.checkpoint([], advance)
        %>      fc.checkpoint([], [])
        %>      fc.checkpoint(msg)
        %>      fc.checkpoint([])
        %>      fc.checkpoint()
        %>
        %>  \endcode
        %>
        %>  \final{checkpoint}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 6:07 PM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function checkpoint(self, msg, advance)
            if ~self.silent
                if  nargin < 3
                    advance = [];
                end
                if  nargin < 2
                    msg = [];
                end
                if  isempty(advance)
                    advance = true;
                end
                if  isempty(msg)
                    msg = "done in " + sprintf("%.6f", string(self.timer.del())) + " seconds.";
                end
                if  advance
                    disp(msg);
                else
                    fprintf(msg);
                end
            end
        end

        %%>  \brief
        %%>  Return a copy of the specified ``field`` (component)
        %%>  of the parent object of class [pm.io.FileContents](@ref FileContents).
        %%>
        %%>  \details
        %%>  This method is an unfortunate result of the lack references in MATLAB.<br>
        %%>  The output of this method is used by the visualization methods of
        %%>  this class to repeatedly sync the internal copy of ``df`` with
        %%>  the original ``df`` component of the parent object.
        %%>
        %%>  \param[in]  field   :   The input scalar MATLAB string containing the
        %%>                          name of a field (component/attribute) of the parent
        %%>                          object whose value will have to be returned.
        %%>
        %%>  \return
        %%>  `val`               :   The output object containing the value of the
        %%>                          specified ``field`` of the parent object.
        %%>
        %%>  \interface{getVal}
        %%>  \code{.m}
        %%>
        %%>      fc = pm.io.FileContents(field)
        %%>      val = fc.getVal(field)
        %%>
        %%>  \endcode
        %%>
        %%>  \final{getVal}
        %%>
        %%>  \author
        %%>  \JoshuaOsborne, May 21 2024, 6:07 PM, University of Texas at Arlington<br>
        %%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        %function val = getVal(self, field)
        %    val = self.(field);
        %end

    end

end