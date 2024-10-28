%>  \brief
%>  Generate and return the relevant Post-Processing Message (ppm) for the current ParaMonte
%>  sampler to be displayed on the MATLAB command line after the sampling is complete.<br>
%>
%>  \param[in]  self    :   The input parent object of class [pm.sampling.Sampler](@ref Sampler)
%>                          which is **implicitly** passed to this dynamic method (not by the user).<br>
%>
%>  \note
%>  This is an internal method of the class [pm.sampling.Sampler](@ref Sampler).
%>
%>  \final{getppm}
%>
%>  \author
%>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin%>
function ppm = getppm(self)
    ppm = "Use the following object methods to read the generated basic output files: " + newline ...
        + newline ...
        + pm.io.tab + self.name + ".readChain()    % Return a list of the contents of the output chain file(s)." + newline ...
        + pm.io.tab + self.name + ".readSample()   % Return a list of i.i.d. sample(s) from the output sample file(s)." + newline ...
        + pm.io.tab + self.name + ".readReport()   % Return a list of summary report(s) from the output report file(s)." + newline ...
        + pm.io.tab + self.name + ".readRestart()  % Return a list of the contents of the ASCII output restart file(s)." + newline ...
        + pm.io.tab + self.name + ".readProgress() % Return a list of the contents of the output progress file(s)." + newline ...
        ;
end