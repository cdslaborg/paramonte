%>  \brief
%>  Finalize the ParaMonte MATLAB sampler simulation run and return nothing.
%>
%>  \param[in]  self    :   The input parent object of class [pm.sampling.Sampler](@ref Sampler)
%>                          which is **implicitly** passed to this dynamic method (not by the user).<br>
%>
%>  \note
%>  This is an internal method of the class [pm.sampling.Sampler](@ref Sampler).
%>
%>  \final{finalize}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:10 AM, University of Texas at Arlington<br>
%>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin%>
function finalize(self)
    if  self.partype == "openmp"
        if ~self.silent
            delete(gcp("nocreate"));
        else
            evalc('delete(gcp("nocreate")');
        end
    end
    munlock(self.mexname);
    clear(self.mexname);
    path(self.matpath);
    if  self.clstype == "gnu"
        setenv('GFORTRAN_STDIN_UNIT' , '-1');
        setenv('GFORTRAN_STDOUT_UNIT', '-1');
        setenv('GFORTRAN_STDERR_UNIT', '-1');
    end
end