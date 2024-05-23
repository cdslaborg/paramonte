%>  \brief
%>  Finalize the ParaMonte MATLAB sampler simulation run.
%>
%>  \param[in]  `None`
%>
%>  \return
%>  `None`
%>  \final{getBorder}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:10 AM, University of Texas at Arlington<br>
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