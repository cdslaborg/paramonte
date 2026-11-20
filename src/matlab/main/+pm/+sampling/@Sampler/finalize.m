%>  \brief
%>  Release run-scoped sampler resources and return nothing.<br>
%>
%>  \details
%>  This is a ``private`` and ``Hidden`` dynamic method of the class [pm.sampling.Sampler](@ref Sampler).<br>
%>  The run cleanup guard calls this method after normal completion or when
%>  MATLAB unwinds the run due to an exception or interrupt. It stops resident
%>  objective tasks without deleting the thread pool, unlocks and clears the
%>  MEX function, restores the saved MATLAB path, and restores the saved GNU
%>  Fortran environment variables when the GNU toolchain was selected.<br>
%>
%>  \param[in]  self    :   The input parent object of class [pm.sampling.Sampler](@ref Sampler)
%>                          which is **implicitly** passed to this dynamic method (not by the user).<br>
%>
%>  \param[in]  dispatcher  :   The run-scoped structure of function handles returned by the
%>                               private ``getObjectiveDispatcher()`` factory for an OpenMP run,
%>                               or empty for a serial or MPI run. May be omitted.<br>
%>  \param[in]  originalEnv :   A three-element cell array of the original values of
%>                               ``GFORTRAN_STDIN_UNIT``, ``GFORTRAN_STDOUT_UNIT``, and
%>                               ``GFORTRAN_STDERR_UNIT``, in that order. Empty or omitted
%>                               values disable environment restoration.<br>
%>
%>  \note
%>  Worker shutdown is idempotent. Cleanup exceptions from worker shutdown and
%>  MEX unlocking/clearing are reported as warnings; cleanup of uninterruptible
%>  operations or a fatal process crash is not guaranteed.<br>
%>  \final{finalize}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:10 AM, University of Texas at Arlington<br>
%>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin<br>
function finalize(self, dispatcher, originalEnv)
    % Stop only this simulation's resident tasks, never a user-owned pool.
    if nargin > 1 && ~isempty(dispatcher)
        try
            dispatcher.stop();
        catch me
            warning('pm:sampling:Sampler:CleanupFailed', '%s', me.message);
        end
    end
    % Restore the path and GNU environment even if worker cleanup failed.
    try
        munlock(self.mexname);
        clear(self.mexname);
    catch me
        warning('pm:sampling:Sampler:MexCleanupFailed', '%s', me.message);
    end
    path(self.matpath);
    if self.clstype == "gnu" && nargin > 2 && ~isempty(originalEnv)
        setenv('GFORTRAN_STDIN_UNIT',  originalEnv{1});
        setenv('GFORTRAN_STDOUT_UNIT', originalEnv{2});
        setenv('GFORTRAN_STDERR_UNIT', originalEnv{3});
    end
end
