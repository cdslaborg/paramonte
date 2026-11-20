%>  \brief Evaluate a batch using the legacy thread-pool parfor implementation.<br>
%>
%>  \details
%>  This private helper runs on an already initialized compatible thread pool.
%>  Each loop iteration evaluates a single state column and contributes its
%>  elapsed time to a reduction. The caller converts total time to the mean
%>  time and measures the overall callback wall time. No persistent queue API
%>  is needed, so this path is usable before MATLAB R2025a.<br>
%>  \note
%>  This helper is internal to Sampler; it is not a public alternative backend.
%>  MATLAB's scheduler determines how iteration work is partitioned.<br>
%>  \param[in] getLogFunc : Scalar objective function handle.
%>  \param[in] state : Candidate states, one per column.
%>  \param[out] logFunc : Ordered objective values.
%>  \param[out] totalTime : Sum of per-proposal objective times, in seconds.
function [logFunc, totalTime] = evaluateObjectiveParfor(getLogFunc, state)
    njob = size(state, 2);
    logFunc = zeros(njob, 1);
    totalTime = 0;
    parfor ijob = 1 : njob
        timer = tic();
        logFunc(ijob) = getLogFunc(state(:, ijob));
        totalTime = totalTime + toc(timer);
    end
end
