%>  \brief
%>  Serve objective-function requests in one long-lived MATLAB thread worker.<br>
%>
%>  \details
%>  The client launches this private function once per resident ``parfeval``
%>  future. Send a readiness message, then block until a proposal chunk is
%>  available or the job queue closes. The objective handle and captured data
%>  remain available throughout the loop. No values are returned as future
%>  outputs; all results travel through the result queue.<br>
%>
%>  \param[in] jobQueue : A shared-destination ``parallel.pool.PollableDataQueue``.
%>                       Each message is a two-element cell array whose first
%>                       element is a vector of one-based proposal indices and
%>                       whose second element is the matching state matrix.
%>  \param[in] resultQueue : A client-owned pollable queue. Each response is a
%>                          3-by-``numJob`` double matrix: row 1 contains indices,
%>                          row 2 contains log densities, and row 3 contains the
%>                          elapsed objective evaluation times in seconds.
%>  \param[in] readyQueue : A client-owned pollable queue to which a logical
%>                         ``true`` readiness message is sent once at startup.
%>  \param[in] getLogFunc : A thread-compatible function handle accepting a
%>                          single state column and returning a real numeric scalar.
%>
%>  \note
%>  Real scalar NaN, positive infinity and negative infinity are passed through
%>  unchanged for native sampler handling. Invalid values cause an exception.
%>  The client detects failed futures and discards the incomplete batch.
%>  MATLAB R2025a or later and Parallel Computing Toolbox are required.<br>
%>
%>  \warning
%>  Do not call this blocking service directly from the MATLAB client. It is
%>  intended for the resident futures managed by ``getObjectiveDispatcher()``.<br>
function serveObjectiveJobs(jobQueue, resultQueue, readyQueue, getLogFunc)
% Internal worker service. One parfeval invocation per worker per run.
% The objective function (including captured data) stays resident here.
% Closing the job queue wakes idle workers without requiring poison messages.
    send(readyQueue, true);
    while true
        [job, received] = poll(jobQueue, Inf);
        if ~received
            return;
        end
        indices = job{1};
        state = job{2};
        numJob = numel(indices);
        result = zeros(3, numJob);
        result(1, :) = indices;
        for ijob = 1 : numJob
            start = tic();
            value = getLogFunc(state(:, ijob));
            elapsed = toc(start);
            if ~isnumeric(value) || ~isscalar(value) || ~isreal(value)
                error('pm:sampling:Dispatcher:InvalidObjectiveValue', ...
                      'The objective must return one real numeric log-density per proposal.');
            end
            % Preserve +/-Inf and NaN for the native sampler's existing policy.
            result(2, ijob) = double(value);
            result(3, ijob) = elapsed;
        end
        % An exception terminates this future; the client detects it and aborts
        % the batch instead of returning partial or fabricated log densities.
        send(resultQueue, result);
    end
end
