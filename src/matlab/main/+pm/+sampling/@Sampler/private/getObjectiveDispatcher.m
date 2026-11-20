%>  \brief
%>  Dispatch batches of objective evaluations for a thread-parallel ParaMonte run.<br>
%>
%>  \details
%>  This private factory creates a run-scoped client dispatcher using nested
%>  function handles and a shared closure workspace.
%>  Call ``start()`` before ``evaluate()`` and ``stop()`` on completion or failure.
%>  The dispatcher may be started again after stopping; it is not passed to workers.
%>
%>  The persistent backend uses a shared-destination job queue, client-owned
%>  result and readiness queues, and one resident future per thread worker.
%>  On R2025a and later the persistent backend is selected automatically.
%>  Earlier thread-pool-capable releases use a timed ``parfor`` for each batch.
%>  No caller option overrides this release-based selection.
%>  Only one synchronous evaluation batch may be in flight. A result is never
%>  returned until all its original proposal indices have been received.<br>
%>
%>  \note
%>  Requires Parallel Computing Toolbox. Persistent queues require MATLAB R2025a
%>  or later. The dispatcher does not delete or resize a user's parallel pool.
%>  It does not implement the native sampling algorithm or invoke MATLAB from
%>  native OpenMP threads.<br>
%>
%>  \warning
%>  Internal implementation detail, not a public sampler specification. The
%>  objective must be thread-safe. No concurrent or nested use of the occupied
%>  pool is supported. Cancellation of uninterruptible code is best effort.<br>
%>
%>  \param[in] serveJobs : Named handle to the private resident service function,
%>                         resolved by the calling Sampler method.
%>  \param[in] evaluateLegacy : Named handle to the private timed parfor helper,
%>                              resolved by the calling Sampler method.
%>  \param[out] dispatcher : Structure of private function handles: ``start``,
%>                           ``evaluate``, ``stop`` and ``getBackend``. Copies of
%>                           this structure refer to the same closure workspace.
%>  \note
%>  The calling method must own an onCleanup guard and call ``stop()`` on every
%>  exit path. A structure has no class destructor. No public factory or class
%>  is registered in the pm.sampling namespace.<br>
function dispatcher = getObjectiveDispatcher(serveJobs, evaluateLegacy)
    % Each factory invocation has its own state; no global or persistent state.
    dispatcherState.backend = "";
    dispatcherState.numWorker = 0;
    dispatcherState.getLogFunc = [];
    dispatcherState.jobQueue = [];
    dispatcherState.resultQueue = [];
    dispatcherState.readyQueue = [];
    dispatcherState.futureList = {};
    dispatcherState.active = false;
    dispatcherState.busy = false;

    dispatcher.start = @start;
    dispatcher.evaluate = @evaluate;
    dispatcher.stop = @stop;
    dispatcher.getBackend = @getBackend;

    %>  \brief Return the resolved backend for this run.<br>
    %>  \param[out] backend : Resolved backend name, or empty before startup.
    function backend = getBackend()
        backend = dispatcherState.backend;
    end

    %>  \brief
    %>  Select a compatible thread pool and initialize the release-dependent dispatcher.<br>
    %>
    %>  \param[in] numThread : Positive finite integer; requested thread count.
    %>  \param[in] getLogFunc : Function handle accepting one state column and returning
    %>                         one real numeric scalar log density.
    %>  \param[in] silent : Logical indicating whether to suppress pool startup output.
    %>
    %>  \details
    %>  Resolve the release-dependent backend and reuse a compatible pool or create
    %>  one if absent. R2022b and later require an exact worker-count match; older
    %>  releases may use a different default/existing count with a warning.
    %>  Persistent mode rejects outstanding pool futures and waits up to 60 seconds
    %>  for all resident workers to report readiness. Failed partial startup is
    %>  aborted, and its original exception is rethrown.<br>
    function start(numThread, getLogFunc, silent)
        if dispatcherState.active || ~isempty(dispatcherState.futureList)
            error('pm:sampling:Dispatcher:AlreadyStarted', ...
                  'A dispatcher can only serve one active run.');
        end
        validateattributes(numThread, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
        if ~isa(getLogFunc, 'function_handle')
            error('pm:sampling:Dispatcher:InvalidObjective', 'Expected an objective function handle.');
        end
        if verLessThan('matlab', '9.8')
            error('pm:sampling:Dispatcher:UnsupportedRelease', ...
                  'MATLAB thread pools require R2020a or newer.');
        end
        % R2025a introduced both Destination="any" and queue close().
        % Do not silently fall back on newer releases if persistent startup fails.
        if verLessThan('matlab', '25.1')
            backend = "parfor";
        else
            backend = "persistent";
        end
        dispatcherState.backend = backend;
        dispatcherState.getLogFunc = getLogFunc;

        pool = getThreadPool(numThread, silent);
        if ~isa(pool, 'parallel.ThreadPool')
            error('pm:sampling:Dispatcher:IncompatiblePool', ...
                  'A non-thread pool is already open. Close it explicitly before running ParaMonte.');
        end
        if pool.NumWorkers ~= numThread
            if verLessThan('matlab', '9.13')
                warning('pm:sampling:Dispatcher:LegacyPoolSize', ...
                        'This MATLAB release cannot size thread pools; using %d workers instead of %d.', ...
                        pool.NumWorkers, numThread);
            else
                error('pm:sampling:Dispatcher:IncompatiblePool', ...
                      'Existing pool has %d workers; requested %d. Resize it explicitly.', ...
                      pool.NumWorkers, numThread);
            end
        end
        dispatcherState.numWorker = pool.NumWorkers;
        if backend == "parfor"
            dispatcherState.active = true;
            return;
        end

        % Do not monopolize a busy user pool with long-lived tasks.
        if ~isempty(pool.FevalQueue.RunningFutures) || ~isempty(pool.FevalQueue.QueuedFutures)
            error('pm:sampling:Dispatcher:BusyPool', ...
                  'The thread pool has outstanding futures; finish them before starting ParaMonte.');
        end
        try
            dispatcherState.jobQueue = parallel.pool.PollableDataQueue('Destination', 'any');
            dispatcherState.resultQueue = parallel.pool.PollableDataQueue();
            dispatcherState.readyQueue = parallel.pool.PollableDataQueue();
            dispatcherState.futureList = cell(dispatcherState.numWorker, 1);
            for iworker = 1 : dispatcherState.numWorker
                dispatcherState.futureList{iworker} = parfeval(pool, serveJobs, 0, ...
                    dispatcherState.jobQueue, dispatcherState.resultQueue, dispatcherState.readyQueue, getLogFunc);
            end
            % A readiness barrier prevents startup cost entering the first batch.
            % Bound the barrier, not the user's objective evaluation time.
            timer = tic();
            for iworker = 1 : dispatcherState.numWorker
                while true
                    remaining = 60 - toc(timer);
                    if remaining <= 0
                        error('pm:sampling:Dispatcher:StartupTimeout', ...
                              'Resident workers did not start within 60 seconds.');
                    end
                    [~, received] = poll(dispatcherState.readyQueue, min(1, remaining));
                    if received
                        break;
                    end
                    checkWorkers();
                end
            end
            checkWorkers();
            close(dispatcherState.readyQueue);
            dispatcherState.readyQueue = [];
            dispatcherState.active = true;
        catch me
            dispatcherState.busy = true; % Abort, rather than drain, partial startup.
            stop();
            rethrow(me);
        end
    end

    %>  \brief
    %>  Evaluate one synchronous batch and restore the original proposal ordering.<br>
    %>
    %>  \param[in] state : An ``ndim``-by-``njob`` numeric matrix of candidate states,
    %>                     one proposal per column.
    %>  \param[out] logFunc : An ``njob``-by-1 double vector of log densities in
    %>                        the order of the input columns.
    %>  \param[out] avgTimePerFunCall : Arithmetic mean of the measured objective
    %>                                 evaluation times, in seconds.
    %>  \param[out] avgCommPerFunCall : Batch wall time minus the preceding mean,
    %>                                 clamped at zero, in seconds. Includes dispatch,
    %>                                 synchronization and potentially load imbalance;
    %>                                 not an isolated communication measurement.
    %>
    %>  \details
    %>  For persistent workers, send chunks of at most
    %>  ``max(1, ceil(njob / (4 * numWorker)))`` consecutive proposals. Poll the
    %>  result queue and validate indices, uniqueness and completeness before
    %>  returning. Check worker health after each one-second poll timeout; no
    %>  arbitrary deadline is imposed on the objective itself. Empty batches
    %>  return an empty column vector and zero timings. Invalid or failed batches
    %>  abort the dispatcher and rethrow the error; partial results are not returned.
    %>  Real scalar NaN and infinite values are preserved in persistent mode for
    %>  handling by the native sampler.<br>
    function [logFunc, avgTimePerFunCall, avgCommPerFunCall] = evaluate(state)
        if ~dispatcherState.active || dispatcherState.busy
            error('pm:sampling:Dispatcher:NotReady', 'Objective dispatcher is not ready.');
        end
        njob = size(state, 2);
        logFunc = zeros(njob, 1);
        avgTimePerFunCall = 0;
        avgCommPerFunCall = 0;
        if njob == 0
            return;
        end
        dispatcherState.busy = true;
        try
            timer = tic();
            if dispatcherState.backend == "parfor"
                % The named private helper does not capture closure state.
                [logFunc, avgTimePerFunCall] = evaluateLegacy(dispatcherState.getLogFunc, state);
            else
                % Typical ParaMonte batches contain one proposal per worker.
                % Coalesce larger batches to bound message count while allowing
                % dynamic assignment from one common queue.
                chunkSize = max(1, ceil(njob / (4 * dispatcherState.numWorker)));
                numChunk = ceil(njob / chunkSize);
                for firstJob = 1 : chunkSize : njob
                    indices = firstJob : min(njob, firstJob + chunkSize - 1);
                    send(dispatcherState.jobQueue, {indices, state(:, indices)});
                end
                seen = false(njob, 1);
                for ichunk = 1 : numChunk
                    while true
                        [result, received] = poll(dispatcherState.resultQueue, 1);
                        if received
                            break;
                        end
                        % No arbitrary objective-function deadline: detect exits
                        % and cancellation without rejecting costly objectives.
                        checkWorkers();
                    end
                    if ~isnumeric(result) || size(result, 1) ~= 3
                        error('pm:sampling:Dispatcher:BadResponse', 'Malformed objective response.');
                    end
                    indices = result(1, :);
                    if isempty(indices) || any(indices < 1 | indices > njob | indices ~= fix(indices) | ~isfinite(indices)) ...
                            || any(diff(indices) ~= 1) || any(seen(indices))
                        error('pm:sampling:Dispatcher:BadResponse', 'Invalid or duplicated objective response.');
                    end
                    seen(indices) = true;
                    logFunc(indices) = result(2, :).';
                    avgTimePerFunCall = avgTimePerFunCall + sum(result(3, :));
                end
                if ~all(seen)
                    error('pm:sampling:Dispatcher:MissingResponse', 'Missing objective response.');
                end
            end
            avgTimePerFunCall = avgTimePerFunCall / njob;
            % Preserve the existing native timing convention (batch wall time
            % minus average per-proposal compute time). This is not pure IPC time.
            avgCommPerFunCall = max(0, toc(timer) - avgTimePerFunCall);
            dispatcherState.busy = false;
        catch me
            % No retries: a failed batch must never reach the native sampler.
            stop();
            rethrow(me);
        end
    end

    %>  \brief
    %>  Stop this dispatcher's resident futures without deleting the thread pool.<br>
    %>  \details
    %>  Close the job queue to release blocked polls. On error or partial startup,
    %>  also close the result queue and cancel owned futures. Allow a single shared
    %>  five-second graceful wait budget, then request cancellation of remaining
    %>  futures, close queues and release client references. Cleanup warnings do
    %>  not deliberately replace the original worker exception. This method is
    %>  idempotent; uninterruptible operations may not stop promptly.<br>
    function stop()
        % Idempotent. Normal completion closes an empty job queue, allowing
        % each blocked poll to exit. Error/interrupt cancels only our futures.
        abort = dispatcherState.busy || ~dispatcherState.active;
        dispatcherState.active = false;
        % Wake blocked polls first. On failure also close the result queue,
        % so a worker that finishes during cancellation cannot keep draining
        % pending proposals indefinitely.
        if abort
            closeQueue(dispatcherState.resultQueue);
        end
        closeQueue(dispatcherState.jobQueue);
        if abort
            cancelWorkers();
        end
        % One deadline for the entire pool, rather than one per worker.
        timer = tic();
        for iworker = 1 : numel(dispatcherState.futureList)
            future = dispatcherState.futureList{iworker};
            if isempty(future), continue; end
            try
                if ~strcmp(future.State, 'finished')
                    remaining = max(0, 5 - toc(timer));
                    if remaining <= 0 || ~wait(future, 'finished', remaining)
                        cancel(future);
                    end
                end
            catch me
                warning('pm:sampling:Dispatcher:CleanupFailed', '%s', me.message);
            end
        end
        cancelWorkers(); % Also covers a worker that missed the graceful deadline.
        closeQueue(dispatcherState.readyQueue);
        closeQueue(dispatcherState.resultQueue);
        dispatcherState.futureList = {};
        dispatcherState.jobQueue = [];
        dispatcherState.readyQueue = [];
        dispatcherState.resultQueue = [];
        dispatcherState.getLogFunc = [];
        dispatcherState.busy = false;
    end

    %>  \brief
    %>  Detect unexpectedly finished or unavailable resident futures.<br>
    %>  \details
    %>  Propagate a worker exception as a ``WorkerFailed`` error with the original
    %>  exception as its cause; otherwise signal ``WorkerExited``.<br>
    function checkWorkers()
        for iworker = 1 : numel(dispatcherState.futureList)
            future = dispatcherState.futureList{iworker};
            if isempty(future), continue; end
            if strcmp(future.State, 'finished') || strcmp(future.State, 'unavailable')
                workerError = future.Error;
                if iscell(workerError) && ~isempty(workerError)
                    workerError = workerError{1};
                end
                if isa(workerError, 'MException')
                    failure = MException('pm:sampling:Dispatcher:WorkerFailed', ...
                                         'Resident objective worker %d failed: %s', iworker, workerError.message);
                    failure = addCause(failure, workerError);
                    throw(failure);
                end
                error('pm:sampling:Dispatcher:WorkerExited', ...
                      'Resident objective worker %d exited unexpectedly.', iworker);
            end
        end
    end
    %>  \brief
    %>  Request cancellation of all unfinished futures owned by this dispatcher.<br>
    function cancelWorkers()
        for iworker = 1 : numel(dispatcherState.futureList)
            future = dispatcherState.futureList{iworker};
            if isempty(future), continue; end
            try
                if ~strcmp(future.State, 'finished')
                    cancel(future);
                end
            catch me
                warning('pm:sampling:Dispatcher:CleanupFailed', '%s', me.message);
            end
        end
    end
end

%>  \brief Obtain or create a thread pool without replacing an existing pool.<br>
%>  \param[in] numThread : Requested worker count.
%>  \param[in] silent : Suppress pool startup messages if true.
%>  \param[out] pool : Existing or newly created pool (validated by start).
function pool = getThreadPool(numThread, silent)
    pool = gcp('nocreate');
    if isempty(pool)
        if verLessThan('matlab', '9.13') % R2022b introduced sized thread pools.
            if silent
                evalc('pool = parpool("threads");');
            else
                pool = parpool('threads');
            end
        elseif silent
            evalc('pool = parpool("threads", numThread);');
        else
            pool = parpool('threads', numThread);
        end
    end
end

%>  \brief
%>  Close a nonempty open queue, reporting cleanup failures as warnings.<br>
%>  \param[in] queue : A pollable data queue or an empty value.
function closeQueue(queue)
    if isempty(queue), return; end
    try
        if ~queue.IsClosed
            close(queue);
        end
    catch me
        warning('pm:sampling:Dispatcher:CleanupFailed', '%s', me.message);
    end
end
