%>  \brief
%>  Display the components of an input MATLAB variable on MATLAB Console recursively.
%>
%>  \details
%>  This function is particularly useful for displaying the
%>  hierarchical contents of a MATLAB ``struct`` or ``cell`` object.
%>
%>  \param[in]  obj     :   The input variable whose contents are to be displayed on MATLAB Console.
%>                          <ol>
%>                              <li>    If the input ``obj`` is a MATLAB ``struct``, the function will
%>                                      recursively show the ``obj`` name, its fieldnames, and their contents.
%>                              <li>    If the input ``obj`` is a cell array, the contents of each cell are displayed.
%>                          </ol>
%>                          (**optional**. default = ``[]``, effectively adding a new line to the command line.)
%>  \param[in]  name    :   The input scalar MATLAB string, representing the actual name of the input ``obj``.<br>
%>                          (**optional**. If missing, the ``obj`` variable name from the caller workspace is used.)
%>  \param[in]  hidden  :   The input scalar MATLAB ``logical``.<br>
%>                          If ``true``, then the contents of the input ``struct`` will not be shown, only the field names.<br>
%>                          (**optional**, default = ``false``)
%>
%>  \interface{show}
%>  \code{.m}
%>
%>      pm.matlab.show();
%>      pm.matlab.show([]);
%>      pm.matlab.show(obj);
%>      pm.matlab.show(obj, name);
%>      pm.matlab.show(obj, [], hidden);
%>      pm.matlab.show(obj, name, hidden);
%>
%>  \endcode
%>
%>  \example{show}
%>  \include{lineno} example/matlab/show/main.m
%>  \output{show}
%>  \include{lineno} example/matlab/show/main.out.m
%>
%>  \final{show}
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function show(obj, name, hidden)

    if  nargin < 3
        hidden = [];
    end
    if  nargin < 2
        name = [];
    end
    if  nargin < 1
        fprintf('\n');
        return;
    end

    if  isempty(hidden)
        hidden = false;
    end
    if  isempty(name)
        name = inputname(1);
    end
    name = convertStringsToChars(name);

    if isstruct(obj)

        %%%% The number of elements to be displayed.

        objsize = numel(obj);
        if ~hidden
            hmax = objsize;
        else
            hmax = min(1, objsize);
        end

        %%%% Recursively display structure including fieldnames.

        for h = 1 : hmax
            fnames = fieldnames(obj(h));
            fnamesCount = length(fnames);
            for i = 1 : fnamesCount
                if objsize > 1
                    siz = size(obj);
                    if ~hidden
                        namei = [name '(' ind2str(siz, h) ').' fnames{i}];
                    else
                        namei = [name '(' ind2str(siz, objsize) ').'  fnames{i}];
                    end
                else
                    namei = [name '.' fnames{i}];
                end
                if isstruct(obj(h).(fnames{i}))
                    pm.matlab.show(obj(h).(fnames{i}), namei, hidden);
                else
                    if iscell(obj(h).(fnames{i}))
                        siz = size(obj(h).(fnames{i}));
                        NC = numel(obj(h).(fnames{i}));
                        if ~hidden
                            jmax = NC;
                        else
                            jmax = 1;
                        end
                        for j = 1 : jmax
                            if ~hidden
                                Namej = [namei '{' ind2str(siz,j) '}'];
                            else
                                Namej = [namei '{' ind2str(siz,NC) '}'];
                            end
                            disp(Namej);
                            if ~hidden
                                disp(obj(h).(fnames{i}){j});
                            end
                        end
                    else
                        disp(namei);
                        if ~hidden
                            disp(obj(h).(fnames{i}));
                        end
                    end
                end
            end
        end

    elseif iscell(obj)

        %%%% Recursively display cell.

        siz = size(obj);
        for i = 1 : numel(obj)
            namei = [name '{' ind2str(siz,i) '}'];
            pm.matlab.show(obj{i}, namei, hidden);
        end

    else

        disp(name);
        if ~hidden
            disp(obj)
        end

    end

end

%>  \cond excluded

function str = ind2str(siz, ndx)
    %%%% Treat vectors and scalars correctly.
    n = length(siz);
    if n == 2
        if siz(1) == 1
            siz = siz(2);
            n = 1;
        elseif siz(2) == 1
            siz = siz(1);
            n = 1;
        end
    end
    k = [1 cumprod(siz(1 : end - 1))];
    ndx = ndx - 1;
    str = '';
    for i = n : -1 : 1
        v = floor(ndx / k(i)) + 1;
        if i == n
            str = num2str(v);
        else
            str = [num2str(v) ',' str];
        end
        ndx = rem(ndx, k(i));
    end
end

%>  \endcond excluded