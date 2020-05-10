function isYes = getUserResponse(varargin)
    if nargin==0 || (nargin==1 && isempty(varargin{1})); msg = ""; end
    while true
        answer = lower(string(input(msg,"s")));
        if strcmp(answer,"y")
            isYes = true;
            break;
        elseif strcmp(answer,"n")
            isYes = false;
            break;
        else
            disp("Invalid answer. Please enter either y or n, then press enter.");
            continue;
        end
    end
end
