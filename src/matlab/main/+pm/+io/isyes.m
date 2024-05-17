function itis = isyes(question)
    %
    %   Return a scalar MATLAB logical that is ``true`` if the
    %   user response to a question on the MATLAB command prompt
    %   is ``Y`` or ``y``, representing YES, otherwise ``false``
    %   if the user response is ``N`` or ``n``.
    %
    %   Continue asking for as long as the user response
    %   is not among the accepted answers mentioned above.
    %
    %   Parameters
    %   ----------
    %
    %       question
    %
    %           The input scalar MATLAB string containing the
    %           question to be displayed on the MATLAB prompt.
    %
    %   Returns
    %   -------
    %
    %       itis
    %
    %           The output scalar MATLAB logical that is ``true`` if the
    %           user response to a question on the MATLAB command prompt
    %           is ``Y`` or ``y``, representing YES, otherwise `false`
    %           if the user response is ``N`` or ``n``.
    %
    %   Interface
    %   ---------
    %
    %       pm.io.isyes()
    %       pm.io.isyes(question)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if nargin == 0
        msg = "";
    elseif nargin == 1
        msg = string(question);
    end
    while true
        % WARNING: MATLAB 2018a compatibility: do NOT convert the arguments to `input()` to string. Arguments must be chars.
        answer = lower(string(input(msg{:}, 's')));
        if strcmp(answer, "y")
            itis = true;
            break;
        elseif strcmp(answer, "n")
            itis = false;
            break;
        else
            disp("Invalid answer. Please enter either y or n, then press enter.");
            continue;
        end
    end
end