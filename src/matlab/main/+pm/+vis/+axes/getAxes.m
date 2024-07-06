%>  \brief
%>  Create axes in the specified input tiled position and return its handle.<br>
%>
%>  \param[in]  varargin    :   Any attributes of the axes, which can be one of the following positional arguments:<br>
%>                              <ol>
%>                                  <li>    ``rows ``   :   The input scalar MATLAB whole number,
%>                                                          representing the number of rows of subplots in the figure.<br>
%>                                  <li>    ``cols ``   :   The input scalar MATLAB whole number,
%>                                                          representing the number of columns of subplots in the figure.<br>
%>                                  <li>    ``cellx``   :   The input scalar MATLAB whole number,
%>                                                          representing the row at which the new axes must be placed in the figure.<br>
%>                                  <li>    ``celly``   :   The input scalar MATLAB whole number,
%>                                                          representing the column at which the new axes must be placed in the figure.<br>
%>                                  <li>    ``spanx``   :   The input scalar MATLAB whole number in the range ``[0, 1]``,
%>                                                          representing the X-span of the current axes in the figure.<br>
%>                                  <li>    ``spany``   :   The input scalar MATLAB whole number in the range ``[0, 1]``,
%>                                                          representing the Y-span of the current axes in the figure.<br>
%>                              </ol>
%>                              or specified as pairs of scalar ``char``,
%>                              followed by the specifications value in the range ``[0, 1]``.<br>
%>                              The following attributes are acceptable:<br>
%>                              <ol>
%>                                  <li>    ``'Spacing'``       , also ``'s'``
%>                                  <li>    ``'SpacingHoriz'``  , also ``'sh'``
%>                                  <li>    ``'SpacingVert'``   , also ``'sv'``
%>                                  <li>    ``'Padding'``       , also ``'p'``
%>                                  <li>    ``'PaddingRight'``  , also ``'pr'``
%>                                  <li>    ``'PaddingLeft'``   , also ``'pl'``
%>                                  <li>    ``'PaddingTop'``    , also ``'pt'``
%>                                  <li>    ``'PaddingBottom'`` , also ``'pb'``
%>                                  <li>    ``'Margin'``        , also ``'m'``
%>                                  <li>    ``'MarginRight'``   , also ``'mr'``
%>                                  <li>    ``'MarginLeft'``    , also ``'ml'``
%>                                  <li>    ``'MarginTop'``     , also ``'mt'``
%>                                  <li>    ``'MarginBottom'``  , also ``'mb'``
%>                                  <li>    ``'Holdaxis'``      , also ``'h'`` : If ``true``, any axes below the current will not be deleted.
%>                              </ol>
%>                              (**optional**. If missing, a fedault value will be used.)
%>
%>  \interface{getAxes}
%>  \code{.m}
%>
%>      h = getAxes(rows, cols, cellno, varargin)
%>      h = getAxes(rows, cols, cellx, celly, varargin)
%>      h = getAxes(rows, cols, cellx, celly, spanx, spany, varargin)
%>
%>  \endcode
%>
%>  \example{getAxes}
%>  \include{lineno} example/vis/axes/getAxes/main.m
%>  \vis{getAxes}
%>  \image html example/vis/axes/getAxes/getAxes.png width=700
%>
%>  \final{getAxes}
%>
%>  This code builds upon the works of 2001-2014 / Aslak Grinsted.<br>
%>
%>  \author
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function h = getAxes(varargin)
    f = gcf;
    userDataArgsOK = 0;
    args = get(f, 'UserData');
    if isstruct(args)
        userDataArgsOK = isfield(args, 'SpacingHorizontal') & isfield(args, 'Holdaxis') & isfield(args, 'rows') & isfield(args, 'cols');
    end
    OKToStoreArgs = isempty(args) | userDataArgsOK;
    if isempty(args) && (~userDataArgsOK)
        args = struct('Holdaxis',0, ...
            'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
            'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
            'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
            'rows',[],'cols',[]);
    end
    args = parseArgs(varargin,args,{'Holdaxis'},{'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
    if (length(args.NumericArguments)>2)
        args.rows=args.NumericArguments{1};
        args.cols=args.NumericArguments{2};
    %remove these 2 numerical arguments
        args.NumericArguments={args.NumericArguments{3:end}};
    end
    if OKToStoreArgs
        set(f,'UserData',args);
    end
    switch length(args.NumericArguments)
       case 0
           return % no arguments but rows/cols....
       case 1
           if numel(args.NumericArguments{1}) > 1 % restore subplot(m,n,[x y]) behaviour
               [x1 y1] = ind2sub([args.cols args.rows],args.NumericArguments{1}(1)); % subplot and ind2sub count differently (column instead of row first) --> switch cols/rows
               [x2 y2] = ind2sub([args.cols args.rows],args.NumericArguments{1}(end));
           else
               x1=mod((args.NumericArguments{1}-1),args.cols)+1; x2=x1;
               y1=floor((args.NumericArguments{1}-1)/args.cols)+1; y2=y1;
           end
    %       x1=mod((args.NumericArguments{1}-1),args.cols)+1; x2=x1;
    %       y1=floor((args.NumericArguments{1}-1)/args.cols)+1; y2=y1;
       case 2
          x1=args.NumericArguments{1};x2=x1;
          y1=args.NumericArguments{2};y2=y1;
       case 4
          x1=args.NumericArguments{1};x2=x1+args.NumericArguments{3}-1;
          y1=args.NumericArguments{2};y2=y1+args.NumericArguments{4}-1;
       otherwise
          error('getAxes argument error')
    end

    cellwidth=((1-args.MarginLeft-args.MarginRight)-(args.cols-1)*args.SpacingHorizontal)/args.cols;
    cellheight=((1-args.MarginTop-args.MarginBottom)-(args.rows-1)*args.SpacingVertical)/args.rows;
    xpos1=args.MarginLeft+args.PaddingLeft+cellwidth*(x1-1)+args.SpacingHorizontal*(x1-1);
    xpos2=args.MarginLeft-args.PaddingRight+cellwidth*x2+args.SpacingHorizontal*(x2-1);
    ypos1=args.MarginTop+args.PaddingTop+cellheight*(y1-1)+args.SpacingVertical*(y1-1);
    ypos2=args.MarginTop-args.PaddingBottom+cellheight*y2+args.SpacingVertical*(y2-1);
    if args.Holdaxis
        h=axes('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
    else
        h=subplot('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
    end
    set(h,'box','on');
    %h=axes('position',[x1 1-y2 x2-x1 y2-y1]);
    set(h,'units',get(gcf,'defaultaxesunits'));
    set(h,'tag','getAxes');
    if (nargout==0), clear h; end;
end

%>  \cond excluded

function argStruct = parseArgs(args, argStruct, varargin)
    % Helper function for parsing varargin.
    %
    % argStruct=parseArgs(varargin,argStruct[,FlagtypeParams[,Aliases]])
    %
    % * argStruct is the structure full of named arguments with default values.
    % * Flagtype params is params that don't require a value. (the value will be set to 1 if it is present)
    % * Aliases can be used to map one argument-name to several argstruct fields
    %
    % example usage:
    % --------------
    % function parseargtest(varargin)
    %
    % %define the acceptable named arguments and assign default values
    % args=struct('Holdaxis',0, ...
    %        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
    %        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
    %        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
    %        'rows',[],'cols',[]);
    %
    % %The capital letters define abrreviations.
    % %  Eg. parseargtest('spacingvertical',0) is equivalent to  parseargtest('sv',0)
    %
    % args=parseArgs(varargin,args, ... % fill the arg-struct with values entered by the user
    %           {'Holdaxis'}, ... %this argument has no value (flag-type)
    %           {'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
    %
    % disp(args)
    %
    %
    %
    %
    % Aslak Grinsted 2004
    % -------------------------------------------------------------------------
    %   Copyright (C) 2002-2004, Aslak Grinsted
    %   This software may be used, copied, or redistributed as long as it is not
    %   sold and this copyright notice is reproduced on each copy made.  This
    %   routine is provided as is without any express or implied warranties
    %   whatsoever.
    persistent matlabver
    if isempty(matlabver)
        matlabver=ver('MATLAB');
        matlabver=str2double(matlabver.Version);
    end
    Aliases={};
    FlagTypeParams='';
    if (length(varargin)>0)
        FlagTypeParams=lower(strvcat(varargin{1}));  %#ok
        if length(varargin)>1
            Aliases=varargin{2};
        end
    end

    %---------------Get "numeric" arguments
    NumArgCount=1;
    while (NumArgCount<=size(args,2))&&(~ischar(args{NumArgCount}))
        NumArgCount=NumArgCount+1;
    end
    NumArgCount=NumArgCount-1;
    if (NumArgCount>0)
        argStruct.NumericArguments={args{1:NumArgCount}};
    else
        argStruct.NumericArguments={};
    end
    %--------------Make an accepted fieldname matrix (case insensitive)
    Fnames=fieldnames(argStruct);
    for i=1:length(Fnames)
        name=lower(Fnames{i,1});
        Fnames{i,2}=name; %col2=lower
        Fnames{i,3}=[name(Fnames{i,1}~=name) ' ']; %col3=abreviation letters (those that are uppercase in the argStruct) e.g. SpacingHoriz->sh
        %the space prevents strvcat from removing empty lines
        Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams)); %Does this parameter have a value?
    end
    FnamesFull=strvcat(Fnames{:,2}); %#ok
    FnamesAbbr=strvcat(Fnames{:,3}); %#ok
    if length(Aliases)>0
        for i=1:length(Aliases)
            name=lower(Aliases{i,1});
            FieldIdx=strmatch(name,FnamesAbbr,'exact'); %try abbreviations (must be exact)
            if isempty(FieldIdx)
                FieldIdx=strmatch(name,FnamesFull); %&??????? exact or not?
            end
            Aliases{i,2}=FieldIdx;
            Aliases{i,3}=[name(Aliases{i,1}~=name) ' ']; %the space prevents strvcat from removing empty lines
            Aliases{i,1}=name; %dont need the name in uppercase anymore for aliases
        end
        %Append aliases to the end of FnamesFull and FnamesAbbr
        FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); %#ok
        FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3})); %#ok
    end
    %--------------get parameters--------------------
    l=NumArgCount+1;
    while (l<=length(args))
        a=args{l};
        if ischar(a)
            paramHasValue=1; % assume that the parameter has is of type 'param',value
            a=lower(a);
            FieldIdx=strmatch(a,FnamesAbbr,'exact'); %try abbreviations (must be exact)
            if isempty(FieldIdx)
                FieldIdx=strmatch(a,FnamesFull);
            end
            if (length(FieldIdx)>1) %shortest fieldname should win
                [mx,mxi]=max(sum(FnamesFull(FieldIdx,:)==' ',2));%#ok
                FieldIdx=FieldIdx(mxi);
            end
            if FieldIdx>length(Fnames) %then it's an alias type.
                FieldIdx=Aliases{FieldIdx-length(Fnames),2};
            end

            if isempty(FieldIdx)
                error(['Unknown named parameter: ' a])
            end
            for curField=FieldIdx' %if it is an alias it could be more than one.
                if (Fnames{curField,4})
                    if (l+1>length(args))
                        error(['Expected a value for parameter: ' Fnames{curField,1}])
                    end
                    val=args{l+1};
                else %FLAG PARAMETER
                    if (l<length(args)) %there might be a explicitly specified value for the flag
                        val=args{l+1};
                        if isnumeric(val)
                            if (numel(val)==1)
                                val=logical(val);
                            else
                                error(['Invalid value for flag-parameter: ' Fnames{curField,1}])
                            end
                        else
                            val=true;
                            paramHasValue=0;
                        end
                    else
                        val=true;
                        paramHasValue=0;
                    end
                end
                if matlabver>=6
                    argStruct.(Fnames{curField,1})=val; %try the line below if you get an error here
                else
                    argStruct=setfield(argStruct,Fnames{curField,1},val); %#ok <-works in old matlab versions
                end
            end
            l=l+1+paramHasValue; %if a wildcard matches more than one
        else
            error(['Expected a named parameter: ' num2str(a)])
        end
    end
end

%>  \endcond