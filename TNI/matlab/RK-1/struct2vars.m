function varargout = struct2vars(varargin)
%STRUCT2VARS Convert scalar structure into variables.
%   STRUCT2VARS(S) unpacks the input structure, instantiating a variable in
%   the caller workspace for each fieldname and its corresponding value. 
%
%   [A,B,C,___] = STRUCT2VARS(S) unpacks each field of the input structure
%   into the specified outputs, i.e. the first field of S will be stored in
%   output variable A, the second field of S will be stored in output
%   variable B, and so forth. The number of outputs must be less than or
%   equal to the number of fieldnames in the structure.
%
%   [___] = STRUCT2VARS(___,NAMES) allows the user to specify which
%   fieldnames are unpacked. The order provided by the user is the order in
%   which fields are unpacked.
%
%   Notes:
%   1) This code is a modified version of the following FEX submission:
%   http://www.mathworks.com/matlabcentral/fileexchange/31532-pack---unpack-variables-to---from-structures-with-enhanced-functionality
%
%   2) For the case when no output arguments are supplied: if a variable 
%   name already exists in the caller workspace, it will be overwritten by 
%   this function.
%
%   See also VARS2STRUCT.
%   Copyright 2016 Matthew R. Eicholtz
% Parse inputs
[s,names] = parseinputs(varargin{:});
n = length(names);
% Unpack structure
nargoutchk(0,n);
if nargout==0 %assign in caller
    % Check to see if any variables are being overwritten
    caller = evalin('caller','whos');
    mask = ismember(names,{caller(:).name}); %overlap with existing variables
    if any(mask)
        str = sprintf('\t%s\n',names{mask});
        warning('The following variables already exist in the caller workspace and will be overwritten:\n%s',str);
    end
    for ii=1:n
        assignin('caller',names{ii},s.(names{ii}));
    end
else %dump into variables
    c = cell(size(names));
    for ii=1:n
        c{ii} = s.(names{ii});
    end
    varargout = c;
end
end
%% Helper functions
function varargout = parseinputs(varargin)
%PARSEINPUTS Custom input parsing function.
    p = inputParser;
    
    p.addRequired('s',@(x) isstruct(x) & isscalar(x));
    p.addOptional('f','all',@iscell);
    
    p.parse(varargin{:});
    
    % How ironic that I would want to use this function right here!
    s = p.Results.s;
    f = p.Results.f;
    
    if strcmp(f,'all')
        f = fieldnames(s);
    end
    varargout = {s,f};
end
