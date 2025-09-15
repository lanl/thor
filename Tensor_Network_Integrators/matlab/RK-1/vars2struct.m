function s = vars2struct(varargin)
    s = struct();
    for i = 1:length(varargin)
        varName = varargin{i};
        if evalin('base', sprintf("exist('%s', 'var')", varName)) % Check if variable exists
            s.(varName) = evalin('base', varName); % Assign variable to struct
        else
            warning('Variable "%s" does not exist in the workspace.', varName);
        end
    end
end