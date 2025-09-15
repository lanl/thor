function data = prepareOutput(data)
    % return only the interior cell results
    data.Q = data.Q(:,4:end-3,4:end-3,4:end-3);  
    data.X = data.X(4:end-3,4:end-3,4:end-3);  
    data.Y = data.Y(4:end-3,4:end-3,4:end-3);  
    data.Z = data.Z(4:end-3,4:end-3,4:end-3);  
end
