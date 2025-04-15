function data = prepareOutput(data)
    % return only the interior cell results
    data.Q = data.Q(:,4:end-3,4:end-3);
    data.X = data.X(4:end-3,4:end-3);
    data.Y = data.Y(4:end-3,4:end-3);
    %
    % dimensionalize variables
    %
    data.Q(1,:,:) = data.Q(1,:,:)*data.Qref(1);
    data.Q(2,:,:) = data.Q(2,:,:)*data.Qref(2);
    data.Q(3,:,:) = data.Q(3,:,:)*data.Qref(3);
    %
    data.X(:) = data.X*data.ref.L;
    data.Y(:) = data.Y*data.ref.L;
    %
end