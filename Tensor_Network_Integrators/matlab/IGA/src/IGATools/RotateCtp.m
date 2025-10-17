function [Xr, Yr, Zr] = RotateCtp(ctpx, ctpy, ctpz, angle_deg, axis)
% RotateCtp - Rotates 3D grid points around a specified axis
%
% Inputs:
%   ctpx, ctpy, ctpz : [n1 x n2 x n3] position arrays
%   angle_deg        : Rotation angle in degrees
%   axis             : Rotation axis: 'x', 'y', or 'z'
%
% Outputs:
%   Xr, Yr, Zr       : Rotated coordinate arrays (same size)

% Convert angle to radians
theta = deg2rad(angle_deg);

% Define rotation matrix
switch lower(axis)
    case 'x'
        R = [1 0 0;
             0 cos(theta) -sin(theta);
             0 sin(theta)  cos(theta)];
    case 'y'
        R = [ cos(theta) 0 sin(theta);
              0          1 0;
             -sin(theta) 0 cos(theta)];
    case 'z'
        R = [cos(theta) -sin(theta) 0;
             sin(theta)  cos(theta) 0;
             0           0          1];
    otherwise
        error('Axis must be ''x'', ''y'', or ''z''.');
end

% Flatten the coordinate arrays
X = ctpx(:);
Y = ctpy(:);
Z = ctpz(:);

% Apply rotation
rotated = R * [X'; Y'; Z'];

% Reshape back to original size
sz = size(ctpx);
Xr = reshape(rotated(1, :), sz);
Yr = reshape(rotated(2, :), sz);
Zr = reshape(rotated(3, :), sz);

end
