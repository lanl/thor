function [INX, m3D] = N3DIndex(knot1, knot2, knot3)

N3DDomain = GetNURBS3DDomain(knot1, knot2, knot3);
m3D       = Getmacro3D(knot1, knot2, knot3);

[NOE, ~, ~]  = size(m3D);
[L, ~, ~]    = size(N3DDomain);

[~, p1]   = EvaluateKnot(knot1);
[~, p2]   = EvaluateKnot(knot2);
[~, p3]   = EvaluateKnot(knot3);

nConn = (p1+1)*(p2+1)*(p3+1);
INX   = zeros(NOE, nConn);

disp('Compute connectivity');
tic
% For each element, vectorized domain overlap check
for k=1:NOE
    m = m3D(k, :); % [x1 x2 y1 y2 z1 z2] for the k-th element
    % Vectorized check for all L domains
    mask = ~(N3DDomain(:,2) <= m(1) | N3DDomain(:,1) >= m(2) | ...
             N3DDomain(:,4) <= m(3) | N3DDomain(:,3) >= m(4) | ...
             N3DDomain(:,6) <= m(5) | N3DDomain(:,5) >= m(6));
    % Get indices of intersecting domains
    idx = find(mask);
    % Only take up to nConn, pad with zeros if fewer
    nn = min(numel(idx), nConn);
    INX(k,1:nn) = idx(1:nn);
    % If fewer than nConn, leave the rest as zeros (MATLAB default)
end
toc
end
