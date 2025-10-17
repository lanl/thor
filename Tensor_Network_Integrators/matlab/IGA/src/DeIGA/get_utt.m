% function utt = get_utt(bc_type, u_in, u_out, u1tt)
% %GET_UTT   Apply 3D Dirichlet BC in TT‐format without ever forming a full array.
% %  bc_type = '1' → utt(1,:,:)   = ua; utt(end,:,:)   = ub
% %  bc_type = '2' → utt(:,1,:)   = ua; utt(:,end,:)   = ub
% %  bc_type = '3' → utt(:,:,1)   = ua; utt(:,:,end)   = ub
% 
%   % 1) unpack the three cores of u1tt
%   n  = u1tt.n;    % [n1; n2; n3]
%   r  = u1tt.r;    % [1; r2; r3; 1]
%   ps = u1tt.ps;
%   C  = u1tt.core;
% 
%   G1 = reshape(C(ps(1):ps(2)-1), [1,    n(1), r(2)]);
%   G2 = reshape(C(ps(2):ps(3)-1), [r(2), n(2), r(3)]);
%   G3 = reshape(C(ps(3):ps(4)-1), [r(3), n(3), 1   ]);
% 
%   % 2) splice in the two boundary channels
%   switch bc_type
%     case "1"  % boundary on dim‐1
%       % new sizes
%       n1 = n(1)+2;  n2 = n(2);    n3 = n(3);
%       R2 = r(2);    R3 = r(3);
%       % core‐1: 1×n1×(R2+2)
%       G1n = zeros(1,n1,R2+2);
%       G1n(1,2:end-1,1:R2) = G1;      % interior
%       G1n(1,1,   R2+1)    = u_in;      % face i1=1
%       G1n(1,end, R2+2)    = u_out;      % face i1=n1
%       % core‐2: (R2+2)×n2×(R3+2)
%       G2n = zeros(R2+2,n2,R3+2);
%       G2n(1:R2,:,1:R3)   = G2;        % interior
%       G2n(R2+1,:,R3+1)   = 1;         % propagate ua
%       G2n(R2+2,:,R3+2)   = 1;         % propagate ub
%       % core‐3: (R3+2)×n3×1
%       G3n = zeros(R3+2,n3,1);
%       G3n(1:R3,:,1)      = G3;        % interior
%       G3n(R3+1,:,1)      = 1;         % emit ua
%       G3n(R3+2,:,1)      = 1;         % emit ub
% 
%     case "2"  % boundary on dim‐2
%       n1 = n(1);    n2 = n(2)+2;  n3 = n(3);
%       R2 = r(2);    R3 = r(3);
%       % core‐1: 1×n1×(R2+2)
%       G1n = cat(3, ...
%         G1, ...
%         ones(1,n1,1), ...
%         ones(1,n1,1) );
%       % core‐2: (R2+2)×n2×(R3+2)
%       G2n = zeros(R2+2,n2,R3+2);
%       G2n(   1:R2, 2:end-1, 1:R3) = G2;   % interior
%       G2n(R2+1, 1,       R3+1)   = u_in;    % pick up ua
%       G2n(R2+2, end,     R3+2)   = u_out;    % pick up ub
%       % core‐3: (R3+2)×n3×1
%       G3n = cat(1, ...
%         G3, ...
%         ones(1,n3,1), ...
%         ones(1,n3,1) );
% 
%     case "3"  % boundary on dim‐3
%       n1 = n(1);    n2 = n(2);    n3 = n(3)+2;
%       R2 = r(2);    R3 = r(3);
%       % core‐1: 1×n1×(R2+2)
%       G1n = cat(3, ...
%         G1, ...
%         ones(1,n1,1), ...
%         ones(1,n1,1) );
%       % core‐2: (R2+2)×n2×(R3+2)
%       G2n = zeros(R2+2,n2,R3+2);
%       G2n(1:R2,:,1:R3) = G2;            % interior
%       G2n(R2+1,:,R3+1)=1;               % pass ua
%       G2n(R2+2,:,R3+2)=1;               % pass ub
%       % core‐3: (R3+2)×n3×1
%       G3n = zeros(R3+2,n3,1);
%       G3n(1:R3,    2:end-1,1) = G3;     % interior
%       G3n(R3+1,    1,      1) = u_in;     % face i3=1
%       G3n(R3+2,    end,    1) = u_out;     % face i3=n3
% 
%     otherwise
%       error('bc_type must be ''1'', ''2'' or ''3''.');
%   end
% 
%   % 3) pack into a struct and call the TT‐constructor
%   n_new    = [ size(G1n,2); size(G2n,2); size(G3n,2) ];
%   r_new    = [ 1;           size(G1n,3); size(G2n,3); 1 ];
%   core_vec = [G1n(:); G2n(:); G3n(:)];
% 
%   S.d    = 3;
%   S.n    = n_new;
%   S.r    = r_new;
%   S.core = core_vec;
%   S.over = 0;         % no overestimate
%   utt    = tt_tensor(S);  % <— uses the “struct” branch in tt_tensor
% 
%   % optionally round to clean up any tiny numerical noise
%   utt = round(utt, 1e-12);
% end
function utt = get_utt(bc_type, u_in, u_out, u1tt)
%GET_UTT   Apply 3D Dirichlet BC in TT‐format without ever forming a full array.
%  bc_type = '1' → utt(1,:,:)   = u_in;  utt(end,:,:)   = u_out
%  bc_type = '2' → utt(:,1,:)   = u_in;  utt(:,end,:)   = u_out
%  bc_type = '3' → utt(:,:,1)   = u_in;  utt(:,:,end)   = u_out

  % unpack the three cores of u1tt
  n  = u1tt.n;    % [n1; n2; n3]
  r  = u1tt.r;    % [1; r2; r3; 1]
  ps = u1tt.ps;
  C  = u1tt.core;

  G1 = reshape(C(ps(1):ps(2)-1), [1,    n(1), r(2)]);
  G2 = reshape(C(ps(2):ps(3)-1), [r(2), n(2), r(3)]);
  G3 = reshape(C(ps(3):ps(4)-1), [r(3), n(3), 1   ]);

  % splice in the two boundary channels
  switch bc_type
    case "1"  % boundary on dim-1
      % new sizes
      n1 = n(1)+2;  n2 = n(2);    n3 = n(3);
      R2 = r(2);    R3 = r(3);

      % core-1: 1×n1×(R2+2)
      G1n = zeros(1, n1, R2+2);
      G1n(1,2:end-1,1:R2) = G1;         % interior
      G1n(1,1,   R2+1)    = u_in;       % face i1=1
      G1n(1,end, R2+2)    = u_out;      % face i1=n1

      % core-2: (R2+2)×n2×(R3+2)
      G2n = zeros(R2+2, n2, R3+2);
      G2n(1:R2,       :, 1:R3)   = G2;  % interior
      G2n(R2+1, :,    R3+1)       = 1;   % propagate u_in
      G2n(R2+2, :,    R3+2)       = 1;   % propagate u_out

      % core-3: (R3+2)×n3×1
      G3n = zeros(R3+2, n3, 1);
      G3n(1:R3, :, 1)      = G3;        % interior
      G3n(R3+1, :, 1)      = 1;         % emit u_in
      G3n(R3+2, :, 1)      = 1;         % emit u_out

    case "2"  % boundary on dim-2
      % new sizes
      n1 = n(1);    n2 = n(2)+2;  n3 = n(3);
      R2 = r(2);    R3 = r(3);

      % core-1: 1×n1×(R2+2)
      G1n = zeros(1, n1, R2+2);
      G1n(1, :,       1:R2) = G1;       % interior
      G1n(1, :,       R2+1) = 1;        % propagate u_in
      G1n(1, :,       R2+2) = 1;        % propagate u_out

      % core-2: (R2+2)×n2×(R3+2)
      G2n = zeros(R2+2, n2, R3+2);
      G2n(1:R2, 2:end-1, 1:R3) = G2;     % interior
      G2n(R2+1, 1,       R3+1) = u_in;   % face j=1
      G2n(R2+2, end,     R3+2) = u_out;  % face j=n2

      % core-3: (R3+2)×n3×1
      G3n = zeros(R3+2, n3, 1);
      G3n(1:R3, :, 1)      = G3;         % interior
      G3n(R3+1, :, 1)      = 1;          % emit u_in
      G3n(R3+2, :, 1)      = 1;          % emit u_out

    case "3"  % boundary on dim-3
      % new sizes
      n1 = n(1);    n2 = n(2);    n3 = n(3)+2;
      R2 = r(2);    R3 = r(3);

      % core-1: 1×n1×(R2+2)
      G1n = zeros(1, n1, R2+2);
      G1n(1, :, 1:R2) = G1;            % interior
      G1n(1, :, R2+1) = 1;             % propagate u_in
      G1n(1, :, R2+2) = 1;             % propagate u_out

      % core-2: (R2+2)×n2×(R3+2)
      G2n = zeros(R2+2, n2, R3+2);
      G2n(1:R2, :, 1:R3) = G2;         % interior
      G2n(R2+1, :, R3+1) = 1;          % pass u_in
      G2n(R2+2, :, R3+2) = 1;          % pass u_out

      % core-3: (R3+2)×n3×1
      G3n = zeros(R3+2, n3, 1);
      G3n(1:R3,    2:end-1, 1) = G3;    % interior
      G3n(R3+1,    1,       1) = u_in;  % face k=1
      G3n(R3+2,    end,     1) = u_out; % face k=n3

    otherwise
      error('bc_type must be ''1'', ''2'', or ''3''.');
  end

  % pack into struct and construct the new TT tensor
  n_new    = [ size(G1n,2); size(G2n,2); size(G3n,2) ];
  r_new    = [ 1;           size(G1n,3); size(G2n,3); 1 ];
  core_vec = [G1n(:); G2n(:); G3n(:)];

  S.d    = 3;
  S.n    = n_new;
  S.r    = r_new;
  S.core = core_vec;
  S.over = 0;            % no overestimate
  utt    = tt_tensor(S); % construct via struct

  % optionally round to clean up numerical noise
  utt = round(utt, 1e-12);
end
