function tt = cross2d_vec(f, n, m, eps, varargin)
  %
  r0 = 2;
  rmax = max(n,m);
  for idx = 1:2:numel(varargin)-1
    switch lower(varargin{idx})
      case 'r0'
        r0 = varargin{idx+1};
      case 'rmax'
        rmax = varargin{idx+1};
      otherwise
        error('Unrecognized option: %s', varargin{idx});
    end
  end
  
  % Initial random bases
  u = randn(n, r0);
  v = randn(m, r0);
  er = 2 * eps;
  
  % Orthogonalize
  [u, ~] = qr(u, 0);
  [v, ~] = qr(v, 0);
  
  % Initial pivots via maxvol on u, v
  indu = maxvol2(u)';
  indv = maxvol2(v)';
  
  % Scale u, v to pivot identity
  uu = u / u(indu, :);
  vv = v / v(indv, :);
  
  indu_add = indu;
  indv_add = indv;
  
  % Build initial cross matrix Phi via indexing
  idx = flatten_index(indu,indv);
  Phi = reshape(f(idx),[length(indu),length(indv)]);
  
  % Iterative refinement
  while er > eps
    % Build uadd and its approximation
    idx       = flatten_index((1:n)',indv_add);
    uadd      = reshape(f(idx),[n,length(indv_add)]); % n x rv_add
    uadd_appr = uu * Phi * vv(indv_add,:).';
  
    % Build vadd and its approximation
    idx       = flatten_index(indu_add,(1:m)');
    vadd      = reshape(f(idx),[length(indu_add),m]).';               % m x ru_add
    vadd_appr = (uu(indu_add, :) * Phi * vv')';   % fixed to match m x ru_add
  
    % Compute error
    er1 = norm(uadd_appr - uadd, 'fro') / norm(uadd, 'fro');
    er2 = norm(vadd_appr - vadd, 'fro') / norm(vadd, 'fro');
    er  = max(er1, er2);
  
    % Update u basis
    comp_u      = uadd - uu * uadd(indu, :);
    [comp_u, ~] = qr(comp_u, 0);
    indu_add    = maxvol2(comp_u)';
    u2          = comp_u / comp_u(indu_add, :);
    u1          = uu - u2 * uu(indu_add, :);
    uu          = [u1, u2];
    fu          = max(abs(uu(:)));
    indu        = [indu; indu_add];
  
    % Update v basis
    comp_v      = vadd - vv * vadd(indv, :);
    [comp_v, ~] = qr(comp_v, 0);
    indv_add    = maxvol2(comp_v)';
    v2          = comp_v / comp_v(indv_add, :);
    v1          = vv - v2 * vv(indv_add, :);
    vv          = [v1, v2];
    fv          = max(abs(vv(:)));
    er          = er * max(fu, fv);
    indv        = [indv; indv_add];
    %
    % Rebuild Phi
    %
    idx = flatten_index(indu,indv);
    Phi = reshape(f(idx),[length(indu), length(indv)]);
    %
    if(length(indu)>rmax)
      break;
    end
    %
    %fprintf("Err = %e  rank = %d\n",er,length(indu));
    %
  end
  
  % Final low-rank factors
  u = uu * Phi;
  v = vv;
  %
  tt=cell2core(tt_tensor,{reshape(u,[1,size(u)]),v'});
  %
end
%
function idx = flatten_index(i,j)
  [i1,j1] = ndgrid(i,j);
  idx = [i1(:) j1(:)];
end