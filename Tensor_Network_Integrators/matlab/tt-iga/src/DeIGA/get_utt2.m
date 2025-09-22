function utt = get_utt2(bc_type, u_in, u_out, u1tt)
%GET_UTT   Apply 3D Dirichlet BC in TT‐format without ever forming a full array.
%  bc_type = '1'  → utt(1,:,:)   = u_in;  utt(end,:,:)   = u_out
%  bc_type = '2'  → utt(:,1,:)   = u_in;  utt(:,end,:)   = u_out
%  bc_type = '3'  → utt(:,:,1)   = u_in;  utt(:,:,end)   = u_out
%  bc_type = '12' → first apply '1', then apply '2'
%  bc_type = '13' → first apply '1', then apply '3'
%  bc_type = '23' → first apply '2', then apply '3'

  % 1) If asked for a 2‐D BC, apply the two single‐dim cases in sequence
  if strlength(bc_type) == 2
    % extract the two digits in the order given
    dims_char = char(bc_type);
    dimA = dims_char(1) - '0';
    dimB = dims_char(2) - '0';
    % sanity checks
    assert(ismember(dimA,1:3) && ismember(dimB,1:3) && dimA~=dimB, ...
      'bc_type must be one of ''12'',''13'',''23''.');
    % apply first BC to the original TT
    utt1 = get_utt(string(dimA), u_in, u_out, u1tt);
    % then apply second BC to the result
    utt  = get_utt(string(dimB), u_in, u_out, utt1);
    return
  end

  % 2) Unpack the three TT‐cores from u1tt
  n  = u1tt.n;    % [n1; n2; n3]
  r  = u1tt.r;    % [1; r2; r3; 1]
  ps = u1tt.ps;
  C  = u1tt.core;

  G1 = reshape(C(ps(1):ps(2)-1), [1,    n(1), r(2)]);
  G2 = reshape(C(ps(2):ps(3)-1), [r(2), n(2), r(3)]);
  G3 = reshape(C(ps(3):ps(4)-1), [r(3), n(3), 1   ]);

  % 3) Splice in boundary channels
  switch bc_type
    case "1"  % boundary on dim‐1
      n1 = n(1)+2;  n2 = n(2);    n3 = n(3);
      R2 = r(2);    R3 = r(3);
      G1n = zeros(1,n1,R2+2);
      G1n(1,2:end-1,1:R2) = G1;           
      G1n(1,1,   R2+1)    = u_in;         
      G1n(1,end, R2+2)    = u_out;        
      G2n = zeros(R2+2,n2,R3+2);
      G2n(1:R2,:,1:R3)   = G2;            
      G2n(R2+1,:,R3+1)   = 1;             
      G2n(R2+2,:,R3+2)   = 1;             
      G3n = zeros(R3+2,n3,1);
      G3n(1:R3,:,1)      = G3;            
      G3n(R3+1,:,1)      = 1;             
      G3n(R3+2,:,1)      = 1;             

    case "2"  % boundary on dim‐2
      n1 = n(1);    n2 = n(2)+2;  n3 = n(3);
      R2 = r(2);    R3 = r(3);
      G1n = cat(3, ...
        G1, ...
        u_in*ones(1,n1,1), ...
        u_out*ones(1,n1,1) );
      G2n = zeros(R2+2,n2,R3+2);
      G2n(1:R2,  2:end-1, 1:R3) = G2;     
      G2n(R2+1,  1,       R3+1)   = 1;    
      G2n(R2+2,  end,     R3+2)   = 1;    
      G3n = cat(1, ...
        G3, ...
        ones(1,n3,1), ...
        ones(1,n3,1) );

    case "3"  % boundary on dim‐3
      n1 = n(1);    n2 = n(2);    n3 = n(3)+2;
      R2 = r(2);    R3 = r(3);
      G1n = cat(3, ...
        G1, ...
        ones(1,n1,1), ...
        ones(1,n1,1) );
      G2n = zeros(R2+2,n2,R3+2);
      G2n(1:R2,:,1:R3) = G2;            
      G2n(R2+1,:,R3+1)=1;               
      G2n(R2+2,:,R3+2)=1;               
      G3n = zeros(R3+2,n3,1);
      G3n(1:R3,    2:end-1,1) = G3;     
      G3n(R3+1,    1,      1) = u_in;   
      G3n(R3+2,    end,    1) = u_out;  

    otherwise
      error('bc_type must be ''1'', ''2'', ''3'', ''12'', ''13'' or ''23''.');
  end

  % 4) Pack the new cores into a TT‐structure and construct
  n_new    = [ size(G1n,2); size(G2n,2); size(G3n,2) ];
  r_new    = [ 1;           size(G1n,3); size(G2n,3); 1 ];
  core_vec = [G1n(:); G2n(:); G3n(:)];

  S.d    = 3;
  S.n    = n_new;
  S.r    = r_new;
  S.core = core_vec;
  S.over = 0;         
  utt    = tt_tensor(S);  

  % 5) Clean up any numerical noise
  utt = round(utt, 1e-12);
end
