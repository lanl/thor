function L2 = CalH1L2_Lshape_cross(utt, knot1n, knot2n, knot3n, ...
                                   ctpxn, ctpyn, ctpzn, ...
                                   tt_tol)

  % Gauss data
  [Gp1, ~, ~] = getIntegralData1D(knot1n);
  [Gp2, ~, ~] = getIntegralData1D(knot2n);
  [Gp3, ~, ~] = getIntegralData1D(knot3n);

  nn1 = length(Gp1);
  nn2 = length(Gp2);
  nn3 = length(Gp3);

  % --- define element evaluator ---
  function [nval, dval, err] = elem(i1, i2, i3)
    [nval, dval, err] = e_elem(utt, i1, i2, i3,...
                               knot1n, knot2n, knot3n,...
                               ctpxn, ctpyn, ctpzn,...
                               Gp1, Gp2, Gp3);
  end

  % --- build two TT‐tensors via amen_cross ---
  L2n = amen_cross([nn1,nn2,nn3], ...
                   @(idx) cross_fun_nD(idx, @elem, 1), ...
                   tt_tol);
  L2d = amen_cross([nn1,nn2,nn3], ...
                   @(idx) cross_fun_nD(idx, @elem, 2), ...
                   tt_tol);
  % note: cross_fun_nD calls elem with 3 subscripts and picks 
  % output #1 or #2 depending on the 3rd argument.

  % global L2
  Nfull = full(L2n);
  Dfull = full(L2d);

  num = sum(Nfull, 'all');
  den = sum(Dfull, 'all');
  L2 = num / den;
  
%   du = zeros(nn1,nn2,nn3);
%   for i1 = 1:nn1
%     for i2 = 1:nn2
%       for i3 = 1:nn3
%         [~,~,du(i1,i2,i3)] = elem(i1,i2,i3);
%       end
%     end
%   end

end

%% subroutines
function [L2n, L2d, du] = e_elem(utt, i1, i2, i3,...
                                 knot1n, knot2n, knot3n,...
                                 ctpxn, ctpyn, ctpzn,...
                                 Gp1, Gp2, Gp3)

[NN, DBS_Dz1, DBS_Dz2, DBS_Dz3, Inon1, Inon2, Inon3] = ...
                 GetBspline3Dn(knot1n, knot2n, knot3n, ...
                               Gp1(i1), Gp2(i2), Gp3(i3));

ctpx = ctpxn(Inon1, Inon2, Inon3);
ctpy = ctpyn(Inon1, Inon2, Inon3);
ctpz = ctpzn(Inon1, Inon2, Inon3);
ul = full(utt(Inon1, Inon2, Inon3));
%%
xe = dot(ctpx(:), NN(:));
ye = dot(ctpy(:), NN(:));
ze = dot(ctpz(:), NN(:));
ue = dot(ul(:), NN(:));

%%
DX1_Dzeta1 = ctpx(:)' * DBS_Dz1(:);
DX2_Dzeta1 = ctpy(:)' * DBS_Dz1(:);
DX3_Dzeta1 = ctpz(:)' * DBS_Dz1(:);

DX1_Dzeta2 = ctpx(:)' * DBS_Dz2(:);
DX2_Dzeta2 = ctpy(:)' * DBS_Dz2(:);
DX3_Dzeta2 = ctpz(:)' * DBS_Dz2(:);

DX1_Dzeta3 = ctpx(:)' * DBS_Dz3(:);
DX2_Dzeta3 = ctpy(:)' * DBS_Dz3(:);
DX3_Dzeta3 = ctpz(:)' * DBS_Dz3(:);
%%
DX1Dz = [DX1_Dzeta1, DX1_Dzeta2, DX1_Dzeta3];
DX2Dz = [DX2_Dzeta1, DX2_Dzeta2, DX2_Dzeta3];
DX3Dz = [DX3_Dzeta1, DX3_Dzeta2, DX3_Dzeta3];

%% Compute Jacobian determinant |J|
j = [DX1Dz; DX2Dz; DX3Dz];

%Analytical solution
ua = (1/(2*pi^2))*sin(pi*xe)*sin(pi*ye);
du = abs(ue - ua);

L2n = abs(det(j))*sqrt(du*du);
L2d = abs(det(j))*sqrt(abs(ua)); 

end
