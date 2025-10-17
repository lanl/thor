function ftt = Decomposef_3D_tt_cross(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn, ffunction, tt_tol)
%Gauss data
[Gp1, ~, ~] = getIntegralData1D(knot1n);
[Gp2, ~, ~] = getIntegralData1D(knot2n);
[Gp3, ~, ~] = getIntegralData1D(knot3n);

nn1 = length(Gp1);
nn2 = length(Gp2);
nn3 = length(Gp3);

% ffg = zeros(nn1, nn2, nn3);
% for i1 = 1:nn1
%   for i2 = 1:nn2
%     for i3 =1:nn3
% 
%       % ffg(i1,i2,i3) = f_elem(i1,i2,i3,ffunction,knot1n,knot2n,knot3n,ctpxn,ctpyn,ctpzn,Gp1,Gp2,Gp3);
%     end
%   end
% end

tempfun = @(i1,i2,i3) f_elem(i1,i2,i3,ffunction,knot1n,knot2n,knot3n,ctpxn,ctpyn,ctpzn,Gp1,Gp2,Gp3);
ftt = amen_cross([nn1,nn2,nn3],@(x) cross_fun_nD(x,tempfun,1),tt_tol);

% %% decompose omega
% ftt = tt_tensor(fJ, tt_tol);

end %function

%% subroutines
function fval = f_elem(i1,i2,i3,ffunction,knot1n,knot2n,knot3n,ctpxn,ctpyn,ctpzn,Gp1,Gp2,Gp3)

[NN, DBS_Dz1, DBS_Dz2, DBS_Dz3, Inon1, Inon2, Inon3] = GetBspline3Dn(knot1n, knot2n, knot3n, ...
  Gp1(i1), Gp2(i2), Gp3(i3));
ctpx = ctpxn(Inon1, Inon2, Inon3);
ctpy = ctpyn(Inon1, Inon2, Inon3);
ctpz = ctpzn(Inon1, Inon2, Inon3);
%%
x = dot(ctpx(:), NN(:));
y = dot(ctpy(:), NN(:));
z = dot(ctpz(:), NN(:));
if ffunction == "1"
  fe = 1;
elseif ffunction == "cosxsinx"
  fe = cos(x)*sin(x);
elseif ffunction == "cosxsiny"
  fe = cos(x)*sin(y);
elseif ffunction == "cosysinx"
  fe = cos(y)*sin(x);
elseif ffunction == "cosysiny"
  fe = cos(y)*sin(y);
  %%%
elseif ffunction == "cosxcosx"
  fe = cos(x)*cos(x);
elseif ffunction == "cosycosy"
  fe = cos(y)*cos(y);
elseif ffunction == "cosxcosy"
  fe = cos(x)*cos(y);
  %%%
elseif ffunction == "sinxsinx"
  fe = sin(x)*sin(x);
elseif ffunction == "sinysiny"
  fe = sin(y)*sin(y);
elseif ffunction == "sinxsiny"
  fe = sin(x)*sin(y);
  %%%
elseif ffunction == "sinxcosy"
  fe = sin(x)*cos(y);
elseif ffunction == "sinpxsinpy"
  fe = sin(pi*x)*sin(pi*y);
elseif ffunction == "sinpxsinpysinpz"
  fe = sin(pi*x)*sin(pi*y)*sin(pi*z);
end
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
fval = abs(det(j))*fe;

end
