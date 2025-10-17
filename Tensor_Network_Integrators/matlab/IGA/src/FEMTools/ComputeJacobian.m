function J = ComputeJacobian(ctpx, ctpy, ctpz, DBS_Dz1, DBS_Dz2, DBS_Dz3)
DX1_Dzeta1 = ctpx(:)' * DBS_Dz1(:);
DX2_Dzeta1 = ctpy(:)' * DBS_Dz1(:);
DX3_Dzeta1 = ctpz(:)' * DBS_Dz1(:);

DX1_Dzeta2 = ctpx(:)' * DBS_Dz2(:);
DX2_Dzeta2 = ctpy(:)' * DBS_Dz2(:);
DX3_Dzeta2 = ctpz(:)' * DBS_Dz2(:);

DX1_Dzeta3 = ctpx(:)' * DBS_Dz3(:);
DX2_Dzeta3 = ctpy(:)' * DBS_Dz3(:);
DX3_Dzeta3 = ctpz(:)' * DBS_Dz3(:);

DX1Dz = [DX1_Dzeta1, DX1_Dzeta2, DX1_Dzeta3];
DX2Dz = [DX2_Dzeta1, DX2_Dzeta2, DX2_Dzeta3];
DX3Dz = [DX3_Dzeta1, DX3_Dzeta2, DX3_Dzeta3];

%% Compute Jacobian determinant |J|
J = [DX1Dz; DX2Dz; DX3Dz];