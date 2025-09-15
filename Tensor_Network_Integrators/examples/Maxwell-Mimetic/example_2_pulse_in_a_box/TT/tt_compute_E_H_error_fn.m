function [err_E, nrm_E, err_H, nrm_H] = tt_compute_E_H_error_fn(D,R)

% compute error in E
%Ex
Il = 1:2^(log2(R.nx-1)-log2(D.nx-1)):R.nx;
Is = 1:2^(log2(R.nx-1)-log2(D.nx-1)):R.nx-1;
%sampling
Sx = tt_sampling_fn(R.tEx,{Is,Il,Il});
Sx = tt_remove_boundary_fn(Sx,[2,3]);
Sy = tt_sampling_fn(R.tEy,{Il,Is,Il});
Sy = tt_remove_boundary_fn(Sy,[1,3]);
Sz = tt_sampling_fn(R.tEz,{Il,Il,Is});
Sz = tt_remove_boundary_fn(Sz,[1,2]);
% [err_E, nrm_E] = error_diagnostics_E(D.tEx,D.tEy,D.tEz,Rsolx,Rsoly,Rsolz,nx,nx,nx);
%
D.tEx = tt_remove_boundary_fn(D.tEx,[2,3]);
err_Ex = sqrt(mean(  (full(D.tEx)-full(Sx)).^2  ));
nrm_Ex = sqrt(mean(  (full(Sx)).^2  ));
% err_Ex1 = sqrt(sum(round(D.tEx-Rsolx,tol).^2)/((nx-1)*nx*nx));
% nrm_Ex = sqrt(sum(round((Rsolx).^2,tol))/((nx-1)*nx*nx));
% %Ey
D.tEy = tt_remove_boundary_fn(D.tEy,[1,3]);
err_Ey = sqrt(mean(  (full(D.tEy)-full(Sy)).^2  ));
nrm_Ey = sqrt(mean(  (full(Sy)).^2  ));
% err_Ey = sqrt(sum(round(D.tEy-Rsoly,tol).^2)/((nx-1)*nx*nx));
% nrm_Ey = sqrt(sum((Rsoly).^2)/((nx-1)*nx*nx));
% %Ez
D.tEz = tt_remove_boundary_fn(D.tEz,[1,2]);
err_Ez = sqrt(mean(  (full(D.tEz)-full(Sz)).^2  ));
nrm_Ez = sqrt(mean(  (full(Sz)).^2  ));
% err_Ez = sqrt(sum(round(D.tEz-Rsolz,tol).^2)/((nx-1)*nx*nx));
% nrm_Ez = sqrt(sum((Rsolz).^2)/((nx-1)*nx*nx));
%final
err_E = err_Ex + err_Ey + err_Ez;
nrm_E = nrm_Ex + nrm_Ey + nrm_Ez;
% err_E/nrm_E

%compute error in H

%sampling
Sx = tt_sampling_fn(R.tHx,{Il,Is,Is});
Sx = tt_remove_boundary_fn(Sx,[1]);

Sy = tt_sampling_fn(R.tHy,{Is,Il,Is});
Sy = tt_remove_boundary_fn(Sy,[2]);

Sz = tt_sampling_fn(R.tHz,{Is,Is,Il});
Sz = tt_remove_boundary_fn(Sz,[3]);
% [err_E, nrm_E] = error_diagnostics_E(D.tEx,D.tEy,D.tEz,Rsolx,Rsoly,Rsolz,nx,nx,nx);
%
D.tHx = tt_remove_boundary_fn(D.tHx,[1]);
err_Hx = sqrt(mean(  (full(D.tHx)-full(Sx)).^2  ));
nrm_Hx = sqrt(mean(  (full(Sx)).^2  ));

% err_Hx = sqrt(sum(round(D.tHx-Rsolx,tol).^2)/((nx-1)^2*nx));
% nrm_Hx = sqrt(sum(round((Rsolx).^2,tol))/((nx-1)^2*nx));
% %Ey
D.tHy = tt_remove_boundary_fn(D.tHy,[2]);
err_Hy = sqrt(mean(  (full(D.tHy)-full(Sy)).^2  ));
nrm_Hy = sqrt(mean(  (full(Sy)).^2  ));

% err_Hy = sqrt(sum(round(D.tHy-Rsoly,tol).^2)/((nx-1)^2*nx));
% nrm_Hy = sqrt(sum((Rsoly).^2)/((nx-1)^2*nx));
% Ez
D.tHz = tt_remove_boundary_fn(D.tHz,[3]);
err_Hz = sqrt(mean(  (full(D.tHz)-full(Sz)).^2  ));
nrm_Hz = sqrt(mean(  (full(Sz)).^2  ));

% err_Hz = sqrt(sum(round(D.tHz-Rsolz,tol).^2)/((nx-1)^2*nx));
% nrm_Hz = sqrt(sum((Rsolz).^2)/((nx-1)^2*nx));
%final
err_H = err_Hx + err_Hy + err_Hz;
nrm_H = nrm_Hx + nrm_Hy + nrm_Hz;
% err_H/nrm_H
fprintf('Err E + H = %.5f \n',err_E + err_H)

