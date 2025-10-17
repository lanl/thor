function Kftt = GetKf_3D_tt_cross(knot1n, knot2n, knot3n, ctpxn, ctpyn, ctpzn,tt_tol)
%Gauss data
[Gp1, ~, ~] = getIntegralData1D(knot1n);
[Gp2, ~, ~] = getIntegralData1D(knot2n);
[Gp3, ~, ~] = getIntegralData1D(knot3n);

nn1 = length(Gp1);
nn2 = length(Gp2);
nn3 = length(Gp3);


% Kf = zeros(nn1, nn2, nn3, 3, 3);
% for i1 = 1:nn1
%   for i2 = 1:nn2
%     for i3 =1:nn3
%       Kf_elem = Kf_elem_fn(i1,i2,i3,ctpxn,ctpyn,ctpzn,knot1n,knot2n,knot3n, Gp1, Gp2, Gp3);
%       Kf(i1,i2,i3,:,:) = Kf_elem;
%     end
%   end
% end

%% tt-cross
tempfun = @(i1,i2,i3) Kf_elem_fn(i1,i2,i3,ctpxn,ctpyn,ctpzn,knot1n,knot2n,knot3n, Gp1, Gp2, Gp3);

Kftt = amen_cross([nn1,nn2,nn3],@(x) cross_fun_nD(x,tempfun,9),tt_tol,'verb',1);
%need to specify tt_tol for tt_reshape
Kftt = tt_reshape(Kftt,[nn1,nn2,nn3,3,3],tt_tol);
Kftt = round(Kftt,tt_tol);
end %function

