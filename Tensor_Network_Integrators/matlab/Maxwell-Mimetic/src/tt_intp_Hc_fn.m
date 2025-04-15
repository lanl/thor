function tHxc = tt_intp_Hc_fn(tHx, idim, tol)
% compute the cell center value for H
G = core2cell(tHx);
G{idim} = 0.5*(G{idim}(:,2:end,:) + G{idim}(:,1:end-1,:));
tHxc = round(cell2core(tt_tensor,G),tol);

end



% 
% 
% % get cores
% n = tHx.n;
% r = tHx.r;
% G = cell(1,numel(n));
% for i = 1:numel(n)
%   G{i} = zeros(r(i),n(i)-(i==idim)*1 + 2, r(i+1));
%   % interpolate the mid-point
%   if i == idim
%     G{i}(:,2:end-1,:) = 0.5*(tHx{i}(:,2:end,:) + tHx{i}(:,1:end-1,:));
%   else
%     G{i}(:,2:end-1,:) = tHx{i};
%   end
% end
% 
% tHxc = round(cell2core(tt_tensor,G),tol);