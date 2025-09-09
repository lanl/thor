function dFtt = tt_oderhs(t,Ftt,a,b,Vfun,et)
% right hand side for ttRK
s = Ftt.n;
N = s(1);
G = core2cell(Ftt); G1 = G; G2 = G;
G1{1}(:,2:N-1,:) = a*G{1}(:,3:N,:) + b*G{1}(:,2:N-1,:) + a*G{1}(:,1:N-2,:);
G1{1}(:,1,:) = b*G{1}(:,1,:)+ a*G{1}(:,2,:);
G1{1}(:,N,:) = b*G{1}(:,N,:)+ a*G{1}(:,N-1,:);
dFtt1 = cell2core(tt_tensor,G1);

N = s(2);
G2{2}(:,2:N-1,:) = a*G{2}(:,3:N,:) + b*G{2}(:,2:N-1,:) + a*G{2}(:,1:N-2,:);
G2{2}(:,1,:) = b*G{2}(:,1,:)+ a*G{2}(:,2,:);
G2{2}(:,N,:) = b*G{2}(:,N,:)+ a*G{2}(:,N-1,:);
dFtt2 = cell2core(tt_tensor,G2);

if t>=5 && t<=15
    Vtt = Vfun{2}; 
else
    Vtt = Vfun{1};
end

dFtt = round(dFtt1 + dFtt2 + Vtt,et);
end