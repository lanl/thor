function dF = fg_oderhs(t,F,a,b,Vfun)
[N,~]= size(F);
dF1 = zeros(N); dF2=dF1;

dF1(2:N-1,:) = a*F(3:N,:) + b*F(2:N-1,:) + a*F(1:N-2,:);
dF1(1,:) = b*F(1,:) + a*F(2,:);
dF1(N,:) = b*F(N,:) + a*F(N-1,:);

dF2(:,2:N-1) = a*F(:,3:N) + b*F(:,2:N-1) + a*F(:,1:N-2);
dF2(:,1) = b*F(:,1) + a*F(:,2);
dF2(:,N) = b*F(:,N) + a*F(:,N-1);

if t>=5 && t<=15
    Vtt = Vfun{2}; 
else
    Vtt = Vfun{1};
end

dF = dF1 + dF2 + Vtt;
end