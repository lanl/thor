function [A_I,A_G] = Space_FDFD_MAT(SN)
%  Compute finite difference matrices for time and space 
% variable
% SN number of pertition in space direction


%%%%  Matrix Corresponding to Space Variable

A_I=sparse(SN-1,SN-1);

for tt=1:SN-1
    if tt == 1
       A_I(tt,tt)=2;
       A_I(tt,tt+1)=-1;
    elseif tt==SN-1
        
        A_I(tt,tt)=2;
        A_I(tt,tt-1)=-1;
    else
        A_I(tt,tt)=2;
        A_I(tt,tt-1)=-1;
        A_I(tt,tt+1)=-1;
    end
end

A_G=sparse(SN+1,SN+1);

for tt=1:SN+1
    if tt == 1
       A_G(tt,tt)=2;
       A_G(tt,tt+1)=-1;
    elseif tt==SN+1
        A_G(tt,tt)=2;
        A_G(tt,tt-1)=-1;
    else
        A_G(tt,tt)=2;
        A_G(tt,tt-1)=-1;
        A_G(tt,tt+1)=-1;
    end
end

end

