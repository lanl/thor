function [A1,A2,M1,M2, B1, B2] = nonlinear_Mat(X)
% A is Nonlinear Stiffness Matrix
% M is Nonlinear Mass Matrix
% B is nonlinear matrix for the stiffness term
% x is X coordinate; y is Y coordinate
x1=X(1);x2=X(2); 
hx=x2-x1;
x_md=(x1+x2)/2;
% Basis function Phi1_x=(x2-x)/hx; Phi2_x=(x-x1)/hx;
% Phi1_y=(y2-y)/hy; Phi2_y=(y-y1)/hy;
% [Phi1_x*Phi1_y; Phi2_x*Phi1_y; Phi1_x*Phi2_y; Phi2_x*Phi2_y]

Phix1 = @(x) PPhix1(x,x2,hx);
Phix2 = @(x) PPhix2(x,x1,hx);

DPhi1=(-1/hx); 
DPhi2=(1/hx);
M1=zeros(2,2); 
M2=zeros(2,2); 
A1=zeros(2,2); 
A2=zeros(2,2);

% Construction of the Mass Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M1(1,1)=(hx/6)*(Phix1(x1)*Phix1(x1)*Phix1(x1) + ...
  4*Phix1(x_md)*Phix1(x_md)*Phix1(x_md)+...
  Phix1(x2)*Phix1(x2)*Phix1(x2));

M1(1,2)=(hx/6)*(Phix1(x1)*Phix1(x1)*Phix2(x1) + ...
  4*Phix1(x_md)*Phix1(x_md)*Phix2(x_md)+...
  Phix1(x2)*Phix1(x2)*Phix2(x2));


M1(2,1)=(hx/6)*(Phix1(x1)*Phix2(x1)*Phix1(x1) + ...
  4*Phix1(x_md)*Phix2(x_md)*Phix1(x_md)+...
  Phix1(x2)*Phix2(x2)*Phix1(x2));

M1(2,2)= (hx/6)*(Phix1(x1)*Phix2(x1)*Phix2(x1) + ...
  4*Phix1(x_md)*Phix2(x_md)*Phix2(x_md)+...
  Phix1(x2)*Phix2(x2)*Phix2(x2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M2(1,1)=(hx/6)*(Phix2(x1)*Phix1(x1)*Phix1(x1) + ...
  4*Phix2(x_md)*Phix1(x_md)*Phix1(x_md)+...
  Phix2(x2)*Phix1(x2)*Phix1(x2));

M2(1,2)=(hx/6)*(Phix2(x1)*Phix1(x1)*Phix2(x1) + ...
  4*Phix2(x_md)*Phix1(x_md)*Phix2(x_md)+...
  Phix2(x2)*Phix1(x2)*Phix2(x2));


M2(2,1)=(hx/6)*(Phix2(x1)*Phix2(x1)*Phix1(x1) + ...
  4*Phix2(x_md)*Phix2(x_md)*Phix1(x_md)+...
  Phix2(x2)*Phix2(x2)*Phix1(x2));

M2(2,2)= (hx/6)*(Phix2(x1)*Phix2(x1)*Phix2(x1) + ...
  4*Phix2(x_md)*Phix2(x_md)*Phix2(x_md)+...
  Phix2(x2)*Phix2(x2)*Phix2(x2));

% Stiffness Matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1(1,1)=(hx/6)*(Phix1(x1)*DPhi1*DPhi1+4*Phix1(x_md)*DPhi1*DPhi1+...
  Phix1(x2)*DPhi1*DPhi1);

A1(1,2)=(hx/6)*(Phix1(x1)*DPhi1*DPhi2+4*Phix1(x_md)*DPhi1*DPhi2+...
  Phix1(x2)*DPhi1*DPhi2);

A1(2,1)=(hx/6)*(Phix1(x1)*DPhi2*DPhi1+4*Phix1(x_md)*DPhi2*DPhi1+...
  Phix1(x2)*DPhi2*DPhi1);

A1(2,2)=(hx/6)*(Phix1(x1)*DPhi2*DPhi2+4*Phix1(x_md)*DPhi2*DPhi2+...
  Phix1(x2)*DPhi2*DPhi2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A2(1,1)=(hx/6)*(Phix2(x1)*DPhi1*DPhi1+4*Phix2(x_md)*DPhi1*DPhi1+...
  Phix2(x2)*DPhi1*DPhi1);

A2(1,2)=(hx/6)*(Phix2(x1)*DPhi1*DPhi2+4*Phix2(x_md)*DPhi1*DPhi2+...
  Phix2(x2)*DPhi1*DPhi2);

A2(2,1)=(hx/6)*(Phix2(x1)*DPhi2*DPhi1+4*Phix2(x_md)*DPhi2*DPhi1+...
  Phix2(x2)*DPhi2*DPhi1);

A2(2,2)=(hx/6)*(Phix2(x1)*DPhi2*DPhi2+4*Phix2(x_md)*DPhi2*DPhi2+...
  Phix2(x2)*DPhi2*DPhi2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix corresponding to the nonlinear term


B1(1,1)=(hx/6)*(Phix1(x1)*DPhi1*Phix1(x1) + ...
  4*Phix1(x_md)*DPhi1*Phix1(x_md)+...
  Phix1(x2)*DPhi1*Phix1(x2));

B1(1,2)=(hx/6)*(Phix1(x1)*DPhi2*Phix1(x1) + ...
  4*Phix1(x_md)*DPhi2*Phix1(x_md)+...
  Phix1(x2)*DPhi2*Phix1(x2));


B1(2,1)=(hx/6)*(Phix1(x1)*DPhi1*Phix2(x1) + ...
  4*Phix1(x_md)*DPhi1*Phix2(x_md)+...
  Phix1(x2)*DPhi1*Phix2(x2));

B1(2,2)= (hx/6)*(Phix1(x1)*DPhi2*Phix2(x1) + ...
  4*Phix1(x_md)*DPhi2*Phix2(x_md)+...
  Phix1(x2)*DPhi2*Phix2(x2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B2(1,1)=(hx/6)*(Phix2(x1)*DPhi1*Phix1(x1) + ...
  4*Phix2(x_md)*DPhi1*Phix1(x_md)+...
  Phix2(x2)*DPhi1*Phix1(x2));

B2(1,2)=(hx/6)*(Phix2(x1)*DPhi2*Phix1(x1) + ...
  4*Phix2(x_md)*DPhi2*Phix1(x_md)+...
  Phix2(x2)*DPhi2*Phix1(x2));


B2(2,1)=(hx/6)*(Phix2(x1)*DPhi1*Phix2(x1) + ...
  4*Phix2(x_md)*DPhi1*Phix2(x_md)+...
  Phix2(x2)*DPhi1*Phix2(x2));

B2(2,2)= (hx/6)*(Phix2(x1)*DPhi2*Phix2(x1) + ...
  4*Phix2(x_md)*DPhi2*Phix2(x_md)+...
  Phix2(x2)*DPhi2*Phix2(x2));


end

%%subroutines
function fy= PPhix1(x,x2,hx)
fy= (x2-x)/hx;
end

function fy= PPhix2(x,x1,hx)
fy= (x-x1)/hx;
end

