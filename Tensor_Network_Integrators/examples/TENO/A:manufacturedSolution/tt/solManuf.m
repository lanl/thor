function [S,Q] = solManuf(xx,i,t,gam)
   %
   % xx : xyz in tt format
   % i  : equation index
   % t  : current time 
   % gam: specific heat ratio
   %
   x=xx(:,1); y=xx(:,2); z=xx(:,3);
   %
   % simple case
   %
   rho_0 = 1.0;     U_0  = 1.0;     V_0  = 1.0;     W_0  = 1.0;     e_0  = 1.0; 
   A_1   = 0.1;     A_2  = 0.1;     A_3  = 0.1;     A_4  = 0.1;     A_5  = 0.1; 
   n_11  = 2.0*pi;  n_21 = 2.0*pi;  n_31 = 2.0*pi;  n_41 = 2.0*pi;  n_51 = 2.0*pi; 
   n_12  = 2.0*pi;  n_22 = 2.0*pi;  n_32 = 2.0*pi;  n_42 = 2.0*pi;  n_52 = 2.0*pi; 
   n_13  = 2.0*pi;  n_23 = 2.0*pi;  n_33 = 2.0*pi;  n_43 = 2.0*pi;  n_53 = 2.0*pi; 
   w_1   = 2.0*pi;  w_2  = 2.0*pi;  w_3  = 2.0*pi;  w_4  = 2.0*pi;  w_5  = 2.0*pi;

   % set trig function arguments
   arg_1 = x*n_11 + y*n_12 +z*n_13 - w_1*t;
   arg_2 = x*n_21 + y*n_22 +z*n_23 - w_2*t;
   arg_3 = x*n_31 + y*n_32 +z*n_33 - w_3*t;
   arg_4 = x*n_41 + y*n_42 +z*n_43 - w_4*t;
   arg_5 = x*n_51 + y*n_52 +z*n_53 - w_5*t;

   % calculate primary variables
   rho = rho_0 + A_1*sin(arg_1);
   U   = U_0   + A_2*sin(arg_2);
   V   = V_0   + A_3*cos(arg_3);
   W   = W_0   + A_4*cos(arg_4);
   e   = e_0   + A_5*cos(arg_5);

   % calculate rho derivatives
   rho_t  = -A_1*w_1*cos(arg_1);
   rho_x  =  A_1*n_11*cos(arg_1);
   rho_y  =  A_1*n_12*cos(arg_1);
   rho_z  =  A_1*n_13*cos(arg_1);

   % calculate U derivatives
   U_t  = -A_2*w_2*cos(arg_2);
   U_x  =  A_2*n_21*cos(arg_2);
   U_y  =  A_2*n_22*cos(arg_2);
   U_z  =  A_2*n_23*cos(arg_2);

   % calculate V derivatives
   V_t  =  A_3*w_3*sin(arg_3);
   V_x  = -A_3*n_31*sin(arg_3);
   V_y  = -A_3*n_32*sin(arg_3);
   V_z  = -A_3*n_33*sin(arg_3);

   % calculate W derivatives
   W_t  =  A_4*w_4*sin(arg_4);
   W_x  = -A_4*n_41*sin(arg_4);
   W_y  = -A_4*n_42*sin(arg_4);
   W_z  = -A_4*n_43*sin(arg_4);

   % calculate e derivatives
   e_t  =  A_5*w_5*sin(arg_5);
   e_x  = -A_5*n_51*sin(arg_5);
   e_y  = -A_5*n_52*sin(arg_5);
   e_z  = -A_5*n_53*sin(arg_5);

   %
   rhoU = rho.*U;
   rhoV = rho.*V;
   rhoW = rho.*W;
   rhoe = rho.*e;
   %
   rhoU_x = rho_x.*U + rho.*U_x;
   rhoU_y = rho_y.*U + rho.*U_y;
   rhoU_z = rho_z.*U + rho.*U_z;
   %
   rhoV_x = rho_x.*V + rho.*V_x;
   rhoV_y = rho_y.*V + rho.*V_y;
   rhoV_z = rho_z.*V + rho.*V_z;
   %
   rhoW_x = rho_x.*W + rho.*W_x;
   rhoW_y = rho_y.*W + rho.*W_y;
   rhoW_z = rho_z.*W + rho.*W_z;
   %
   rhoe_t = rho_t.*e + rho.*e_t;
   rhoe_x = rho_x.*e + rho.*e_x;
   rhoe_y = rho_y.*e + rho.*e_y;
   rhoe_z = rho_z.*e + rho.*e_z;
   %
   u2pv2   = U.^2 + V.^2 + W.^2;
   u2pv2_t = 2*(U.*U_t + V.*V_t + W.*W_t);
   u2pv2_x = 2*(U.*U_x + V.*V_x + W.*W_x);
   u2pv2_y = 2*(U.*U_y + V.*V_y + W.*W_y);
   u2pv2_z = 2*(U.*U_z + V.*V_z + W.*W_z);

   % calculate pressure and its derivatives
   p   = (gam - 1.0)*rho.*e;
   p_x = (gam - 1.0)*rhoe_x;
   p_y = (gam - 1.0)*rhoe_y;
   p_z = (gam - 1.0)*rhoe_z;
   
   % calculate total energy and its derivatives
   Et    = rhoe     + 0.5*rho.*u2pv2;
   Et_t  = rhoe_t   + 0.5*rho_t.*u2pv2 + 0.5*rho.*u2pv2_t;
   Et_x  = rhoe_x   + 0.5*rho_x.*u2pv2 + 0.5*rho.*u2pv2_x;
   Et_y  = rhoe_y   + 0.5*rho_y.*u2pv2 + 0.5*rho.*u2pv2_y;
   Et_z  = rhoe_z   + 0.5*rho_z.*u2pv2 + 0.5*rho.*u2pv2_z;
   
   % calculate the time derivative of the conserved variable vector
   if(i==1)
      %
      S = rho_t  ...
        + rhoU_x ...
        + rhoV_y...
        + rhoW_z; 
      % 
      Q = rho; 
      %
   elseif(i==2)
      %
      S = rho_t.*U + rho.*U_t ...
        + rhoU_x.*U + rhoU.*U_x + p_x ...
        + rhoV_y.*U + rhoV.*U_y ...
        + rhoW_z.*U + rhoW.*U_z; 
      %
      Q = rhoU;  
      %
   elseif(i==3)
      %
      S = rho_t.*V + rho.*V_t ...
        + rhoU_x.*V + rhoU.*V_x ...
        + rhoV_y.*V + rhoV.*V_y + p_y ...
        + rhoW_z.*V + rhoW.*V_z; 
      %
      Q = rhoV;
      %
   elseif(i==4)
      %
      S = rho_t.*W + rho.*W_t ...
        + rhoU_x.*W + rhoU.*W_x ...
        + rhoV_y.*W + rhoV.*W_y ...
        + rhoW_z.*W + rhoW.*W_z + p_z;
      %
      Q = rhoW;
      %
   elseif(i==5)
      %
      S = Et_t ...
        + U_x.*(Et + p) + U.*(Et_x + p_x) ...
        + V_y.*(Et + p) + V.*(Et_y + p_y)...
        + W_z.*(Et + p) + W.*(Et_z + p_z);
      %
      Q = Et;
      %
   end

   
    