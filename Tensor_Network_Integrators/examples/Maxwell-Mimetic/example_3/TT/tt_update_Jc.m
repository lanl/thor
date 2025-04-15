function[ tJx, tJy, tJz ] = tt_update_Jc(dx, dy, dz, dt, nx, ny, nz, time, tol)

% set time
tm=time;
te=time+dt/2;
%% update tJx
xgrid = dx*([1:nx-1]-1);
ygrid = dy*([1:ny]-1);
zgrid = dz*([1:nz]-1);

temp_1 = tt_Ex_value(xgrid+dx/2,ygrid,zgrid,te); %compute Ex_1
temp_0 = tt_Ex_value(xgrid+dx/2,ygrid,zgrid,te-dt); %compute Ex_0
tJx = round(temp_1-temp_0, tol)/dt;

temp_1 = tt_Hz_value(xgrid+dx/2, ygrid+dy/2, zgrid, tm); %compute Hz_1 
temp_0 = tt_Hz_value(xgrid+dx/2, ygrid+dy/2-dy, zgrid, tm); %compute Hz_0 
tJx = round(tJx - (temp_1-temp_0)/dy, tol);

temp_1 = tt_Hy_value(xgrid+dx/2, ygrid, zgrid+dz/2, tm, tol); %compute Hy_1 
temp_0 = tt_Hy_value(xgrid+dx/2, ygrid, zgrid+dz/2-dz, tm, tol); %compute Hy_0 
tJx = round(tJx + (temp_1-temp_0)/dz, tol);

%% update tJy
xgrid = dx*([1:nx]-1);
ygrid = dy*([1:ny-1]-1);
zgrid = dz*([1:nz]-1);

temp_1 = tt_Ey_value(xgrid,ygrid+dy/2,zgrid,te); %compute Ex_1
temp_0 = tt_Ey_value(xgrid,ygrid+dy/2,zgrid,te-dt); %compute Ex_0
tJy = round(temp_1-temp_0, tol)/dt;

temp_1 = tt_Hx_value(xgrid, ygrid+dy/2, zgrid+dz/2, tm, tol); %compute Hz_1 
temp_0 = tt_Hx_value(xgrid, ygrid+dy/2, zgrid+dz/2-dz, tm, tol); %compute Hz_0 
tJy = round(tJy - (temp_1-temp_0)/dz, tol);

temp_1 = tt_Hz_value(xgrid+dx/2, ygrid+dy/2, zgrid, tm); %compute Hy_1 
temp_0 = tt_Hz_value(xgrid+dx/2-dx, ygrid+dy/2, zgrid, tm); %compute Hy_0 
tJy = round(tJy + (temp_1-temp_0)/dx, tol);

%% update tJz
xgrid = dx*([1:nx]-1);
ygrid = dy*([1:ny]-1);
zgrid = dz*([1:nz-1]-1);

temp_1 = tt_Ez_value(xgrid,ygrid,zgrid+dz/2,te, tol); %compute Ez_1
temp_0 = tt_Ez_value(xgrid,ygrid,zgrid+dz/2,te-dt, tol); %compute Ez_0
tJz = round(temp_1-temp_0, tol)/dt;

temp_1 = tt_Hy_value(xgrid+dx/2, ygrid, zgrid+dz/2, tm, tol); %compute Hy_1 
temp_0 = tt_Hy_value(xgrid+dx/2-dx, ygrid, zgrid+dz/2, tm, tol); %compute Hy_0 
tJz = round(tJz - (temp_1-temp_0)/dx, tol);

temp_1 = tt_Hx_value(xgrid, ygrid+dy/2, zgrid+dz/2, tm, tol); %compute Hx_1 
temp_0 = tt_Hx_value(xgrid, ygrid+dy/2-dy, zgrid+dz/2, tm, tol); %compute Hx_0 
tJz = round(tJz + (temp_1-temp_0)/dy, tol);

end