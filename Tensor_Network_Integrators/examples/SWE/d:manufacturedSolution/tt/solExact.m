function Q=solExact(x,y,t,data)
    %
    %  Q: nondimensional conserved variables
    %  x,y,t: nondimensional x/y coordinates and time 
    %
    data=manuf(x,y,t,data);
    %
    Q = cell(1,3);
    %
    for i=1:data.tt.Neq
        Q{i}    = applyQuadrature(data.tt.Q0{i},data.quad,data.h);
        Q{i}{1} = Q{i}{1}/data.vol;
        Q{i}    = round(Q{i},data.tt.eps_rk(i));
    end
    %
end