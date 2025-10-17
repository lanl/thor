function obj = tt_cell_centric_update_materials_after_advection_3d(obj, dir)
% Update material quantities after advection in 3D
% This function updates the material volume fractions, densities, and energies
% in the computational domain after advection is performed in a given direction (1: x, 2: y, 3: z).

global ttmesh tt_tol;

sizes_edge = ttmesh.ncell_prob + 2 * ttmesh.nbdry_prob; % Number of cells including boundary layers
nmat = size(obj.vol_for_cell,4);
% For 1D (x-direction)
if dir == 1
  dimidx = 3;
  Iy = 1:sizes_edge(2); % Row index range
  Iz = 1:sizes_edge(3); % Depth index range
  Ix = 1:sizes_edge(1)-1; % Column index range (one less to access the right side)
  ops = {'none','none','diff','none'};
  CIdx = {{Iz,0},{Iy,0},{Ix+1,Ix},{1:nmat,0}};
  factor = [1,1,1,1];
  % % Update the cell properties for volume, mass, and energy across the edge
  % obj.vol_for_cell(Iz, Iy, Ix, :) = obj.vol_for_cell(Iz, Iy, Ix, :) + obj.dvol_for_edge(Iz, Iy, Ix, :) ...
  %   - obj.dvol_for_edge(Iz, Iy, Ix+1, :);
  %
  % obj.mass_for_cell(Iz, Iy, Ix, :) = obj.mass_for_cell(Iz, Iy, Ix, :) + obj.dmass_for_edge(Iz, Iy, Ix, :) ...
  %   - obj.dmass_for_edge(Iz, Iy, Ix+1, :);
  %
  % obj.ener_for_cell(Iz, Iy, Ix, :) = obj.ener_for_cell(Iz, Iy, Ix, :) + obj.dener_for_edge(Iz, Iy, Ix, :) ...
  %   - obj.dener_for_edge(Iz, Iy, Ix+1, :);

  % For 2D (y-direction)
elseif dir == 2
  dimidx = 2;
  Iy = 1:sizes_edge(2)-1; % Row index range (one less to access the lower side)
  Iz = 1:sizes_edge(3);   % Depth index range
  Ix = 1:sizes_edge(1);   % Column index range
  ops = {'none','diff','none','none'};
  CIdx = {{Iz,0},{Iy+1,Iy},{Ix,0},{1:nmat,0}};
  factor = [1,1,1,1];
  % Update the cell properties for volume, mass, and energy across the edge
  % obj.vol_for_cell(Iz, Iy, Ix, :) = obj.vol_for_cell(Iz, Iy, Ix, :) + obj.dvol_for_edge(Iz, Iy, Ix, :) ...
  %   - obj.dvol_for_edge(Iz, Iy+1, Ix, :);
  %
  % obj.mass_for_cell(Iz, Iy, Ix, :) = obj.mass_for_cell(Iz, Iy, Ix, :) + obj.dmass_for_edge(Iz, Iy, Ix, :) ...
  %   - obj.dmass_for_edge(Iz, Iy+1, Ix, :);
  %
  % obj.ener_for_cell(Iz, Iy, Ix, :) = obj.ener_for_cell(Iz, Iy, Ix, :) + obj.dener_for_edge(Iz, Iy, Ix, :) ...
  %   - obj.dener_for_edge(Iz, Iy+1, Ix, :);

  % For 3D (z-direction)
elseif dir == 3
  dimidx = 1;
  Iy = 1:sizes_edge(2); % Row index range
  Iz = 1:sizes_edge(3)-1; % Depth index range (one less to access the lower side)
  Ix = 1:sizes_edge(1);   % Column index range
  ops = {'diff','none','none','none'};
  CIdx = {{Iz+1,Iz},{Iy,0},{Ix,0},{1:nmat,0}};
  factor = [1,1,1,1];
  % Update the cell properties for volume, mass, and energy across the edge
  % obj.vol_for_cell(Iz, Iy, Ix, :) = obj.vol_for_cell(Iz, Iy, Ix, :) + obj.dvol_for_edge(Iz, Iy, Ix, :) ...
  %   - obj.dvol_for_edge(Iz+1, Iy, Ix, :);
  %
  % obj.mass_for_cell(Iz, Iy, Ix, :) = obj.mass_for_cell(Iz, Iy, Ix, :) + obj.dmass_for_edge(Iz, Iy, Ix, :) ...
  %   - obj.dmass_for_edge(Iz+1, Iy, Ix, :);
  %
  % obj.ener_for_cell(Iz, Iy, Ix, :) = obj.ener_for_cell(Iz, Iy, Ix, :) + obj.dener_for_edge(Iz, Iy, Ix, :) ...
  %   - obj.dener_for_edge(Iz+1, Iy, Ix, :);
end

dvoltt = FD_tt_ops(obj.dvol_for_edge, ops, CIdx, factor);
dvoltt{dimidx} = cat(2,dvoltt{dimidx},zeros(dvoltt.r(dimidx),1,dvoltt.r(dimidx+1)));
obj.vol_for_cell = round(obj.vol_for_cell + dvoltt,tt_tol);

dmasstt = FD_tt_ops(obj.dmass_for_edge, ops, CIdx, factor);
dmasstt{dimidx} = cat(2,dmasstt{dimidx},zeros(dmasstt.r(dimidx),1,dmasstt.r(dimidx+1)));
obj.mass_for_cell = round(obj.mass_for_cell + dmasstt,tt_tol);

denertt = FD_tt_ops(obj.dener_for_edge, ops, CIdx, factor);
denertt{dimidx} = cat(2,denertt{dimidx},zeros(denertt.r(dimidx),1,denertt.r(dimidx+1)));
obj.ener_for_cell = round(obj.ener_for_cell + denertt,tt_tol);