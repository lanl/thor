function [Index1E, Index2E, Index3E, m3D] = N3DIndex_grid(knot1, knot2, knot3)
  % same up through el_* and nx,ny,nz
  m3D   = Getmacro3D(knot1, knot2, knot3);
  [~,p1] = EvaluateKnot(knot1);
  [~,p2] = EvaluateKnot(knot2);
  [~,p3] = EvaluateKnot(knot3);
  %n1 = length(knot1)-p1-1;  n2 = length(knot2)-p2-1;  n3 = length(knot3)-p3-1;
  el_x = find(diff(knot1)>0);
  el_y = find(diff(knot2)>0);
  el_z = find(diff(knot3)>0);
  nx = numel(el_x);  ny = numel(el_y);  nz = numel(el_z);

  % preallocate
  Index1E = zeros(nx, ny, nz, p1+1);
  Index2E = zeros(nx, ny, nz, p2+1);
  Index3E = zeros(nx, ny, nz, p3+1);

  for i = 1:nx
    ix = (el_x(i)-p1) : el_x(i);
    for j = 1:ny
      iy = (el_y(j)-p2) : el_y(j);
      for k = 1:nz
        iz = (el_z(k)-p3) : el_z(k);

        % assign the (p+1)-vector into that last dim
        Index1E(i,j,k,:) = ix;
        Index2E(i,j,k,:) = iy;
        Index3E(i,j,k,:) = iz;
      end
    end
  end
end
