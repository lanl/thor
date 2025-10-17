function write2vtk_structured(dims, X, Q, Qname, filename, timestep)
  %
  [~,~,ext] = fileparts(filename);
  if ~strcmpi(ext,'.vtk')
    filename = [filename, '.vtk'];
  end
  %
  if numel(dims)==2
    nx = dims(1); ny = dims(2); nz = 1;
  elseif numel(dims)==3
    nx = dims(1); ny = dims(2); nz = dims(3);
  else
    error('dims must be length-2 or length-3 vector.');
  end
  Npts = nx*ny*nz;
  [N, Dim] = size(X);
  assert(N==Npts, 'Size(X,1) must equal prod(dims).');
  assert(Dim==2 || Dim==3, 'X must have 2 or 3 columns.');
  %
  pts = [X, zeros(Npts,3-Dim)]';
  %
  fid = fopen(filename, 'w+');
  if fid<0
    error('Could not open %s for writing.', filename);
  end
  %
  fprintf(fid, '# vtk DataFile Version 2.0\n');
  if exist('timestep','var')
    if isnumeric(timestep)
      fprintf(fid, 'Time Step: %g\n', timestep);
    else
      fprintf(fid, 'Time Step: %s\n', timestep);
    end
  else
    fprintf(fid, 'VTK output from MATLAB\n');
  end
  fprintf(fid, 'BINARY\n');
  fprintf(fid, 'DATASET STRUCTURED_GRID\n');
  %
  fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
  fprintf(fid, 'POINTS %d float\n', Npts);
  fwrite(fid, pts, 'float32', 'b');
  %
  [Np, nVar] = size(Q);
  assert(Np==Npts, 'Q must have one row per grid point.');
  fprintf(fid, '\nPOINT_DATA %d\n', Npts);
  for j = 1:nVar
    name = Qname{j};
    fprintf(fid, '\nSCALARS %s float 1\n', name);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fwrite(fid, single(Q(:,j)), 'float32', 'b');
  end
  fclose(fid);
end
  