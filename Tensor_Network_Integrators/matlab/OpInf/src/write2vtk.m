function write2vtk(Connectivity, X, Q, Qname, filename, timestep)
  %WRITE2VTK  Write a binary legacy VTK UnstructuredGrid with cell-centered data.
  %
  %  Connectivity : Nelems×Nc  (1-based node indices per cell; Nc=3 for triangles, 4 for tetras)
  %  X            : Nnodes×Dim (Dim=2 or 3)
  %  Q            : Nelems×nVar  cell-centered solution variables
  %  Qname        : 1×nVar cell array of variable names (strings)
  %  filename     : output file path (will get “.vtk” appended if needed)
  %  timestep     : (opt) numeric or string, written as the title line
  
    %--- ensure .vtk extension
    [~,~,ext] = fileparts(filename);
    if ~strcmpi(ext, '.vtk')
      filename = [filename, '.vtk'];
    end
  
    %--- sizes & checks
    [Nnodes, Dim] = size(X);
    [Nelems, Nc]  = size(Connectivity);
    assert((Dim==2 && Nc==3) || (Dim==3 && Nc==4), ...
           'Connectivity must be Nx3 (triangles) or Nx4 (tetras) for 2D/3D.')
  
    %--- choose VTK cell type code
    if Dim==2
      cellType = 5;    % VTK_TRIANGLE
    else
      cellType = 10;   % VTK_TETRA
    end
  
    %--- pad to 3D points
    pts = [X, zeros(Nnodes, 3-Dim)]';
  
    %--- open in binary mode
    fid = fopen(filename, 'w+');
    if fid<0
      error('Could not open %s for writing.', filename);
    end
  
    %--- HEADER / TITLE / FORMAT
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
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
  
    %--- POINTS block
    fprintf(fid, 'POINTS %d float\n', Nnodes);
    fwrite(fid, pts,    'float32', 'b');  % big-endian floats
  
    %--- CELLS block
    totalIdx = Nelems*(Nc+1);
    fprintf(fid, '\nCELLS %d %d\n', Nelems, totalIdx);
    cellBlock = [ repmat(Nc, Nelems,1), Connectivity-1 ]';  % zero-based
    fwrite(fid, cellBlock, 'int32', 'b');
  
    %--- CELL_TYPES block
    fprintf(fid, '\nCELL_TYPES %d\n', Nelems);
    fwrite(fid, repmat(cellType, Nelems,1), 'int32', 'b');
  
    %--- CELL_DATA block (one SCALARS array per Q-field)
    fprintf(fid, '\nCELL_DATA %d\n', Nelems);
    for j = 1:numel(Qname)
      fprintf(fid, '\nSCALARS %s float 1\n', Qname{j});
      fprintf(fid, 'LOOKUP_TABLE default\n');
      fwrite(fid, single(Q(:,j)), 'float32', 'b');
    end
  
    fclose(fid);
  end
  