classdef xdmfio < handle
  % xdmfio: XDMF I/O class
  properties
    fid_h5 = -1;   % file ID for the HDF5 file
    gid_h5 = -1;   % group ID for the HDF5
    fid_xdmf = -1; % XDMF file ID

    fname_h5   = 'matlab_output.h5';   % HDF5 file name (default: output.h5)
    fname_xdmf = 'matlab_output.xdmf'; % XDMF file name (default: output.xdmf)
  end

  methods

    function obj = xdmfio(varargin)
      % xdmfio Constructor for the xdmfio class.
      %
      % This constructor initializes an xdmfio object with specified HDF5 and
      % XDMF filenames. If no filenames are provided, default values are used.
      %
      % Syntax:
      %   obj = xdmfio()
      %   obj = xdmfio('fname_h5', h5_filename, 'fname_xdmf', xdmf_filename)
      %
      % Inputs (Name-Value Pairs):
      %   'fname_h5'   - (char) Name of the HDF5 file. Default is 'output.h5'.
      %   'fname_xdmf' - (char) Name of the XDMF file. Default is 'output.xdmf'.
      %
      % Outputs:
      %   obj          - (xdmfio) Instance of the xdmfio class.
      %
      % Example:
      %   % Create an xdmfio object with default filenames
      %   obj = xdmfio();
      %
      %   % Create an xdmfio object with custom filenames
      %   obj = xdmfio('fname_h5', 'mydata.h5', 'fname_xdmf', 'mydata.xdmf');
      %
      p = inputParser;
      addParameter(p, 'fname_h5', 'output.h5', @ischar);
      addParameter(p, 'fname_xdmf', 'output.xdmf', @ischar);
      parse(p, varargin{:});

      obj.fname_h5 = p.Results.fname_h5;
      obj.fname_xdmf = p.Results.fname_xdmf;
    end


    function create_new_output_files(obj)
      % create_new_output_files Creates new XDMF and HDF5 output files.
      %
      % This method initializes new XDMF and HDF5 files associated with the
      % `xdmfio` object. The existing files (if any) are overwritten. The XDMF
      % file is populated with the basic structure required for temporal data.
      %
      % Syntax:
      %   obj.create_new_output_files()
      %
      % Outputs: none
      %
      % Uses:
      %   obj.fname_xdmf - (char) The name of the XDMF file to create.
      %   obj.fname_h5   - (char) The name of the HDF5 file to create.
      %
      % Example:
      %   % Create a new instance of xdmfio and create new output files
      %   xdmf = xdmfio('fname_xdmf', 'output.xdmf', 'fname_h5', 'output.h5');
      %   xdmf.create_new_output_files();
      %
      % Notes:
      %   - This method overwrites any existing files with the same names as
      %     `obj.fname_xdmf` and `obj.fname_h5`.
      %   - The XDMF file is initialized with a basic structure suitable for
      %     temporal data collections.
      obj.fid_xdmf = fopen(obj.fname_xdmf, 'w');
      obj.fid_h5 = H5F.create(obj.fname_h5, ...
        'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
      H5F.close(obj.fid_h5);
      fprintf(obj.fid_xdmf, [...
        '<Xdmf Version=\"3.0\">\n',...
        '  <Domain><Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n',...
        '  </Grid></Domain>\n', ...
        '</Xdmf>']);
      fclose(obj.fid_xdmf);
    end


    function line_count = count_lines_in_file(~,filename)
      % This function counts the number of lines in the specified file.
      %
      % Inputs:
      %   filename - (char) The name of the file to count lines from.
      %
      % Outputs:
      %   line_count - (int) The number of lines in the file. If an error
      %                occurs during file opening, it returns -1.

      fileID = fopen(filename, 'r');
      if fileID == -1
        fprintf('count_lines_in_file: error opening file\n');
        line_count = -1;
        return;
      end

      line_count = 0;
      while ~feof(fileID)
        ch = fgetl(fileID);
        if ischar(ch)
          line_count = line_count + 1;
        end
      end
      fclose(fileID);
    end


    function h5_open_existing_rdwr(obj)
      % This function opens an existing HDF5 file in read/write mode.
      try
        obj.fid_h5 = H5F.open(obj.fname_h5, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
      catch
        fprintf('h5_open_existing_rdwr: error - Unable to open file\n');
        obj.fid_h5 = -1;  % Failure
      end
      % Leave the file open
    end


    function xdmf_open_file_append(obj)
      % This function opens an XDMF file for appending by first reading the
      % file, removing it, recreating it, and then writing all lines except
      % the last two.
      %
      obj.fid_xdmf = -1;

      % Try to read the file into lines
      try
        fileContent = fileread(obj.fname_xdmf);
        lines = splitlines(fileContent);
      catch
        fprintf('xdmf_open_file_append: error reading lines\n');
        return;
      end

      % Remove the file
      if exist(obj.fname_xdmf, 'file') == 2
        delete(obj.fname_xdmf);
      else
        fprintf('xdmf_open_file_append: error deleting file\n');
        return;
      end

      % Create the file again
      obj.fid_xdmf = fopen(obj.fname_xdmf, 'w');
      if obj.fid_xdmf == -1
        fprintf('xdmf_open_file_append: error creating file\n');
        return;
      end

      % Write all lines except the last two (closing lines)
      for i = 1:length(lines)-2
        fprintf(obj.fid_xdmf, '%s\n', lines{i});
      end

      % Leave the file open
    end


    function h5_open_group(obj, group_name)
      % h5_open_group Creates or opens an HDF5 group.
      %
      % This function attempts to create a group within an HDF5 file. If the
      % group already exists, it opens the group instead. If neither
      % operation is successful, it returns -1.
      %
      % Syntax:
      %   obj.h5_open_group(group_name)
      %
      % Inputs:
      %   group_name  - (char) Name of the group to create or open.
      %
      % Uses:
      %   obj.fid_h5  - (hid_t) Identifier of the HDF5 file.
      %
      % Modifies:
      %   obj.gid_h5  - (hid_t) Identifier of the HDF5 group.
      %
      % Attempt to create the group
      try
        obj.gid_h5 = H5G.create(obj.fid_h5, group_name, ...
          'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
      catch
        % If creation fails, clear the error state and try to open the group
        H5E.clear();
        try
          obj.gid_h5 = H5G.open(obj.fid_h5, group_name, 'H5P_DEFAULT');
        catch
          fprintf('h5_open_group: error creating or opening group\n');
          obj.gid_h5 = -1;
        end
      end
    end


    function open_group(obj, group_name)
      % open_group   Creates/opens an HDF5 group; opens XDMF in append
      % mode
      %
      % This function creates or appends a group in an HDF5 file.
      % Opens XDMF file and removes two closing tags lines.
      %
      % Syntax:
      %   obj.open_group(group_name)
      %
      % Inputs:
      %   group_name  - (char) Name of the group to create or open.
      %
      obj.h5_open_existing_rdwr();
      obj.h5_open_group(group_name);
      obj.xdmf_open_file_append();
    end


    function h5_close(obj)
      % h5_close   Closes the current HDF5 group and the HDF5 file.
      %
      if obj.gid_h5 ~= -1
        H5G.close(obj.gid_h5);
      end
      obj.gid_h5 = -1;

      if obj.fid_h5 ~= -1
        H5F.close(obj.fid_h5);
      end
      obj.fid_h5 = -1;
    end


    function h5_write_1d(obj, dataset_name, data)
      % h5_write_1d - Writes a 1D array to an HDF5 file or group.
      %
      % This function writes a 1D array to an HDF5 file or group. It creates
      % the necessary dataspace and dataset within the specified file or
      % group, then writes the provided data to the dataset.
      %
      % Syntax:
      %   obj.h5_write_1d(dataset_name, data)
      %
      % Inputs:
      %   dataset_name - (char) Name of the dataset to create and write to.
      %   data         - (1D array of double) The data to write to the dataset.
      %
      % Uses:
      %   gid_h5       - (H5ML.id) group ID where to write the data to.
      %
      % Example:
      %   xio = xmdfio();
      %   xio.create_new_output_files();
      %
      fogr_id = obj.gid_h5;
      len = numel(data);

      % Define the dimensions of the 1D array (dataspace)
      dataspace_id = H5S.create_simple(1, len, []);
      if dataspace_id < 0
        fprintf('h5_write_1d: Error creating dataspace\n');
        H5G.close(fogr_id);
        return;
      end

      % Create the dataset within the file or group
      dataset_id = H5D.create(fogr_id, dataset_name, ...
        'H5T_NATIVE_DOUBLE', dataspace_id, ...
        'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
      if dataset_id < 0
        fprintf('h5_write_1d: Error creating dataset\n');
        H5S.close(dataspace_id);
        H5G.close(fogr_id);
        return;
      end

      % Write the 2D array data to the dataset
      try
        H5D.write(dataset_id, 'H5T_NATIVE_DOUBLE', ...
          'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data(:));
      catch
        fprintf('h5_write_1d: Error writing data to dataset\n');
      end

      % Close the dataset and dataspace
      H5D.close(dataset_id);
      H5S.close(dataspace_id);

    end


    function h5_write_2d(obj, dataset_name, data)
      % h5_write_2d Writes a 2D array to an HDF5 file or group.
      %
      % This function writes a 2D array to an HDF5 file or group. It creates
      % the necessary dataspace and dataset within the specified file or
      % group, then writes the provided data to the dataset.
      %
      % Syntax:
      %   obj.h5_write_2d(dataset_name, data)
      %
      % Inputs:
      %   dataset_name - (char) Name of the dataset to create and write to.
      %   data         - (2D array of double) The data to write to the dataset.
      %
      % Uses:
      %   gid_h5       - (H5ML.id) group ID where to write the data to.
      %
      % Example:
      %   xio = xmdfio();
      %   xio.create_new_output_files();
      %
      fogr_id = obj.gid_h5;
      dims = flip(size(data));

      % Define the dimensions of the 2D array (dataspace)
      dataspace_id = H5S.create_simple(2, dims, []);
      if dataspace_id < 0
        fprintf('h5_write_2d: Error creating dataspace\n');
        H5G.close(fogr_id);
        return;
      end

      % Create the dataset within the file or group
      dataset_id = H5D.create(fogr_id, dataset_name, ...
        'H5T_NATIVE_DOUBLE', dataspace_id, ...
        'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
      if dataset_id < 0
        fprintf('h5_write_2d: Error creating dataset\n');
        H5S.close(dataspace_id);
        H5G.close(fogr_id);
        return;
      end

      % Write the 2D array data to the dataset
      try
        H5D.write(dataset_id, 'H5T_NATIVE_DOUBLE', ...
          'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data');
      catch
        fprintf('h5_write_2d: Error writing data to dataset\n');
      end

      % Close the dataset and dataspace
      H5D.close(dataset_id);
      H5S.close(dataspace_id);

    end


    function h5_write_3d(obj, dataset_name, data)
      % h5_write_3d  Writes a 3D array to an HDF5 file or group.
      %
      % This function writes a 3D array to an HDF5 file or group. It creates
      % the necessary dataspace and dataset within the specified file or
      % group, then writes the provided data to the dataset.
      %
      % Syntax:
      %   obj.h5_write_3d(dataset_name, data)
      %
      % Inputs:
      %   dataset_name - (char) Name of the dataset to create and write to.
      %   data         - (3D array of double) The data to write to the dataset.
      %
      % Uses:
      %   gid_h5       - (H5ML.id) group ID where to write the data to.
      %
      % Example:
      %   xio = xmdfio();
      %   xio.create_new_output_files();
      %
      fogr_id = obj.gid_h5;
      dims = flip(size(data));

      % Define the dimensions of the 2D array (dataspace)
      dataspace_id = H5S.create_simple(3, dims, []);
      if dataspace_id < 0
        fprintf('h5_write_3d: Error creating dataspace\n');
        H5G.close(fogr_id);
        return;
      end

      % Create the dataset within the file or group
      dataset_id = H5D.create(fogr_id, dataset_name, ...
        'H5T_NATIVE_DOUBLE', dataspace_id, ...
        'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
      if dataset_id < 0
        fprintf('h5_write_3d: Error creating dataset\n');
        H5S.close(dataspace_id);
        H5G.close(fogr_id);
        return;
      end

      % Write the 3D array data to the dataset
      % TODO: the data may need some transposing
      data_tsp = permute(data, [3, 2, 1]);
      try
        H5D.write(dataset_id, 'H5T_NATIVE_DOUBLE', ...
          'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data_tsp);
      catch
        fprintf('h5_write_3d: Error writing data to dataset\n');
      end

      % Close the dataset and dataspace
      H5D.close(dataset_id);
      H5S.close(dataspace_id);

    end


    function xdmf_write_mesh(obj, mesh)
      % xdmf_write_mesh Writes mesh information to the XDMF file.
      %
      % This method writes mesh data to the XDMF file associated with the
      % `xdmfio` object. It takes the properties from a `c_mesh` object and
      % writes the necessary grid, topology, and geometry information
      % based on the dimensionality of the mesh.
      %
      % Syntax:
      %   obj.xdmf_write_mesh(mesh)
      %
      % Inputs:
      %   mesh - (c_mesh) Object containing mesh properties such as dimension,
      %          domain boundaries, number of cells, and number of boundary
      %          cells.
      %
      % Outputs:
      %   None. The function modifies the XDMF file associated with the
      %   `xdmfio` object.
      %
      % Example:
      %   % Assuming xdmf is an instance of xdmfio and mesh is an instance of c_mesh
      %   xio.xdmf_write_mesh(mesh);
      %
      % See also: xdmfio, fopen, fclose

      % Extract properties from the mesh object
      dim = mesh.dim_prob;
      xl = mesh.xl_prob;
      xr = mesh.xr_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;

      % Calculate the grid spacing
      dx = (xr - xl) ./ ncell;

      % Check the dimensionality and write appropriate XDMF data
      if dim == 2
        % Write 2D grid information
        fprintf(obj.fid_xdmf, "\n" + ...
          "<Grid Name=""StructuredGrid"">\n" + ...
          " <Topology TopologyType=""2DCoRectMesh"" Dimensions=""%d %d""/>\n" + ...
          " <Geometry GeometryType=""ORIGIN_DXDY"">\n" + ...
          "  <DataItem DataType=""Float"" Dimensions=""2"" Format=""XML"">\n" + ...
          "    %14.7e %14.7e\n" + ...
          "  </DataItem>\n" + ...
          "  <DataItem DataType=""Float"" Dimensions=""2"" Format=""XML"">\n" + ...
          "    %14.7e %14.7e\n" + ...
          "  </DataItem>\n" + ...
          " </Geometry>\n", ...
          ncell(2) + 2*nbdry + 1, ncell(1) + 2*nbdry + 1, ...
          xl(2) - nbdry*dx(2), xl(1) - nbdry*dx(1), dx(2), dx(1));
      elseif dim == 3
        % Write 3D grid information
        fprintf(obj.fid_xdmf, "\n" + ...
          "<Grid Name=""StructuredGrid"">\n" + ...
          " <Topology TopologyType=""3DCoRectMesh"" Dimensions=""%d %d %d""/>\n" + ...
          " <Geometry GeometryType=""ORIGIN_DXDYDZ"">\n" + ...
          "  <DataItem DataType=""Float"" Dimensions=""3"" Format=""XML"">\n" + ...
          "    %14.7e %14.7e %14.7e\n" + ...
          "  </DataItem>\n" + ...
          "  <DataItem DataType=""Float"" Dimensions=""3"" Format=""XML"">\n" + ...
          "    %14.7e %14.7e %14.7e\n" + ...
          "  </DataItem>\n" + ...
          " </Geometry>\n", ...
          ncell(3) + 2*nbdry + 1, ncell(2) + 2*nbdry + 1, ncell(1) + 2*nbdry + 1, ...
          xl(3) - nbdry*dx(3), xl(2) - nbdry*dx(2), xl(1) - nbdry*dx(1), ...
          dx(3), dx(2), dx(1));
      end
    end


    function xdmf_write_dataset(obj, mesh, timestamp_group, varname, var_type)
      % xdmf_write_dataset Writes dataset information to the XDMF file.
      %
      % This method writes dataset information to the XDMF file associated with the
      % `xdmfio` object. The mesh properties such as dimension, domain boundaries,
      % grid spacing, and number of cells are obtained from the `mesh` object, while
      % the dataset name and type are passed as arguments.
      %
      % Syntax:
      %   obj.xdmf_write_dataset(mesh, varname, var_type)
      %
      % Inputs:
      %   mesh      - (c_mesh) Object containing mesh properties such as dimension,
      %               domain boundaries, and number of cells.
      %   timestamp_group - (char) Name of the time group.
      %   varname   - (char) Name of the variable (dataset) to be written.
      %   var_type  - (variable_type) Enum specifying whether the variable is a
      %               CELL_VARIABLE or NODE_VARIABLE.
      %
      % Outputs:
      %   None. The function modifies the XDMF file associated with the
      %   `xdmfio` object.
      %
      % Example:
      %   % Assuming xdmf is an instance of xdmfio and mesh is an instance of c_mesh
      %   xdmf.xdmf_write_dataset(mesh, 'temperature', variable_type.CELL_VARIABLE);
      %
      % See also: xdmfio, c_mesh, fopen, fclose

      % Extract properties from the mesh object
      dim = mesh.dim_prob;
      xl = mesh.xl_prob;
      ncell = mesh.ncell_prob + 2*mesh.nbdry_prob;

      % Check the dimensionality and write appropriate XDMF data
      if dim == 2
        if var_type == variable_type.CELL_VARIABLE
          % Write 2D cell-centered data information
          fprintf(obj.fid_xdmf, "\n" + ...
            " <Attribute Name=""%s"" AttributeType=""Scalar"" Center=""Cell"">\n" + ...
            "   <DataItem Dimensions=""%d %d""\n" + ...
            "     NumberType=""Float"" Precision=""4"" Format=""HDF"">\n" + ...
            "     %s:/%s/%s\n" + ...
            "   </DataItem>\n" + ...
            " </Attribute>\n", varname, ...
            ncell(2), ncell(1), ...
            obj.fname_h5, timestamp_group, varname);
        elseif var_type == variable_type.NODE_VARIABLE
          % Write 2D node-centered data information
          fprintf(obj.fid_xdmf, "\n" + ...
            " <Attribute Name=""%s"" AttributeType=""Scalar"" Center=""Node"">\n" + ...
            "   <DataItem Dimensions=""%d %d""\n" + ...
            "     NumberType=""Float"" Precision=""4"" Format=""HDF"">\n" + ...
            "     %s:/%s/%s\n" + ...
            "   </DataItem>\n" + ...
            " </Attribute>\n", varname, ...
            ncell(2) + 1, ncell(1) + 1, ...
            obj.fname_h5, timestamp_group, varname);
        end
      elseif dim == 3
        if var_type == variable_type.CELL_VARIABLE
          % Write 3D cell-centered data information
          fprintf(obj.fid_xdmf, "\n" + ...
            " <Attribute Name=""%s"" AttributeType=""Scalar"" Center=""Cell"">\n" + ...
            "   <DataItem Dimensions=""%d %d %d""\n" + ...
            "     NumberType=""Float"" Precision=""4"" Format=""HDF"">\n" + ...
            "     %s:/%s/%s\n" + ...
            "   </DataItem>\n" + ...
            " </Attribute>\n", varname, ...
            ncell(3), ncell(2), ncell(1), ...
            obj.fname_h5, timestamp_group, varname);
        elseif var_type == variable_type.NODE_VARIABLE
          % Write 3D node-centered data information
          fprintf(obj.fid_xdmf, "\n" + ...
            " <Attribute Name=""%s"" AttributeType=""Scalar"" Center=""Node"">\n" + ...
            "   <DataItem Dimensions=""%d %d %d""\n" + ...
            "     NumberType=""Float"" Precision=""4"" Format=""HDF"">\n" + ...
            "     %s:/%s/%s\n" + ...
            "   </DataItem>\n" + ...
            " </Attribute>\n", varname, ...
            ncell(3) + 1, ncell(2) + 1, ncell(1) + 1, ...
            obj.fname_h5, timestamp_group, varname);
        end
      end
    end


    function xdmf_write_timestamp(obj, timestamp)
      % xdmf_write_timestamp    Writes a timestamp to the XDMF file.
      %
      % Inputs:
      %   timestamp - (double) The timestamp value to write in the XDMF file.
      %
      fprintf(obj.fid_xdmf, "\n <Time Value=""%14.7e"" />\n", timestamp);
    end


    function xdmf_close_group(obj)
      % xdmf_close_group   Closes the current XDMF group in the file.
      %
      % This method closes the current `<Grid>` group in the XDMF file
      % associated with the `xdmfio` object by writing the closing tag.
      fprintf(obj.fid_xdmf, "</Grid>\n\n");
    end


    function xdmf_close(obj)
      % xdmf_close_file Closes the XDMF file after writing closing tags.
      %
      % This method writes the necessary closing tags in XDMF file
      %
      % Example:
      %   % Assuming the object 'xio' that was used to create and open files,
      %   % close the XDMF file properly:
      %   xio.xdmf_close_file();
      %
      fprintf(obj.fid_xdmf, [...
        '  </Grid></Domain>\n' ...
        '</Xdmf>']);
      fclose(obj.fid_xdmf);
      obj.fid_xdmf = -1;
    end


    function close(obj)
      % close   Closes the currently open group, then HDF5 and XDMF
      %
      obj.h5_close();
      obj.xdmf_close_group();
      obj.xdmf_close();
    end


    function write_dump(obj, mesh, mat, ncycle, timestamp)
      stepstr = sprintf('Step#%d', ncycle);
      obj.open_group(stepstr);
      obj.xdmf_write_mesh(mesh);

      if mesh.dim_prob == 2
        obj.h5_write_2d('rho', mesh.rho_for_2dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'rho', variable_type.CELL_VARIABLE);

        obj.h5_write_2d('vx', mesh.vel_for_2dnode(:,:,1));
        obj.xdmf_write_dataset(mesh, stepstr, 'vx', variable_type.NODE_VARIABLE);

        obj.h5_write_2d('vy', mesh.vel_for_2dnode(:,:,2));
        obj.xdmf_write_dataset(mesh, stepstr, 'vy', variable_type.NODE_VARIABLE);

        obj.h5_write_2d('ei', mesh.ei_for_2dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'ei', variable_type.CELL_VARIABLE);

        obj.h5_write_2d('pres', mesh.pres_for_2dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'pres', variable_type.CELL_VARIABLE);

        obj.h5_write_2d('cs', mesh.cs_for_2dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'cs', variable_type.CELL_VARIABLE);

        obj.h5_write_2d('nmat', mesh.nmat_for_2dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'nmat', variable_type.CELL_VARIABLE);

        obj.h5_write_2d('matid', mesh.matid_for_2dcell - 1);
        obj.xdmf_write_dataset(mesh, stepstr, 'matid', variable_type.CELL_VARIABLE);

        obj.h5_write_1d('nmixcell', mat.nmixcell);
        obj.h5_write_1d('nmixcell_int', mat.nmixcell_int);
        obj.h5_write_1d('nmixcell_mpoly', mat.nmixcell_mpoly);
        obj.h5_write_1d('lsize', sum(mat.nmat_in_mixcell));
      else
        obj.h5_write_3d('rho', mesh.rho_for_3dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'rho', variable_type.CELL_VARIABLE);

        obj.h5_write_3d('vx', mesh.vel_for_3dnode(:,:,:,1));
        obj.xdmf_write_dataset(mesh, stepstr, 'vx', variable_type.NODE_VARIABLE);

        obj.h5_write_3d('vy', mesh.vel_for_3dnode(:,:,:,2));
        obj.xdmf_write_dataset(mesh, stepstr, 'vy', variable_type.NODE_VARIABLE);

        obj.h5_write_3d('vz', mesh.vel_for_3dnode(:,:,:,3));
        obj.xdmf_write_dataset(mesh, stepstr, 'vz', variable_type.NODE_VARIABLE);

        obj.h5_write_3d('ei', mesh.ei_for_3dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'ei', variable_type.CELL_VARIABLE);

        obj.h5_write_3d('pres', mesh.pres_for_3dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'pres', variable_type.CELL_VARIABLE);

        obj.h5_write_3d('cs', mesh.cs_for_3dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'cs', variable_type.CELL_VARIABLE);

        obj.h5_write_3d('nmat', mesh.nmat_for_3dcell);
        obj.xdmf_write_dataset(mesh, stepstr, 'nmat', variable_type.CELL_VARIABLE);

        obj.h5_write_3d('matid', mesh.matid_for_3dcell - 1);
        obj.xdmf_write_dataset(mesh, stepstr, 'matid', variable_type.CELL_VARIABLE);

        obj.h5_write_1d('nmixcell', mat.nmixcell);
        obj.h5_write_1d('nmixcell_int', mat.nmixcell_int);
        obj.h5_write_1d('nmixcell_mpoly', mat.nmixcell_mpoly);
        obj.h5_write_1d('lsize', sum(mat.nmat_in_mixcell));
      end

      if mat.nmixcell > 0
        lsize = sum(mat.nmat_in_mixcell);
        obj.h5_write_1d('nmat_in_mixcell', mat.nmat_in_mixcell);
        ijk_tsp = mat.ijk_in_mixcell';
        obj.h5_write_1d('ijk_in_mixcell', ijk_tsp(:)-1);
        obj.h5_write_1d('matids_in_mixcell', obj.flatten_cell1d(mat.matids_in_mixcell)-1);
        obj.h5_write_1d('vf_in_mixcell', obj.flatten_cell1d(mat.vf_in_mixcell));
        obj.h5_write_1d('rho_in_mixcell', obj.flatten_cell1d(mat.rho_in_mixcell));
        obj.h5_write_1d('pres_in_mixcell', obj.flatten_cell1d(mat.pres_in_mixcell));
        obj.h5_write_1d('ei_in_mixcell', obj.flatten_cell1d(mat.ei_in_mixcell));
      end

      % obj.xdmf_write_timestamp(timestamp);
      obj.xdmf_write_timestamp(timestamp + 1e-5*ncycle); % hack
      obj.close();
    end

    function flattened_cell_array = flatten_cell1d(~, carr)
      lsize = 0;
      for i = 1:numel(carr)
        lsize = lsize + numel(carr{i});
      end
      flattened_cell_array = zeros(1, lsize);
      offset = 1;
      for i = 1:numel(carr)
        flattened_cell_array(offset:offset + numel(carr{i}) - 1) ...
          = carr{i}(1:end);
        offset = offset + numel(carr{i});
      end
    end


  end
end

