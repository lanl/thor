classdef unittest < handle

    properties
        testdir = '../testdata';
    end

    methods
        function var = h5read_transposed(~,filename, varname)
            var = h5read(filename, varname);
            dims = size(var);
            var = reshape(var, [dims(2), dims(1)])';
        end

        function [total, passed] = test_rec_rec(obj, varargin)
            % test_rec_rec  Runs a test for geom.rec_rec
            %
            % Test Matlab function geom.rec_rec by comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/test_rec_rec.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_rec_rec();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/test_rec_rec.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;

            % The geom object
            geom = c_geom();

            ifinquiry = h5read(filename, '/ifinquiry');
            geop = h5read(filename, '/geop');
            vcell = h5read(filename, '/vcell');
            dim = h5read(filename, '/dim');
            xxl = h5read(filename, '/xxl');
            xxr = h5read(filename, '/xxr');
            xl = h5read(filename, '/xl');
            dx = h5read(filename, '/dx');

            ifmixed = h5read(filename, '/ifmixed');
            vol = h5read(filename, '/vol');

            total = length(ifinquiry);
            passed = 0;

            for i = 1:total
                [ifmixed_, vol_] = geom.rec_rec(ifinquiry(i), geop(i), vcell(i), dim(i), xxl(:,i), xxr(:,i), xl(:,i), dx(:,i));

                if (abs(vol_ - vol(i)) > 1e-14) || (ifmixed_ ~= ifmixed(i))
                    fprintf('failed test %4d: expected {mixed,vol}={%d, %f}, got {%d, %f}\n',...
                        i, ifmixed(i), vol(i), ifmixed_, vol_);
                else
                    passed = passed + 1;
                    if verbose
                        fprintf('passed test %4d: expected {mixed,vol}={%d, %f}, got {%d, %f}\n',...
                        i, ifmixed(i), vol(i), ifmixed_, vol_);
                    end

                end
            end
        end

        function [total, passed] = test_sound_speed_solid(obj, varargin)
            % test_sound_speed_solid  Runs a test for eos.sound_speed_solid
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/test_rec_rec.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_sound_speed_solid();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/sound_speed_solid.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;

            % The geom object
            eos = c_eos();

            rho = h5read(filename, '/rho');
            cs = h5read(filename, '/cs');

            total = length(rho);
            passed = 0;

            for i = 1:total
                cs_ = eos.sound_speed_solid(rho(i));

                if (abs(cs_ - cs(i)) > 1e-14)
                    fprintf('failed test %4d: expected {cs}={%f}, got {%f}\n',...
                        i, cs(i), cs_);
                else
                    passed = passed + 1;
                    if verbose
                        fprintf('passed test %4d: expected {cs}={%f}, got {%f}\n',...
                        i, cs(i), cs_);
                    end

                end
            end
        end

        function [meshes, num_read] = read_meshes_from_hdf5(obj, varargin)
            % read_meshes_from_hdf5 Reads the HDF5 file containing mesh data and creates c_mesh objects.
            %
            % This function reads an HDF5 file that contains multiple mesh groups (e.g., mesh_0, mesh_1, ...)
            % and creates a corresponding `c_mesh` object for each mesh group.
            %
            % Optional Inputs (via varargin):
            %   filename - (string, optional) The path to the HDF5 file containing the meshes.
            %
            % Outputs:
            %   meshes   - (cell array of c_mesh objects) The created mesh objects.
            %   num_read - (integer) Number of meshes that have been read
            %              successfully.
            %

            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/test_meshes.h5'], @isstring);
            parse(p, varargin{:});
            filename = p.Results.filename;

            info = h5info(filename);  % Get HDF5 file info
            num_meshes = length(info.Groups);  % Number of groups (meshes) in the file

            % Initialize a cell array to store the mesh objects
            meshes = cell(1, num_meshes);
            bdry = c_bdry();
            bdry_arg = [bdry_type.bdry_transmitted, bdry_type.bdry_transmitted];

            % Loop through each mesh group
            failed = 0;
            for i = 1:num_meshes
                group_name = info.Groups(i).Name;  % e.g., '/mesh_0', '/mesh_1', etc.

                % Read the mesh attributes
                dim = h5readatt(filename, group_name, 'dim');
                ncell = h5readatt(filename, group_name, 'ncell');
                nbdry = h5readatt(filename, group_name, 'nbdry');
                xl = h5readatt(filename, group_name, 'xl');
                xr = h5readatt(filename, group_name, 'xr');
                dx = h5readatt(filename, group_name, 'dx');

                % Check
                if norm(dx(1:dim) - (xr(1:dim) - xl(1:dim))./double(ncell(1:dim))) > 1e-14
                    failed = failed + 1;
                    continue;
                end

                % Create a new c_mesh object with the attributes from the HDF5 file
                mesh = c_mesh();
                mesh.dim_prob = dim;
                mesh.ncell_prob = ncell;
                mesh.nbdry_prob = nbdry;
                mesh.xl_prob = xl;
                mesh.xr_prob = xr;
                mesh.dx = dx;

                mesh.rho_for_2dcell   = obj.h5read_transposed(filename, [group_name, '/rho']);
                mesh.ei_for_2dcell    = obj.h5read_transposed(filename, [group_name, '/ei']);
                mesh.pres_for_2dcell  = obj.h5read_transposed(filename, [group_name, '/pres']);
                mesh.nmat_for_2dcell  = obj.h5read_transposed(filename, [group_name, '/nmat']);
                mesh.cs_for_2dcell    = obj.h5read_transposed(filename, [group_name, '/cs']);

                % in Matlab, unlike C, material IDs start with 1, so we
                % need to add 1 to the test data
                mesh.matid_for_2dcell = obj.h5read_transposed(filename, [group_name, '/matid']) + 1;

                nnode_ext = ncell; nnode_ext = nnode_ext + 2*nbdry + 1;
                mesh.vel_for_2dnode = zeros(nnode_ext(2), nnode_ext(1), dim);
                mesh.vel_for_2dnode(:,:,1) = obj.h5read_transposed(filename, [group_name, '/vel_x']);
                mesh.vel_for_2dnode(:,:,2) = obj.h5read_transposed(filename, [group_name, '/vel_y']);

                mesh.vav_for_2dnode = zeros(nnode_ext(2), nnode_ext(1), dim);
                mesh.vav_for_2dnode(:,:,1) = obj.h5read_transposed(filename, [group_name, '/vav_x']);
                mesh.vav_for_2dnode(:,:,2) = obj.h5read_transposed(filename, [group_name, '/vav_y']);

                % Store the mesh object in the cell array
                meshes{i} = mesh;
            end
            num_read = num_meshes - failed;
        end

        function [total, passed] = test_bdry_cell_2d(obj, meshes, varargin)
            % test_bdry_cell_2d  Runs a test for bdry.bdry_cell_2d
            %
            % Test Matlab function bdry_cell_2d by comparing results to C output, saved
            % in HDF5 file. The input meshes are provided by the `meshes` argument.
            % The output meshes are read from the optional `filename` file.
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/bdry_cell_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_bdry_cell_2d(meshes);
            %   passed = ut.test_bdry_cell_2d(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/bdry_cell_2d.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;

            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;

            % Initialize a cell array to store the mesh objects
            bdry = c_bdry();
            bdry_arg = [bdry_type.bdry_transmitted, bdry_type.bdry_transmitted];

            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1;
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end


            for i = imn:imx
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                ncell_ext = ncell + 2 * nbdry;  % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;


                % Function call
                bdry.bdry_cell_2d(meshes{i}, bdry_arg, bdry_arg);

                % Read fields after function call
                %group_name = info.Groups(i).Name;  % e.g., '/mesh_0', '/mesh_1', etc.
                test_group_name = sprintf('/mesh_%d', i - 1);

                rho_after   = obj.h5read_transposed(filename, [test_group_name, '/rho']);
                ei_after    = obj.h5read_transposed(filename, [test_group_name, '/ei']);
                pres_after  = obj.h5read_transposed(filename, [test_group_name, '/pres']);
                nmat_after  = obj.h5read_transposed(filename, [test_group_name, '/nmat']);
                matid_after = obj.h5read_transposed(filename, [test_group_name, '/matid']) + 1;

                % Check
                nfailed = failed + 1;
                if ~isequal(meshes{i}.nmat_for_2dcell, nmat_after)
                    if verbose, fprintf("%s: nmat_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end

                if ~isequal(meshes{i}.matid_for_2dcell, matid_after)
                    if verbose, fprintf("%s: matid_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end

                if ~all(abs(meshes{i}.rho_for_2dcell - rho_after) < tol)
                    if verbose, fprintf("%s: rho_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end

                if ~all(abs(meshes{i}.ei_for_2dcell - ei_after) < tol)
                    if verbose, fprintf("%s: ei_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end

                if ~all(abs(meshes{i}.pres_for_2dcell - pres_after) < tol)
                    if verbose, fprintf("%s: pres_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if failed < nfailed
                    passed = passed + 1;
                end
            end
        end

        function [total, passed] = test_bdry_node_2d(obj, meshes, varargin)
            % test_bdry_node_2d  Runs a test for bdry.bdry_node_2d
            %
            % Test Matlab function bdry_node_2d by comparing results to C output, saved
            % in HDF5 file. The input meshes are provided by the `meshes` argument.
            % The output meshes are read from the optional `filename` file.
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/bdry_node_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_bdry_node_2d(meshes);
            %   passed = ut.test_bdry_node_2d(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/bdry_node_2d.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;

            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;

            % Initialize a cell array to store the mesh objects
            bdry = c_bdry();
            bdry_arg = [bdry_type.bdry_transmitted, bdry_type.bdry_transmitted];

            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1;
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                ncell_ext = ncell + 2 * nbdry;  % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;

                % Function call
                bdry.bdry_node_2d(meshes{i}, bdry_arg, bdry_arg);

                % Read fields after function call
                test_group_name = sprintf('/mesh_%d', i - 1);

                vel_after = zeros(nnode_ext(2), nnode_ext(1), dim);
                vel_after(:,:,1) = obj.h5read_transposed(filename, [test_group_name, '/vel_x']);
                vel_after(:,:,2) = obj.h5read_transposed(filename, [test_group_name, '/vel_y']);

                vav_after = zeros(nnode_ext(2), nnode_ext(1), dim);
                vav_after(:,:,1) = obj.h5read_transposed(filename, [test_group_name, '/vav_x']);
                vav_after(:,:,2) = obj.h5read_transposed(filename, [test_group_name, '/vav_y']);

                % Check
                nfailed = failed + 1;
                if ~isequal(meshes{i}.vel_for_2dnode, vel_after)
                    if verbose, fprintf("%s: vel_for_2dnode differs\n", test_group_name), end
                    failed = nfailed;
                end

                if ~isequal(meshes{i}.vav_for_2dnode, vav_after)
                    if verbose, fprintf("%s: vav_for_2dnode differs\n", test_group_name), end
                    failed = nfailed;
                end
                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_bdry_cell_1var_2d(obj, meshes, varargin)
            % test_bdry_cell_1var_2d  Runs a test for bdry.bdry_cell_1var_2d
            %
            % Test Matlab function bdry_cell_1var_2d by comparing results to C output, saved
            % in HDF5 file. The input meshes are provided by the `meshes` argument.
            % The output meshes are read from the optional `filename` file.
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/bdry_cell_1var_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_bdry_cell_2d(meshes);
            %   passed = ut.test_bdry_cell_2d(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/bdry_cell_1var_2d.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;

            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;

            % Initialize a cell array to store the mesh objects
            bdry = c_bdry();
            bdry_arg = [bdry_type.bdry_transmitted, bdry_type.bdry_transmitted];

            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1;
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end


            for i = imn:imx
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                ncell_ext = ncell + 2 * nbdry;  % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;


                % Function call
                bdry.nmat_for_cell = meshes{i}.nmat_for_2dcell;
                bdry.bdry_cell_1var_2d(meshes{i}, bdry_arg, bdry_arg);

                % Read fields after function call
                test_group_name = sprintf('/mesh_%d', i - 1);
                nmat_after  = obj.h5read_transposed(filename, [test_group_name, '/nmat']);

                % Check
                nfailed = failed + 1;
                if ~isequal(bdry.nmat_for_cell, nmat_after)
                    if verbose, fprintf("%s: nmat_for_cell differs\n", test_group_name), end
                    failed = nfailed;
                end

                if failed < nfailed
                    passed = passed + 1;
                end
            end
        end

        function [total, passed] = test_bdry_cell_ragged_2d(obj, meshes, varargin)
            % test_bdry_cell_ragged_2d   Runs a test for bdry.bdry_cell_ragged_2d
            %
            % Test Matlab function bdry_cell_ragged_2d by comparing results to C output, saved
            % in HDF5 file. The input meshes are provided by the `meshes` argument.
            % The output meshes are read from the optional `filename` file.
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/bdry_cell_ragged_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_bdry_cell_2d(meshes);
            %   passed = ut.test_bdry_cell_2d(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/bdry_cell_ragged_2d.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;

            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;

            % Initialize a cell array to store the mesh objects
            bdry = c_bdry();
            bdry_arg = [bdry_type.bdry_transmitted, bdry_type.bdry_transmitted];

            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1;
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end


            for i = imn:imx
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                ncell_ext = ncell + 2 * nbdry;  % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;
                test_group_name = sprintf('/mesh_%d', i - 1);

                % Function call
                bdry.nmat_for_cell = meshes{i}.nmat_for_2dcell;
                bdry.bdry_cell_1var_2d(meshes{i}, bdry_arg, bdry_arg);

                bdry.matid_for_cell = cell(ncell_ext(2), ncell_ext(1));
                bdry.vol_for_cell  = cell(ncell_ext(2), ncell_ext(1));
                bdry.mass_for_cell = cell(ncell_ext(2), ncell_ext(1));
                bdry.ener_for_cell = cell(ncell_ext(2), ncell_ext(1));

                matid  = h5read(filename, [test_group_name, '/matid_before']);
                vol    = h5read(filename, [test_group_name, '/vol_before']);
                mass   = h5read(filename, [test_group_name, '/mass_before']);
                ener   = h5read(filename, [test_group_name, '/ener_before']);
                clen = size(matid, 1);
                
                expected = {};
                expected.matid  = h5read(filename, [test_group_name, '/matid_after']);
                expected.vol    = h5read(filename, [test_group_name, '/vol_after']);
                expected.mass   = h5read(filename, [test_group_name, '/mass_after']);
                expected.ener   = h5read(filename, [test_group_name, '/ener_after']);
                
                offset = 0;
                for j=1:ncell_ext(2)
                    for k=1:ncell_ext(1)
                        nm = bdry.nmat_for_cell(j, k);
                        bdry.matid_for_cell{j, k} = matid(offset+1:offset+nm);
                        bdry.vol_for_cell{j, k}  = vol(offset+1:offset+nm);
                        bdry.ener_for_cell{j, k} = ener(offset+1:offset+nm);
                        bdry.mass_for_cell{j, k} = mass(offset+1:offset+nm);
                        offset = offset + nm;
                    end
                end

                % Function call
                bdry.bdry_cell_ragged_2d(meshes{i}, bdry_arg, bdry_arg);

                % Reformat outputs
                actual = {};
                actual.matid = zeros(clen, 1);
                actual.vol   = zeros(clen, 1);
                actual.mass  = zeros(clen, 1);
                actual.ener  = zeros(clen, 1);
                                
                offset = 0;
                for j=1:ncell_ext(2)
                    for k=1:ncell_ext(1)
                        nm = bdry.nmat_for_cell(j, k);
                        actual.matid(offset+1:offset+nm) = bdry.matid_for_cell{j, k}(1:nm);
                        actual.vol(offset+1:offset+nm)   = bdry.vol_for_cell{j, k}(1:nm);
                        actual.mass(offset+1:offset+nm)  = bdry.mass_for_cell{j, k}(1:nm);
                        actual.ener(offset+1:offset+nm)  = bdry.ener_for_cell{j, k}(1:nm);
                        offset = offset + nm;
                    end
                end

                % Check
                nfailed = failed + 1;
                if ~isequal(expected.matid, actual.matid)
                    if verbose, fprintf("%s: matid_for_cell differs\n", test_group_name), end
                    failed = nfailed;
                end

                if ~all(abs(expected.vol - actual.vol) < tol)
                    if verbose, fprintf("%s: vol_for_cell differs\n", test_group_name), end
                    failed = nfailed;
                end

                if ~all(abs(expected.mass - actual.mass) < tol)
                    if verbose, fprintf("%s: mass_for_cell differs\n", test_group_name), end
                    failed = nfailed;
                end

                if ~all(abs(expected.ener - actual.ener) < tol)
                    if verbose, fprintf("%s: ener_for_cell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if failed < nfailed
                    passed = passed + 1;
                end
            end
        end

        function [total, passed] = test_find_interface2d(obj, varargin)
            % test_find_interface2d  Runs a test for vof2d.find_interface2d
            %
            % Test Matlab function vof2d.find_interface2d comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/find_interface2d.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_find_interface2d();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/find_interface2d.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;

            vof2d = c_vof2d();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            for i = 1:total
                group_name = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                coords   = h5read(filename, [group_name, '/coords']);
                nodelist = h5read(filename, [group_name, '/nodelist']);
                nodelist = nodelist + 1;
                norm     = h5read(filename, [group_name, '/norm']);
                node_loc = h5read(filename, [group_name, '/node_loc']);
                distance = h5read(filename, [group_name, '/distance']);
                nnode = length(nodelist);

                expected = {};
                expected.coords_new = h5read(filename, [group_name, '/coords_new']);
                expected.nodes_interface = h5read(filename, [group_name, '/nodes_interface']);
                expected.nodelist_upper = h5read(filename, [group_name, '/nodelist_upper']);
                expected.nodelist_lower = h5read(filename, [group_name, '/nodelist_lower']);

                % adjust index and transpose the expected nodelists
                expected.nodes_interface = expected.nodes_interface + 1;
                expected.nodelist_upper = expected.nodelist_upper + 1;
                expected.nodelist_lower = expected.nodelist_lower + 1;

                actual = {};
                [~, actual.coords_new, ~, actual.nodes_interface, ~, ...
                    actual.nodelist_lower, ~, actual.nodelist_upper] = ...
                     vof2d.find_interface2d(nnode, coords, nodelist, norm, distance, node_loc);

                % Check
                nfailed = failed + 1;
                if length(actual.coords_new) ~= length(expected.coords_new)
                    if verbose, fprintf("%s: coords_new length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~all(abs(actual.coords_new - expected.coords_new) < tol)
                        if verbose, fprintf("%s: coords_new differs\n", group_name), end
                        failed = nfailed;
                    end
                end

                if length(actual.nodes_interface) ~= length(expected.nodes_interface)
                    if verbose, fprintf("%s: nodes_interface length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~isempty(actual.nodes_interface) ...
                    && ~isequal(actual.nodes_interface, expected.nodes_interface)
                        if verbose, fprintf("%s: nodes_interface differ\n", group_name), end
                        failed = nfailed;
                    end
                end

                if length(actual.nodelist_upper) ~= length(expected.nodelist_upper)
                    if verbose, fprintf("%s: nodelist_upper length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~isempty(actual.nodelist_upper) ...
                    && ~isequal(actual.nodelist_upper, expected.nodelist_upper)
                        if verbose, fprintf("%s: nodelist_upper differ\n", group_name), end
                        failed = nfailed;
                    end
                end

                if length(actual.nodelist_lower) ~= length(expected.nodelist_lower)
                    if verbose, fprintf("%s: nodelist_lower length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~isempty(actual.nodelist_lower) ...
                    && ~isequal(actual.nodelist_lower, expected.nodelist_lower)
                        if verbose, fprintf("%s: nodelist_lower differ\n", group_name), end
                        failed = nfailed;
                    end
                end
                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_cal_poly_area(obj, varargin)
            % test_cal_poly_area  Runs a test for util.cal_poly_area
            %
            % Test Matlab function util.cal_poly_area by comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/find_interface2d.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_cal_poly_area();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/cal_poly_area.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;

            util = c_util();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            for i = 1:total
                group_name = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                coords   = h5read(filename, [group_name, '/coords']);
                nodelist = h5read(filename, [group_name, '/nodelist']);
                nodelist = nodelist + 1;
                nnode = length(nodelist);

                expected = {};
                expected.area = h5read(filename, [group_name, '/area']);

                actual = {};
                actual.area = util.cal_poly_area(nnode, coords, nnode, nodelist);

                % Check
                nfailed = failed + 1;
                if ~all(abs(actual.area - expected.area) < tol)
                    if verbose, fprintf("%s: area differs\n", group_name), end
                    failed = nfailed;
                end

                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_rz_area(obj, varargin)
            % test_rz_area  Runs a test for util.rz_area
            %
            % Test Matlab function util.rz_area by comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/find_interface2d.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_rz_area();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/rz_area.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;

            util = c_util();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            for i = 1:total
                group_name = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                rz   = h5read(filename, [group_name, '/rz']);
                nn = length(rz)/2;

                expected = {};
                expected.vol = h5read(filename, [group_name, '/vol']);

                actual = {};
                actual.vol = util.rz_area(nn, rz);

                % Check
                nfailed = failed + 1;
                if ~(abs(actual.vol - expected.vol) < tol)
                    if verbose, fprintf("%s: vol differs\n", group_name), end
                    failed = nfailed;
                end

                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_bounds_2d(obj, varargin)
            % test_bounds_2d  Runs a test for vof2d.bounds_2d
            %
            % Test Matlab function vof2d.bounds_2d comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/bounds_2d.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_bounds_2d();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/bounds_2d.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'testnum', 0, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            testnum = p.Results.testnum;

            global util vof2d;
            vof2d = c_vof2d();
            util = c_util();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            if testnum == 0
                imn = 1;
                imx = total;
            else
                imn = testnum;
                imx = testnum;
                total = 1;
            end

            for i = imn:imx
                group_name = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                geop        = h5read(filename, [group_name, '/geop']);
                vf_to_match = h5read(filename, [group_name, '/vf_to_match']);
                volume      = h5read(filename, [group_name, '/volume']);
                norm        = h5read(filename, [group_name, '/norm']);
                coords   = h5read(filename, [group_name, '/coords']);
                nodelist = h5read(filename, [group_name, '/nodelist']);
                nodelist = nodelist + 1;
                nnode = length(nodelist);
                node_order_for_ds = h5read(filename, [group_name, '/node_order_for_ds']);
                node_order_for_ds = node_order_for_ds + 1;
                ds_ea_node = h5read(filename, [group_name, '/ds_ea_node']);

                expected = {};
                expected.ds_lower = h5read(filename, [group_name, '/ds_lower']);
                expected.ds_upper = h5read(filename, [group_name, '/ds_upper']);
                expected.vf_lower = h5read(filename, [group_name, '/vf_lower']);
                expected.vf_upper = h5read(filename, [group_name, '/vf_upper']);
                expected.coords_new = h5read(filename, [group_name, '/coords_new']);
                expected.vol_matched = h5read(filename, [group_name, '/vol_matched']);
                expected.nodes_interface = h5read(filename, [group_name, '/nodes_interface']);
                expected.nodelist_upper = h5read(filename, [group_name, '/nodelist_upper']);
                expected.nodelist_lower = h5read(filename, [group_name, '/nodelist_lower']);

                % adjust index and transpose the expected nodelists
                expected.nodes_interface = expected.nodes_interface + 1;
                expected.nodelist_upper = expected.nodelist_upper + 1;
                expected.nodelist_lower = expected.nodelist_lower + 1;
                
                actual = {};
                [actual.ds_lower, actual.ds_upper, actual.vf_lower, actual.vf_upper, ...
                  actual.nnode_new, actual.coords_new, actual.vol_matched, ...
                  actual.nnode_interface, actual.nodes_interface, actual.nnode_lower, ...
                  actual.nodelist_lower, actual.nnode_upper, actual.nodelist_upper] = ...
                  vof2d.bounds_2d(geop, vf_to_match, volume, norm, nnode, ...
                            coords, nodelist, node_order_for_ds, ds_ea_node);

                % Check
                nfailed = failed + 1;
                if ~(abs(actual.ds_lower - expected.ds_lower) < tol)
                    if verbose, fprintf("%s: ds_lower differs\n", group_name), end
                    failed = nfailed;
                end

                if ~(abs(actual.ds_upper - expected.ds_upper) < tol)
                    if verbose, fprintf("%s: ds_upper differs\n", group_name), end
                    failed = nfailed;
                end

                if ~(abs(actual.vf_lower - expected.vf_lower) < tol)
                    if verbose, fprintf("%s: vf_lower differs\n", group_name), end
                    failed = nfailed;
                end

                if ~(abs(actual.vf_upper - expected.vf_upper) < tol)
                    if verbose, fprintf("%s: vf_upper differs\n", group_name), end
                    failed = nfailed;
                end

                if length(actual.coords_new) ~= length(expected.coords_new)
                    if verbose, fprintf("%s: coords_new length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~all(abs(actual.coords_new - expected.coords_new) < tol)
                        if verbose, fprintf("%s: coords_new differs\n", group_name), end
                        failed = nfailed;
                    end
                end

                if ~isequal(actual.vol_matched, expected.vol_matched)
                    if verbose, fprintf("%s: vol_matched differs\n", group_name), end
                    failed = nfailed;
                end

                if length(actual.nodes_interface) ~= length(expected.nodes_interface)
                    if verbose, fprintf("%s: nodes_interface length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~isequal(actual.nodes_interface, expected.nodes_interface)
                        if verbose, fprintf("%s: nodes_interface differ\n", group_name), end
                        failed = nfailed;
                    end
                end

                if length(actual.nodelist_upper) ~= length(expected.nodelist_upper)
                    if verbose, fprintf("%s: nodelist_upper length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~isequal(actual.nodelist_upper, expected.nodelist_upper)
                        if verbose, fprintf("%s: nodelist_upper differ\n", group_name), end
                        failed = nfailed;
                    end
                end

                if length(actual.nodelist_lower) ~= length(expected.nodelist_lower)
                    if verbose, fprintf("%s: nodelist_lower length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~isequal(actual.nodelist_lower, expected.nodelist_lower)
                        if verbose, fprintf("%s: nodelist_lower differ\n", group_name), end
                        failed = nfailed;
                    end
                end

                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_order_nodes_along_norm(obj, varargin)
            % test_order_nodes_along_norm  Runs a test for util.order_nodes_along_norm
            %
            % Test Matlab function util.order_nodes_along_norm comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/order_nodes_along_norm.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_order_nodes_along_norm();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/order_nodes_along_norm.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'testnum', 0, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            testnum = p.Results.testnum;

            global util;
            util = c_util();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            if testnum == 0
                imn = 1;
                imx = total;
            else
                imn = testnum;
                imx = testnum;
                total = 1;
            end

            for i = imn:imx
                group_name = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                dim        = h5read(filename, [group_name, '/dim']);
                norm        = h5read(filename, [group_name, '/norm']);
                nnode = h5read(filename, [group_name, '/nnode']);
                coords   = h5read(filename, [group_name, '/coords']);

                expected = {};
                expected.node_order = h5read(filename, [group_name, '/node_order']);
                expected.node_order = expected.node_order + 1;
                expected.ds_ea_node = h5read(filename, [group_name, '/ds_ea_node']);

                actual = {};
                [actual.node_order, actual.ds_ea_node] = ...
                    util.order_nodes_along_norm(dim, norm, nnode, coords);

                % Check
                nfailed = failed + 1;
                if ~isequal(actual.node_order, expected.node_order)
                    if verbose, fprintf("%s: node_order differs\n", group_name), end
                    failed = nfailed;
                end

                if ~(abs(actual.ds_ea_node - expected.ds_ea_node) < tol)
                    if verbose, fprintf("%s: ds_ea_node differs\n", group_name), end
                    failed = nfailed;
                end

                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_cal_distance2d(obj, varargin)
            % test_bounds_2d  Runs a test for vof2d.cal_distance2d
            %
            % Test Matlab function vof2d.cal_distance2d comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/cal_distance2d.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_cal_distance2d();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/cal_distance2d.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'testnum', 0, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            testnum = p.Results.testnum;

            global util vof2d;
            vof2d = c_vof2d();
            util = c_util();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            if testnum == 0
                imn = 1;
                imx = total;
            else
                imn = testnum;
                imx = testnum;
                total = 1;
            end

            for i = imn:imx
                group_name = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                geop        = h5read(filename, [group_name, '/geop']);
                vf_to_match = h5read(filename, [group_name, '/vf_to_match']);
                volume      = h5read(filename, [group_name, '/volume']);
                norm        = h5read(filename, [group_name, '/norm']);
                coords   = h5read(filename, [group_name, '/coords']);
                nodelist = h5read(filename, [group_name, '/nodelist']);
                nodelist = nodelist + 1;
                nnode = length(nodelist);

                expected = {};
                expected.coords_new = h5read(filename, [group_name, '/coords_new']);
                expected.distance = h5read(filename, [group_name, '/distance']);
                expected.nodes_interface = h5read(filename, [group_name, '/nodes_interface']);
                expected.nodelist_upper = h5read(filename, [group_name, '/nodelist_upper']);
                expected.nodelist_lower = h5read(filename, [group_name, '/nodelist_lower']);

                % adjust index and transpose the expected nodelists
                expected.nodes_interface = expected.nodes_interface + 1;
                expected.nodelist_upper = expected.nodelist_upper + 1;
                expected.nodelist_lower = expected.nodelist_lower + 1;
                
                actual = {};
                [actual.distance, actual.nnode_new, actual.coords_new, ...
                    actual.nnode_interface, actual.nodes_interface, actual.nnode_lower, ...
                    actual.nodelist_lower, actual.nnode_upper, actual.nodelist_upper] = ...
                    vof2d.cal_distance2d(geop, vf_to_match, volume, norm, ...
                            nnode, coords, nodelist);
                  

                % Check
                nfailed = failed + 1;
                if length(actual.coords_new) ~= length(expected.coords_new)
                    if verbose, fprintf("%s: coords_new length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~all(abs(actual.coords_new - expected.coords_new) < tol)
                        if verbose, fprintf("%s: coords_new differs\n", group_name), end
                        failed = nfailed;
                    end
                end

                if ~(abs(actual.distance - expected.distance) < tol)
                    if verbose, fprintf("%s: ds_lower differs\n", group_name), end
                    failed = nfailed;
                end


                if length(actual.nodes_interface) ~= length(expected.nodes_interface)
                    if verbose, fprintf("%s: nodes_interface length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~isequal(actual.nodes_interface, expected.nodes_interface)
                        if verbose, fprintf("%s: nodes_interface differ\n", group_name), end
                        failed = nfailed;
                    end
                end

                if length(actual.nodelist_upper) ~= length(expected.nodelist_upper)
                    if verbose, fprintf("%s: nodelist_upper length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~isequal(actual.nodelist_upper, expected.nodelist_upper)
                        if verbose, fprintf("%s: nodelist_upper differ\n", group_name), end
                        failed = nfailed;
                    end
                end

                if length(actual.nodelist_lower) ~= length(expected.nodelist_lower)
                    if verbose, fprintf("%s: nodelist_lower length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~isequal(actual.nodelist_lower, expected.nodelist_lower)
                        if verbose, fprintf("%s: nodelist_lower differ\n", group_name), end
                        failed = nfailed;
                    end
                end

                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_cal_cell_zgrad2d(obj, varargin)
            % test_cal_cell_zgrad2d  Runs a test for mesh.cal_cell_zgrad2d
            %
            % Test Matlab function util.rz_area by comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/find_interface2d.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   [tot, passd] = ut.test_cal_cell_zgrad2d();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/cal_cell_zgrad2d.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;

            mesh = c_mesh();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            for i = 1:total
                group_name = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                dx   = h5read(filename, [group_name, '/dx']);
                var9 = h5read(filename, [group_name, '/var']);
                var = reshape(var9, [3, 3]);
                
                expected = {};
                expected.grad = h5read(filename, [group_name, '/grad'])';

                actual = {};
                actual.grad = mesh.cal_cell_zgrad2d(dx, var);

                % Check
                nfailed = failed + 1;
                if ~all(abs(actual.grad - expected.grad) < tol)
                    if verbose, fprintf("%s: grad differs\n", group_name), end
                    failed = nfailed;
                end

                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end
        
        function [total, passed] = test_compute_divu(obj, meshes, varargin)
            % test_compute_divu  Runs a test for update.compute_divu
            %
            % Test Matlab function update.compute_divu by comparing results
            % to C output, saved in HDF5 file.
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/bdry_node_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_compute_divu(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/compute_divu.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;

            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;

            % Initialize a cell array to store the mesh objects
            upd = c_update();

            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1;
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                upd.divu_for_2dcell = [];        % Reset divu_for_2dcell
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                ncell_ext = ncell + 2 * nbdry;  % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;

                % Function call
                upd.compute_divu(meshes{i});

                % Read fields after function call
                test_group_name = sprintf('/mesh_%d', i - 1);

                vel_after = zeros(nnode_ext(2), nnode_ext(1), dim);
                vel_after(:,:,1) = obj.h5read_transposed(filename, [test_group_name, '/vel_x']);
                vel_after(:,:,2) = obj.h5read_transposed(filename, [test_group_name, '/vel_y']);

                divu_after = obj.h5read_transposed(filename, [test_group_name, '/divu']);

                % Check
                nfailed = failed + 1;
                if ~isequal(meshes{i}.vel_for_2dnode, vel_after)
                    if verbose, fprintf("%s: vel_for_2dnode differs\n", test_group_name), end
                    failed = nfailed;
                end

                if ~all(abs(upd.divu_for_2dcell - divu_after) < tol)
                    if verbose, fprintf("%s: divu_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_compute_qvis(obj, meshes, varargin)
            % test_compute_qvis  Runs a test for update.compute_qvis
            %
            % Test Matlab function update.compute_qvis by comparing results 
            % to C output, saved in HDF5 file. 
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to 
            %                (default: '[obj.testdir]/bdry_node_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all 
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_compute_qvis(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/compute_qvis.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;
        
            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;
        
            % Initialize a cell array to store the mesh objects
            upd = c_update();
       
            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1; 
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                upd.divu_for_2dcell = [];        % Reset divu_for_2dcell
                upd.qvis_for_2dcell = [];        % Reset qvis_for_2dcell
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                ncell_ext = ncell + 2 * nbdry;  % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;

                % Function call
                upd.compute_divu(meshes{i});
                upd.compute_qvis(meshes{i});

                % Read fields after function call
                test_group_name = sprintf('/mesh_%d', i - 1);

                vel_after = zeros(nnode_ext(2), nnode_ext(1), dim);
                vel_after(:,:,1) = obj.h5read_transposed(filename, [test_group_name, '/vel_x']);
                vel_after(:,:,2) = obj.h5read_transposed(filename, [test_group_name, '/vel_y']);
                divu_after = obj.h5read_transposed(filename, [test_group_name, '/divu']);
                qvis_after = obj.h5read_transposed(filename, [test_group_name, '/qvis']);
                

                % Check
                nfailed = failed + 1;
                if ~isequal(meshes{i}.vel_for_2dnode, vel_after)
                    fprintf("%s: vel_for_2dnode differs\n", test_group_name);
                    failed = nfailed;
                end
                if ~all(abs(upd.divu_for_2dcell - divu_after) < tol)
                    if verbose, fprintf("%s: divu_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(upd.qvis_for_2dcell - qvis_after) < tol)
                    if verbose, fprintf("%s: qvis_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_compute_force(obj, meshes, varargin)
            % test_compute_force  Runs a test for update.compute_force
            %
            % Test Matlab function update.compute_force by comparing results 
            % to C output, saved in HDF5 file. 
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to 
            %                (default: '[obj.testdir]/bdry_node_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all 
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_compute_force(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %

            global tiny;
            tiny = 1.0e-30;

            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/compute_force.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;
        
            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;
        
            % Initialize a cell array to store the mesh objects
            upd = c_update();
       
            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1; 
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                % skip cases with nbdry < 1
                if nbdry < 1
                    total = total - 1;
                    continue;
                end
                upd.divu_for_2dcell = [];        % Reset divu_for_2dcell
                upd.qvis_for_2dcell = [];        % Reset qvis_for_2dcell
                upd.force_for_2dnode = [];       % Reset force_for_2dnode
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                ncell_ext = ncell + 2 * nbdry;   % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;


                % Read in direction 
                test_group_name = sprintf('/mesh_%d', i - 1);
                direction = obj.h5read_transposed(filename, ...
                    [test_group_name, '/direction']) + 1; % +1 to account for Matlab 1 indexing 

                % Function call
                upd.compute_divu(meshes{i});
                upd.compute_qvis(meshes{i});
                upd.compute_force(meshes{i}, direction);

                % Read fields after function call
                vel_after = zeros(nnode_ext(2), nnode_ext(1), dim);
                vel_after(:,:,1) = obj.h5read_transposed(filename, [test_group_name, '/vel_x']);
                vel_after(:,:,2) = obj.h5read_transposed(filename, [test_group_name, '/vel_y']);
                                
                divu_after = obj.h5read_transposed(filename, [test_group_name, '/divu']);
                qvis_after = obj.h5read_transposed(filename, [test_group_name, '/qvis']);
                force_after = obj.h5read_transposed(filename, [test_group_name, '/force']);      

                % Check
                nfailed = failed + 1;
                if ~isequal(meshes{i}.vel_for_2dnode, vel_after)
                    fprintf("%s: vel_for_2dnode differs\n", test_group_name);
                    failed = nfailed;
                end
                if ~all(abs(upd.divu_for_2dcell - divu_after) < tol*norm(divu_after))
                    if verbose, fprintf("%s: divu_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(upd.qvis_for_2dcell - qvis_after) < tol*norm(qvis_after))
                    if verbose, fprintf("%s: qvis_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end

                % for the force, only check the inner nodes; boundary nodes
                % are not computed and may contain junk data
                if ~all(abs(upd.force_for_2dnode(nbdry+1:end-nbdry,nbdry+1:end-nbdry)...
                                   - force_after(nbdry+1:end-nbdry,nbdry+1:end-nbdry))...
                                   < tol*norm(force_after(nbdry+1:end-nbdry,nbdry+1:end-nbdry)))
                    if verbose, fprintf("%s: force_for_2dnode differs\n", test_group_name), end
                    failed = nfailed;
                end
                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_update_vel_comp(obj, meshes, varargin)
            % test_update_vel_comp  Runs a test for update.update_vel_comp
            %
            % Test Matlab function update.update_vel_comp by comparing results 
            % to C output, saved in HDF5 file. 
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to 
            %                (default: '[obj.testdir]/bdry_node_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all 
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_update_vel_comp(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %

            global tiny;
            tiny = 1.0e-30;

            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/update_vel_comp.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;
        
            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;
        
            % Initialize a cell array to store the mesh objects
            upd = c_update();
       
            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1; 
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                % skip cases with nbdry < 1
                if nbdry < 1
                    total = total - 1;
                    continue;
                end
                upd.divu_for_2dcell = [];        % Reset divu_for_2dcell
                upd.qvis_for_2dcell = [];        % Reset qvis_for_2dcell
                upd.force_for_2dnode = [];       % Reset force_for_2dnode
                
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                ncell_ext = ncell + 2 * nbdry;   % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;


                % Read in direction 
                test_group_name = sprintf('/mesh_%d', i - 1);
                direction = obj.h5read_transposed(filename, ...
                    [test_group_name, '/direction']) + 1; % +1 to account for Matlab 1 indexing 
                dt = obj.h5read_transposed(filename, ...
                    [test_group_name, '/dt']); 

                % Function call
                upd.compute_divu(meshes{i});
                upd.compute_qvis(meshes{i});
                upd.compute_force(meshes{i}, direction);
                upd.update_vel_comp(meshes{i}, dt, direction)


                % Read fields after function call
                vel_after = zeros(nnode_ext(2), nnode_ext(1), dim);
                vel_after(:,:,1) = obj.h5read_transposed(filename, [test_group_name, '/vel_x']);
                vel_after(:,:,2) = obj.h5read_transposed(filename, [test_group_name, '/vel_y']);

                vav_after = zeros(nnode_ext(2), nnode_ext(1), dim);
                vav_after(:,:,1) = obj.h5read_transposed(filename, [test_group_name, '/vav_x']);
                vav_after(:,:,2) = obj.h5read_transposed(filename, [test_group_name, '/vav_y']);
                                
                divu_after = obj.h5read_transposed(filename, [test_group_name, '/divu']);
                qvis_after = obj.h5read_transposed(filename, [test_group_name, '/qvis']);
                force_after = obj.h5read_transposed(filename, [test_group_name, '/force']);      

                % Check
                nfailed = failed + 1;
                if ~all(abs(meshes{i}.vel_for_2dnode - vel_after) < tol*norm(vel_after, 'fro'))
                    if verbose, fprintf("%s: vel_for_2dnode differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(meshes{i}.vav_for_2dnode - vav_after) < tol*norm(vav_after, 'fro'))
                    if verbose, fprintf("%s: vav_for_2dnode differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(upd.divu_for_2dcell - divu_after) < tol*norm(divu_after))
                    if verbose, fprintf("%s: divu_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(upd.qvis_for_2dcell - qvis_after) < tol*norm(qvis_after))
                    if verbose, fprintf("%s: qvis_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end


                % for the force, only check the inner nodes; boundary nodes
                % are not computed and may contain junk data
                if ~all(abs(upd.force_for_2dnode(nbdry+1:end-nbdry,nbdry+1:end-nbdry)...
                                   - force_after(nbdry+1:end-nbdry,nbdry+1:end-nbdry))...
                                   < tol*norm(force_after(nbdry+1:end-nbdry,nbdry+1:end-nbdry)))
                    if verbose, fprintf("%s: force_for_2dnode differs\n", test_group_name), end
                    failed = nfailed;
                end

                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end
        
        function [total, passed] = test_update_density(obj, meshes, varargin)
            % test_update_density  Runs a test for update.update_density
            %
            % Test Matlab function update.update_vel_comp by comparing results 
            % to C output, saved in HDF5 file. 
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to 
            %                (default: '[obj.testdir]/bdry_node_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all 
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_update_density(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            
            global tiny;
            tiny = 1.0e-30;
    
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/update_density.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;
        
            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;
        
            % Initialize a cell array to store the mesh objects
            upd = c_update();
       
            % Loop through meshes
            failed = 0;
            passed = 0;
    
            if mesh_index < 1
                imn = 1; 
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                mat = c_mat();  % reset mat for each mesh 
                upd.divu_for_2dcell = [];        % Reset divu_for_2dcell
    
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                ncell_ext = ncell + 2 * nbdry;   % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;
        
                % Read in input data from hdf5  
                test_group_name = sprintf('/mesh_%d', i - 1);
                dt = obj.h5read_transposed(filename, ...
                    [test_group_name, '/dt']); 
    
    
                mat.nmat_in_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/nmat_mixcell']); % 1d array, nmixcell size
                mat.nmixcell = length(mat.nmat_in_mixcell);
                ijk_in_mixcell = double(obj.h5read_transposed(filename, ...
                    [test_group_name, '/ijk_mixcell'])) + 1; % 1d array 
                mat.ijk_in_mixcell = reshape(ijk_in_mixcell, dim, mat.nmixcell)'; 
                rho_mixcell_old = obj.h5read_transposed(filename, ...
                    [test_group_name, '/rho_mixcell_old']); % is input to c update_density routine
                rho_mixcell_new = obj.h5read_transposed(filename, ...
                    [test_group_name, '/rho_mixcell_new']); % is output from c update_density routine
                vf_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/vf_mixcell']);
    
                % Unpack mixcell variables into cell arrays 
                % and store in mat
                mat.vf_in_mixcell = cell(mat.nmixcell, 1); 
                mat.rho_in_mixcell = cell(mat.nmixcell, 1);
                rho_mix_after = cell(mat.nmixcell, 1);
                offset = 1;             
                for mx = 1:mat.nmixcell
                    mat.vf_in_mixcell{mx} = vf_mixcell(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    mat.rho_in_mixcell{mx} = rho_mixcell_old(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    rho_mix_after{mx} = rho_mixcell_new(offset:offset+mat.nmat_in_mixcell(mx)-1); % for comparison 
                    offset = offset+mat.nmat_in_mixcell(mx);
                end            
        
                % Function call
                upd.compute_divu(meshes{i});
                cournt = upd.update_density(meshes{i}, mat, dt); 
        
                % Read fields after function call
                vel_after = zeros(nnode_ext(2), nnode_ext(1), dim);
                vel_after(:,:,1) = obj.h5read_transposed(filename, [test_group_name, '/vel_x']);
                vel_after(:,:,2) = obj.h5read_transposed(filename, [test_group_name, '/vel_y']);                            
                divu_after = obj.h5read_transposed(filename, [test_group_name, '/divu']);
                rho_after = obj.h5read_transposed(filename, [test_group_name, '/rho']);
                cournt_after = obj.h5read_transposed(filename, [test_group_name, '/courant']);
      
        
                % Check
                nfailed = failed + 1;
                if ~all(abs(meshes{i}.vel_for_2dnode - vel_after) < tol*norm(vel_after, 'fro'))
                    if verbose, fprintf("%s: vel_for_2dnode differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(upd.divu_for_2dcell - divu_after) < tol*norm(divu_after))
                    if verbose, fprintf("%s: divu_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(meshes{i}.rho_for_2dcell - rho_after) < tol*norm(rho_after))
                    if verbose, fprintf("%s: rho_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(cell2mat(cellfun(@(x,y) abs(x-y)< tol*norm(x), rho_mix_after, mat.rho_in_mixcell, 'UniformOutput', false)))
                    if verbose, fprintf("%s: rho_in_mixcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if (abs(cournt - cournt_after) > tol)
                    if verbose, fprintf("%s: Courant number differs \n", test_group_name), end
                    failed = nfailed; 
                end        
                if failed < nfailed
                    passed = passed + 1;
                end
            end
        end
        
        function [total, passed] = test_update_energy(obj, meshes, varargin)
            % test_update_energy  Runs a test for update.update_energy
            %
            % Test Matlab function update.update_vel_comp by comparing results 
            % to C output, saved in HDF5 file. 
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to 
            %                (default: '[obj.testdir]/bdry_node_2d.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all 
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_update_energy(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            
            global tiny;
            tiny = 1.0e-30;
    
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/update_energy.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;
        
            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;
        
            % Initialize a cell array to store the mesh objects
            upd = c_update();
       
            % Loop through meshes
            failed = 0;
            passed = 0;
    
            if mesh_index < 1
                imn = 1; 
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                mat = c_mat();  % reset mat for each mesh 
                upd.divu_for_2dcell = [];        % Reset divu_for_2dcell
    
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                ncell_ext = ncell + 2 * nbdry;   % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;
        
                % Read in input data from hdf5  
                test_group_name = sprintf('/mesh_%d', i - 1);
                dt = obj.h5read_transposed(filename, ...
                    [test_group_name, '/dt']); 
    
    
                mat.nmat_in_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/nmat_mixcell']); % 1d array, nmixcell size
                mat.nmixcell = length(mat.nmat_in_mixcell);
                ijk_in_mixcell = double(obj.h5read_transposed(filename, ...
                    [test_group_name, '/ijk_mixcell'])) + 1; % 1d array 
                mat.ijk_in_mixcell = reshape(ijk_in_mixcell, dim, mat.nmixcell)'; 
                rho_mixcell_old = obj.h5read_transposed(filename, ...
                    [test_group_name, '/rho_mixcell_old']);
                rho_mixcell_new = obj.h5read_transposed(filename, ...
                    [test_group_name, '/rho_mixcell_new']);
                vf_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/vf_mixcell']);
                es_mixcell_old = obj.h5read_transposed(filename, ...
                    [test_group_name, '/es_mixcell_old']);
                pres_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/pres_mixcell']);
                ei_mixcell_old = obj.h5read_transposed(filename, ...
                    [test_group_name, '/ei_mixcell_old']);
                ei_mixcell_new = obj.h5read_transposed(filename, ...
                    [test_group_name, '/ei_mixcell']);
                meshes{i}.es_for_2dcell_old = obj.h5read_transposed(filename, ...
                    [test_group_name, '/es']);
                meshes{i}.rho_for_2dcell_old = meshes{i}.rho_for_2dcell;
    
    
                % Unpack mixcell variables into cell arrays 
                % and store in mat
                mat.vf_in_mixcell = cell(mat.nmixcell, 1); 
                mat.rho_in_mixcell = cell(mat.nmixcell, 1);
                rho_in_mixcell_old = cell(mat.nmixcell, 1);
                es_in_mixcell_old = cell(mat.nmixcell, 1);
                mat.pres_in_mixcell = cell(mat.nmixcell, 1);
                mat.ei_in_mixcell = cell(mat.nmixcell, 1);
                ei_mix_after = cell(mat.nmixcell, 1);
                rho_mix_after = cell(mat.nmixcell, 1);
                offset = 1;             
                for mx = 1:mat.nmixcell
                    mat.vf_in_mixcell{mx} = vf_mixcell(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    mat.rho_in_mixcell{mx} = rho_mixcell_old(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    mat.pres_in_mixcell{mx} = pres_mixcell(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    mat.ei_in_mixcell{mx} = ei_mixcell_old(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    rho_in_mixcell_old{mx} = rho_mixcell_old(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    es_in_mixcell_old{mx} = es_mixcell_old(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    rho_mix_after{mx} = rho_mixcell_new(offset:offset+mat.nmat_in_mixcell(mx)-1); % for comparison 
                    ei_mix_after{mx} = ei_mixcell_new(offset:offset+mat.nmat_in_mixcell(mx)-1); % for comparison 
                    offset = offset+mat.nmat_in_mixcell(mx);
                end            
        
                % Function call
                upd.compute_divu(meshes{i});
                cournt = upd.update_density(meshes{i}, mat, dt); 
                upd.update_energy(meshes{i}, mat, dt, rho_in_mixcell_old, es_in_mixcell_old);
        
                % Read fields after function call
                vel_after = zeros(nnode_ext(2), nnode_ext(1), dim);
                vel_after(:,:,1) = obj.h5read_transposed(filename, [test_group_name, '/vel_x']);
                vel_after(:,:,2) = obj.h5read_transposed(filename, [test_group_name, '/vel_y']);                            
                divu_after = obj.h5read_transposed(filename, [test_group_name, '/divu']);
                rho_after = obj.h5read_transposed(filename, [test_group_name, '/rho']);
                cournt_after = obj.h5read_transposed(filename, [test_group_name, '/courant']);
                ei_after = obj.h5read_transposed(filename, [test_group_name, '/ei']);
      
        
                % Check
                nfailed = failed + 1;
                if ~all(abs(meshes{i}.vel_for_2dnode - vel_after) < tol*norm(vel_after, 'fro'))
                    if verbose, fprintf("%s: vel_for_2dnode differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(upd.divu_for_2dcell - divu_after) < tol*norm(divu_after))
                    if verbose, fprintf("%s: divu_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(meshes{i}.rho_for_2dcell - rho_after) < tol*norm(rho_after))
                    if verbose, fprintf("%s: rho_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(cell2mat(cellfun(@(x,y) abs(x-y)< tol*norm(x), rho_mix_after, mat.rho_in_mixcell, 'UniformOutput', false)))
                    if verbose, fprintf("%s: rho_in_mixcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(abs(meshes{i}.ei_for_2dcell - ei_after) < tol*norm(ei_after))
                    if verbose, fprintf("%s: ei_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if ~all(cell2mat(cellfun(@(x,y) abs(x-y)< tol*norm(x), ei_mix_after, mat.ei_in_mixcell, 'UniformOutput', false)))
                    if verbose, fprintf("%s: ei_in_mixcell differs\n", test_group_name), end
                    failed = nfailed;
                end
                if (abs(cournt - cournt_after) > tol)
                    if verbose, fprintf("%s: Courant number differs \n", test_group_name), end
                    failed = nfailed; 
                end        
        
                if failed < nfailed
                    passed = passed + 1;
                end
        
            end

        end 

        function [total, passed] = test_reconstruct2d_nmat_pagosa(obj, varargin)
            % test_reconstruct2d_nmat_pagosa  Runs a test for vof2d::reconstruct2d_nmat_pagosa
            %
            % Test Matlab function vof2d.test_reconstruct2d_nmat_pagosa comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/reconstruct2d_nmat_pagosa.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   [t, p] = ut.test_reconstruct2d_nmat_pagosa();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/reconstruct2d_nmat_pagosa.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'testnum', 0, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            testnum = p.Results.testnum;

            global util vof2d xdm mesh;
            vof2d = c_vof2d();
            util = c_util();
            xdm = xdmfio();
            mesh = c_mesh();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            if testnum == 0
                imn = 1;
                imx = total;
            else
                imn = testnum;
                imx = testnum;
                total = 1;
            end

            for i = imn:imx
                group_name    = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                geop          = h5read(filename, [group_name, '/geop']);
                xl            = h5read(filename, [group_name, '/xl'])';
                dx            = h5read(filename, [group_name, '/dx']);
                nmat_2dsmesh  = obj.h5read_transposed(filename, [group_name, '/nmat_2dsmesh']);

                matid_raw = h5read(filename, [group_name, '/matid_2dsmesh']);
                matid_raw = reshape(matid_raw, [2,3,3]);
                matid_2dsmesh = cell([3, 3]);
                for i = 1:3
                    for j = 1:3
                        matid_2dsmesh{j, i} = matid_raw(1:nmat_2dsmesh(j,i), i, j);
                    end
                end

                vf_raw = h5read(filename, [group_name, '/vf_2dsmesh']);
                vf_raw = reshape(vf_raw, [2,3,3]);
                vf_2dsmesh = cell([3, 3]);
                for i = 1:3
                    for j = 1:3
                        vf_2dsmesh{j, i} = vf_raw(1:nmat_2dsmesh(j,i), i, j);
                    end
                end
                
                expected = {};
                expected.coords_tot = h5read(filename, [group_name, '/coords_tot']);
                expected.nnode_for_minterface = h5read(filename, [group_name, '/nnode_for_minterface']);
                expected.nodes_for_minterface = h5read(filename, [group_name, '/nodes_for_minterface']);
                expected.nnode_for_mpoly = h5read(filename, [group_name, '/nnode_for_mpoly']);
                expected.nodes_for_mpoly = h5read(filename, [group_name, '/nodes_for_mpoly']);

                % adjust index and transpose the expected nodelists
                expected.nodes_for_mpoly = expected.nodes_for_mpoly + 1;
                expected.nodes_for_minterface = expected.nodes_for_minterface + 1;
                
                actual = {};
                [actual.coords_tot, actual.nnode_tot, actual.nodes_for_minterface, ...
                  actual.nodes_for_mpoly, actual.nnode_for_minterface, actual.nnode_for_mpoly] = ...
                    vof2d.reconstruct2d_nmat_pagosa(geop, xl, dx, nmat_2dsmesh, matid_2dsmesh, vf_2dsmesh);
                actual.coords_tot = actual.coords_tot(1:actual.nnode_tot*2);
                  

                % Check
                nfailed = failed + 1;

                if length(actual.coords_tot) ~= actual.nnode_tot*2
                    if verbose, fprintf("%s: len(coords_tot) =/= nnode_tot\n", group_name), end
                    failed = nfailed;
                elseif length(actual.coords_tot) ~= length(expected.coords_tot)
                    if verbose, fprintf("%s: coords_tot length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~all(abs(actual.coords_tot - expected.coords_tot) < tol)
                        if verbose, fprintf("%s: coords_tot differs\n", group_name), end
                        failed = nfailed;
                    end
                end

                if length(actual.nodes_for_minterface) ~= actual.nnode_for_minterface
                    if verbose, fprintf("%s: len(nodes_for_minterface) =/= nnode_for_minterface\n", group_name), end
                    failed = nfailed;
                elseif length(actual.nodes_for_minterface) ~= length(expected.nodes_for_minterface)
                    if verbose, fprintf("%s: nodes_for_minterface length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~all(actual.nodes_for_minterface' == expected.nodes_for_minterface)
                        if verbose, fprintf("%s: nodes_for_minterface differs\n", group_name), end
                        failed = nfailed;
                    end
                end

                nodes_for_mpoly_arr = xdm.flatten_cell1d(actual.nodes_for_mpoly);
                if numel(nodes_for_mpoly_arr) ~= sum(actual.nnode_for_mpoly)
                    if verbose, fprintf("%s: len(nodes_for_mpoly) =/= sum(nnode_for_mpoly)\n", group_name), end
                    failed = nfailed;
                elseif length(nodes_for_mpoly_arr) ~= length(expected.nodes_for_mpoly)
                    if verbose, fprintf("%s: nodes_for_mpoly length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~all(nodes_for_mpoly_arr' == expected.nodes_for_mpoly)
                        if verbose, fprintf("%s: nodes_for_mpoly differs\n", group_name), end
                        failed = nfailed;
                    end
                end

                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_cal_mixcell_zgrad2d(obj, varargin)
            % test_cal_mixcell_zgrad2d  Runs a test for mat::cal_mixcell_zgrad2d
            %
            % Test Matlab function mat.cal_mixcell_zgrad2d comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/cal_mixcell_zgrad2d.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   [t, p] = ut.test_cal_mixcell_zgrad2d();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/cal_mixcell_zgrad2d.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'testnum', 0, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            testnum = p.Results.testnum;

            global util vof2d xdm;
            vof2d = c_vof2d();
            util = c_util();
            xdm = xdmfio();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            if testnum == 0
                imn = 1;
                imx = total;
            else
                imn = testnum;
                imx = testnum;
                total = 1;
            end

            for i = imn:imx
                group_name    = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                ncell = obj.h5read_transposed(filename, [group_name, '/ncell']);
                nbdry = h5read(filename, [group_name, '/nbdry']);
                xl = obj.h5read_transposed(filename, [group_name, '/xl'])';
                dx = obj.h5read_transposed(filename, [group_name, '/dx']);
                xr = xl + dx.*double(ncell);
                mesh = c_mesh();
                mesh.set_mesh(2, xl, xr, ncell, double(nbdry), 2, [1, 2]);
                mesh.nmat_for_2dcell  = obj.h5read_transposed(filename, [group_name, '/nmat_for_2dcell']);
                mesh.matid_for_2dcell  = obj.h5read_transposed(filename, [group_name, '/matid_for_2dcell']);
                mesh.mixcell_for_2dcell  = obj.h5read_transposed(filename, [group_name, '/mixcell_for_2dcell']);

                mesh.matid_for_2dcell = mesh.matid_for_2dcell + 1;
                mesh.mixcell_for_2dcell = mesh.mixcell_for_2dcell + 1;
                
                mat = c_mat();
                mat.nmixcell = h5read(filename, [group_name, '/nmixcell']);
                mat.nmat_in_mixcell = h5read(filename, [group_name, '/nmat_in_mixcell']);
                matids = h5read(filename, [group_name, '/matids_in_mixcell']);
                vfs = h5read(filename, [group_name, '/vf_in_mixcell']);
                ijk = h5read(filename, [group_name, '/ijk_in_mixcell']);
                vfgrad = h5read(filename, [group_name, '/vfgrad_in_mixcell']);
                
                mat.ijk_in_mixcell    = ijk' + 1;
                mat.matids_in_mixcell = cell(1, mat.nmixcell);
                mat.vf_in_mixcell     = cell(1, mat.nmixcell);
                mat.vfgrad_mixcell    = cell(1, mat.nmixcell);

                expected = {};
                expected.vfgrad_mixcell = cell(1, mat.nmixcell_int);
                for mx = 1:mat.nmixcell
                    mat.matids_in_mixcell{mx} = matids(:,mx) + 1;
                    mat.vf_in_mixcell{mx}     = vfs(:,mx);
                    mat.vfgrad_mixcell{mx}    = zeros(2);
                    expected.vfgrad_mixcell{mx} = reshape(vfgrad(4*mx-3:4*mx), [2, 2])';
                end
                
                mat.cal_mixcell_zgrad2d(mesh);
                 
                % Check
                nfailed = failed + 1;
                 
                err = 0.0;
                for mx = 1:mat.nmixcell
                    err = err + norm(mat.vfgrad_mixcell{mx} - expected.vfgrad_mixcell{mx});
                end
                if err < tol
                    if verbose, fprintf("%s: coords_tot differs\n", group_name), end
                    failed = nfailed;
                end
                
                if failed < nfailed
                    passed = passed + 1;
                end
            end
        end

        function [total, passed] = test_get_mpoly(obj, varargin)
            % test_get_mpoly  Runs a test for mat::get_mpoly
            %
            % Test Matlab function mat.cal_mixcell_zgrad2d comparing results to C output,
            % saved in HDF5 file.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/cal_mixcell_zgrad2d.h5').
            %   verbose    - (logical, default: false).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   [t, p] = ut.test_get_mpoly();
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/get_mpoly.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'testnum', 0, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            testnum = p.Results.testnum;

            global util vof2d xdm;
            vof2d = c_vof2d();
            util = c_util();
            xdm = xdmfio();
            info = h5info(filename);
            total = length(info.Groups);  % Number of groups (tests) in the file
            passed = 0;
            failed = 0;
            tol = 1e-12;

            if testnum == 0
                imn = 1;
                imx = total;
            else
                imn = testnum;
                imx = testnum;
                total = 1;
            end

            for i = imn:imx
                group_name    = info.Groups(i).Name;  % e.g., '/test_0', '/test_1', etc.
                ncell = obj.h5read_transposed(filename, [group_name, '/ncell']);
                nbdry = h5read(filename, [group_name, '/nbdry']);
                xl = obj.h5read_transposed(filename, [group_name, '/xl']);
                dx = obj.h5read_transposed(filename, [group_name, '/dx']);
                xr = xl + dx.*double(ncell);
                mesh = c_mesh();
                mesh.set_mesh(2, xl, xr, ncell, double(nbdry), 2, [1, 2]);
                mesh.dx = dx';
                mesh.nmat_for_2dcell  = obj.h5read_transposed(filename, [group_name, '/nmat_for_2dcell']);
                mesh.matid_for_2dcell  = obj.h5read_transposed(filename, [group_name, '/matid_for_2dcell']);
                mesh.mixcell_for_2dcell  = obj.h5read_transposed(filename, [group_name, '/mixcell_for_2dcell']);

                mesh.matid_for_2dcell = mesh.matid_for_2dcell + 1;
                mesh.mixcell_for_2dcell = mesh.mixcell_for_2dcell + 1;
                
                mat = c_mat();
                mat.nmixcell = h5read(filename, [group_name, '/nmixcell']);
                mat.nmat_in_mixcell = h5read(filename, [group_name, '/nmat_in_mixcell']);
                matids = h5read(filename, [group_name, '/matids_in_mixcell']);
                vfs = h5read(filename, [group_name, '/vf_in_mixcell']);
                ijk = h5read(filename, [group_name, '/ijk_in_mixcell']);
                
                mat.ijk_in_mixcell    = ijk' + 1;
                mat.matids_in_mixcell = cell(1, mat.nmixcell);
                mat.vf_in_mixcell     = cell(1, mat.nmixcell);

                for mx = 1:mat.nmixcell
                    mat.matids_in_mixcell{mx} = matids(:,mx) + 1;
                    mat.vf_in_mixcell{mx}     = vfs(:,mx);
                    mat.vfgrad_mixcell{mx}    = zeros(2);
                end
                
                mat.get_mpoly(mesh);

                expected = {};
                expected.nmixcell_mpoly = h5read(filename, [group_name, '/nmixcell_mpoly']);
                expected.mix_mpoly_to_mix = h5read(filename, [group_name, '/mix_mpoly_to_mix']);
                expected.mix_to_mpoly_mix = h5read(filename, [group_name, '/mix_to_mpoly_mix']);
                
                expected.nnode_for_mixcell = h5read(filename, [group_name, '/nnode_for_mixcell']);
                expected.coords_for_mixcell = h5read(filename, [group_name, '/coords_for_mixcell']);
                
                expected.nnode_for_mpoly_in_mixcell = h5read(filename, [group_name, '/nnode_for_mpoly_in_mixcell']);
                expected.nodes_for_mpoly_in_mixcell = h5read(filename, [group_name, '/nodes_for_mpoly_in_mixcell']);
                expected.nodes_for_minterface_in_mixcell = h5read(filename, [group_name, '/nodes_for_minterface_in_mixcell']);
                
                expected.nnode_for_mixcell_total = h5read(filename, [group_name, '/nnode_for_mixcell_total']);
                expected.nnode_for_mix_mpoly_total = h5read(filename, [group_name, '/nnode_for_mix_mpoly_total']);
                
                 
                % Check
                nfailed = failed + 1;
                 
                err = 0.0;
                if mat.nmixcell_mpoly ~= expected.nmixcell_mpoly
                    if verbose, fprintf("%s: nmixcell_mpoly differs\n", group_name), end
                    failed = nfailed;
                end
                mm = mat.nmixcell_mpoly;
                if ~all(mat.mix_mpoly_to_mix(1:mm) == expected.mix_mpoly_to_mix(1:mm) + 1)
                    if verbose, fprintf("%s: mix_mpoly_to_mix differs\n", group_name), end
                    failed = nfailed;
                end
                if ~all(mat.mix_to_mpoly_mix == expected.mix_to_mpoly_mix + 1)
                    if verbose, fprintf("%s: mix_to_mpoly_mix differs\n", group_name), end
                    failed = nfailed;
                end
                if length(mat.nnode_in_mixcell(1:mm)) ~= length(expected.nnode_for_mixcell)
                    if verbose, fprintf("%s: nnode_in_mixcell length differ\n", group_name), end
                    failed = nfailed;
                else
                    if ~all(mat.nnode_in_mixcell(1:mm) == expected.nnode_for_mixcell)
                        if verbose, fprintf("%s: nnode_in_mixcell differs\n", group_name), end
                        failed = nfailed;
                    end
                end

                offset = 0;
                for mxp = 1:mat.nmixcell_mpoly
                    if length(mat.coords_in_mixcell{mxp}) ~= 2*mat.nnode_in_mixcell(mxp)
                        if verbose, fprintf("%s: coords_in_mixcell have wrong length\n", group_name), end
                        failed = nfailed;
                        break;
                    else
                        i1 = offset + 1; 
                        i2 = offset + 2*mat.nnode_in_mixcell(mxp);
                        if abs(mat.coords_in_mixcell{mxp} - expected.coords_for_mixcell(i1:i2)) > tol
                            if verbose, fprintf("%s: nnode_in_mixcell differs\n", group_name), end
                            failed = nfailed;
                            break;
                        end
                        offset = i2;                       
                    end
                end

                offset = 0;
                for mxp = 1:mat.nmixcell_mpoly
                    mx = mat.mix_mpoly_to_mix(mxp);
                    if length(mat.nnode_for_mpoly_in_mixcell{mxp}) ~= mat.nmat_in_mixcell(mx)
                        if verbose, fprintf("%s: nnode_for_mpoly_in_mixcell has wrong length\n", group_name), end
                        failed = nfailed;
                        break;
                    else
                        i1 = offset + 1; 
                        i2 = offset + mat.nmat_in_mixcell(mx);
                        if ~all(mat.nnode_for_mpoly_in_mixcell{mxp} == expected.nnode_for_mpoly_in_mixcell(i1:i2))
                            if verbose, fprintf("%s: nnode_for_mpoly_in_mixcell differ\n", group_name), end
                            failed = nfailed;
                            break;
                        end
                        offset = i2;                       
                    end
                end

                % if there is a fail at this point, continue
                if failed == nfailed, continue, end

                offset = 0;
                for mxp = 1:mat.nmixcell_mpoly
                    mx = mat.mix_mpoly_to_mix(mxp);
                    for nm = 1:mat.nmat_in_mixcell(mx)
                        if length(mat.nodes_for_mpoly_in_mixcell{mxp}{nm}) ~= ...
                                  mat.nnode_for_mpoly_in_mixcell{mxp}(nm)
                            if verbose, fprintf("%s: nodes_for_mpoly_in_mixcell has wrong length\n", group_name), end
                            failed = nfailed;
                            break;
                        else
                            i1 = offset + 1; 
                            i2 = offset + mat.nnode_for_mpoly_in_mixcell{mxp}(nm);
                            if ~all(mat.nodes_for_mpoly_in_mixcell{mxp}{nm} == ...
                                    expected.nodes_for_mpoly_in_mixcell(i1:i2)' + 1)
                                if verbose, fprintf("%s: nodes_for_mpoly_in_mixcell differ\n", group_name), end
                                failed = nfailed;
                                break;
                            end
                            offset = i2;                       
                        end
                    end
                    if failed == nfailed, break, end
                end

                offset = 0;
                for mxp = 1:mat.nmixcell_mpoly
                    mx = mat.mix_mpoly_to_mix(mxp);
                    for nm = 1:mat.nmat_in_mixcell(mx) - 1
                        if length(mat.nodes_for_minterface_in_mixcell{mxp}) ~= 2
                            if verbose, fprintf("%s: nodes_for_minterface_in_mixcell length =/= 2\n", group_name), end
                            failed = nfailed;
                            break;
                        else
                            i1 = offset + 1;
                            i2 = offset + 2;
                            if ~all(mat.nodes_for_minterface_in_mixcell{mxp} == ...
                                    expected.nodes_for_minterface_in_mixcell(i1:i2)' + 1)
                                if verbose, fprintf("%s: nodes_for_minterface_in_mixcell differ\n", group_name), end
                                failed = nfailed;
                                break;
                            end
                            offset = offset + 2;                       
                        end
                    end
                    if failed == nfailed, break, end
                end

                if failed < nfailed
                    passed = passed + 1;
                end
            end
        end

        function [total, passed] = test_courant_from_cs(obj, meshes, varargin)
            % test_courant_from_cs  Runs a test for mesh.courant_from_cs
            %
            % Test Matlab function mesh.courant_from_cs by comparing results
            % to C output, saved in HDF5 file.
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/courant_from_cs.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_courant_from_cs(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/courant_from_cs.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;

            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;

            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1;
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                ncell_ext = ncell + 2 * nbdry;  % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;
                
                % Read in input data from hdf5  
                test_group_name = sprintf('/mesh_%d', i - 1);
                dt = obj.h5read_transposed(filename, ...
                    [test_group_name, '/dt']); 
                    
                % Function call
                cournt = meshes{i}.courant_from_cs(dt);

                % Read fields after function call
                test_group_name = sprintf('/mesh_%d', i - 1);

                cournt_after = obj.h5read_transposed(filename, [test_group_name, '/courant']);


                % Check
                nfailed = failed + 1;
                if (abs(cournt - cournt_after) > tol)
                    if verbose, fprintf("%s: Courant number differs \n", test_group_name), end
                    failed = nfailed; 
                end    
                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_mesh_sound_speed(obj, meshes, varargin)
            % test_mesh_sound_speed  Runs a test for mesh.mesh_sound_speed
            %
            % Test Matlab function mesh.mesh_sound_speed by comparing results
            % to C output, saved in HDF5 file.
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/mesh_sound_speed.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_mesh_sound_speed(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/mesh_sound_speed.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;

            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;
            
            % The geom object
            global eos tiny 
            eos = c_eos();
            tiny = 1.0e-30;

            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1;
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                ncell_ext = ncell + 2 * nbdry;   % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;
                
                % Read in input data from hdf5  
                test_group_name = sprintf('/mesh_%d', i - 1);
                nmat = obj.h5read_transposed(filename, [test_group_name, '/nmat']); 
                matids = obj.h5read_transposed(filename, [test_group_name, '/matids']) + 1; 
                is_solid = obj.h5read_transposed(filename, [test_group_name, '/issolid']); 
                gamma_ea_mat = obj.h5read_transposed(filename, [test_group_name, '/gamma']);
                
                meshes{i}.nmat_mesh = nmat; 
                meshes{i}.matids_mesh = matids;

                % Function call
                meshes{i}.mesh_sound_speed(nmat, matids, is_solid, gamma_ea_mat);

                % Read fields after function call
                test_group_name = sprintf('/mesh_%d', i - 1);
                cs_after = obj.h5read_transposed(filename, [test_group_name, '/cs']);

                % Check
                nfailed = failed + 1;
                if ~all(abs(meshes{i}.cs_for_2dcell - cs_after) < tol*norm(cs_after))
                    if verbose, fprintf("%s: cs_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end 
                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end

        function [total, passed] = test_mat_sound_speed(obj, meshes, varargin)
            % test_mat_sound_speed  Runs a test for mat.mat_sound_speed
            %
            % Test Matlab function mat.mat_sound_speed by comparing results
            % to C output, saved in HDF5 file.
            %
            % Inputs:
            %   meshes     - (cell array of c_mesh) Input meshes.
            %
            % Optional Inputs (via varargin):
            %   filename   - (string) The filename to save the test results to
            %                (default: '[obj.testdir]/mat_sound_speed.h5').
            %   verbose    - (logical, default: false).
            %   mesh_index - (numeric) Index of a specific mesh to process. If -1, all
            %                meshes will be processed (default: -1).
            %
            % Outputs:
            %   total      - (numeric) Total number of tests.
            %   passed     - (numeric) The number of tests that passed.
            %
            % Example Usage:
            %   ut = unittest();
            %   passed = ut.test_mat_sound_speed(meshes, 'mesh_index', 3, 'filename', 'testdata/output.h5');
            %
            % parse optional arguments
            p = inputParser;
            addOptional(p, 'filename', [obj.testdir, '/mat_sound_speed.h5'], @isstring);
            addOptional(p, 'verbose', false, @islogical);
            addOptional(p, 'mesh_index', -1, @isnumeric);
            parse(p, varargin{:});
            filename = p.Results.filename;
            verbose = p.Results.verbose;
            mesh_index = p.Results.mesh_index;

            info = h5info(filename);  % Get HDF5 file info
            total = length(info.Groups);  % Number of groups (meshes) in the file
            assert(total == length(meshes));
            tol = 1e-12;
            
            % The geom object
            global eos tiny 
            eos = c_eos();
            tiny = 1.0e-30;

            % Loop through meshes
            failed = 0;
            passed = 0;

            if mesh_index < 1
                imn = 1;
                imx = total;
            else
                imn = mesh_index;
                imx = mesh_index;
                total = 1;
            end
            for i = imn:imx
                mat = c_mat();  % reset mat for each mesh 

                dim = meshes{i}.dim_prob;
                ncell = meshes{i}.ncell_prob;    % Number of cells in each spatial dimension
                nbdry = meshes{i}.nbdry_prob;    % Number of boundary cells
                ncell_ext = ncell + 2 * nbdry;   % Number of cells including boundaries
                nnode_ext = ncell_ext; nnode_ext = nnode_ext + 1;
                
                % Read in input data from hdf5  
                test_group_name = sprintf('/mesh_%d', i - 1);
                nmat = obj.h5read_transposed(filename, [test_group_name, '/nmat']); 
                matids = obj.h5read_transposed(filename, [test_group_name, '/matids']) + 1; 
                is_solid = obj.h5read_transposed(filename, [test_group_name, '/issolid']); 
                gamma_ea_mat = obj.h5read_transposed(filename, [test_group_name, '/gamma']);
                
                meshes{i}.nmat_mesh = nmat; 
                meshes{i}.matids_mesh = matids; 

                mat.nmat_prob = nmat;
                mat.matids_prob = matids;

                mat.nmat_in_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/nmat_mixcell']); % 1d array, nmixcell size
                mat.nmixcell = length(mat.nmat_in_mixcell);
                ijk_in_mixcell = double(obj.h5read_transposed(filename, ...
                    [test_group_name, '/ijk_mixcell'])) + 1; % 1d array 
                mat.ijk_in_mixcell = reshape(ijk_in_mixcell, dim, mat.nmixcell)'; 
                rho_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/rho_mixcell']);
                vf_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/vf_mixcell']);
                pres_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/pres_mixcell']);
                matids_mixcell = obj.h5read_transposed(filename, ...
                    [test_group_name, '/matids_mixcell']);

                % Unpack mixcell variables into cell arrays 
                % and store in mat
                mat.vf_in_mixcell = cell(mat.nmixcell, 1); 
                mat.rho_in_mixcell = cell(mat.nmixcell, 1);
                mat.pres_in_mixcell = cell(mat.nmixcell, 1);
                mat.matids_in_mixcell = cell(mat.nmixcell, 1);     

                offset = 1;             
                for mx = 1:mat.nmixcell
                    mat.vf_in_mixcell{mx} = vf_mixcell(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    mat.rho_in_mixcell{mx} = rho_mixcell(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    mat.pres_in_mixcell{mx} = pres_mixcell(offset:offset+mat.nmat_in_mixcell(mx)-1);
                    mat.matids_in_mixcell{mx} = matids_mixcell(offset:offset+mat.nmat_in_mixcell(mx)-1) + 1;                 
                    offset = offset+mat.nmat_in_mixcell(mx);
                end            
    

                % Function call
                mat.mat_sound_speed(meshes{i}, nmat, matids, is_solid, gamma_ea_mat);

                % Read fields after function call
                test_group_name = sprintf('/mesh_%d', i - 1);
                cs_after = obj.h5read_transposed(filename, [test_group_name, '/cs']);

                % Check
                nfailed = failed + 1;
                if ~all(abs(meshes{i}.cs_for_2dcell - cs_after) < tol*norm(cs_after))
                    if verbose, fprintf("%s: cs_for_2dcell differs\n", test_group_name), end
                    failed = nfailed;
                end 
                if failed < nfailed
                    passed = passed + 1;
                end

            end
        end
        
        function success = run_all(obj)
            % run_all  Runs all the tests
            %
            % Outputs:
            %  - passed : (logical) True is all tests were run successfully
            total = 0;
            passed = 0;
            fmtstr = '%-25s: [%-3d / %-3d]: %s\n';

            fprintf ('Running unit tests:\n');
            [t, p] = obj.test_rec_rec();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'rec_rec', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_sound_speed_solid();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'sound_speed_solid', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_find_interface2d();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'find_interface2d', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_cal_poly_area();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'cal_poly_area', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_rz_area();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'rz_area', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_bounds_2d();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'bounds_2d', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_order_nodes_along_norm();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'order_nodes_along_norm', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_cal_distance2d();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'cal_distance2d', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_cal_cell_zgrad2d();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'cal_cell_zgrad2d', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_reconstruct2d_nmat_pagosa();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'reconstruct2d_nmat_pagosa', p, t, 'pass'*q + 'fail'*(~q));

            [t, p] = obj.test_cal_mixcell_zgrad2d();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'cal_mixcell_zgrad2d', p, t, 'pass'*q + 'fail'*(~q));
            
            [t, p] = obj.test_get_mpoly();
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'get_mpoly', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_bdry_cell_2d(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'bdry_cell_2d', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_bdry_node_2d(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'bdry_node_2d', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_bdry_cell_1var_2d(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'bdry_cell_1var_2d', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_bdry_cell_ragged_2d(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'bdry_cell_ragged_2d', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_compute_divu(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'compute_divu', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_compute_qvis(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'compute_qvis', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_compute_force(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'compute_force', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_update_vel_comp(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'update_vel_comp', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_update_density(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'update_density', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_update_energy(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'update_energy', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_courant_from_cs(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'courant_from_cs', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_mesh_sound_speed(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'mesh_sound_speed', p, t, 'pass'*q + 'fail'*(~q));

            [meshes, ~] = obj.read_meshes_from_hdf5();
            [t, p] = obj.test_mat_sound_speed(meshes);
            q = (t == p); passed = passed + q; total = total + 1;
            fprintf(fmtstr, 'mat_sound_speed', p, t, 'pass'*q + 'fail'*(~q));

            success = (passed == total);
            fprintf('-----\n%-4s: %d out of %d passed', ...
                'OK  '*success + 'FAIL'*(~success), passed, total)

        end
    end
end

