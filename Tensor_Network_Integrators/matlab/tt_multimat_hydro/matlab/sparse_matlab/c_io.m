classdef c_io < handle
    %INIT_STRUCT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % todo
    end
    methods
        function fileid = create_file(obj, filename, t, ncycle)
            % todo
            fileid = 0;
        end

        function meshid = write_mesh_mat(obj, fileid, meshname, mesh)
            % mesh is an object of class c_mesh
            dim = mesh.dim_prob;
            xl_prob = mesh.xl_prob;
            xr_prob = mesh.xr_prob;
            ncell_prob = mesh.ncell_prob;
            nbdry_prob = mesh.nbdry_prob;
            nmat_mesh = mesh.nmat_mesh;
            matids_mesh = mesh.matids_mesh;
            nmat_for_2dcell = mesh.nmat_for_2dcell;
            matid_for_2dcell = mesh.matid_for_2dcell;
            nmat_for_3dcell = mesh.nmat_for_3dcell;
            matid_for_3dcell = mesh.matid_for_3dcell;

            % Original logic here using these extracted properties
            meshid = 0;
        end

        function write_mesh_vars(obj, fileid, meshid, mesh)
            % mesh is an object of class c_mesh
            dim = mesh.dim_prob;
            ncell_prob = mesh.ncell_prob;
            nbdry_prob = mesh.nbdry_prob;
            nmat_for_2dcell = mesh.nmat_for_2dcell;
            rho_for_2dcell = mesh.rho_for_2dcell;
            ei_for_2dcell = mesh.ei_for_2dcell;
            pres_for_2dcell = mesh.pres_for_2dcell;
            vel_for_2dnode = mesh.vel_for_2dnode;
            nmat_for_3dcell = mesh.nmat_for_3dcell;
            rho_for_3dcell = mesh.rho_for_3dcell;
            ei_for_3dcell = mesh.ei_for_3dcell;
            pres_for_3dcell = mesh.pres_for_3dcell;
            vel_for_3dnode = mesh.vel_for_3dnode;

            % Original logic here using these extracted properties
        end
        function close_file(obj, fileid)
            % todo
        end

    end

end