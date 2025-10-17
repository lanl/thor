classdef c_vof < handle
    %INIT_STRUCT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vfmin double = 1.0e-06;   % smallest volume fraction, otherwise, considered zero. 
    end
    methods
        function vof_init(obj, dim, xl_prob, xr_prob, ncell_prob, nbdry_prob, nmat, t, ncycle, vf1d, rho1d, pres1d, ei1d)
            % Initialize variables
            filename = '';
            ncell_ext = zeros(1, 3);
            sizes = zeros(1, 4);
            llsize = 1;
            matids = [];
            dx = zeros(1, 3);
            
            % Pass mesh data
            [xl, xr, ncell, nbdry, nmat_for_2dcell, matid_for_2dcell, mixcell_for_2dcell, ...
                nmat_for_3dcell, matid_for_3dcell, mixcell_for_3dcell] = obj.mesh_pass_mat();
        
            ifmesh_initialized_for_mat = true;
        
            % Compute extended cell sizes
            for i = 1:dim
                ncell_ext(i) = ncell_prob(i) + nbdry_prob + nbdry_prob;
                sizes(i+1) = ncell_ext(i);
                llsize = llsize * ncell_ext(i);
            end
        
            % Initialize matids
            matids = 0:(nmat-1);
        
            sizes(1) = nmat;
            if dim == 2
                vf2d     = reshape(vf1d, sizes(2), sizes(1));
                rho2d    = reshape(rho1d, sizes(2), sizes(1));
                pres2d   = reshape(pres1d, sizes(2), sizes(1));
                ei2d     = reshape(ei1d, sizes(2), sizes(1));
                
                if isempty(nmat_for_2dcell)
                    ifmesh_initialized_for_mat = false;
                    nmat_for_2dcell  = zeros(ncell_ext(2), ncell_ext(1));
                    matid_for_2dcell = zeros(ncell_ext(2), ncell_ext(1));
                    mixcell_for_2dcell = -ones(ncell_ext(2), ncell_ext(1));
                end
                
                for j = 1:ncell_ext(2)
                    for i = 1:ncell_ext(1)
                        nmat_for_2dcell(j, i) = 1;
                        matid_for_2dcell(j, i) = 0;
                        mixcell_for_2dcell(j, i) = -1; % no mixed cell  
                    end
                end
            elseif dim == 3
                vf3d     = reshape(vf1d, sizes(3), sizes(2), sizes(1));
                rho3d    = reshape(rho1d, sizes(3), sizes(2), sizes(1));
                pres3d   = reshape(pres1d, sizes(3), sizes(2), sizes(1));
                ei3d     = reshape(ei1d, sizes(3), sizes(2), sizes(1));
                
                if isempty(nmat_for_3dcell)
                    nmat_for_3dcell  = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));
                    matid_for_3dcell = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));
                    mixcell_for_3dcell = -ones(ncell_ext(3), ncell_ext(2), ncell_ext(1));
                end
                
                for k = 1:ncell_ext(3)
                    for j = 1:ncell_ext(2)
                        for i = 1:ncell_ext(1)
                            nmat_for_3dcell(k, j, i) = 1;
                            matid_for_3dcell(k, j, i) = 0;
                            mixcell_for_3dcell(k, j, i) = -1;
                        end
                    end
                end
            end
        
            nmix = 0;
            lsize = 0;
            if dim == 2
                for j = 1:ncell_ext(2)
                    for i = 1:ncell_ext(1)
                        nm = 0;
                        vfs = vf2d{j, i};
                        for m = 1:nmat
                            if vfs(m) >= obj.vfmin
                                nm = nm + 1;
                                matid_for_2dcell(j, i) = m - 1;
                            end
                        end
                        nmat_for_2dcell(j, i) = nm;
                        if nm > 1
                            mixcell_for_2dcell(j, i) = nmix;
                            nmix = nmix + 1;
                            lsize = lsize + nm;
                        end
                    end
                end
            elseif dim == 3
                for k = 1:ncell_ext(3)
                    for j = 1:ncell_ext(2)
                        for i = 1:ncell_ext(1)
                            nm = 0;
                            vfs = vf3d{k, j, i};
                            for m = 1:nmat
                                if vfs(m) >= obj.vfmin
                                    nm = nm + 1;
                                    matid_for_3dcell(k, j, i) = m - 1;
                                end
                            end
                            nmat_for_3dcell(k, j, i) = nm;
                            if nm > 1
                                mixcell_for_3dcell(k, j, i) = nmix;
                                nmix = nmix + 1;
                                lsize = lsize + nm;
                            end
                        end
                    end
                end
            end
        
            ijk_for_mixcell = zeros(nmix, dim);
            nmat_for_mixcell = zeros(1, nmix);
            vf_for_mixcell = zeros(nmix, lsize);
            rho_for_mixcell = zeros(nmix, lsize);
            pres_for_mixcell = zeros(nmix, lsize);
            ei_for_mixcell = zeros(nmix, lsize);
            matids_for_mixcell = zeros(nmix, lsize);
        
            offset = 0;
            nmix = 0;
            if dim == 2
                for j = 1:ncell_ext(2)
                    for i = 1:ncell_ext(1)
                        if nmat_for_2dcell(j, i) < 2, continue; end
                        ijk_for_mixcell(nmix, :) = [i, j];
                        nm = 0;
                        vfs = vf2d{j, i};
                        matids_for_mixcell(nmix, offset+1:offset+nmat_for_2dcell(j,i)) = find(vfs >= obj.vfmin) - 1;
                        vf_for_mixcell(nmix, offset+1:offset+nmat_for_2dcell(j,i)) = vfs(vfs >= obj.vfmin);
                        rho_for_mixcell(nmix, offset+1:offset+nmat_for_2dcell(j,i)) = rho2d{j, i}(vfs >= obj.vfmin);
                        pres_for_mixcell(nmix, offset+1:offset+nmat_for_2dcell(j,i)) = pres2d{j, i}(vfs >= obj.vfmin);
                        ei_for_mixcell(nmix, offset+1:offset+nmat_for_2dcell(j,i)) = ei2d{j, i}(vfs >= obj.vfmin);
                        nmat_for_mixcell(nmix) = sum(vfs >= obj.vfmin);
                        offset = offset + nmat_for_2dcell(j, i);
                        nmix = nmix + 1;
                    end
                end
            elseif dim == 3
                for k = 1:ncell_ext(3)
                    for j = 1:ncell_ext(2)
                        for i = 1:ncell_ext(1)
                            if nmat_for_3dcell(k, j, i) < 2, continue; end
                            ijk_for_mixcell(nmix, :) = [i, j, k];
                            nm = 0;
                            vfs = vf3d{k, j, i};
                            matids_for_mixcell(nmix, offset+1:offset+nmat_for_3dcell(k,j,i)) = find(vfs >= obj.vfmin) - 1;
                            vf_for_mixcell(nmix, offset+1:offset+nmat_for_3dcell(k,j,i)) = vfs(vfs >= obj.vfmin);
                            rho_for_mixcell(nmix, offset+1:offset+nmat_for_3dcell(k,j,i)) = rho3d{k, j, i}(vfs >= obj.vfmin);
                            pres_for_mixcell(nmix, offset+1:offset+nmat_for_3dcell(k,j,i)) = pres3d{k, j, i}(vfs >= obj.vfmin);
                            ei_for_mixcell(nmix, offset+1:offset+nmat_for_3dcell(k,j,i)) = ei3d{k, j, i}(vfs >= obj.vfmin);
                            nmat_for_mixcell(nmix) = sum(vfs >= obj.vfmin);
                            offset = offset + nmat_for_3dcell(k, j, i);
                            nmix = nmix + 1;
                        end
                    end
                end
            end
        
            if ~ifmesh_initialized_for_mat
                obj.mesh_get_mat(dim, xl_prob, xr_prob, ncell_prob, nbdry_prob, ...
                    nmat_for_2dcell, matid_for_2dcell, mixcell_for_2dcell, ...
                    nmat_for_3dcell, matid_for_3dcell, mixcell_for_3dcell);
            end
        
            obj.mat_get_mix(dim, nmat, nmix, ijk_for_mixcell, nmat_for_mixcell, matids_for_mixcell, ...
                vf_for_mixcell, rho_for_mixcell, pres_for_mixcell, ei_for_mixcell);
        
            % Calculate dx
            for i = 1:dim
                dx(i) = (xr_prob(i) - xl_prob(i)) / ncell_prob(i);
            end
        
            % Get material polygons/polyhedrons
            obj.get_mpoly(dim, ncell_prob, nbdry_prob, xl_prob, dx, nmix, ...
                nmat_for_2dcell, matid_for_2dcell, mixcell_for_2dcell, ...
                nmat_for_3dcell, matid_for_3dcell, mixcell_for_3dcell);
        
            % Free memory
            if dim == 2
                clear vf2d rho2d pres2d ei2d;
            elseif dim == 3
                clear vf3d rho3d pres3d ei3d;
            end
        
            % Write data to file
            filename = sprintf('file_%08d', ncycle);
            [fileid, meshid] = obj.create_file(filename, t, 0);
            obj.write_mesh_mat(fileid, 'mesh', dim, xl_prob, xr_prob, ncell_prob, nbdry_prob, ...
                nmat, matids, nmat_for_2dcell, matid_for_2dcell, ...
                nmat_for_3dcell, matid_for_3dcell, meshid);
        
            obj.write_mat(fileid, dim);
            obj.close_file(fileid);
        
            % Free matids
            clear matids;

    end
end