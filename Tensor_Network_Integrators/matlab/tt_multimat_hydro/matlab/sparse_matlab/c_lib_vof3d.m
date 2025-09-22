classdef c_lib_vof3d
    %LIB_VOF3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % libdir = "/home/korobkin/num/william_c_code";
        libdir = "../";
        libname = "libvof3d"; 
        header_file = "vof3d_lib.h";
    end
    
    methods
        function obj = c_lib_vof3d(obj)
            % constructor method
            loadlibrary(obj.libdir + "/" + obj.libname + ".so", ...
                        obj.libdir + "/" + obj.header_file);
        end

        function ls(obj)
            libfunctions(obj.libname,'-full');
        end

        function [coords_new, nodelist_for_interface, ...
                faceindexlow_of_interface, faceindexhgh_of_interface, ... 
                nface_lower, nnode_for_face_lower, nodelist_for_face_lower, ... 
                nface_upper, nnode_for_face_upper, nodelist_for_face_upper] = ...
                interface3d(obj, dx_scale, nface, nnode, coords, ...
                            nnode_for_face, nodelist_for_face, cell_vol, ...
                            vf_to_match, normal)
            nodelist_for_face = nodelist_for_face  - 1;
            verbose = false;
            
            if verbose
                % inputs
                dx_scale = 0.01;
                nface = 6;
                nnode = 8;
                coords = [0.0, 0.0, 0.0, ...
                1.0, 0.0, 0.0, ...
                1.0, 1.0, 0.0, ...
                0.0, 1.0, 0.0, ...
                0.0, 0.0, 1.0, ...
                1.0, 0.0, 1.0, ...
                1.0, 1.0, 1.0, ...
                0.0, 1.0, 1.0];
                nnode_for_face = [4, 4, 4, 4, 4, 4];
                nodelist_for_face = [0, 4, 7, 3, ...
                          1, 2, 6, 5, ...
                          0, 1, 5, 4, ...
                          2, 3, 7, 6, ...
                          0, 3, 2, 1, ...
                          4, 5, 6, 7];
                cell_vol = 1.0;
                vf_to_match = 0.96551724137931039; %0.5; 
                normal = [-1, 0, 0]; % [-1, 1, 1]/sqrt(3);
            end

            

            % outputs
            out_nnode_new = 0;
            out_coords_new = zeros(100, 1);
            out_nnode_for_interface = 0; 
            out_nodelist_for_interface = zeros(20, 1);
            out_faceindexlow_of_interface = 0;
            out_faceindexhgh_of_interface = 0;
            out_nface_lower = 0;
            out_nnode_for_face_lower = zeros(nface + 1, 1);
            out_nodelist_for_face_lower = zeros(8*(nface + 1), 1);
            out_nface_upper = 0;
            out_nnode_for_face_upper = zeros(nface + 1, 1);
            out_nodelist_for_face_upper = zeros(8*(nface + 1), 1);
            

        [~, ~, ~, ~, out_nnode_new, ... 
                  out_coords_new, ... 
                  out_nnode_for_interface, ... 
                  out_nodelist_for_interface, ... 
                  faceindexlow_of_interface, ... 
                  faceindexhgh_of_interface, ... 
                  out_nface_lower, ... 
                  out_nnode_for_face_lower, ... 
                  out_nodelist_for_face_lower, ... 
                  out_nface_upper, ... 
                  out_nnode_for_face_upper, ... 
                  out_nodelist_for_face_upper] = calllib(obj.libname, "matlab_interface3d", ...
                  dx_scale, nface, nnode, coords', nnode_for_face, nodelist_for_face, ...
                  cell_vol, vf_to_match, normal, ... 
                  out_nnode_new, ... 
                  out_coords_new, ... 
                  out_nnode_for_interface, ... 
                  out_nodelist_for_interface, ... 
                  out_faceindexlow_of_interface, ... 
                  out_faceindexhgh_of_interface, ... 
                  out_nface_lower, ... 
                  out_nnode_for_face_lower, ... 
                  out_nodelist_for_face_lower, ... 
                  out_nface_upper, ... 
                  out_nnode_for_face_upper, ... 
                  out_nodelist_for_face_upper);
            
            coords_new = out_coords_new(1:3*out_nnode_new);
            nodelist_for_interface = out_nodelist_for_interface(1:out_nnode_for_interface) + 1;
            nface_lower = out_nface_lower;
            nnode_for_face_lower = out_nnode_for_face_lower(1:out_nface_lower);
            nodelist_for_face_lower = out_nodelist_for_face_lower(1:sum(nnode_for_face_lower)) + 1;
            nface_upper = out_nface_upper;
            nnode_for_face_upper = out_nnode_for_face_upper(1:out_nface_upper);
            nodelist_for_face_upper = out_nodelist_for_face_upper(1:sum(nnode_for_face_upper)) + 1;

            if verbose
                fprintf("faceindexlow_of_interface = %d\n", faceindexlow_of_interface);
                fprintf("faceindexhgh_of_interface = %d\n", faceindexhgh_of_interface);
    
                offset = 1;
                fprintf("nodelist_for_face_lower = {\n");
                for fc = 1:out_nface_lower
                    for nd = 1:nnode_for_face_lower(fc)
                        fprintf("%4d ", nodelist_for_face_lower(offset));
                        offset = offset + 1;
                    end
                    fprintf("\n");
                end
                fprintf("};\n");
                
                offset = 1;
                fprintf("nodelist_for_face_upper = {\n");
                for fc = 1:out_nface_upper
                    for nd = 1:nnode_for_face_upper(fc)
                        fprintf("%4d ", nodelist_for_face_upper(offset));
                        offset = offset + 1;
                    end
                    fprintf("\n");
                end
                fprintf("};\n");
            end
            
        end

        function [nface_lower, nnode_lower, coords_lower, ...
                  nnode_ea_face_lower, nodelist_for_face_lower, ...
                  nface_upper, nnode_upper, coords_upper, ...
                  nnode_ea_face_upper, nodelist_for_face_upper] = ...
                  polyhedron_plane(obj, nface, nnode, coords, nnode_ea_face, ...
                  nodelist_for_face, norm_plane, ds_plane)
            nodelist_for_face = nodelist_for_face  - 1;
            verbose = false;
            
            if verbose
                % inputs
                nface = 6;
                nnode = 8;
                coords = [0.0, 0.0, 0.0, ...
                1.0, 0.0, 0.0, ...
                1.0, 1.0, 0.0, ...
                0.0, 1.0, 0.0, ...
                0.0, 0.0, 1.0, ...
                1.0, 0.0, 1.0, ...
                1.0, 1.0, 1.0, ...
                0.0, 1.0, 1.0];
                nnode_ea_face = [4, 4, 4, 4, 4, 4];
                nodelist_for_face = [0, 4, 7, 3, ...
                          1, 2, 6, 5, ...
                          0, 1, 5, 4, ...
                          2, 3, 7, 6, ...
                          0, 3, 2, 1, ...
                          4, 5, 6, 7];
                norm_plane = [-1, 0, 0]; % [-1, 1, 1]/sqrt(3);
                ds_plane = 0.005;
            end

            % outputs
            out_nface_lower = 0;
            out_nnode_lower = 0;
            out_coords_lower = zeros(50,1);
            out_nnode_ea_face_lower = zeros(10,1);
            out_nodelist_for_face_lower = zeros(50,1);

            out_nface_upper = 0;
            out_nnode_upper = 0;
            out_coords_upper = zeros(50,1);
            out_nnode_ea_face_upper = zeros(10,1);
            out_nodelist_for_face_upper = zeros(50,1);
            
            [~, ~, ~, ~, out_nface_lower, out_nnode_lower, out_coords_lower, ...
                     out_nnode_ea_face_lower, out_nodelist_for_face_lower, ...
                     out_nface_upper, out_nnode_upper, out_coords_upper, ...
                     out_nnode_ea_face_upper, out_nodelist_for_face_upper] = ...
            calllib(obj.libname, "matlab_polyhedron_plane", nface, nnode, coords', ...
                     nnode_ea_face, nodelist_for_face, norm_plane, ds_plane, ...
                     out_nface_lower, out_nnode_lower, out_coords_lower, ...
                     out_nnode_ea_face_lower, out_nodelist_for_face_lower, ...
                     out_nface_upper, out_nnode_upper, out_coords_upper, ...
                     out_nnode_ea_face_upper, out_nodelist_for_face_upper);
            nface_lower = out_nface_lower;
            nnode_lower = out_nnode_lower;
            coords_lower = out_coords_lower(1:3*out_nnode_lower);
            nnode_ea_face_lower = out_nnode_ea_face_lower(1:out_nface_lower);
            nodelist_for_face_lower = out_nodelist_for_face_lower(1:sum(nnode_ea_face_lower)) + 1;
            nface_upper = out_nface_upper;
            nnode_upper = out_nnode_upper;
            coords_upper = out_coords_upper(1:3*out_nnode_upper);
            nnode_ea_face_upper = out_nnode_ea_face_upper(1:out_nface_upper);
            nodelist_for_face_upper = out_nodelist_for_face_upper(1:sum(nnode_ea_face_upper)) + 1;
            
            if verbose    
                offset = 1;
                fprintf("nodelist_for_face_lower = {\n");
                for fc = 1:out_nface_lower
                    for nd = 1:nnode_ea_face_lower(fc)
                        fprintf("%4d ", nodelist_for_face_lower(offset));
                        offset = offset + 1;
                    end
                    fprintf("\n");
                end
                fprintf("};\n");
                
                offset = 1;
                fprintf("nodelist_for_face_upper = {\n");
                for fc = 1:out_nface_upper
                    for nd = 1:nnode_ea_face_upper(fc)
                        fprintf("%4d ", nodelist_for_face_upper(offset));
                        offset = offset + 1;
                    end
                    fprintf("\n");
                end
                fprintf("};\n");
            end
        end

        function delete(obj)
            unloadlibrary(obj.libname);
        end
        
    end
end

