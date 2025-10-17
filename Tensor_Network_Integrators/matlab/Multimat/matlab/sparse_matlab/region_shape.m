classdef region_shape < handle
    properties
        type shape_type        % The type of shape, using the Shape_Type enum
        parameters = []        % An array of parameters describing the shape
        % // for shape_universe, parameters = NULL;
        % // for shape_sphere,   parameters[0] = radius, parameters[1:dim] = center of radius 
        % // for shape_cylinder, parameters[0:dim] = radius and center of the one circle 
        % //                     parameters[dim+1:dim+dim] = radius and center of another circle
        % // for shape_quad,     parameters[0:7] = (x0,y0; x1,y1; x2,y2; x3,y3) of vertices
        % // for shape_triangle, parameters[0:5] = (x0,y0; x1,y1; x2,y2) of three vertices
        % // for shape_hex,      parameters[0:15] = (x0,y0; x1,y1; ...; x15,y15) in the following order
        % // 
        % //             7__________6  
        % //            /.         /|   
        % //           / .        / |
        % //          4----------5  |  
        % //          |  3.......|..2   
        % //          | .        |  /       
        % //          |.         | /     
        % //          0----------1/       
        % //                      
        % //     z         
        % //     |   y  
        % //     |  /  
        % //     | / 
        % //      -------> x   

    end
    
    methods
        function obj = region_shape(mtype, parameters)
            % Constructor to initialize the region_shape
            obj.type = mtype;
            obj.parameters = parameters;
        end
    end
end