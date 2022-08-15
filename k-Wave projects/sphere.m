
function sphere_indicies = sphere(x_len, y_len, z_len, x_center, y_center, z_center, inner_rad, outer_rad, default_val, boundary_val, interior_val)

% given inputs for dimensions of gridspace from xyz_len, 

%returns a matrix sphere_indicies where desired values given the depth of
%the radius is assigned a value of 1 while other values are assigned a
%value of zero.  


sphere_indicies = default_val * ones(x_len,y_len,z_len);

%function of this specfic sphere

for xiter = 1:x_len
    for yiter = 1:y_len
        for ziter = 1:z_len
            rad = rad_val(xiter, yiter, ziter, x_center, y_center, z_center);
            if rad <= outer_rad
                if inner_rad <= rad
                    
                    
                sphere_indicies(xiter, yiter, ziter) = boundary_val;
                end
            end
            if rad < inner_rad
                sphere_indicies(xiter, yiter, ziter) = interior_val;
        end
    end
end



end 

            
            
            
            
            
            