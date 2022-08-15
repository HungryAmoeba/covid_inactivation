function rad = rad_val(xiter, yiter, ziter, x_center, y_center, z_center)
%calculates radius values given the inputs for center and iteration

rad = ((xiter-x_center)^2+(yiter-y_center)^2 + (ziter-z_center)^2)^(1/2);

end

