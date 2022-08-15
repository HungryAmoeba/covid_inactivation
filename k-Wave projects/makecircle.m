function circle = makecircle(Nx,Ny,inner_val, outer_val, x_center, y_center, rad)

circle = outer_val * ones(Nx,Ny);

for xiter = 1:Nx
    for yiter = 1:Ny
        if ((xiter-x_center)^2 + (yiter-y_center)^2 < rad^2)
            circle(xiter, yiter) = inner_val;
        end
    end
end


end