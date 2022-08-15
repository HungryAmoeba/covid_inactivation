% Simulation in three dimensions using k-Wave
% 
% Assumes that the surrounding medium is water at 20C
% 
% Has speed 1481 m/s
% 
% Medium density is 1000 kg/m^3
% 
% Need to specify and find membrane coefficients and specficy exazt size. 
% For model, assume that the lipid membrane 100 nm in diameter. Lipid bilayer is 4 nm thick


% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
%note that f_max = medium_sound_speed/(2*dx)

Nx = 64;            % number of grid points in the x direction
Ny = 64;            % number of grid points in the y direction
Nz = 64;            % number of grid points in the z direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
dz = 0.1e-3;        % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
% needs adaptation to circular function. 

%produces logical array where the values inside the sphere
xyz_len = 64;
x_center = 50;
y_center = 32;
z_center = 32;
inner_rad = 5;
outer_rad = 6;


% medium.sound_speed = 1500 * ones(Nx, Ny, Nz);	% [m/s]
% medium.sound_speed(sphere_indicies) = 1800;        % [m/s]
% medium.density = 1000 * ones(Nx, Ny, Nz);       % [kg/m^3]
% medium.density(sphere_indicies) = 1200;          % [kg/m^3]

%1800 and 1200 are currently arbitrary
medium.sound_speed = sphere(xyz_len, xyz_len, xyz_len, x_center, y_center, z_center, inner_rad, outer_rad, 1500, 1800, 800);
medium.density = sphere(xyz_len, xyz_len, xyz_len, x_center, y_center, z_center, inner_rad, outer_rad, 1000, 1200, 700);

% create initial pressure distribution using makeBall
% later can this pressure distribution be adapted to simulate the freq? 
ball_magnitude = 8;    % [Pa]
ball_x_pos = 10;        % [grid points]
ball_y_pos = 32;        % [grid points]
ball_z_pos = 32;        % [grid points]
ball_radius = 5;        % [grid points]
ball_1 = ball_magnitude * makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

% Test to represent viral particle in voxelplot 
ball_magnitude = 0.00000001;    % [Pa]
ball_x_pos = 50;        % [grid points]
ball_y_pos = 32;        % [grid points]
ball_z_pos = 32;        % [grid points]
ball_radius = 6;        % [grid points]
ball_2 = ball_magnitude * makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

source.p0 = ball_1 + ball_2;

% define a series of Cartesian points to collect the data
% later fill this in to the layers of the membrane

% x = (-22:2:22) * dx;            % [m]
% y = (-22:2:22) * dy;    % [m]
% z = (-22:2:22) * dz;            % [m]

x = (-11:1:1) * dx;
y = (-11:1:1) * dy;
z = (-11:1:1) * dz;
sensor.mask = [x; y; z];

% input arguments
input_args = {'PlotLayout', true, 'PlotPML', false, ...
    'DataCast', 'single', 'CartInterp', 'nearest'}; %, ... 
    %'RecordMovie', true, 'MovieName', 'example_movie_1'}; %delete this part if buggy

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the simulation layout using voxelplot
voxelPlot(double(source.p0 | cart2grid(kgrid, sensor.mask)));
view([50, 20]);

% plot the simulated sensor data
figure;
imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;

