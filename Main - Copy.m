clc;
clear all;

%% GRID
% 2D grid implicitly storing node positions

% NxN dimension of grid
N = 10;

% # of particles
P = N * N * 16;

% grid spacing
h = 1;

% masses
grid_m = zeros(N, N);

% velocities
grid_v = zeros(N, 2, N); % need to make matrix this way for right vector outputs

% forces
grid_int_f = zeros(N, 2, N);
grid_ext_f = zeros(N, 2, N);

% accelerations
grid_a = zeros(N, 2, N);

% basis function outputs (store them for faster computation)
basis_ij = zeros(N, N, P);
grad_basis_ij = zeros(N, 2, N, P);

%% PARTICLES

% positions
ptcl_pos = zeros(P, 2);

% velocities
ptcl_v = zeros(P, 2);

% velocity gradients (matrix)
ptcl_vgrad = zeros(2, 2, P);% matrices need to be stored this way for right vector outputs

% masses
ptcl_m = zeros(P);

% volumes
ptcl_volume = zeros(P);

% deformation gradients (matrix)
ptcl_def = zeros(2, 2, P);

% stresses (matrix)
ptcl_stress = zeros(2, 2, P);

% external forces
ptcl_ext_f = zeros(P, 2);

%% OTHER INITIALIZATIONS

% modulus of elasticity
E = 10.8; %gold

% time-step
dt = 0.1;

% intialize particle volumes, masses, deformation gradient
for ind_p = 1:P
    ptcl_volume(ind_p) = 1;
    ptcl_m(ind_p) = 1;
    ptcl_def(:, :, ind_p) = eye(2);
    
    % gravity
    ptcl_ext_f(ind_p, 2) = -9.8;
end

% initialize particle positions
ind_i = 1;
ind_j = 1;
for ind_p = 1:P
    ptcl_pos(ind_p, :) = [ind_i * h/4, ind_j * h/4];
    ind_i = ind_i + 1;
    if ind_i > 40
        ind_j = ind_j + 1;
        ind_i = 1;
    end
end
% save the non-deformed volumes
ptcl_init_volume = ptcl_volume;


%% ALGORITHM
figure
plot(ptcl_pos(:, 1), ptcl_pos(:, 2), 'b*');
M(1) = getframe;
% LOOP
for k = 1:50
    
    ptcl_v(1, :)
    % Step 0: Pre-computations
    for ind_i = 1:N
        for ind_j = 1:N
            
            for ind_p = 1:P
                % compute basis function at (i, j) for position of particle p [1]
                % and store it
                basis_ij(ind_i, ind_j, ind_p) = basis(ind_i, h, ptcl_pos(ind_p, 1)) ...
                    * basis(ind_j, h, ptcl_pos(ind_p, 2));

                % create gradient of basis function vector [2] and store it
                grad_basis_ij(ind_i, :, ind_j, ind_p) = [grad_basis(ind_i, h, ptcl_pos(ind_p, 1)) ... 
                    * basis(ind_j, h, ptcl_pos(ind_p, 2)), ...
                    grad_basis(ind_j, h, ptcl_pos(ind_p, 2)) ...
                    * basis(ind_i, h, ptcl_pos(ind_p, 1))];
            end
        end
    end



    % Step 1: Projecting particle data on to the grid
    for ind_i = 1:N
        for ind_j = 1:N

            % mass at particular node (ind_x, ind_y)
            mass = 0;

            % mass summed via basis function
            % m_i = sum(p = 1 to P) { (phi_i(p) * m_p) }
            for ind_p = 1:P

                % sum new mass [3]
                mass = mass + basis_ij(ind_i, ind_j, ind_p) * ptcl_m(ind_p);
            end

            % store final mass at grid node
            grid_m(ind_i, ind_j) = mass;


            % velocity at particular node (ind_x, ind_y)
            vel = [0 0];

            % velocity summed via basis function
            for ind_p = 1:P

                % sum new velocity [4]
                vel = vel + basis_ij(ind_i, ind_j, ind_p) ...
                    * ptcl_m(ind_p) * ptcl_v(ind_p, :);
            end

            % store final velocity at grid node
            grid_v(ind_i, :, ind_j) = vel;

        end
    end

    % Step 2: Evaluating field variables at material points to be used in the
    % constitutive model
    for ind_p = 1:P

        % velocity gradient at particle at (ind_p) is a 2x2 matrix
        vgrad = [0 0; 0 0];

        % velocity gradient summed via basis function and nodal velocities
        for ind_i = 1:N
            for ind_j = 1:N

                % sum velocity gradient [5]
                vgrad = vgrad + grad_basis_ij(ind_i, :, ind_j, ind_p).' ...
                    * grid_v(ind_i, :, ind_j);
            end
        end

        % store final velocity gradient at particle point
        ptcl_vgrad(:, :, ind_p) = vgrad;
    end

    % Step 3: Constitutive model evaluation
    % Hyperelastic constitutive model
    for ind_p = 1:P
        % create deformation gradient "~" [6]
        def_til = eye(2) + ptcl_vgrad(:, :, ind_p) * dt;

        % evaluate new particle deformation gradient [7]
        ptcl_def(:, :, ind_p) = def_til * ptcl_def(:, :, ind_p);

        % evaluate particle stress [8]
        ptcl_stress(:, :, ind_p) = E * (ptcl_def(:, :, ind_p) - eye(2));
    end

    % Step 4: Evaluation of internal forces (and external) on the grid

    % evaluate particle volumes
    for ind_p = 1:P
        ptcl_volume(ind_p) = det(ptcl_def(:, :, ind_p)) * ptcl_init_volume(ind_p);
    end

    % evaluate particle forces
    for ind_i = 1:N
        for ind_j = 1:N

            f_ext = [0 0];
            f_int = [0; 0];

            for ind_p = 1:P

                % sum external forces on nodes [9]
                f_ext = f_ext + basis_ij(ind_i, ind_j, ind_p) * ptcl_ext_f(ind_p, :);

                % sum internal forces on nodes [11]
                f_int = f_int - ptcl_stress(:, :, ind_p) ...
                    * grad_basis_ij(ind_i, :, ind_j, ind_p).' ...
                    * ptcl_volume(ind_p);
            end

            % transpose f_int to make it match other vectors
            f_int = f_int.';

            % store forces on grid nodes)
            grid_int_f(ind_i, :, ind_j) = f_int;
            grid_ext_f(ind_i, :, ind_j) = f_ext;
        end
    end

    % Step 5: Advecting of the Lagrangian particles

    % Calculating grid accelerations and (future) velocities via Newton's Method
    for ind_i = 1:N
        for ind_j = 1:N

            % compute grid accelerations [12]
            grid_a(ind_i, :, ind_j) = (grid_int_f(ind_i, :, ind_j) + ...
                grid_ext_f(ind_i, :, ind_j)) / grid_m(ind_i, ind_j);

            % compute predicted future grid velocities [13].
            % storing them in the grid_v array because why not.
            % these values will not actually be used in the next time step, it
            % is just for computation right now. The next time step grid
            % velocities will be calculated normally in Step 1, which uses
            % particle velocities that are calculated from the grid velocities
            % here.
            grid_v(ind_i, :, ind_j) = grid_v(ind_i, :, ind_j) + ... 
                grid_a(ind_i, :, ind_j) * dt;
        end
    end

    % Calculate new particle velocities and positions
    for ind_p = 1:P

        % create variable for particle velocities and positions (to avoid the 
        % indexing time using arrays)
        vel = ptcl_v(ind_p, :);
        pos = ptcl_pos(ind_p, :);

        for ind_i = 1:N
            for ind_j = 1:N
                % sum new vel [14]
                vel = vel + basis_ij(ind_i, ind_j, ind_p) * grid_a(ind_i, :, ind_j) * dt;

                % sum new pos [15]
                pos = pos + basis_ij(ind_i, ind_j, ind_p) * grid_v(ind_i, :, ind_j) * dt;
            end
        end
        
        % store new particle velocities and positions
        ptcl_v(ind_p, :) = vel;
        ptcl_pos(ind_p, :) = pos;
    end
    
    
    plot(ptcl_pos(:, 1), ptcl_pos(:, 2), 'b*');
    M(k+1) = getframe;
    %movie2avi(M, 'mpm_movie');
end



