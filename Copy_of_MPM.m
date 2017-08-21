clc;
clear all;

%% Initilization

mass = 0.001; % all particles have same mass
h = 1; % spacing
dt = 0.01; % time differential
T = 100; % # of time steps

% Field constants
gravity = -9.8;
MoE = 88; % Modulus of elasticity
PoR = 0.48; % Poisson's ratio

% mxm grid nodes (m+1 x m+1 grid nodes)(stored in 2D matrix)
m = 20;

% N particles (stored in 1D matrix)
N = 100;

% Un - 2D position array
Xp = zeros(N, 2);
dXp = zeros(N, 2);
Xp_prev = zeros(N, 2);
dXp_prev = zeros(N, 2);
Xp_init = zeros(N, 2);

% Mn - Mass array
Mp = ones(N, 1);
Mp = Mp./100;

% Vn - Velocity array
Vp = zeros(N, 2);

% Voln - Volume array
Volp = zeros(N, 1);
for p = 1:N
    Volp(p) = (h/10)^2;
end

% Def - Deformation gradient
Defp = zeros(2, 2, N);
% Piola Kirchoff stress tensor
Piola = zeros(2, 2, N);
for p = 1:N
    Defp(:,:,p) = eye(2, 2);
    Piola(:,:,p) = eye(2);
end

% Potential energy, Neohookian Elasticity
energy = zeros(N, 1);
Energy_total = 0;
Energy_total_prev = 0;
dEnergy_total = 0;
dEnergy_total_prev = 0;

% Vector that tells you the row and column index of a grid cell
GridCells = zeros(m^2, 2);
gindex = 1;
for gj = 0:(m-1)
    for gi = 0:(m-1)
        GridCells(gindex, :) = [gi, gj];
        gindex = gindex + 1;
    end
end

%% Initializing positons
% initialize particle positions
ind_i = 0;
ind_j = 0;
for ind_p = 1:N
    Xp(ind_p, :) = [0 + ind_i * h/10, 0 + ind_j * h/10];
    ind_i = ind_i + 1;
    if ind_i > 9
        ind_j = ind_j + 1;
        ind_i = 0;
    end
end
Xp(:, 2) = Xp(:, 2) + 5;
Xp(:, 1) = Xp(:, 1) + 5;
% Keep initial position
Xp_init = Xp;
Xp_prev = Xp;
Xp_init(:,1) = Xp_init(:, 1) + 0.01;
% Up(:, 2) = (Up(:,2) + 5);
% Up(:, 1) = (Up(:,1) + 3);

%Un(1, :) = [2, 2.2];

%% LOOP


% For visualization
figure
scatter(Xp(:, 1), Xp(:, 2), 15, energy, 'filled');
colorbar;
grid on;
set(gca,'xtick',[0:1:m]);
set(gca,'ytick',[0:1:m]);
axis([-1 21 -1 21]);
title('Object simulation');
M(1) = getframe;

for k = 1:T
    %% STEP 0 - Compute interpolation weights
    Wpg = zeros(N, m^2);
    WpgGrad = zeros(N, 2, m^2);
    for p = 1:N
        for ind = 1:m^2
            Wpg(p, ind) = BilinearInterpolation(Xp(p,1), Xp(p,2), GridCells(ind, 1), ...
                GridCells(ind, 2), h);
            WpgGrad(p, :, ind) = BiLinGrad(Xp(p,1), Xp(p,2), GridCells(ind, 1), ...
                GridCells(ind, 2), h);
        end 
    end
    
    %% STEP 1 - Project data from particles to grid (mass and velocity)
    
    % Mg - grid mass vector. The position of these nodes is implicitly
    % understood because of the uniform grid spacing h.
    Mg = zeros(m^2, 1);
    
    % Vg - grid velocity vector
    Vg = zeros(m^2, 2);
    
    % Ug - grid displacement vector
    UDispg = zeros(m^2, 2);
    % Approximate UDispp
    UDispp = Xp - Xp_init;
    
    
    for ind = 1:m^2
        % Project mass from material points to grid nodes
        for p = 1:N
            Mg(ind) = Mg(ind) + Wpg(p, ind) .* Mp(p);
        end
        
        % Project velocity from material points to grid nodes
        for p = 1:N
            Vg(ind,:) = Vg(ind,:) + Wpg(p, ind) * repmat(Mp(p), [1 2]) .* Vp(p,:);
        end
        if Mg(ind) == 0
            Mg(ind) = 1;
            
        end
        Vg(ind,:) = Vg(ind,:) ./ repmat(Mg(ind), [1 2]);
        
%         % Project displacement from material points to grid nodes
%         for p = 1:N
%             UDispg(ind, :) = UDispg(ind, :) + BilinearInterpolation(Xp(p,1), Xp(p,2), GridCells(ind, 1), ...
%                 GridCells(ind, 2), h) * UDispp(p,:);
%         end
    end

    
    %% STEP 2 - Evaluate field variables at particle points to be used in consitutive model
    
    % Hyper Elasticity Model
    
%     % Create matrix for velocity gradients of particles
%     VGradp = zeros(2, 2, N);
%     
%     % Create matrix for deformation gradients
%     
%     for p = 1:N
%         for ind = 1:m^2
%             VGradp(:, :, p) = VGradp(:, :, p) + WpgGrad(p, :, ind).' * Vg(ind, :);
%         end
%     end
%     
%     % Create a vector that stores the grid velocities w.r.t. particles
%     % within a grid, one for each corner of a cell
%     Vg1 = Vg((GridCellj)*m + GridCelli+1, :);
%     Vg2 = Vg((GridCellj+1)*m + GridCelli+1, :);
%     Vg3 = Vg((GridCellj+1)*m + GridCelli+1+1, :);
%     Vg4 = Vg((GridCellj)*m + GridCelli+1+1, :);
%     
%     % Compute velocity gradients
%     for p = 1:N
%         VGradp(:, :, p) = VGradp(:, :, p) + BasisGrad1(p,:).' * Vg1(p,:);
%         VGradp(:, :, p) = VGradp(:, :, p) + BasisGrad2(p,:).' * Vg2(p,:);
%         VGradp(:, :, p) = VGradp(:, :, p) + BasisGrad3(p,:).' * Vg3(p,:);
%         VGradp(:, :, p) = VGradp(:, :, p) + BasisGrad4(p,:).' * Vg4(p,:);
%     end
%     
    % Create matrix for displacement gradients of particles
%     UDispGradp = ParticleGradient(UDispg, GridCelli, GridCellj, ...
%         GridCellParticleCount, gridCell, BasisGrad1, BasisGrad2, BasisGrad3, BasisGrad4, m, N);
%      UDispGradp = zeros(2, 2, N);
%      for p = 1:N
%         for ind = 1:m^2
%            	UDispGradp(:, :, p) = UDispGradp(:, :, p) + BiLinGrad(Xp(p, 1), Xp(p, 2), ...
%                 GridCells(ind, 1), GridCells(ind, 2), h).' * UDispg(ind, :);
%         end
%     end
%     % Compute deformation gradients
%     for p = 1:N
%         Defp(:, :, p) = inv(eye(2) - UDispGradp(:, :, p));
%     end
    
    
    
    %% STEP 4 - Evaluation of forces on the grid
    
    % External forces (gravity)
    % Calculate gravity directly on grid
    Feg = zeros(m^2, 2);
    
    % F = mg
    Feg(:, 2) = Mg * gravity;
    
    % internal forces
    Fig = zeros(m^2, 2);
    
    % Piola Kirchoff Stress tensor
    Piola = PiolaKirchoff(Defp, MoE, PoR, N);
    
    % Compute internal forces
    for ind = 1:m^2
        
        for p = 1:N
          Fig(ind, :) = Fig(ind, :) - (Piola(:, :, p) * (Defp(:, :, p).') * (WpgGrad(p, :, ind).')).' * Volp(p);
            %Fig(ind, :) = Fig(ind, :) - (Piola(:, :, p) * (Defp(:, :, p).') * (WpgGrad(p, :, ind).')).' * Volp(p);
%         Fig(GridCellj(p)*m + GridCelli(p)+1, :) = Fig(GridCellj(p)*m + GridCelli(p)+1, :) - ...
%             (Strp(:, :, p) * BasisGrad1(p, :).').' * VolDefp(p);
%         Fig((GridCellj(p)+1)*m + GridCelli(p)+1, :) = Fig((GridCellj(p)+1)*m + GridCelli(p)+1, :) - ...
%             (Strp(:, :, p) * BasisGrad2(p, :).').' * VolDefp(p);
%         Fig((GridCellj(p)+1)*m + GridCelli(p)+1+1, :) = Fig((GridCellj(p)+1)*m + GridCelli(p)+1+1, :) - ...
%             (Strp(:, :, p) * BasisGrad3(p, :).').' * VolDefp(p);
%         Fig(GridCellj(p)*m + GridCelli(p)+1+1, :) = Fig(GridCellj(p)*m + GridCelli(p)+1+1, :) - ...
%             (Strp(:, :, p) * BasisGrad4(p, :).').' * VolDefp(p);
        end
    end

%     
    %% STEP 5 - Advecting of the Lagrangian Particles
    
    % Update acceleration and velocity at each grid node
    
    % Change all 0's in mass matrix to 1 to avoid dividing by 0
    for ind = 1:m^2
        if Mg(ind) == 0
            Mg(ind) = 1;
        end
    end
    Ag = (Fig + Feg) ./ repmat(Mg, [1 2]);
    
    % Change NaN to 0, since NaN is caused by diving by 0 mass. Just change
    % acceleration to 0 in that case
    % Am(isnan(Am)) = 0;
    
    Vg = Vg + Ag * dt;
    
    % Save previous positions
    Xp_prev = Xp;
    dXp_prev = dXp;
    
%     % Semi-Implicit update equation
%      fun = @(v_next) -M*v_next - dt^2 * K * v_next ...
%          + M * v_curr + dt * (f_ext - f_int);
%     
    % Update velocity each particle
    for p = 1:N
        posx = Xp(p, 1);
        posy = Xp(p, 2);
        for ind = 1:m^2
            Vp(p, :) = Vp(p, :) + Wpg(p, ind) * Ag(ind, :) * dt;
            Xp(p, :) = Xp(p, :) + Wpg(p, ind) * Vg(ind, :) * dt;
        end
%         Xp(p, 2) - posy
        dXp = Xp - Xp_prev;
    end
    
    % Check if they hit the 2nd bottom of the grid
    for p = 1:N
        if Xp(p, 2) <= 0
            Xp(p, 2) = 0;
            Vp(p, 2) = 0;
        elseif Xp(p, 2) >= m-1
            Xp(p, 2) = m-1;
            Vp(p, 2) = 0;
        end
        if Xp(p, 1) <= 0
            Xp(p, 1) = 0;
            Vp(p, 1) = 0;
        elseif Xp(p, 1) >= m - 1
            Xp(p, 1) = m - 1;
            Vp(p, 1) = 0;
        end
    end
    
    %% STEP 3 - Evaluate the constitutive model
    % Hyper elasticity
    
    % Allocate stress matrix
    Strp = zeros(2, 2, N);
    
    % Update deformation gradient and get stress matrix using the
    % constitutive model
    Deftp = zeros(2, 2, N);
    for p = 1:N
        Deftp(:, :, p) = eye(2);
        for ind = 1:m^2
            Deftp(:, :, p) = Deftp(:, :, p) + dt * (Vg(ind, :)).' * WpgGrad(p,:,ind);
        end
    end

    % Compute deformation gradient of each particle directly
    for p = 1:N
        Defp(:, :, p) = Deftp(:, :, p) * Defp(:, :, p);
        Strp(:, :, p) = MoE * (Defp(:,:,p) - eye(2));
    end
    
    % Get deformed volumes
    VolDefp = zeros(N, 1);
    for p = 1:N
        VolDefp(p) = det(Defp(:, :, p)) * Volp(p);
    end
    
    % Save previous energy_total
    Energy_total_prev = Energy_total;
    % Compute energy at each particle
    energy = Neohookian(Defp, MoE, PoR, N);
    Energy_total = sum(energy .* Volp);
    % Save previous dEnergy_total
    dEnergy_total_prev = dEnergy_total;
    dEnergy_total = Energy_total - Energy_total_prev;
    
    % stiffness matrix
    d2E = (dEnergy_total - dEnergy_total_prev) / dt;
    
    K = zeros(2, 2, N);
    for p = 1:N
        dx2 = (dXp(p, 1) - dXp_prev(p, 1)) / dt;
        dx = dXp(p, 1);
        dy = dXp(p, 2);
        dy2 = (dXp(p, 2) - dXp_prev(p, 2)) / dt;
        K(1, 1, p) = d2E/dx2;
        K(1, 2, p) = d2E/(dx*dy);
        K(2, 1, p) = K(1, 2, p);
        K(2, 2, p) = d2E/dy2;
    end
    
    
    
    %k*dt;
    
    %% Update figure

    
%     % Create stiffness matrix
%     % First compute differentials
%     d_Un = Un - Un_prev;
%     d_energy = energy - energy_prev;
%     % Then stiffness matrix
%     K = zeros(2,2,N);
%     for p = 1:N
%         ensq = d_energy(p)^2;
%         K(1, 1, p) = ensq / d_Un(p, 1)^2;
%         K(1, 2, p) = ensq / (d_Un(p, 1) * d_Un(p, 2));
%         K(2, 1, p) = K(1, 2, p);
%         K(2, 2, p) = ensq / d_Un(p, 2)^2;
%     end
%     VGradp
%     Deftp
%     Defp
%     Fig
%     Ag
    % Plot particles on figure, update movie frame
%     subplot(2,1,1);
    scatter(Xp(:, 1), Xp(:, 2), 15, energy, 'filled');
    colorbar;
    grid on;
    set(gca,'xtick',[0:1:m]);
    set(gca,'ytick',[0:1:m]);
    axis([-1 21 -1 21]);
    title('Object simulation');
%     subplot(2,1,2);
%     plot(particle(:, 1), energy(:,1), 'r*');
    
    M(k+1) = getframe;
    %movie2avi(M, 'mpm_movie');

end


