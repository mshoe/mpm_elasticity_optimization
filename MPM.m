clc;
clear all;

%% Initilization

mass = 0.001; % all particles have same mass
h = 1; % spacing
D = h/16; % particle diameter
dt = 0.01; % time differential
T = 200; % # of time steps

% Field constants
gravity = -9.8;
MoE = 1e6; % Modulus of elasticity
PoR = 0.47; % Poisson's ratio

% mxm grid nodes (m+1 x m+1 grid nodes)(stored in 2D matrix)
m = 2;

% N particles (stored in 1D matrix)
N = 64;
P = linspace(1, N, N); % array that has particle indexes

% Un - 2D position array
Xp = zeros(N, 2);
dXp = zeros(N, 2);
Xp_prev = zeros(N, 2);
dXp_prev = zeros(N, 2);
Xp_init = zeros(N, 2);

% Mn - Mass array
Mp = ones(N, 1);
Mp = Mp./100;

% Fep - external forces
Fep = zeros(N, 2);
Fep(:,2) = Mp * gravity;

% Vn - Velocity array
Vp = zeros(N, 2);

% Voln - Volume array
Volp = zeros(N, 1);
for p = 1:N
    Volp(p) = pi * (D/2)^2;
end

% Def - Deformation gradient
Defp = zeros(2, 2, N);
for p = 1:N
    Defp(:,:,p) = eye(2, 2);
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
    Xp(ind_p, :) = [0.1 + ind_i * D, 0.1 + ind_j * D];
    ind_i = ind_i + 1;
    if ind_i > 7
        ind_j = ind_j + 1;
        ind_i = 0;
    end
end
% for ind_p = 1:N
%     Xp(ind_p, :) = [-0.2 * cos(2*pi/16 * ind_p), -0.2 * sin(2*pi/16 * ind_p)];
% end

% Xp(7, 2) = Xp(7, 2) - 0.05;
% Xp(6, 2) = Xp(6, 2) - 0.1;
% Xp(5, 2) = Xp(5, 2) - 0.15;
% Xp(4, 2) = Xp(4, 2) - 0.15;
% Xp(3, 2) = Xp(3, 2) - 0.1;
% Xp(2, 2) = Xp(2, 2) - 0.05;



Xp(:, 1) = Xp(:, 1) + 0.18 ;%+ 0.5;
Xp(:, 2) = Xp(:, 2) + 0.3;

% Xp(1, :) = [-0.25, -0.25];
% Xp(2, :) = [-0.25, 0.25];
% Xp(3, :) = [0.25, 0.25];
% Xp(4, :) = [0.25, -0.25];
%Xp(5, :) = [0.27, -0.27];

% Xp(:, 2) = Xp(:, 2) + 0.2;
% Xp(:, 1) = Xp(:, 1) + 5;
% Keep initial position
%Xp(:,1) = Xp(:,1) * 0.75;
%Xp(:,2) = Xp(:,2) * 0.75;
% Xp = Xp + 1.5;
Xp_init = Xp;% * 1.1;


%% LOOP


% For visualization
figure
scatter(Xp(:, 1), Xp(:, 2), 15, energy, 'filled');
colorbar;
grid on;
set(gca,'xtick',[0:1:m]);
set(gca,'ytick',[0:1:m]);
axis([0 m-1 0 m-1]);
% axis([3 8 -1 4]);
% axis([-1 11 -1 11]);
title('MPM');
M(1) = getframe;

%% STEP -1 - Calculate CoM
    CoM = [0, 0];
    for p = 1:N
        CoM = CoM + Xp_init(p, :);
    end
    CoM = CoM / N;

for k = 1:T
    %% STEP 0 - Compute interpolation weights
    Wpg = Weights(Xp, GridCells, N, m, h);
    WpgGrad = WeightsGrad(Xp, GridCells, N, m, h);
    
    
    
    %% STEP 1 - Project data from particles to grid (mass and velocity)
    
    % Mg - grid mass vector. The position of these nodes is implicitly
    % understood because of the uniform grid spacing h.
    Mg = zeros(m^2, 1);
    
    % Vg - grid velocity vector
    Vg = zeros(m^2, 2);
    
    % Ug - grid displacement vector
    Ug = zeros(m^2, 2);
    % Approximate UDispp
    Up = Xp - Xp_init;
    %Up = Xp - repmat(CoM, [N 1]);
    %Up = Up * 0.5;
    
    % Xg
    Xg = zeros(m^2, 2);
    
    % External force
    
    
%    Feg = zeros(m^2, 2);
    
    % OPTIMIZATION PROJECTION TO GRID
    
    
    % Optimization projection of particle displacement
%     funUgx = @(Ugx) 0.5 * sum((sum(repmat(Ugx, [1 N]) .* Wpg.', 1) - Up(:, 1).').^2, 2);
%     funUgy = @(Ugy) 0.5 * sum((sum(repmat(Ugy, [1 N]) .* Wpg.', 1) - Up(:, 2).').^2, 2);
% 
%     Ug(:, 1) = fminunc(funUgx, Ug(:,1));
%     Ug(:, 2) = fminunc(funUgy, Ug(:,1));

    Ug(:, 1) = ProjPG(Up(:, 1), Wpg, N, m);
    Ug(:, 2) = ProjPG(Up(:, 2), Wpg, N, m);
    
    
    % Optimization projection of mass
%     funMg = @(Mg) 0.5 * sum((sum(repmat(Mg, [1 N]) .* Wpg.', 1) - Mp.').^2, 2);
%     
%     Mg = fminunc(funMg, Mg);
    Mg = ProjPG(Mp, Wpg, N, m);
    
    % Optimization projection of velocity
%     funVgx = @(Vgx) 0.5 * sum((sum(repmat(Vgx, [1 N]) .* Wpg.', 1) - Vp(:, 1).').^2, 2);
%     funVgy = @(Vgy) 0.5 * sum((sum(repmat(Vgy, [1 N]) .* Wpg.', 1) - Vp(:, 2).').^2, 2);
%     
%     Vg(:, 1) = fminunc(funVgx, Vg(:,1));
%     Vg(:, 2) = fminunc(funVgy, Vg(:,2));
    Vg(:, 1) = ProjPG(Vp(:, 1), Wpg, N, m);
    Vg(:, 2) = ProjPG(Vp(:, 2), Wpg, N, m);
    
    Xg(:, 2) = ProjPG(Xp(:, 2), Wpg, N, m);
    
    % Optimization projection of exteral force
    
%     funFegy = @(Fegy) 0.5 * sum((sum(repmat(Fegy, [1 N]) .* Wpg.', 1) - Fep(:, 2).').^2, 2);
%     Feg(:, 2) = fminunc(funFegy, Feg(:, 2));

%    Feg(:, 2) = ProjPG(Fep(:, 2), Wpg, N, m);
%     for ind = 1:m^2
%         % Project mass from material points to grid nodes
%         for p = 1:N
%             Mg(ind) = Mg(ind) + Wpg(p, ind) .* Mp(p);
%         end
%         
%         % Project velocity from material points to grid nodes
%         for p = 1:N
%             Vg(ind,:) = Vg(ind,:) + Wpg(p, ind) * repmat(Mp(p), [1 2]) .* Vp(p,:);
%         end
%         if Mg(ind) == 0
%             Mg(ind) = 1;
%             
%         end
%         Vg(ind,:) = Vg(ind,:) ./ repmat(Mg(ind), [1 2]);
%         
%         % Project displacement from material points to grid nodes using
%         % optmization
%         
%         
% %         for p = 1:N
% %             Ug(ind, :) = Ug(ind, :) + Wpg(p, ind)  * repmat(Mp(p), [1 2]) .* Up(p,:);
% %         end
% %         
% %         Ug(ind,:) = Ug(ind,:) ./ repmat(Mg(ind), [1 2]);
% 
%     end

    
    %% STEP 2 - Evaluate field variables at particle points to be used in consitutive model
    
    % Hyper Elasticity Model
    

    
    % Create deformation gradient
    Defp = DeformationGrad(Ug, WpgGrad, N, m);

    % Compute Mass Inertia matrix (create virtual projection of mass from
    % grid to particles)
    MMg = zeros(2*m^2, 2*m^2);
    
%     for indx = 1:2*m^2
%         for indy = 1:2*m^2
%             for p = 1:N
%                 MMg(indx, indy) = MMg(indx, indy) + Wpg(p, ceil(indx/2)) * Wpg(p, ceil(indy/2)) * Mp(p);
%             end
%         end
%     end
    for p = 1:N
        Wp = Wpg(p,1) * eye(2);
        for ind = 2:m^2
            Wp = [Wp Wpg(p, ind) * eye(2)];
        end
        %Wp = [Wpg(p, 1)*eye(2) Wpg(p, 2)*eye(2) Wpg(p, 3)*eye(2) Wpg(p, 4)*eye(2)];
        MMg = MMg + Mp(p) * (Wp.' * Wp);
    end
    
    
    %% STEP 3 - Evaluate the constitutive model
%     % Hyper elasticity
%     
%     % Allocate stress matrix
%     Strp = zeros(2, 2, N);
% % % %     
% % % %     % Update deformation gradient and get stress matrix using the
% % % %     % constitutive model
% % % %     Deftp = VGradp * dt;
% % % % 
%     % Compute stress matrix of each particle directly
%     for p = 1:N
% %         Deftp(:, :, p) = Deftp(:, :, p) + eye(2, 2);
% %         Defp(:, :, p) = Deftp(:, :, p) * Defp(:, :, p);
%         Strp(:, :, p) = MoE * (Defp(:,:,p) - eye(2));
%     end
% % %     
    % Get deformed volumes
    VolDefp = zeros(N, 1);
    for p = 1:N
        VolDefp(p) = det(Defp(:, :, p)) * Volp(p);
    end
%     
% %     % Save previous energy_total
% %     Energy_total_prev = Energy_total;
    % Compute energy at each particle
    energy = Neohookian(Defp, MoE, PoR, N);
    Energy_total = sum(energy .* Volp);
% %     % Save previous dEnergy_total
% %     dEnergy_total_prev = dEnergy_total;
% %     dEnergy_total = Energy_total - Energy_total_prev;
%     
% %     % stiffness matrix
% %     d2E = (dEnergy_total - dEnergy_total_prev) / dt;
%     
% %     K = zeros(2, 2, N);
% %     for p = 1:N
% %         dx2 = (dXp(p, 1) - dXp_prev(p, 1)) / dt;
% %         dx = dXp(p, 1);
% %         dy = dXp(p, 2);
% %         dy2 = (dXp(p, 2) - dXp_prev(p, 2)) / dt;
% %         K(1, 1, p) = d2E/dx2;
% %         K(1, 2, p) = d2E/(dx*dy);
% %         K(2, 1, p) = K(1, 2, p);
% %         K(2, 2, p) = d2E/dy2;
% %     end
%     
%    Piola = PiolaKirchoff(Defp, MoE, PoR, N);
%     
%     k*dt
%     
    %% STEP 4 - Evaluation of forces on the grid
%     
%     % External forces (gravity)
%     % Calculate gravity directly on grid
%     Feg = zeros(m^2, 2);
%     
%     % F = mg
%     Feg(:, 2) = Mg * gravity;
% %     
%     % internal forces
%     Fig = zeros(m^2, 2);
% % %     
% % % 
% % %     
%     % Compute internal forces
%     for ind = 1:m^2
%         
%         for p = 1:N
%           Fig(ind, :) = Fig(ind, :) - WpgGrad(p,:,ind) * Strp(:, :, p) * VolDefp(p);
% %            Fig(ind, :) = Fig(ind, :) - (Piola(:, :, p) * (Defp(:, :, p).') * (WpgGrad(p, :, ind).')).' * Volp(p);
% %          Fig(GridCellj(p)*m + GridCelli(p)+1, :) = Fig(GridCellj(p)*m + GridCelli(p)+1, :) - ...
% %              (Strp(:, :, p) * BasisGrad1(p, :).').' * VolDefp(p);
% %         Fig((GridCellj(p)+1)*m + GridCelli(p)+1, :) = Fig((GridCellj(p)+1)*m + GridCelli(p)+1, :) - ...
% %             (Strp(:, :, p) * BasisGrad2(p, :).').' * VolDefp(p);
% %         Fig((GridCellj(p)+1)*m + GridCelli(p)+1+1, :) = Fig((GridCellj(p)+1)*m + GridCelli(p)+1+1, :) - ...
% %             (Strp(:, :, p) * BasisGrad3(p, :).').' * VolDefp(p);
% %         Fig(GridCellj(p)*m + GridCelli(p)+1+1, :) = Fig(GridCellj(p)*m + GridCelli(p)+1+1, :) - ...
% %             (Strp(:, :, p) * BasisGrad4(p, :).').' * VolDefp(p);
%         end
%     end
% 
% %     
    %% STEP 5 - Advecting of the Lagrangian Particles
%     

    Gp = zeros(N, 1);
    Gp = Gp -9.8;
    
    Gg = zeros(m^2, 2);
    Gg(:,2) = ProjPG(Gp, Wpg, N, m);

% IMPLICIT UPDATE OPTIMIZATION
% Attempt #1
%     funUpdate = @(dVg) 0.5 * dVg.' * MMg * dVg ...
%         - dt * dVg.' * MMg * Lin2Dto1D(Gg) ...
%         + sum(Neohookian(DeformationGrad(Ug  + dt * Vg + dt * Lin1Dto2D(dVg), ...
%         WpgGrad, N, m), MoE, PoR, N) .* ...
%         DetDefGrad3D(DeformationGrad(Ug + dt * Vg + dt * Lin1Dto2D(dVg),WpgGrad, N, m), N) .* Volp, 1);
% Attempt #2
    funUpdate = @(dVg) 0.5 * dVg.' * MMg * dVg ... % Kinetic Energy
        - dt * dVg.' * MMg * Lin2Dto1D(Gg) ... % Potential External Energy
        + sum(Neohookian(DeformationGrad(Ug + dt * Vg + Lin1Dto2D(dt * dVg), ... % Potential Interal Energy
        WeightsGrad(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Wpg, N, m), GridCells, N, m, h), N, m), MoE, PoR, N) .* ...
        DetDefGrad3D(DeformationGrad(Ug + dt * Vg + Lin1Dto2D(dt * dVg), ...
        WeightsGrad(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Wpg, N, m), GridCells, N, m, h), N, m), N) .* Volp, 1);
% Attempt #3
%     funUpdate = @(dVg) 0.5 * dVg.' * MMg * dVg ... % Kinetic Energy
%         - dt * dVg.' * MMg * Lin2Dto1D(Gg) ... % Potential External Energy
%         + sum(Neohookian(DeformationGrad(Ug + dt * Vg + Lin1Dto2D(dt * dVg), ... % Potential Interal Energy
%         WeightsGrad(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Weights(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Wpg, N, m), GridCells, N, m, h), N, m), GridCells, N, m, h), N, m), MoE, PoR, N) .* ...
%         DetDefGrad3D(DeformationGrad(Ug + dt * Vg + Lin1Dto2D(dt * dVg), ...
%         WeightsGrad(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Weights(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Wpg, N, m), GridCells, N, m, h), N, m), GridCells, N, m, h), N, m), N) .* Volp, 1);
% Attempt #4
%     funUpdate = @(dVg) 0.5 * dVg.' * MMg * dVg ... % Kinetic Energy
%         - dt * dVg.' * MMg * Lin2Dto1D(Gg) ... % Potential External Energy
%         + sum(Neohookian(DeformationGrad(Ug + dt * Vg + Lin1Dto2D(dt * dVg), ... % Potential Interal Energy
%         WeightsGrad(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Weights(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Wpg, N, m), GridCells, N, m, h), N, m), GridCells, N, m, h), N, m), MoE, PoR, N) .* ...
%         DetDefGrad3D(DeformationGrad(Ug + dt * Vg + Lin1Dto2D(dt * dVg), ...
%         WeightsGrad(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Weights(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Wpg, N, m), GridCells, N, m, h), N, m), GridCells, N, m, h), N, m), N) .* Volp, 1);
% Attemp #5
%     funUpdate = @(dVg) 0.5 * dVg.' * MMg * dVg ... % Kinetic Energy
%         - dt * dVg.' * MMg * Lin2Dto1D(Gg) ... % Potential External Energy
%         + sum(Neohookian(DeformationGrad(Ug + dt * Vg + Lin1Dto2D(dt * dVg), ... % Potential Interal Energy
%         WeightsGrad(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Weights(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Weights(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Wpg, N, m), GridCells, N, m, h), N, m), GridCells, N, m, h), N, m), GridCells, N, m, h), N, m), MoE, PoR, N) .* ...
%         DetDefGrad3D(DeformationGrad(Ug + dt * Vg + Lin1Dto2D(dt * dVg), ...
%         WeightsGrad(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Weights(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Weights(Xp + dt * ProjGP(Vg + Lin1Dto2D(dVg), Wpg, N, m), GridCells, N, m, h), N, m), GridCells, N, m, h), N, m), GridCells, N, m, h), N, m), N) .* Volp, 1);

    DVg = zeros(2*m^2, 1);
% Delta Velocity estimate:
%     DVg = (Feg + Fig)./repmat(Mg, [1 2]) / dt;

    options = optimoptions('fminunc', 'MaxFunEvals', 100000000);
    DVg2 = fminunc(funUpdate, DVg, options);
    
    DVg2 = Lin1Dto2D(DVg2);
    %det(DeformationGrad(Ug + dt * Vg + dt * dVg, WpgGrad, N, m)) * Volp)

    % Update acceleration and velocity at each grid node
    
%     % Change all 0's in mass matrix to 1 to avoid dividing by 0
%     for ind = 1:m^2
%         if Mg(ind) == 0
%             Mg(ind) = 1;
%         end
%     end
%     Ag = (Fig + Feg) ./ repmat(Mg, [1 2]);
% %     maxAg = 50000;
% %     for ind = 1:m^2
% %         if Ag(ind, 1) > maxAg
% %             Ag(ind, 1) = maxAg;
% %         elseif Ag(ind, 1) < -maxAg
% %             Ag(ind, 1) = -maxAg;
% %         end
% %         
% %         if Ag(ind, 2) > maxAg
% %             Ag(ind, 2) = maxAg;
% %         elseif Ag(ind, 2) < -maxAg
% %             Ag(ind, 2) = -maxAg;
% %         end
% %         
% %     end
%     %Ag = zeros(m^2, 2);%temporary
%     
%     % Change NaN to 0, since NaN is caused by diving by 0 mass. Just change
%     % acceleration to 0 in that case
%     % Am(isnan(Am)) = 0;
%     
%     Vg = Vg + Ag * dt;
    DVp = ProjGP(DVg2, Wpg, N, m);
    
%     if DVp(3,2) > DVp(1,2) + 0.1
%         break;
%     end;
    
        % constrain velocity of bottom particles
%     for p = 1:8
%         DVp(p , 2) = DVp(2, 2);
%     end
    
    Vp = Vp + DVp;
    Xp = Xp + Vp * dt;
%     
%     % Save previous positions
%     Xp_prev = Xp;
%     dXp_prev = dXp;
%     
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
% %     % Update deformation gradient
% %     Deftp = VGradp * dt;
% % % 
% %     % Compute deformation gradient of each particle directly
% %     for p = 1:N
% %         Deftp(:, :, p) = Deftp(:, :, p) + eye(2, 2);
% %         Defp(:, :, p) = Deftp(:, :, p) * Defp(:, :, p);
% %     end
%     
% %     % Semi-Implicit update equation
% %      fun = @(v_next) -M*v_next - dt^2 * K * v_next ...
% %          + M * v_curr + dt * (f_ext - f_int);
% %     
    % Update position and velocity at each particle
%     for p = 1:N
%         posx = Xp(p, 1);
%         posy = Xp(p, 2);
%         for ind = 1:m^2
%             Vp(p, :) = Vp(p, :) + Wpg(p, ind) * Ag(ind, :) * dt;
%             Xp(p, :) = Xp(p, :) + Wpg(p, ind) * Vg(ind, :) * dt;
%         end
% %         Xp(p, 2) - posy
%         dXp = Xp - Xp_prev;
%     end
    
    % Check if they hit the 2nd bottom of the grid
    for p = 1:N
        if Xp(p, 2) <= 0
            Xp(p, 2) = 0;
        elseif Xp(p, 2) >= m-1
            Xp(p, 2) = m-1;
        end
        if Xp(p, 1) <= 0
            Xp(p, 1) = 0;
        elseif Xp(p, 1) >= m - 1
            Xp(p, 1) = m - 1;
        end
    end
    
%     Radp = sqrt(VolDefp/pi); % Radius of particle domain
    
%     % SIMPLE COLLISION CHECK
%     for p1 = 1:N
%         for p2 = 1:N
%             dist = abs(Xp(p1,:) - Xp(p2,:)) - (Radp(p1) + Radp(p2));
%             if dist < 0
%                 Xp(p1,:) = Xp(p1, : ) - dist * norm(Xp(p1, :) - Xp(p2,:));
%             end
%         end
%     end

    
    % Forces
%     Forg = Lin1Dto2D(DVg)/dt;
%     Forg(:, 2) = Forg(:, 2) + Gg(:, 2);
%     Forg = Forg .* repmat(Mg, [1 2]);
%     Forp = ProjGP(Forg, Wpg, N, m);

    
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
    k
    Vp(3,:)
    DVp(3,:)
    scatter(Xp(:, 1), Xp(:, 2), 15 * 322 * VolDefp, energy, 'filled');
    %quiver(Xp(:, 1), Xp(:, 2), DVp(:,1), DVp(:, 2));
    colorbar;
    grid on;
    set(gca,'xtick',[0:1:m]);
    set(gca,'ytick',[0:1:m]);
%     axis([3 8 -1 4]);
%     axis([-1 11 -1 11]);
  axis([0 m-1 0 m-1]);
    title('MPM');
%     subplot(2,1,2);
%     plot(particle(:, 1), energy(:,1), 'r*');
    
    M(k+1) = getframe;
    %movie2avi(M, 'mpm_movie');

end


