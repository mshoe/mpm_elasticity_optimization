clc;
clear all;

%% Initilization
% all particles have same mass
mass = 0.001;
h = 1; % spacing
pinfl = h/10; % material point radius of influence
dt = 0.01; % time differential
T = 1; % # of time steps

% Field constants
gravity = -9.8;
MoE = 18.83; % Modulus of elasticity
PoR = 0.48; % Poisson's ratio

% mxm grid cells (m+1 x m+1 grid nodes)(stored in 2D matrix)
m = 1;

% N particles (stored in 1D matrix)
N = 4;

% Un - 2D position array
Xp = zeros(N, 2);
Xp_Init = zeros(N, 2);

% Mn - Mass array
Mp = ones(N, 1);

% Vn - Velocity array
Vp = zeros(N, 2);

% Voln - Volume array
Volp = zeros(N, 1);
for p = 1:N
    Volp(p) = (h)^2;
end

% Def - Deformation gradient
Defp = zeros(2, 2, N);
for p = 1:N
    Defp(:,:,p) = eye(2, 2);
end

% Potential energy, Neohookian Elasticity
energy = zeros(N, 1);
energy_total_prev = 0;
energy_total = 0;
dEnergy_total = 0;
dEnergy_total_prev = 0;


%% Basis functions

% % First create unitstep function, which will return 1 if the
% % point is in the grid cell, otherwise 0
% ustp = @(x, y, ind_i, ind_j) (x >= (ind_i-1)*h) .* (x <= ind_i*h) .* (y >=(ind_j-1)*h) .* (y <= ind_j*h);

% Vector that tells you the row and column index of a grid cell
GridCells = zeros((m+1)^2, 2);
gindex = 1;
for gj = 0:m
    for gi = 0:m
        GridCells(gindex, :) = [(gi-0.5), (gj-0.5)];
        gindex = gindex + 1;
    end
end

% Then create the basic alphas
A1 = @(x, ind) (1 - (x - ind * h)/h);
A2 = @(x, ind) (1 + (x - ind * h)/h);
% Then create the basis function for the four corners
% bottom left
Shape1 = @(x, y, ind_i, ind_j) A1(x, ind_i) .* A1(y, ind_j);%(1 - (x - ind_i*h)/h) .* (1 - (y - ind_j*h)/h);% .* ustp(x, y, ind_i, ind_j);
% top left
Shape2 = @(x, y, ind_i, ind_j) A1(x, ind_i) .* A2(y, ind_j+1);%(1 - (x - ind_i*h)/h) .* (1 + (y - ind_j*h)/h);% .* ustp(x, y, ind_i, ind_j);
% top right
Shape3 = @(x, y, ind_i, ind_j) A2(x, ind_i+1) .* A2(y, ind_j+1);%(1 + (x - ind_i*h)/h) .* (1 + (y - ind_j*h)/h);% .* ustp(x, y, ind_i, ind_j);
% bottom right
Shape4 = @(x, y, ind_i, ind_j) A2(x, ind_i+1) .* A1(y, ind_j);%(1 + (x - ind_i*h)/h) .* (1 - (y - ind_j*h)/h);% .* ustp(x, y, ind_i, ind_j);

% Then create the gradient of the basis function for the 4 corners
ShapeGrad1 = @(x, y, ind_i, ind_j) [-1/h * A1(y, ind_j), -1/h * A1(x, ind_i)];
ShapeGrad2 = @(x, y, ind_i, ind_j) [-1/h * A2(y, ind_j), 1/h * A1(x, ind_i+1)];
ShapeGrad3 = @(x, y, ind_i, ind_j) [1/h * A2(y, ind_j+1), 1/h * A2(x, ind_i+1)];
ShapeGrad4 = @(x, y, ind_i, ind_j) [1/h * A1(y, ind_j+1), -1/h * A2(x, ind_i)];

%% Initializing positons
% initialize particle positions
ind_i = 0;
ind_j = 0;
% for ind_p = 1:N
%     Xp(ind_p, :) = [0 + ind_i * h/10, 0 + ind_j * h/10];
%     ind_i = ind_i + 1;
%     if ind_i > 9
%         ind_j = ind_j + 1;
%         ind_i = 0;
%     end
% end
% Keep initial position
Xp_Init = Xp;
Xp(1, :) = [-0.25, -0.25];
Xp(2, :) = [-0.25, 0.25];
Xp(3, :) = [0.25, 0.25];
Xp(4, :) = [0.25, -0.25];
% Xp(:, 2) = (Xp(:,2) + 1);
% Xp(:, 1) = (Xp(:,1) + 0);

%Un(1, :) = [2, 2.2];
%% LOOP


% For visualization
figure
scatter(Xp(:, 1), Xp(:, 2), 15, energy, 'filled');
colorbar;
grid on;
set(gca,'xtick',[0:1:m]);
set(gca,'ytick',[0:1:m]);
axis([-0.5 0.5 -0.5 0.5]);
%axis([-1 21 -1 21]);
title('Object simulation');

for k = 1:T
    %% STEP 0 - Assign particle to a cell
    % Create particle position vectors for x and y
    PX = Xp(:, 1);
    PY = Xp(:, 2);

    % Create grid position vector for each cell (one vector for left
    % x-corner position and bottom y-position, aka the bottome left corner)
    GX = zeros((m+1)^2, 1);
    GY = zeros((m+1)^2, 1);
    gindex = 1;

    for gj = 0:m
        for gi = 0:m
            GX(gindex) = (gi-0.5) * h;
            GY(gindex) = (gj-0.5) * h;
            gindex = gindex + 1;
        end
    end

    % Apply transformations to vectors to make them matrices suitable for
    % operations of P - G
    PX = repmat(PX.', [(m+1)^2 1]);
    PY = repmat(PY.', [(m+1)^2 1]);
    GX = repmat(GX, [1 N]);
    GY = repmat(GY, [1 N]);

    % Create new matrix PGX and PGY that calculates the alpha
    % The alpha is (x - x_0) / dx in the bottom left corner case
    PGX = 1 - (PX - GX)/h;
    PGY = 1 - (PY - GY)/h;

    % Now create a vector that tells which elements are in the range [0 1]
    IX = PGX > 0 & PGX <= 1;
    IY = PGY > 0 & PGY <= 1;

    % Now create the vector that stores the values of the alphas in
    % range, and stores 0 for alphas not in range
    PGX = IX .* PGX;
    PGY = IY .* PGY;

    % Final step to determine if a particle is in a cell
    PG = PGX .* PGY;
    
    % Get a vector that contains the gridCell index and particle index of
    % every particle that is in a grid cell
    [gridCell, particle] = find(PG);
    
%     % Create a vector that counts # of particles in each grid cell
%     GridCellParticleCount = zeros((m+1)^2,1);
%     for p = 1:N
%         GridCellParticleCount(gridCell(p)) = GridCellParticleCount(gridCell(p)) + 1;
%     end
     
    %% STEP 1 - Project data from particles to grid (mass and velocity)
    
    % Create vectors describing the cell which a particle is in
    GridCelli = GridCells(gridCell, 1);
    GridCellj = GridCells(gridCell, 2);
    % Create vectors describing the particle positions
    PX = Xp(:, 1);
    PY = Xp(:, 2);
    
    % Compute basis function over these cells
    Basis1 = Shape1(PX, PY, GridCelli, GridCellj);
    Basis2 = Shape2(PX, PY, GridCelli, GridCellj);
    Basis3 = Shape3(PX, PY, GridCelli, GridCellj);
    Basis4 = Shape4(PX, PY, GridCelli, GridCellj);
    
    % Compute the gradient of the basis functions over these cells,
    % different for each corner of cell
    BasisGrad1 = ShapeGrad1(PX, PY, GridCelli, GridCellj);
    BasisGrad2 = ShapeGrad2(PX, PY, GridCelli, GridCellj);
    BasisGrad3 = ShapeGrad3(PX, PY, GridCelli, GridCellj);
    BasisGrad4 = ShapeGrad4(PX, PY, GridCelli, GridCellj);

    % Reset the mass and velocity matrices
    
    % Mg - grid mass vector. The position of these nodes is implicitly
    % understood because of the uniform grid spacing h.
    Mg = zeros((m+1)^2, 1);
    
    % Vg - grid velocity vector
    Vg = zeros((m+1)^2, 2);
    
    % Ug - grid displacement vector
    UDispg = zeros((m+1)^2, 2);
    
    % Project mass from material points to each corner of each grid cell
    % using matrix multiplication
    Mg((GridCellj)*m + GridCelli + 1) = Mg((GridCellj)*m + GridCelli + 1) + Basis1 .* Mp;
    Mg((GridCellj+1)*m + GridCelli + 1) = Mg((GridCellj+1)*m + GridCelli + 1) + Basis2 .* Mp;
    Mg((GridCellj+1)*m + GridCelli + 1 + 1) = Mg((GridCellj+1)*m + GridCelli + 1 + 1) + Basis3 .* Mp;
    Mg((GridCellj)*m + GridCelli + 1 + 1) = Mg((GridCellj)*m + GridCelli + 1 + 1) + Basis4 .* Mp;
    
    % Project velocity from material points to each corner of each grid cell
    Vg((GridCellj)*m + GridCelli+1, :) = Vg((GridCellj)*m + GridCelli+1, :) ...
        + repmat(Basis1, [1 2]) .* repmat(Mp, [1 2]) .* Vp ...
        ./ repmat(Mg((GridCellj)*m + GridCelli+1), [1 2]);
    Vg((GridCellj+1)*m + GridCelli+1, :) = Vg((GridCellj+1)*m + GridCelli+1, :) ...
        + repmat(Basis2, [1 2]) .* repmat(Mp, [1 2]) .* Vp ...
        ./ repmat(Mg((GridCellj+1)*m + GridCelli+1), [1 2]);
    Vg((GridCellj+1)*m + GridCelli+1+1, :) = Vg((GridCellj+1)*m + GridCelli+1+1, :) ...
        + repmat(Basis3, [1 2]) .* repmat(Mp, [1 2]) .* Vp ...
        ./ repmat(Mg((GridCellj+1)*m + GridCelli+1+1), [1 2]);
    Vg((GridCellj)*m + GridCelli+1+1, :) = Vg((GridCellj)*m + GridCelli+1+1, :) ...
        + repmat(Basis4, [1 2]) .* repmat(Mp, [1 2]) .* Vp ...
        ./ repmat(Mg((GridCellj)*m + GridCelli+1+1), [1 2]);
    
    % Project displacement from material points to each corner of each grid
    % cell
    % Position displacement per particle
    UDispp = Xp - Xp_Init;
    UDispg((GridCellj)*m + GridCelli + 1, :) = UDispg((GridCellj)*m + GridCelli + 1, :) ...
        + repmat(Basis1, [1 2]) .* UDispp;
    UDispg((GridCellj+1)*m + GridCelli + 1, :) = UDispg((GridCellj+1)*m + GridCelli + 1, :) ... 
        + repmat(Basis2, [1 2]) .* UDispp;
    UDispg((GridCellj+1)*m + GridCelli + 1 + 1, :) = UDispg((GridCellj+1)*m + GridCelli + 1 + 1 , :) ...
        + repmat(Basis3, [1 2]) .* UDispp;
    UDispg((GridCellj)*m + GridCelli + 1 + 1, :) = UDispg((GridCellj)*m + GridCelli + 1 + 1, :) ...
        + repmat(Basis4, [1 2]) .* UDispp;
    
    %% STEP 2 - Evaluate field variables at particle points to be used in consitutive model
    
    % Hyper Elasticity Model
    
    % Create matrix for velocity gradients of particles
    VGradp = zeros(2, 2, N);
    
    % Create a vector that stores the grid velocities w.r.t. particles
    % within a grid, one for each corner of a cell
    Vg1 = Vg((GridCellj)*m + GridCelli+1, :);
    Vg2 = Vg((GridCellj+1)*m + GridCelli+1, :);
    Vg3 = Vg((GridCellj+1)*m + GridCelli+1+1, :);
    Vg4 = Vg((GridCellj)*m + GridCelli+1+1, :);
    
    % Compute velocity gradients
    for p = 1:N
        VGradp(:, :, p) = VGradp(:, :, p) + BasisGrad1(p,:).' * Vg1(p,:);
        VGradp(:, :, p) = VGradp(:, :, p) + BasisGrad2(p,:).' * Vg2(p,:);
        VGradp(:, :, p) = VGradp(:, :, p) + BasisGrad3(p,:).' * Vg3(p,:);
        VGradp(:, :, p) = VGradp(:, :, p) + BasisGrad4(p,:).' * Vg4(p,:);
    end
%     
    % Create matrix for displacement gradients of particles
    UDispGradp = ParticleGradient(UDispg, GridCelli, GridCellj, ...
        BasisGrad1, BasisGrad2, BasisGrad3, BasisGrad4, m, N);
    
    % Compute deformation gradients
    for p = 1:N
        Defp(:, :, p) = inv(eye(2) - UDispGradp(:, :, p));
    end
    
    %% STEP 3 - Evaluate the constitutive model
    % Hyper elasticity
    
    % Allocate stress matrix
    Strp = zeros(2, 2, N);
    
    % Update deformation gradient and get stress matrix using the
    % constitutive model
    Deftp = VGradp * dt;

    % Compute deformation gradient of each particle directly
    for p = 1:N
%         % Use gridCell of current particle
%         cell = gridCell(p);
%         % Compute distance between current particle and every other
%         % particle in same grid cell
%         for p2 = 1:N
%             if gridCell(p2) == cell
%                 if norm(Un(p2) - Un(p)) > pinfl
%                     Defp(:, :, p)
%                 end
%             end
%         end
        Deftp(:, :, p) = Deftp(:, :, p) + eye(2, 2);
%         Defp(:, :, p) = Deftp(:, :, p) * Defp(:, :, p);
        Strp(:, :, p) = MoE * (Defp(:,:,p) - eye(2));
    end
    
        % Get deformed volumes
    VolDefp = zeros(N, 1);
    for p = 1:N
        VolDefp(p) = det(Defp(:, :, p)) * Volp(p);
    end
    
    % Save previous energy_total
%     energy_total_prev = energy_total;
    % Compute energy at each particle
    energy = Neohookian(Defp, MoE, PoR, N);
    energy_total = sum(energy .* VolDefp);
    % Save previous dEnergy_total
%     dEnergy_total_prev = dEnergy_total;
%     dEnergy_total = energy_total - energy_total_prev;
    k*dt
    
    
    
    %% STEP 4 - Evaluation of forces on the grid
    
    % External forces (gravity)
    % Calculate gravity directly on grid
    Feg = zeros((m+1)^2, 2);
    
    % F = mg
    Feg(:, 2) = Mg * gravity;
    
    % internal forces
    Fig = zeros((m+1)^2, 2);
    

    
    % Compute internal forces
    for p = 1:N
        Fig(GridCellj(p)*m + GridCelli(p)+1, :) = Fig(GridCellj(p)*m + GridCelli(p)+1, :) - ...
            (Strp(:, :, p) * BasisGrad1(p, :).').' * VolDefp(p);
        Fig((GridCellj(p)+1)*m + GridCelli(p)+1, :) = Fig((GridCellj(p)+1)*m + GridCelli(p)+1, :) - ...
            (Strp(:, :, p) * BasisGrad2(p, :).').' * VolDefp(p);
        Fig((GridCellj(p)+1)*m + GridCelli(p)+1+1, :) = Fig((GridCellj(p)+1)*m + GridCelli(p)+1+1, :) - ...
            (Strp(:, :, p) * BasisGrad3(p, :).').' * VolDefp(p);
        Fig(GridCellj(p)*m + GridCelli(p)+1+1, :) = Fig(GridCellj(p)*m + GridCelli(p)+1+1, :) - ...
            (Strp(:, :, p) * BasisGrad4(p, :).').' * VolDefp(p);
    end
    
    %% STEP 5 - Advecting of the Lagrangian Particles
    
    % Update acceleration and velocity at each grid node
    
    % Change all 0's in mass matrix to 1 to avoid dividing by 0
    for ind = 1:(m+1)^2
        if Mg(ind) == 0
            Mg(ind) = 1;
        end
    end
    Ag = (Feg) ./ repmat(Mg, [1 2]);
    
    % Change NaN to 0, since NaN is caused by diving by 0 mass. Just change
    % acceleration to 0 in that case
    % Am(isnan(Am)) = 0;
    
    Vg = Vg + Ag * dt;
    
    % Update velocity each particle
    Vp = Vp + repmat(Basis1, [1 2]) .* Ag((GridCellj)*m + GridCelli+1, :) * dt;
    Vp = Vp + repmat(Basis2, [1 2]) .* Ag((GridCellj+1)*m + GridCelli+1, :) * dt;
    Vp = Vp + repmat(Basis3, [1 2]) .* Ag((GridCellj+1)*m + GridCelli+1+1, :) * dt;
    Vp = Vp + repmat(Basis4, [1 2]) .* Ag((GridCellj)*m + GridCelli+1+1, :) * dt;
    
    % Save previous position
    Un_prev = Xp;
    % Update position at each particle
    Xp = Xp + repmat(Basis1, [1 2]) .* Vg((GridCellj)*m + GridCelli+1, :) * dt;
    Xp = Xp + repmat(Basis2, [1 2]) .* Vg((GridCellj+1)*m + GridCelli+1, :) * dt;
    Xp = Xp + repmat(Basis3, [1 2]) .* Vg((GridCellj+1)*m + GridCelli+1+1, :) * dt;
    Xp = Xp + repmat(Basis4, [1 2]) .* Vg((GridCellj)*m + GridCelli+1+1, :) * dt;
    
    % Check if they hit the bottom of the grid
    for p = 1:N
        if Xp(p, 2) <= 0
            Xp(p, 2) = 0;
            Vp(p, 2) = 0;
        elseif Xp(p, 2) >= m
            Xp(p, 2) = m;
            Vp(p, 2) = 0;
        end
        if Xp(p, 1) <= 0
            Xp(p, 1) = 0;
            Vp(p, 1) = 0;
        elseif Xp(p, 1) >= m
            Xp(p, 1) = m;
            Vp(p, 1) = 0;
        end
    end
    
    
    
    %% scatter plot
    scatter(Xp(:, 1), Xp(:, 2), 15, energy, 'filled');
    colorbar;
    grid on;
    set(gca,'xtick',[0:1:m]);
    set(gca,'ytick',[0:1:m]);
    axis([-0.5 0.5 -0.5 0.5]);
%    axis([-1 21 -1 21]);
    title('Object simulation');
%     subplot(2,1,2);
%     plot(particle(:, 1), energy(:,1), 'r*');
    
    M(k+1) = getframe;
    %movie2avi(M, 'mpm_movie');

end


