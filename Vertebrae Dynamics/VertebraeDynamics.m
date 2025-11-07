clear
clc

%% Constants
constants;

% Initial transformation matrix and twist
T_0 = TMatExponential3([0 1 0]', pi/2, [0 0 0]') * TMatExponential3([0 0 1]', pi, [0 0 0]'); % Transformation defining the equilibrium rotations
T_translation = [1 0 0 0;...
                 0 1 0 1.5777;...
                 0 0 1 0;...
                 0 0 0 1]; % Transformation matrix defining the translation from equilibrium state
T_0 = T_0 * T_translation; % Total initial transformation
Vb_0 = [0 0 0 0 0 0]'; % initial twist

% Time of Simulation
tspan = [0 0.5];
dt = 0.001;
t = tspan(1) : dt : tspan(2);

% Initialization of arrays
Vb = [Vb_0 zeros(6, length(t) - 1)];
Vb_dot = zeros(6, length(t));

Fb = zeros(6, length(t));

T = cat(3, T_0, zeros(4, 4, length(t) - 1));

for i = 1 : (length(t) - 1)
    % Calculation of body wrench at current step's configuration
    Fb(:, i) = getWrenchT(T(:, :, i), params);

    % Adjoint Vb
    adVb(:, :, i) = [vec2so3(Vb(1:3, i)) vec2so3(Vb(4:6, i)); 
                    zeros(3) vec2so3(Vb(1:3, i))]; % Adjoint

    % Calculation of Vb dot and numerical forward step of Vb
    Vb_dot(:, i) = Gb \ (Fb(:, i) - adVb(:, :, i) * (Gb * Vb(:, i)));
    Vb(:, i + 1) = Vb(:, i) + Vb_dot(:, i) * dt;

    % Evolution of T    
    T(:, :, i + 1) = T(:, :, i) * screw2TMat(Vb(:, i + 1) * dt);

    % Energy in the system
    Ek(i) = 0.5 * Vb(:, i)' * Gb * Vb(:, i); % Kinetic Energy
    Es(i) = getSpringEnergy(T(:, :, i), params); % Spring Potential Energy

end

figure(1)
plot(t(1: end - 1), Ek)
xlabel("Time (s)")
ylabel("Kinetic Energy (kg cm^2 / s^2)")
title("Kinetic Energy over Time")

figure(2)
plot(t(1: end - 1), Es)
xlabel("Time (s)")
ylabel("Spring Energy (kg cm^2 / s^2)")
title("Spring Energy over Time")