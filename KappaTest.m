% KappaTest.m: Script of the Kappa Test for testing an optimal reconstruction of the 
% Neumann boundary condition from data in a one dimensional Heat Equation (PM).
% Includes both unconstrained and constrained formulations of the problem, in the 
% latter case with an exact (via retraction) and an approximate (via projection onto 
% the subspace tangent to the constraint manifold) enforcement of the constraint.
%
% Tests with respect to the L2 topology, since the Sobolev gradient is simply a 
% filtered version of the L2 element
%
% Author: Pritpal 'Pip' Matharu
% Numerical Analysis, Department of Mathematics
% KTH Royal Institute of Technology
% Date: 2023/11/30
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function KappaTest
% Initialize workspace, load functions, and set default step sizes
HeatOptFuncs
%% Flags optimization problems
% Whether or not the additional energy constraint is added to the model
% Constr = 0: test "unconstrained" gradient 
% Constr = 1: test the normal element, in the "constrained" gradient
Constr = 1;

% Initial condition for PDE
% IC = 0: use zero initial condition
% IC = 1: use non-zero initial condition
IC = 1;

%% Kappa Test parameters
% Range of epsilon for kappa test
ep = logspace(-12, 1, 14);

% Storage Matrix
kappa  = zeros(length(ep), 1);

% Set varying step sizes, for kappa test
d_x = [0.01, 0.001]; % Spatial step size 
d_t = [0.01, 0.001]; % Temporal step size

% Loop through each step size
for m = 1:length(d_t)
for n = 1:length(d_x)    
    %% Set step sizes and variables for current set of parameters
    dx = d_x(n); % Step size in the spatial (x) domain for current set of parameters
    dt = d_t(m); % Step size in the time (t) domain for current set of parameters
    
    % Initialize variables for Heat Equation and Adjoint solvers
    [x, t, ~, u0, phi, phi_true, A] = heat_init(IC, dx, dt);
    
    % Perturbation function
    phi_prime = 4*(t).*exp(-4*t + pi*t);
    
    %% Target Solution
    % Solve Heat Equation and Gradient system for each time step forwards
    % using the true solution for the cost function
    u_b = heatsolve(phi_true, u0, A, dx, dt);
    
    %% Original Functional
    % Solve Heat Equation and cost function, with no perturbation
    [J0, u_r, u] = costfun(phi, u0, A, dx, dt, u_b, t);

    % Energy in system using initial guess
    E0  = 0.5*trapz( t, trapz( x, u.^2 ) );

    % Solve Adjoint system for Riesz representer (backwards in time)
    % Constr = 1: Determine normal element for "constrained" gradient
    % Constr = 0: Determine L2 gradient for "unconstrained" gradient
    RR_J = gradsolve(u, u_r, u_b, A, dx, dt, Constr);

    % Determine Kappa Test denominator, from Riesz form
    kap  = trapz( t, RR_J.*phi_prime );
    
    %% Perturbation Loop
    for p = 1:length(ep)
        % Solve Heat Equation and cost functional with Perturbation
        [J, ~, u] = costfun( phi+ep(p)*phi_prime, u0, A, dx, dt, u_b, t );
                
        % Determine numerator and evaluate Kappa Test
        switch Constr
            % Perform Kappa Test for normal element for "constrained" gradient
            case 1
                % Energy in system with Perturbation
                knumer = 0.5*trapz( t, trapz( x, u.^2 ) ) - E0;

            % Perform Kappa Test for L2 gradient for "unconstrained" gradient
            otherwise
                % Calculate Cost Functional with Perturbation, determine numerator
                knumer = J - J0;

        end % End of Kappa Test Numerator switch

        % Solve for numerator of Kappa Test
        kappa(p) = knumer/(kap*ep(p));

    end % End of perturbation loop p

    %% Plotting
    % Set linear and log difference vectors for plotting Kappa Test
    kappa_linear = kappa;
    kappa_log = log10(abs(kappa - 1));

    % Initialize plotting
    if n == 1 && m == 1
        close all
        % Kappa Test
        figure(1)
        clf;
        ax = gca;
        semilogx(ep, kappa_linear, '--o', 'LineWidth', 1.2)
        hold on

        % Log-based difference Kappa test
        figure(2)
        clf;
        ax1 = gca;
        semilogx(ep, kappa_log, '--o', 'LineWidth', 1.2)
        hold on

        % First legend entry
        Legend{1} = sprintf('$\\Delta x= %.1d, \\Delta t= %.1d$', dx, dt);
    else
        % Kappa Test
        figure(1)
        Linplt = semilogx(ep, kappa_linear, '--o', 'LineWidth', 1.2);

        % Log-based difference Kappa test
        figure(2)
        Logplt = semilogx(ep, kappa_log, '--o', 'LineWidth', 1.2);
        
        % Update legend entry
        Legend{end+1}     = sprintf('$\\Delta x= %.1d, \\Delta t= %.1d$', dx, dt);
    end % End of plotting initialization

    % Fill in Marker for refined time steps
    if m == 2
        Linplt.MarkerFaceColor = Linplt.Color;
        Logplt.MarkerFaceColor = Logplt.Color;
    end
    
end % End of spatial step size loop n

% Reset colouring scheme for plotting
ax.ColorOrderIndex  = 1;
ax1.ColorOrderIndex = 1;
end % End of temporal step size loop m

%% Finalize Plots
% Kappa test plot
figure(1)
hold off
ylim([0 2])
xlim([ep(1) ep(end)])
set(gca,'XTick',[1e-15 1e-10 1e-5 1e0])
set(gca,'XTickLabel',{'$10^{-15}$','$10^{-10}$','$10^{-5}$','$10^{0}$'})
xlabel('$\epsilon$', 'FontSize', 20)
stry = sprintf('$\\kappa_{%i}(\\epsilon)$', Constr+1);
ylabel(stry, 'FontSize', 20)
legend(Legend, 'location', 'best')
title('Kappa Test', 'FontSize', 20)

%  Log-based difference Kappa Test
figure(2)
hold off
ylim([-inf 1])
xlim([ep(1) ep(end)])
set(gca,'YTick',[-6 -5 -4 -3 -2 -1 0])
set(gca,'YTickLabel',{'$10^{-6}$', '$10^{-5}$', '$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$', '$10^{0}$'})
set(gca,'XTick',[1e-15 1e-10 1e-5 1e0])
set(gca,'XTickLabel',{'$10^{-15}$','$10^{-10}$','$10^{-5}$','$10^{0}$'})
xlabel('$\epsilon$', 'FontSize', 20)
stry = sprintf('$|1 - \\kappa_{%i}(\\epsilon)|$', Constr+1);
ylabel(stry, 'FontSize', 20)
legend(Legend, 'location', 'best')
title('Kappa Test - Log Difference', 'FontSize', 18)

toc
end % End of function
