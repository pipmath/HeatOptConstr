% Heat_Optimization.m: Script performs an optimal reconstruction of the Neumann 
% boundary condition from data in a one dimensional Heat Equation (PM).
% Includes both unconstrained and constrained formulations of the problem, in the 
% latter case with an exact (via retraction) and an approximate (via projection onto 
% the subspace tangent to the constraint manifold) enforcement of the constraint.
%
% Author: Pritpal 'Pip' Matharu
% Numerical Analysis, Department of Mathematics
% KTH Royal Institute of Technology
% Date: 2023/11/30
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Heat_Optimization
% Initialize workspace, load functions, and set default step sizes
HeatOptFuncs
%% Flags optimization problems
% Whether or not the additional energy constraint is added to the model
% Constr = 0: test "unconstrained" gradient
% Constr = 1: test the normal element, in the "constrained" gradient
Constr = [0, 1];

% Initial condition for PDE
% IC = 0: use zero initial condition
% IC = 1: use non-zero initial condition
IC     = 1;

%% Parameters for optimiziation
% Tolerance for optimization
tol        = 1e-7;
% Maximum number of iterations
max_iter   = 1000;
% Length scale for solving Sobolev gradient, which ensures suitable regularization
L          = 0.01;
% Starting step length for line search
tau0       = 0.1;
% Bracketing Factor, to enlarge & shrink the fminbnd interval
brack_fact = 2;
% Clearing frequency for Conjugate Gradients (FQ = 1 mean standard gradient descent)
FQ         = 10;

% Loop through both "unconstrained" and "constrained" problems
for l = 1:length(Constr)
    fprintf('Constraint = %d \n', Constr(l))

    % Whether or not the retraction operator will be used
    % Retrac = 0: no retraction operator
    % Retrac = 1: if constraint is enforced and initial condition is zero, use retraction
    Retrac = Constr(l) == 1 && IC ==0;

    % Frequency clearing for conjugate gradient method momentum term beta
    % NOTE: We only use conjugate gradient method when constraint is not applied
    fq = FQ(Constr(l) ~= 1); fq(Constr(l) == 1) = 1;

    % Initialize variables for Heat Equation and Adjoint solvers
    [x, t, t_len, u0, phi, phi_true, A] = heat_init(IC, dx, dt);

    %% Target Solution
    % Solve Heat Equation and Gradient system for each time step forwards
    % using the true solution for the cost function
    [u_b, u] = heatsolve(phi_true, u0, A, dx, dt);

    % Energy in system
    Etrue    = 0.5*trapz( t, trapz( x, u.^2 ) );

    %% ORIGINAL FUNCTIONAL - Calculate values for initial guess
    % Calculate Cost Functional (no need for projection)
    [J1, u_r, u] = costfun(phi, u0, A, dx, dt, u_b, t); % Value of initial cost function
    J2           = 0;                                   % Dummy value

    % Energy in system defining manifold
    E1 = 0.5*trapz(t, trapz(x, u.^2));

    fprintf('True Energy E0 = %d \n', Etrue)
    fprintf('Energy E1      = %d \n', E1)
    fprintf('Iteration: %i, Cost Function: %d \n', 0, J1)

    %% Storage arrays for iterations
    ur    = zeros(t_len, max_iter); % RHS solution
    Phi   = zeros(t_len, max_iter); % Neumann Boundary Function
    J     = zeros(max_iter, 1);     % Cost Functional
    grad  = zeros(t_len, max_iter); % Standard gradient applied
    CG    = zeros(t_len, max_iter); % Conjugate gradients (if not being used, CG = grad)
    Tau   = zeros(max_iter, 1);     % Step lengths for applying gradient
    E     = zeros(max_iter, 1);     % Energy in system
    Exit  = zeros(max_iter, 1);     % Exit flags from line search method
    Num   = zeros(max_iter, 1);     % Number of function evaluations
    bval  = zeros(max_iter, 1);     % Momentum term from Polak-Ribiere method

    % Initialization iteration counter
    p = 1;

    % Store values for initial guess
    ur(:, p)     = u_r;
    Phi(:, p)    = phi;
    J(p, 1)      = J1;
    Tau(p, 1)    = tau0;
    E(p, 1)      = E1;
    Eiter        = E1;

    % Constrain is enforced but no retraction operator is available, we limit magnitude
    % of the step size in the bracketing method
    if IC && Constr(l) == 1
        shrnkcond = sprintf('tau0 > 1e-6 && abs(Etest/Eiter - 1) > 0.001');
        expndcond = sprintf('tau0 < 1e6  && abs(Etest/Eiter - 1) < 0.001');
    else
        % Only ensure step size does not get too small or large
        shrnkcond = sprintf('tau0 > 1e-6');
        expndcond = sprintf('tau0 < 1e6');
    end

    %% Loop Iteration
    % Iterate until condition is met
    while ( (abs(J1 - J2))/(abs(J1)) > tol && p <= max_iter )

        %% Determine Sobolev Gradient
        % Solve Adjoint system and determine L2 gradient (backwards in time)
        grad_L2 = gradsolve(u, u_r, u_b, A, dx, dt, 0);

        % Determine Sobolev gradient, for additional regularity
        del_J   = Sobolevgrad(grad_L2, dt, L);

        % Sobolev gradient with or without constraint
        switch Constr(l)
            % Add additional energy constraint on system, via additional adjoint solve
            case 1
                % Solve Adjoint system for constraint (backwards in time)
                normal_L2 = gradsolve(u, u_r, u_b, A, dx, dt, 1);

                % Determine Sobolev gradient, for additional regularity
                N_J       = Sobolevgrad(normal_L2, dt, L);

                % Derivative of gradient and normal vector for inner product
                dJ = derivfunc(del_J, dt);
                dN = derivfunc(N_J, dt);
                % Numerator of the projection term
                numer_proj = trapz( t, N_J.*del_J ) + (L^2)*trapz( t, dN.*dJ );
                % Denominator of the projection term
                denom_proj = trapz( t, N_J.*N_J )   + (L^2)*trapz( t, dN.*dN );
                % Enforce energy constraint, and obtain new gradient
                % Apply steepest DESCENT direction
                del_J = -( del_J - (numer_proj/denom_proj) * N_J );

                % Determine regular Sobolev gradient
            otherwise
                % Take negative value, to go in the DESCENT direction
                del_J = -del_J;
        end

        %% Conjugate Gradient method
        % Use Conjugate gradient method, with the Polak-Ribiere method including
        % frequency resetting
        if p >= 2 && mod(p, fq)~=0
            % Using the Polak Ribere Method
            bPR = ( del_J'*(del_J - delk) )/( norm(delk)^2 );

            % Use value to create the conjugate gradient
            PR  = del_J + bPR.*pr;
        else
            % Frequency clearing to reset the conjugate-gradient procedure
            bPR = 0;

            % Use gradient (previous gradient cleared)
            PR  = del_J;
        end % End of conjugate gradient statement

        % Ensure that the cost does not equal zero (only for the first iteration)
        if p ~= 1
            % Set Cost functional equal to current iteration
            J1 = J2;
        end

        % Bracketing the minimum from the right, for the Brent method
        Feval = 1;        % Count number of function calls
        Gtau  = [0.0 J1]; % Store step length and respective value of cost function

        % Calculate the cost function for the given step length and gradient
        if Retrac
            % Calculate cost functional and flux with retraction
            J2      =    costfun( phi + tau0*PR, u0, A, dx, dt, u_b, t, x, E1 );
            phitest = RetracConstr( phi + tau0*PR, u0, A, dx, dt, t, x, E1 );
        else
            % Calculate cost function and flux without retraction
            J2      = costfun(phi + tau0*PR, u0, A, dx, dt, u_b, t);
            phitest = phi + tau0*PR;
        end

        % Solve Heat Equation (forward in time)
        [~, u] = heatsolve(phitest, u0, A, dx, dt);
        % Energy in system
        Etest  = 0.5*trapz( t, trapz(x, u.^2) );


        %% Modifying the interval for step lengths via bracketing
        if J2 > J1
            % If the Cost Functional at tau0 is GREATER than the current iteration of
            % the evaluated Cost Functional, we try to shrink the bracket interval that
            % we evaluate the fminbnd function

            % Shrink bracket interval, until appropriate value is found
            while J2 > J1 && eval(shrnkcond)
                % Shrink by a constant bracketing factor
                tau0 = tau0/brack_fact;

                % Calculate the cost function for the given step length and gradient
                if Retrac
                    % Retraction
                    J2 = costfun(phi + tau0*PR, u0, A, dx, dt, u_b, t, x, E1);
                else
                    J2 = costfun(phi + tau0*PR, u0, A, dx, dt, u_b, t);
                end

                % Update the number of times the function is called
                Feval = Feval + 1;
                % Update store the step length used, and the value of cost function
                Gtau  = [Gtau; [tau0, J2] ];
            end % End of while statement for shrinking interval

            % To ensure that we did not shrink the bracket too much!
            tau0 = tau0*brack_fact;

        else
            % If the Cost Functional at tau0 is less than the current iteration of the
            % evaluated Cost Functional, we try to expand the bracket that we evaluate
            % the fminbnd function at because there may exist an even larger interval
            % that could give us an even smaller (and better) step length for the Cost
            % Functional!

            % Expand bracket interval, until appropriate value is found
            while J2 < J1 && eval(expndcond)
                % Expand by a constant bracketing factor
                tau0 = tau0*brack_fact;

                % Calculate the cost function for the given step length and gradient
                if Retrac
                    % Calculate cost functional and flux with retraction
                    J2      =    costfun( phi + tau0*PR, u0, A, dx, dt, u_b, t, x, E1 );
                    phitest = RetracConstr( phi + tau0*PR, u0, A, dx, dt, t, x, E1 );
                else
                    J2      = costfun( phi + tau0*PR, u0, A, dx, dt, u_b, t );
                    % Calculate cost function and flux without retraction
                    phitest = phi + tau0*PR;
                end
                % Expanding interval may lead to large deviations to constraint, so limit
                % drift by checking the relative difference in the constraint

                % Solve Heat Equation (forward in time)
                [~, u] = heatsolve(phitest, u0, A, dx, dt);
                % Energy in system
                Etest  = 0.5*trapz(t, trapz(x, u.^2));

                % Update the number of times the function is called
                Feval  = Feval + 1;
                % Update store the step length used, and the value of cost function
                Gtau   = [Gtau; [tau0, J2] ];
            end % End of while statement for expanding interval
        end % End of if statement for bracketing

        %% Determine optimal value for steplength
        % Determine step length along the gradient, using line search method
        if Retrac
            % Retraction
            [tau0, J2, Gtau, exitflag, output] = costfunmin(phi, PR, u0, A, dx, dt, u_b, t, tau0, Gtau, x, E1);
        else
            [tau0, J2, Gtau, exitflag, output] = costfunmin(phi, PR, u0, A, dx, dt, u_b, t, tau0, Gtau);
        end
        % Update number of function evaluations
        Feval = Feval + length(Gtau);

        %%  Update Neumann Boundary condition
        phi = phi + tau0*PR;
        % Use retraction operator, if available
        if Retrac
            phi = RetracConstr(phi, u0, A, dx, dt, t, x, E1);
        end

        % Solve Heat Equation (forward in time)
        [u_r, u] = heatsolve(phi, u0, A, dx, dt);
        % Energy in system
        Eiter    = 0.5*trapz(t, trapz(x, u.^2));

        % Display the current iteration (based on frequency clearing)
        fprintf('Iteration: %i, Cost Function: %d \n', p, J1)

        %% Storing Values
        % Current iteration values
        ur(:, p+1)     = u_r;
        Phi(:, p+1)    = phi;
        J(p+1, 1)      = J2;
        Tau(p+1, 1)    = tau0;
        E(p+1, 1)      = Eiter;

        % For diagnostics, save results from line search
        Exit(p, 1)     = exitflag;
        Num(p, 1)      = Feval;
        bval(p, 1)     = bPR;
        CG(:, p)       = PR;
        grad(:, p)     = del_J;
        % Diagnostics for minimization process
        if ( exitflag ~= 1 )
            fprintf('  PROBLEM: exitflag=%d \n', exitflag);
        end

        % Store values for the Polak Ribere method
        delk = del_J;
        pr   = PR;

        % Increment iteration counter
        p = p+1;

    end % End of while loop for minimization

    %% Finalized and store values
    % Storing optimal values
    ur  = ur(:, 1:p);
    Phi = Phi(:, 1:p);
    J   = J(1:p, 1);
    Tau = Tau(1:p, 1);
    E   = E(1:p, 1);

    % Storing diagnostic quanitites
    Exit = Exit(1:p-1, 1);
    Num  = Num(1:p-1, 1);
    bval = bval(1:p-1, 1);
    CG   = CG(:, 1:p-1);
    grad = grad(:, 1:p-1);

    fprintf('Final Energy E1 = %d \n', Eiter)

    % Stop timing
    time = toc;

    % Display Optimization information
    fprintf('\n Number of iterations: %i, Final Value: %i, Energy Ratio: %d, Time: %i \n', p, J(end), E(end)/E(1), time)

    % String for saving variables
    str = sprintf('OptimizationSummary_IC%d_Constr%d.mat', IC, Constr(l));
    % Save information
    save(str, 't', 'u_b', 'ur', 'u_r', 'Phi', 'phi_true', 'phi', 'p', 'time', 'J', 'E', 'Tau', 'Phi', 'ur', 'grad', 'Exit', 'Num', 'CG', 'IC')

    % Store separate variables for ease of plotting
    if Constr(l) ~= 1
        % Unconstrained Problem
        u_r0   = u_r;
        phi0   = phi;
        Jvals0 = J;
        Evals0 = E;
    else
        % Constrained Problem
        u_r1   = u_r;
        phi1   = phi;
        Jvals1 = J;
        Evals1 = E;
    end

end % End of main loop

%% Plot Results from optimization
plt_optvals(t, u_r0, u_r1, phi0, phi1, Jvals0, Jvals1, Evals0, Evals1, u_b, ur, phi_true, Phi, IC);

end % End of function

%% ******************************* PLOTTTING FUNCTION *******************************
function plt_optvals(t, u_r0, u_r1, phi0, phi1, J0, J1, E0, E1, u_b, ur, phi_true, Phi, IC)
% Various linestyles for plotting
line_sty = {'-.', '--', '*', ':', ':'};

% Legend Entries
Legend{1}  = sprintf('$\\bar{\\phi}$, True Function');
Legend{2}  = sprintf('$\\phi_{0}$, Initial Guess');
Legend{3}  = sprintf('$\\check{\\phi} \\in \\mathcal{S}$, Problem %d.A', IC+1);
Legend{4}  = sprintf('$\\check{\\phi} \\in \\mathcal{M}$, Problem %d.B', IC+1);

Legend2{1} = sprintf('$\\phi \\in \\mathcal{S}$, Problem %d.A', IC+1);
Legend2{2} = sprintf('$\\phi \\in \\mathcal{M}$, Problem %d.B', IC+1);

%% Plot Right Time history
figure(1);clf;
ax1 = gca;
hold on
box on
plot(t(1:end), u_b(1:end), '-k')
ax1.ColorOrderIndex = 5;
plot(t(1:25:end), ur(1:25:end, 1), line_sty{5}, 'LineWidth', 2.0, 'MarkerSize', 10)
ax1.ColorOrderIndex = 1;
plot(t(1:end), u_r0(1:end), line_sty{1}, 'LineWidth', 1.2, 'MarkerSize', 10)
ax1.ColorOrderIndex = 2;
plot(t(1:end), u_r1(1:end), line_sty{2}, 'LineWidth', 1.2, 'MarkerSize', 10)
hold off
set(gca,'FontSize', 14)
set(gca,'TickLabelInterpreter','latex')
xlabel('$t$', 'Interpreter','LaTex', 'FontSize', 20)
ylabel('$u(\phi(t))|_{x=b}$', 'Interpreter','LaTex', 'FontSize', 20)
legend(Legend, 'Interpreter', 'LaTex', 'location', 'best')
title('Temperature of the Right Endpoint', 'FontSize', 20)

%% Plot Left Neumann boundary condition
figure(2);clf;
ax2 = gca;
hold on
box on
plot(t, phi_true, '-k')
ax2.ColorOrderIndex = 5;
plot(t(1:25:end), Phi(1:25:end, 1), line_sty{5}, 'LineWidth', 2.0, 'MarkerSize', 10)
ax2.ColorOrderIndex = 1;
plot(t(1:end), phi0(1:end, 1), line_sty{1}, 'LineWidth', 1.2, 'MarkerSize', 10)
ax2.ColorOrderIndex = 2;
plot(t(1:end), phi1(1:end, 1), line_sty{2}, 'LineWidth', 1.2, 'MarkerSize', 10)
hold off
xlabel('$t$', 'FontSize', 20)
ylabel('$\phi(t)$', 'FontSize', 20)
ylim([-20 20])
legend(Legend, 'location', 'best')
title('Left Neumann Boundary Condition', 'FontSize', 20)

%% Plot Cost functional iterations
figure(3);clf;
ax3 = gca;
hold on
box on
ax3.ColorOrderIndex = 1;
semilogy(0:length(J0)-1, J0(1:length(J0))/J0(1), line_sty{1}, 'LineWidth', 1.2, 'MarkerSize', 10)
ax3.ColorOrderIndex = 2;
semilogy(0:length(J1)-1, J1(1:length(J1))/J1(1), line_sty{2}, 'LineWidth', 1.2, 'MarkerSize', 10)
hold off
xlabel('$n$', 'FontSize', 20)
ylabel('$\mathcal{J}(\phi^{(n)})/\mathcal{J}(\phi_0)$', 'FontSize', 20)
set(ax3, 'yscale', 'log')
legend(Legend2, 'location', 'best')
title('Normalized Error Functional', 'FontSize', 20)

%% Plot Energy iterations
Xval = 50;
Xend = length(E0);
TileS = 6;
TileN = 9;
% Create figure and tile structure
hFig = figure(4);
figtile = tiledlayout(1,TileN,'TileSpacing','none');
% Setup axes and ticks
bgAx = axes(figtile,'XTick',[],'YTick',[],'Box','off');
bgAx.Layout.TileSpan = [1 TileN];
color = get(hFig,'Color');
set(gca,'YColor',color,'TickDir','out')
ax41 = axes(figtile);
ax41.Layout.Tile = 1;
ax41.Layout.TileSpan = [1 TileS];
hold on
set(ax41, 'yscale', 'linear')
set(ax41, 'xscale', 'linear')
% Set limits
ax41.Box = 'off';

ax41.ColorOrderIndex = 1;
plot(0:Xval, E0(1:Xval+1)/E0(1), line_sty{1}, 'LineWidth', 1.2, 'MarkerSize', 10)
ax41.ColorOrderIndex = 2;
plot(0:Xval, E1(1:Xval+1)/E1(1), line_sty{2}, 'LineWidth', 1.2, 'MarkerSize', 10)

% Create second plot
ax42 = axes(figtile);
ax42.Layout.Tile = TileS+1;
ax42.Layout.TileSpan = [1 TileN-TileS];
hold on
set(ax42, 'yscale', 'linear')
set(ax42, 'xscale', 'linear')

% Set second plot axis and ticks values
ax42.YAxis.Visible = 'off';
ax42.XAxis.Scale = 'linear';
ax42.Box = 'off';

ax42.ColorOrderIndex = 1;
plot(Xval:Xend-1, E0(Xval+1:Xend)/E0(1), line_sty{1}, 'LineWidth', 1.2, 'MarkerSize', 10)
ax42.ColorOrderIndex = 2;
plot(Xval:Xend-1, E1(Xval+1:Xend)/E0(1), line_sty{2}, 'LineWidth', 1.2, 'MarkerSize', 10)

xline(ax42, Xval, ':', 'color', [.5 .5 .5])
xline(Xend)
yline(ax41, 1.2)
yline(ax42, 1.2)

set(ax41,'FontSize', 14)
set(ax42,'FontSize', 14)
figtile.XLabel.String = '$n$';
figtile.YLabel.String = '$[E(\cdot;\phi^{(n)})]_T/[E(\cdot;\phi_0)]_T$';
figtile.XLabel.Interpreter = 'latex';
figtile.YLabel.Interpreter = 'latex';
figtile.XLabel.FontSize = 20;
figtile.YLabel.FontSize = 20;

% Set labels
linkaxes([ax41 ax42], 'y')
tmpy = ylim;
ylim([tmpy(1) 1.2])
xlim(ax41,[0 Xval])
xlim(ax42,[Xval Xend])
ax42.XTickLabelRotation = 0;
axtick1 = [0 10 20 30 40];
axtick2 = [50 250 500 750 1000];
set(ax41, 'xtick', axtick1) % Remove overlapping value
set(ax42, 'xtick', axtick2) % Remove overlapping value

Lgnd = legend(Legend2, 'location', 'best');
Lgnd.AutoUpdate = 'off';
title('Normalized Energy Constraint', 'FontSize', 20)

end
