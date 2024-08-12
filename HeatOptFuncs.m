% HeatOptFuncs.m: File containing functions for optimal reconstruction of the Neumann
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
%% Loading file sets up the workspace, sets default values, and loads the functions
function HeatOptFuncs
%% Initialize workspace
clear
close all
clc
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',14);
disp(datetime)
tic % Start timing function

%% Setup default domain and discretization for solvers
% Default Step Sizes
dx = 0.01;   % Step size in the spatial (x) domain
dt = 0.001;  % Step size in the time (t) domain

%% Assign function calls and time steps in file
assignin('caller', 'dt',           dt);
assignin('caller', 'dx',           dx);

assignin('caller', 'heatsolve',    @heatsolve);
assignin('caller', 'gradsolve',    @gradsolve);
assignin('caller', 'Sobolevgrad',  @Sobolevgrad);
assignin('caller', 'derivfunc',    @derivfunc);
assignin('caller', 'costfun',      @costfun);
assignin('caller', 'costfunmin',   @costfunmin);
assignin('caller', 'RetracConstr', @RetracConstr);
assignin('caller', 'heat_init',    @heat_init);

%% *********************************** FUNCTIONS ************************************
% -------------------------------------------------------------------------
% FUNCTION: heat_init
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2023/11/30
%
% Initialize the Heat Equation solver
%
% INPUT
% IC ........ Flag to indicate type of initial condition to use
% dx ........ Step size in the spatial domain
% dt ........ Step size in the time domain
%
% OUTPUT
% x ............ Spatial domain, discretized equispaced (0 <= x <= 1)
% t ............ Time domain, discretized equispaced    (0 <= t <= 1)
% t_len ........ Reference length for time array
% u0 ........... Initial Condition Function
% phi .......... Initial guess for Neumann boundary condition at left endpoint (x = a)
% phi_true ..... True Neumann boundary condition for comparison
% A ............ Matrix for solving PDES. Crank-Nicolson method is used. 
%
% FORMAT
% [x, t, t_len, u0, phi, phi_true, A] = heat_init(IC, dx, dt)
%
% -------------------------------------------------------------------------
    function [x, t, t_len, u0, phi, phi_true, A] = heat_init(IC, dx, dt)
        %% Setup discretization domain for solvers
        % Boundaries
        x_a = 0;  % Left Boundary Point  (x = a)
        x_b = 1;  % Right Boundary Point (x = b)
        t_i = 0;  % Initial time         (t = 0)
        t_f = 1;  % Final time           (t = T)
        x   = (x_a:dx:x_b)'; % Create spatial array, with equally spaced intervals
        t   = (t_i:dt:t_f)'; % Create time array, with equally spaced intervals

        % Reference lengths
        x_len = length(x); % Reference for spatial array
        t_len = length(t); % Reference for time array
        h     = dt/(dx^2); % Magnitude reference

        %% Functions for PDE
        % Set Initial Condition
        switch IC
            case 0
                % Use Zero IC
                u0 = zeros(size(x));
            case 1
                % Use Non-zero IC
                u0  = 10*cos(pi*x);
        end

        % Initial Guess for Neumann Boundary Function
        phi = 18*sin((pi/2)*t).*exp(-4*t);

        % True Neumann boundary function, as a reference solution
        phi_true = exp(-4)*t.*(-1000*cos(15*t*(pi/2))+18);

        %% Matrix Setup (Crank-Nicolson method)
        % Determine values in matrix for PDE solvers. The matrix A is constant.
        A_U = sparse(1:x_len-1, 2:x_len, (-0.5*h),x_len, x_len); % Upper Diagonal Entries
        A_D = sparse(1:x_len, 1:x_len, (1 + h),x_len, x_len);    % Main  Diagonal Entries
        A_L = sparse(2:x_len, 1:x_len-1, (-0.5*h),x_len, x_len); % Lower Diagonal Entries
        A   = A_U + A_D + A_L;

        A(1, 1) = -3; % Left Boundary Condition, First Row
        A(1, 2) =  4; % Left Boundary Condition, First Row
        A(1, 3) = -1; % Left Boundary Condition, First Row
        A(x_len, x_len)     = -3; % Right Boundary Condition, Last Row
        A(x_len, x_len - 1) =  4; % Right Boundary Condition, Last Row
        A(x_len, x_len - 2) = -1; % Right Boundary Condition, Last Row

    end % End of function heat_init

% -------------------------------------------------------------------------
% FUNCTION: heatsolve
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2023/11/30
%
% Solves the Heat Equation and the outputs the right side values
%
% INPUT
% phi ....... Current Neumann Boundary condition at left endpoint (x = a)
% u0 ........ Initial Condition Function
% A ......... Matrix for time stepping
% dx ........ Step size in the spatial domain
% dt ........ Step size in the time domain
%
% OUTPUT
% u_r ....... Heat Equation solved at right endpoint (x = b)
% u ......... Solution to Heat Equation
%
% FORMAT
% [u_r, u] = heatsolve(phi, u0, A, dx, dt)
%
% -------------------------------------------------------------------------
    function [u_r, u] = heatsolve(phi, u0, A, dx, dt)

        % Reference lengths
        x_len = length(u0);  % Reference length for spatial dimension
        t_len = length(phi); % Reference length for time dimension
        h     = dt/(dx^2);   % Magnitude reference

        % Storage Matrices
        u     = zeros(x_len, t_len);
        y_PDE = zeros(x_len, 1);

        % Apply Initial Condition at t = 0
        u(:, 1) = u0;

        % Solve Heat equation system for each time step forwards
        for n = 1:t_len-1
            % Solve time step forwards
            y_PDE(1)         = 3*u(1, n) - 4*u(2, n) + u(3, n) + 2*dx*( phi(n) + phi(n+1) );
            y_PDE(2:x_len-1) = 0.5*h*u(1:x_len-2, n) + u(2:x_len-1, n) - h*u(2:x_len-1, n) + 0.5*h*u(3:x_len, n);
            y_PDE(x_len)     = u(x_len-2, n) - 4*u(x_len-1, n) + 3*u(x_len, n);

            % Solve next time step
            u(:, n+1) = A \ y_PDE;
        end

        % Heat equation solved at right limit
        u_r = u(end, :)';

    end % End of function heatsolve

% -------------------------------------------------------------------------
% FUNCTION: gradsolve
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2023/11/30
%
% Solves the Adjoint System and determines the L2 gradient
%
% INPUT
% u ......... Heat Equation solution
% u_r ....... Heat Equation solved at right limit
% u_b ....... Target function
% A ......... Matrix for next time step
% dx ........ Step size in the spatial domain
% dt ........ Step size in the time domain
% Constr .... Flag to indicate whether we solve gradient or normal element
%
% OUTPUT
% del_J ..... L2 Gradient of system
%
% FORMAT
% [grad_L2] = gradsolve(u, u_r, u_b, A, dx, dt, Constr)
%
% -------------------------------------------------------------------------
    function [RR_L2] = gradsolve(u, u_r, u_b, A, dx, dt, Constr)

        % Reference lengths
        t_len = length(u_r); % Reference length for time dimension
        x_len = length(A);   % Reference length for spatial dimension
        h     = dt/(dx^2);   % Magnitude reference

        % Storage Matrices
        v     = zeros(x_len, t_len);
        y_adj = zeros(x_len, 1);

        % Evaluation switch in adjoint solver, depending on gradient or normal element solve
        switch Constr
            % Modify RHS of Heat equation PDE, for normal element adjoint solve
            case 1
                eval1 = sprintf('(dt/2)*(u(2:x_len-1, m)+u(2:x_len-1, m-1))');
                eval2 = sprintf('0');

            % Modify boundary condition at right end point (x=b), for gradient adjoint solve
            otherwise
                eval1 = sprintf('0');
                eval2 = sprintf('2*dx*((u_r(m-1) + u_r(m)) - (u_b(m-1) + u_b(m)))');

        end % End of evaluation switch


        % Solve adjoint system for each time step BACKWARDS
        for m = t_len:-1:2
            % Solve time step backwards
            y_adj(1)         = 3*v(1, m) - 4*v(2, m) + v(3, m);
            y_adj(2:x_len-1) = 0.5*h*v(1:x_len-2, m) + v(2:x_len-1, m) - h*v(2:x_len-1, m) + 0.5*h*v(3:x_len, m) + eval(eval1);
            y_adj(x_len)     = 3*v(x_len, m) - 4*v(x_len-1, m) + v(x_len-2, m) - eval(eval2);

            % Solve previous time step
            v(:, m-1) = A \ y_adj;
        end

        % L2 Gradient of system
        RR_L2 = -v(1, :)';

    end % End of function

% -------------------------------------------------------------------------
% FUNCTION: Sobolevgrad
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2023/11/30
%
% Solves Sobolev (H^1) gradient of the Heat Equation, given the L^2 gradient
%
% INPUT
% grad_L2 ....... L2 gradient in time, solved using adjoint analysis
% dt ............ Step size in the time domain
% L ............. Sobolev parameter, to ensure regularity
%
% OUTPUT
% del_J ..... Sobolev Gradient of system
%
% FORMAT
% [del_J] = Sobolevgrad(grad_L2, dt, L)
%
% -------------------------------------------------------------------------
    function [del_J] = Sobolevgrad(grad_L2, dt, L)

        % Reference length
        t_len = length(grad_L2); % Reference length for time dimension

        % Storage vector
        del_J = zeros(t_len, 1);

        % Create Matrix, in order to solve the PDE resulting from the Sobolev Gradient System
        P = sparse(1:t_len-3, 2:t_len-2, ((-1)*((L/dt)^2)), t_len-2, t_len-2); % Upper Diagonal Entries
        Q = sparse(1:t_len-2, 1:t_len-2, (2*((L/dt)^2) + 1),t_len-2, t_len-2); % Diagonal Entries
        R = sparse(2:t_len-2, 1:t_len-3, ((-1)*((L/dt)^2)), t_len-2, t_len-2); % Lower Diagonal Entries
        A = P+Q+R;

        % Solve the gradient for the middle terms
        grad_H1 = A \ (grad_L2(2:end-1));

        % Sobolev gradient of system
        del_J(2:end-1) = grad_H1;

    end % End of function

% -------------------------------------------------------------------------
% FUNCTION: derivfunc
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2023/11/30
%
% Determine first derivative using a simple symmetric difference quotient
%
% INPUT
% f ............. Function to differentiate with respect to time
% dt ............ Step size in the time domain
% t_len ......... Reference length for time dimension
%
% OUTPUT
% df ..... Derivative of function f with respect to time
%
% FORMAT
% [df] = derivfunc(f, dt)
%
% -------------------------------------------------------------------------
    function df = derivfunc(f, dt)

        % Reference length
        t_len = length(f); % Reference length for time dimension

        % Storage vector
        df = zeros(t_len, 1);

        % Determine derivative
        df(2:end-1) = ( f(3:end) - f(1:end-2) )/(2*dt);

    end % End of function derivfunc

% -------------------------------------------------------------------------
% FUNCTION: costfun
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2023/11/30
%
% Calculates the Cost Functional for the given Heat Equation system
%
% USES THE RIGHT BOUNDARY
%
% INPUT
% phi ....... Current Neumann Boundary condition
% u0 ........ Initial Condition Function
% A ......... Matrix for next time step
% dx ........ Step size in the spatial domain
% dt ........ Step size in the time domain
% u_b ....... Ideal function (Right boundary)
% t ......... Time domain
%
% OUTPUT
% J ......... Cost Functional
% u_r ....... Heat Equation solved at right limit (optional)
% u ......... Solution to Heat Equation (optional)
% FORMAT
% [J, u_r] = costfun(phi, u0, A, dx, dt, u_b, t)
% -------------------------------------------------------------------------
    function [J, u_r, u] = costfun(phi, u0, A, dx, dt, u_b, t, x, E1)

        % Retraction to constraint manifold
        if nargin == 10
            phi = RetracConstr(phi, u0, A, dx, dt, t, x, E1);
        end

        % Heat equation solved at right limit
        [u_r, u] = heatsolve(phi, u0, A, dx, dt);

        % Calculate Cost Functional
        J = 0.5*trapz( t, ( (u_r - u_b).^2 ) );

    end % End of function costfun


% -------------------------------------------------------------------------
% FUNCTION: costfunmin
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2023/11/30
%
% Minimizing the cost functional, to obtain the optimal step length to take
%
% INPUT
% phi ....... Current Neumann Boundary condition
% grad ...... Conjugate gradient
% u0 ........ Initial Condition Function
% A ......... Matrix for next time step
% dx ........ Step size in the spatial domain
% dt ........ Step size in the time domain
% u_b ....... Ideal function (Right boundary)
% t ......... Time domain
%
% MAINTAIN A RECORD OF VALUES
% tau ....... Initial step length
% history ... History, containing step length, and value of cost functional
%
% OUTPUT
% tau ....... Optimal step length
% fval ...... Value of cost functional, corresponding to optimal tau
% history ... History, containing step length, and value of cost functional
% exitflag .. Exitflag, to ensure the minimization performed correctly
% output .... Diagnostic, of the minimization process
%
% FORMAT
% [tau, fval, history, exitflag, output] = costfunmin(phi, grad, u0, A, dx, dt, u_b, t, tau, history)
% -------------------------------------------------------------------------
    function [tau, fval, history, exitflag, output] = costfunmin(phi, grad, u0, A, dx, dt, u_b, t, tau, history, x, E1)

        % Set options, to output the history
        options = optimset('OutputFcn', @myoutput);

        % The cost functional is line-minimized with respect to 's'
        if nargin == 12
            [tau, fval, exitflag, output] = fminbnd(@(s) costfun(phi + s*grad, u0, A, dx, dt, u_b, t, x, E1), 0.0, tau, options);
        else
            [tau, fval, exitflag, output] = fminbnd(@(s) costfun(phi + s*grad, u0, A, dx, dt, u_b, t), 0.0, tau, options);
        end

        % Nested function, to output the history
        function stop = myoutput(x,optimvalues,state)
            % Flag for exiting function
            stop = false;
            % If still iterating, update history vector
            if isequal(state,'iter')
                % Append the step-size and corresponding value to history
                history = [history; [x, optimvalues.fval] ];
            end
        end % End of nested function
    end % End of function costfunmin

% -------------------------------------------------------------------------
% FUNCTION: RetracConstr
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2023/11/30
%
% Use retraction operator to map back to manifold
%
% INPUT
% phi ....... Current Neumann Boundary condition
% u0 ........ Initial Condition Function
% A ......... Matrix for next time step
% dx ........ Step size in the spatial domain
% dt ........ Step size in the time domain
% t ......... Time domain
% x ......... Spatial domain
% E1 ........ Energy value defining manifold
%
% OUTPUT
% phiR ..... Neumann Boundry condition, 
%
% FORMAT
% phiR = RetracConstr(phi, u0, A, dx, dt, t, x, E1)
%
% -------------------------------------------------------------------------
    function phiR = RetracConstr(phi, u0, A, dx, dt, t, x, E1)

        % Solve Heat equation for current value of phi
        [~, utmp] = heatsolve(phi, u0, A, dx, dt);

        % Energy in system for current value of phi
        Ephi = 0.5*trapz(t, trapz(x, utmp.^2));

        % Retraction operator via simple normalization
        phiR = (sqrt(E1)/sqrt(Ephi))*phi;

    end % End of function RetracConstr

end
