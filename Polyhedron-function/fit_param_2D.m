function feval = fit_param_2D(x_sol,vargin)
% Inputs
% x_sol: solution structure
% vargin: parameters/data

% Output
% feval: error

% Assign parameters
theta = vargin{1};
hr = vargin{2};
periodic_N = vargin{3};

% Guess solution
shift_val = x_sol(1);
pre_fact2 = x_sol(2);
pre_fact4 = x_sol(3);
pre_fact6 = x_sol(4);


% Calculate guess values
theta_mod = mod(theta+periodic_N/2-shift_val/2,periodic_N)-periodic_N/2;
hr_guess = pre_fact2*(theta_mod).^2 + pre_fact4*(theta_mod).^4 + pre_fact6*(theta_mod).^6 + 1;

% Calculate error
feval = sum((hr - hr_guess).^2);


end