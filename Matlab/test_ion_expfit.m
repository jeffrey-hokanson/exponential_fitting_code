% Test the exponential fitting routine using an example from Ion
warning off

% Parameters from a test by Ion-Victor
omega = 1i*pi*linspace(0,1,21).';
omega = omega(1:20);

p = length(omega);
a = (1:p)';

n = 60;

% Construct synthetic measurements
y = mkV(omega,n)*a;

% Perturb initial estimates
omega0 = omega + log(0.7) +1e-2*(randn(p,1)+1i*randn(p,1)); % This yields a initial estimates in the right half plane

fprintf('Allowing initial estimates in the right half plane.\n');

% Compute condition number of example
fprintf('Condition number of true parameters : %g\n',cond(mkV(omega,n)));
fprintf('Condition number of initial estimate: %g\n',cond(mkV(omega0,n)));

% Fit
[omega_fit,a_fit] = expfit_varpro(y,omega0);

% Compute residual
r = y - mkV(omega_fit,n)*a_fit;
fprintf('Norm of residual: %g\n',norm(r))


fprintf('\n');
fprintf('Restricting initial estimates in the left half plane.\n');

omega0 = min(real(omega),0) + 1i*imag(omega);	% Restrict initial estimates to left half plane

% Compute condition number of example
fprintf('Condition number of true parameters : %g\n',cond(mkV(omega,n)));
fprintf('Condition number of initial estimate: %g\n',cond(mkV(omega0,n)));

% Fit
[omega_fit,a_fit] = expfit_varpro(y,omega0);



% Compute residual
r = y - mkV(omega_fit,n)*a_fit;
fprintf('Norm of residual: %g\n',norm(r))

