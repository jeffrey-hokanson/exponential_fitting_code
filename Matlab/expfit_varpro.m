function [omega,a] = expfit_varpro(y,omega,options)
% function [omega,a] = expfit_varpro(y,omega,options)
% 
% Solve the exponential fitting problem using data y
% from initial estimate omega passing options to lsqnonlin
% to solve the nonlinear least squares problem.
%
% Copyright 2013 Jeffrey Hokanson
% Distributed under the GPLv2 License: http://www.gnu.org/licenses/gpl.html

n = length(y);
p = length(omega);

% By default, lsqnonlin does not use Jacobians.  Turn these
% on and also allow user setting of tolerances.
default_options = optimset('display','off','Jacobian','on','TolFun',1e-15,'TolX',1e-15);
if nargin == 3
	options = optimset(default_options,options);
else
	options = default_options;
end

% As lsqnonlin only allows real parameters, split the variables
om = [real(omega);imag(omega)];

% Construct upper bounds for the real part of omega
om_ub = [100/n*ones(p,1);inf*ones(p,1)];

% Optimize the Variable Projection Jacobian
om = lsqnonlin(@(x) mkrJ(x),om,[],om_ub,options);

% Reconstruct omega
omega = om(1:p)+1i*om(p+1:end);

% Determine linear parameters (optional)
if nargout >1
	V = mkV(omega,n);
	a = V\y;
end


function [r_real,J_real] = mkrJ(x)
% Compute the residual and Jacobian using Variable Projection

% reform complex omega
omega_local = x(1:p)+1i*x(p+1:end);

% Construct complex residual 
V = mkV(omega_local,n);
[Q,R] = qr(V,0);
b = Q'*y;
r = y-Q*b;

% Construct real residual
r_real = [real(r); imag(r)];


if nargout > 1
	% Construct complex Jacobian
	Vp = mkVp(omega_local,n);
	a = R\b;
	L = (Vp - Q*Q'*Vp)*diag(a);

	% The computation of K can be removed to speed execution
	%K = Q*( (R')\diag(Vp'*r));
	%J = -L-K;
	J = -L;
	
	% Construct real version
	J_real = [real(J), -imag(J);
	     imag(J), real(J)];
end
end % mkrJ


end % expfit_varpro

