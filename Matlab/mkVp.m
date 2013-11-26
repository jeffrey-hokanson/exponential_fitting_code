function V = mkVp(omega,n)
% function Vp = mkVp(omega,n)
%
% Forms the matrix \dot{V}(\omega) where each column is the derivative
% of the corresponding row in V(\omega) (see mkV.m)
%
% [Vp]_{j,k} = j e^{j \omega_k}
% 
% omega  - p parameters, each corresponding to one column
% n      - number of rows to create
% Vp     - n x p output matrix
% Copyright 2013 Jeffrey Hokanson
% Distributed under the GPLv2 License: http://www.gnu.org/licenses/gpl.html


p = length(omega);
omega = reshape(omega,p,1);

% actual formula should read:
% V(2:end,j) = (0:n-1)'.*exp(omega(j)*(1:n-1)'.*sc);
% To avoid possible numerical instability, we pull the scaling 
% into the exponential.
sc = log( (1:n-1)');

V = zeros(n,p);
V(2:end,:) = exp((1:n-1)'*omega.' + sc*ones(1,p));

% TODO: improve timing using BSXFUN
