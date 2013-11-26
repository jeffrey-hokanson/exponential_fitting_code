function V = mkV(omega,n)
% function V = mkV(omega,n)
%
% Forms the matrix [V(omega)]_{j,k} = e^{j omega_k}
% 
% omega  - p parameters, each corresponding to one column
% n      - number of rows to create
% V      - n x p output matrix
%
% Copyright 2013 Jeffrey Hokanson
% Distributed under the GPLv2 License: http://www.gnu.org/licenses/gpl.html

p = length(omega);
omega = reshape(omega,p,1);
V = exp((0:n-1)'*omega.');

% TODO: improve timing using BSXFUN
