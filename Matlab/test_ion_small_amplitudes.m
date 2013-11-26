% Test small amplitude termination criterion
p = 5;
p_fit = 6;
n = 30;

a = (1:p)';
omega = 2i*pi*linspace(0,1,p+1).';
omega = omega(1:p);

omega0 = 2i*pi*linspace(0,1,p_fit).';
omega0 = omega0(1:p_fit) + 1e-2*(randn(p_fit,1)+1i*randn(p_fit,1));
omega0 = min(real(omega0),0) + 1i*imag(omega0);

y = mkV(omega,n)*a;

% Noisy measurements
g = (randn(n,1)+1i*randn(n,1))/sqrt(2);
yt = y + 1e-4*g;


opts = optimset('tolx',1e-20);

[omega_fit,a_fit] = expfit_varpro(yt,omega0,opts);

[abs_a_fit,I] = sort(abs(a_fit),'descend');
omega_fit = omega_fit(I);

fprintf('Im(Omega) \t Amplitude\n')
for j = 1:p_fit
	fprintf('% 5e\t %5e\n',imag(omega_fit(j)),abs_a_fit(j));
end
