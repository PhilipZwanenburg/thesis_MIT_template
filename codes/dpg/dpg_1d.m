format short
clear,clc,
addpath external
addpath external/Hesthaven

% Demkowicz2011:     [alpha,beta,gamma] = [1 0 0]
% Step function v_f:                    = Still possible?

% [alpha,beta,gamma] >= 0; sum(alpha:gamma) > 0.
alpha = 2;
beta  = 3;
gamma = 4;
h = 1/2^4;

basis = 'm'
ps = 2;
dp = 1;
p = ps+dp;

p_c = p+0;
[r,w,~] = cub(p_c,1,'GL');
n = p_c;


if (basis == 'n') % nodal (GLL)
	[tmp,~,~] = cub(p,1,'GLL');
	T = inv(phi(p,tmp));
%	[tmp,~,~] = cub(ps,1,'GLL');
%	Ts = inv(phi(ps,tmp));
elseif (basis == 'm') % modal
	T = eye(p+1);
end
Ts = eye(ps+1);

W = diag(w);
psi_m1 = phi(p,-1)*T;
psi_p1 = phi(p,+1)*T;
grad_psi_rv = grad_phi(p,r)*T;
psi_H1 = 2/h*grad_psi_rv'*W*grad_psi_rv;
n = size(psi_m1',1);
zero = zeros(n,n);

% v_f
b = [psi_p1 -psi_m1]';

K = [
 psi_H1 zero;
 zero psi_H1;
] + alpha*[
 psi_p1'*psi_p1 zero;
 zero zero;
] + beta*[
 psi_p1'*[psi_p1 -psi_m1];
-psi_m1'*[psi_p1 -psi_m1];
] + gamma*[
 zero zero;
 zero psi_m1'*psi_m1;
] + (beta+gamma)*[
 psi_m1'*psi_m1 zero;
 zero zero;
] + (alpha+beta)*[
 zero zero;
 zero psi_p1'*psi_p1;
];

cond(K)
v_f_coef = K\b % OK! Exact step.

% v_s
psi_rv = phi(p,r)*T;
phi_rv = phi(ps,r)*Ts;

b = [-grad_psi_rv'*W*phi_rv];

K = [
 psi_H1
] + (alpha+beta)*[
psi_p1'*psi_p1
] + (beta+gamma)*[
psi_m1'*psi_m1
];

cond(K)
v_s_coef = K\b

% Note: Integrals of Legendre polynomials of degree greater than 1 are equal to zero on __both__ end points. See the
% link below for:
% \int P_n(r) dr = 1/(2n+1)*(P_{n+1}(r)-P_{n-1}(r))
% link: http://functions.wolfram.com/Polynomials/LegendreP/21/ShowAll.html
%
% This implies that only the v_s component corresponding to the p0 basis function adds a term for the test norm!





% plotting
p_p = 20;

[r_p,~,~] = cub(p_p,1,'GLL');

phi_p = phi(p,r_p)*T;
v_s = phi_p*v_s_coef;

clf; hold on;
subplot(3,1,1); hold on; grid on;
plot(r_p,v_s(:,1),'-bo');
if (ps > 0)
subplot(3,1,2); hold on; grid on;
plot(r_p,v_s(:,2),'-bo');
end
if (ps > 1)
subplot(3,1,3); hold on; grid on;
plot(r_p,v_s(:,3),'-bo');
end







% Check
ind = 0;
grad_psi_rv'*W*(2/h*grad_psi_rv*v_s_coef(:,ind+1)+phi_rv(:,ind+1)) ...
+ (alpha+beta)*psi_p1'*psi_p1*v_s_coef(:,ind+1) ...
+ (beta+gamma)*psi_m1'*psi_m1*v_s_coef(:,ind+1)

if (ind == 0) 
	a1 = -2/sqrt(2)*(4/h+alpha+2*beta+gamma-(alpha-gamma)^2/(alpha+2*beta+gamma))^-1;
	a0 = (gamma-alpha)/(alpha+2*beta+gamma)*a1;

	v_ex = @(r) [a0+a1*r];
	gv_ex = @(r) [a1+0*r];

	grad_psi_rv'*W*(2/h*gv_ex(r)+phi_rv(:,ind+1)) ...
	+ (alpha+beta)*psi_p1'*v_ex(1) ...
	+ (beta+gamma)*psi_m1'*v_ex(-1)
elseif (ind == 1)
	v_ex = @(r) [sqrt(3/2)/4*h*(1-r.^2)];
	gv_ex = @(r) [sqrt(3/2)/4*h*(-2*r)];
	grad_psi_rv'*W*(2/h*gv_ex(r)+phi_rv(:,ind+1))
	psi_p1'*v_ex(1)
end

subplot(3,1,ind+1);
plot(r_p,v_ex(r_p),'-rs');
