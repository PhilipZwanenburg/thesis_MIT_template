format shorte
clear,clc,
addpath external
addpath external/Hesthaven

% Demkowicz2011:     [alpha,beta,gamma] = [1 0 0]
% Step function v_f:                    = Still possible?

h = 1/2^2;
% [alpha,beta,gamma] >= 0; alpha+beta+gamma > 0.
%alpha = 2; beta  = 0; gamma = 0;
alpha = 0; beta  = 0; gamma = 2;
alpha = 2; beta  = 3; gamma = 2;

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
v_f_coef = K\b;

n_coef = size(v_f_coef,1)/2;
v_f_coef_l = v_f_coef((1:n_coef)+0)
v_f_coef_r = v_f_coef((1:n_coef)+n_coef)

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
v_fl = phi_p*v_f_coef_l;
v_fr = phi_p*v_f_coef_r;

clf; hold on;
subplot(4,1,1); hold on; grid on;
plot(r_p-1,v_fl,'-bo');
plot(r_p+1,v_fr,'-bo');

subplot(4,1,2); hold on; grid on;
plot(r_p,v_s(:,1),'-bo');
if (ps > 0)
subplot(4,1,3); hold on; grid on;
plot(r_p,v_s(:,2),'-bo');
end
if (ps > 1)
subplot(4,1,4); hold on; grid on;
plot(r_p,v_s(:,3),'-bo');
end







% Checking

% v_s
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

subplot(4,1,ind+2);
plot(r_p,v_ex(r_p),'-rs');

% v_f



A = [
alpha+2*beta+gamma alpha-gamma -beta beta;
alpha-gamma 4/h+alpha+2*beta+gamma -beta beta;
-beta -beta alpha+2*beta+gamma alpha-gamma;
beta beta alpha-gamma 4/h+alpha+2*beta+gamma;
];

b = [1 1 -1 1]';

tmp = A\b
a0_l = tmp(1);
a1_l = tmp(2);
a0_r = tmp(3);
a1_r = tmp(4);

vl_ex = @(r) [a0_l+a1_l*r];
vr_ex = @(r) [a0_r+a1_r*r];

subplot(4,1,1);
plot(r_p-1,vl_ex(r_p),'rs');
plot(r_p+1,vr_ex(r_p),'rs');

a_lr = [
((alpha*gamma-gamma^2)*h+alpha);
0;
-((alpha^2-alpha*gamma)*h+2*alpha)/2;
(alpha^2-alpha*gamma)*h/2;
]/(alpha^2-2*alpha*gamma-(alpha*gamma^2-gamma^3)*h)
