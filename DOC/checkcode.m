function [alpha,dalpha,b0,db0]	= checkcode(X,y,lambda,gam,l)
%%FOLDUP
% [alpha,dalpha,b0,db0]	= checkcode(X,y,lambda,gam,l)
%
% code to compute the preimage of the l-factor WPLS regression coefficient 
% to fit the model y \approx X diag(lambda) diag(lambda) X' alpha
% also computes the jacobian of the preimage alpha with respect to lambda.
% and the offset and its gradient wrt lambda.
%
% input:
%  X           an m by n matrix, column centered.
%  y           an m by 1 vector, column centered.
%  lambda      an n vector of channel weightings.
%  gam         an m vector of response weightings.
%  l           the number of pls factors.
% output:
%  alpha       a m by 1 vector of the regression coefficients.
%  dalpha      a m by n matrix of the Jacobian of beta wrt lambda.
%  b0          the scalar offset (centers data)
%  db0         the gradient of same.
%
% nb: assumes m << n
%
% suggested use:
% m=100;n=80;
% checkcode(center(randn(m,n)),center(rand(m,1)),rand(n,1),ones(m,1),3);
% checkcode((randn(m,n)),(rand(m,1)),rand(n,1),ones(m,1),3);
%
% Author: Steven Pav 
% Created: 2006.04.11
% Copyright 2006
% $Id: checkcode.m 104 2006-05-16 20:55:21Z spav $
%%UNFOLD

%sanity checking%FOLDUP
[m,n]					= size(X);
lambda				= lambda(:);
gam           = gam(:);
y							= y(:);

if ((length(y) != m) || (length(lambda) != n) || (length(gam) != m))
	error('size mismatch.');
end
%UNFOLD

%allocate storage%FOLDUP
V	= zeros(m,l);dV	= zeros(m,n,l);
w = zeros(l,1);dw = zeros(l,n);
q = zeros(l+1,1);dq = zeros(l+1,n);

%used so much we compute once and store.
Xlam = X * diag(lambda);XLXg = Xlam * Xlam' * diag(gam);				%nm^2 hit.
%UNFOLD

T0 = Tk = ones(m,1);dTk = zeros(m,n);Vk = y;dVk = zeros(m,n);
t0 = sum(gam);
for k=0:l
	rk = y' * (gam .* Tk);  drk = dTk' * (gam .* y);
	tk = Tk' * (gam .* Tk); dtk = 2 * dTk' * (gam .* Tk);
	q(k+1) = rk/tk;         dq(k+1,:) = (drk - q(k+1) * dtk)' / tk;
	if (k < l)
		Vk = Vk - q(k+1) * Tk;
		V(:,k+1) = Vk;
		dV(:,:,k+1) = (dVk = dVk - q(k+1)*dTk - Tk * dq(k+1,:));
		Qk = XLXg * Vk;dQk = XLXg * dVk + 2 * Xlam * diag(X' * (gam .* Vk));
		Uk = Qk .- (gam' * Qk / t0);
		dUk = dQk - T0 * (gam' * dQk / t0);
		w(k+1) = Tk' * (gam .* Uk) / tk;
		dw(k+1,:)	= (dTk' * (gam .* Uk) + dUk' * (gam .* Tk) - w(k+1) .* dtk)' / tk;
		dTk = dUk - w(k+1) * dTk - Tk * dw(k+1,:);
		Tk = Uk - w(k+1) * Tk;
	end
end

%now compute M\q and its jacobian%FOLDUP
iMq	= zeros(l,1);diMq = zeros(l,n);
iMq(l) = q(l+1);diMq(l,:)	= dq(l+1,:);

for k=(l-1):-1:1
	iMq(k)		= q(k+1) - w(k+1) * iMq(k+1);
	diMq(k,:)	= dq(k+1,:) - (w(k+1) .* diMq(k+1,:) + iMq(k+1) .* dw(k+1,:));
end
%UNFOLD

%now compute VM\q and its jacobian%FOLDUP
alpha			= V * iMq;
dalpha		= V * diMq;			%+ more stuff:

%ack! no tensor product in octave/Matlab :(
for k=2:l		%dV(:,:,1) is all zeros?
	dalpha	+= iMq(k) .* dV(:,:,k);
endfor
%now the offset
b0				= gam' * (y - XLXg * alpha) ./ sum(gam);
db0				= - (XLXg * dalpha + 2 * Xlam * diag(X' * (gam .* alpha)))' * gam ./ sum(gam);
%UNFOLD

%check against the gold standard, 
Xlc = Xlam - ones(m,1) * gam' * Xlam / sum(gam);
yc = y .- gam' * y / sum(gam);

[W,T,P,qpls,bpls,b0pls,minRMSEP] = plsr(Xlc,yc,l);
%b0pls 
b0pls = gam' * (y - Xlam * bpls) ./ sum(gam);

%qpls - q(2:l+1)
%qpls
%q

beta			= Xlam' * (gam .* alpha);
printf('norm bdiff: %g, abs(b0diff): %g\n',norm(beta - bpls),abs(b0pls - b0));

plot(beta, bpls,'*r;;');
mesh(W - Xlam' * diag(gam) * V);

%T' * T

endfunction
