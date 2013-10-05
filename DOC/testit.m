
m						= 40;
n						= 50;
l						= min(m,n);
k						= 18;
errlvl			= 0.1;
errlvl			= 0;

wch					= 5;

Ts					= randn(m,k);
Xs					= Ts * randn(k,n);
ys					= Ts * randn(k,1);
ys				 += errlvl * randn(size(ys));

W = T = P = q = [];

Xk					= center(Xs);
yk					= center(ys);

for j=1:l%FOLDUP
	W(:,j)			= Xk' * yk;
	T(:,j)			= Xk * W(:,j);
	tweight			= T(:,j)' * T(:,j);
	P(:,j)			= Xk' * T(:,j) ./ tweight;
	q(j,1)			= yk' * T(:,j) ./ tweight;

	if (j > 1)
		printf("||X_%d (P_%d - W_%d P_%d' P_%d)|| = %g\n",j,j-1,j,j,j-1,norm( \
		 Xk * (P(:,j-1) - W(:,j) * P(:,j)' * P(:,j-1)) ));
	endif
	%old way
%		Xk				 -= (T(:,j) * T(:,j)' ./ tweight) * Xk;
	%faster?
	Xk				 -= T(:,j) * P(:,j)';

	printf("max(max(X_%d)) = %g\n",j+1,max(max(Xk)));
	printf("T_%d'y = %g\n",j,yk' * T(:,j));
	printf("q_%d = %g\n",j,q(j,1));

	if (j >= wch)
		printf("||X_%d P_%d|| = %g\n",j+1,wch,norm( Xk * P(:,wch) ));
	endif

endfor%UNFOLD

