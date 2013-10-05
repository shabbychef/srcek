function testpls

m = 100;
n = 200;
l = 6;

X	= randn(m,n);
y = rand(m,1);

%works under these conditions:
%X = center(X);
%y = center(y);
%aber, without is nfg?

X2 = X1 = X;
XX = X*X';

T1 = T2 = ones(m,1);
V2 = y;

for k=0:l
	printf('Tdiff (%d): %g\n',k,norm(T1 - T2));
	t1 = T1' * T1;t2 = T2' * T2;
	printf('tdiff (%d): %g\n',k,abs(t1-t2));
	q1 = y' * T1 / t1;q2 = y' * T2 / t2;
	printf('qdiff (%d): %g\n',k,abs(q1-q2));

	P1 = X1' * T1 / t1;P2 = X' * T2 / t2;
	printf('Pdiff (%d): %g\n',k,norm(P1 - P2));
	
	if (k < l)
		X1 = X1 - T1 * P1';				%X2 = X2 - T2 * T2' * X2 / t2;

		W1 = X1' * y;V2 = V2 - q2 * T2;W2 = X' * V2;
		printf('Wdiff (%d): %g\n',k+1,norm(W1 - W2));

		U2 = center(XX * V2);U1 = center(X * W2);
		printf('Udiff (%d): %g\n',k+1,norm(U1 - U2));

		w1 = P1' * W1;w2 = T2' * U2 / t2;
		printf('wdiff (%d): %g\n',k+1,abs(w1-w2));

		if (k == 0)
			printf('w2 = %g\n',w2);
			T1 = X1 * W1;T2 = U2;
		else
			T1 = X1 * W1;T2 = U2 - w2 * T2;
		endif
	endif
endfor



endfunction


