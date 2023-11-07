clc;
A=randn(4)+1i*randn(4);
A=A*A';

[V,D]=eig(A,'vector');
sort(D)

D_=real(sort(eig(A+3*V(:,1)*V(:,1)')))
