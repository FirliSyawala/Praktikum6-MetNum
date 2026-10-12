clc; clear;
A = [ 5 -1 -1  0  0;
     -1  5 -1 -1  0;
     -1 -1  4 -1 -1;
      0  0  1  4 -2;
      0  1 -1  1  4];
b = [-1; 2; 6; 2; -1];
N = 5;
maxit = 20;
tol = 0.001;
nama = {'p','q','r','s','t'};

fprintf('\n===== METODE GAUSS-SEIDEL =====\n');
fprintf('Iter |     p     |     q     |     r     |     s     |     t     | Selisih\n');
xS = zeros(N,1);
iterGS = 0;
for it = 1:maxit
    lama = xS;
    for i = 1:N
        idx = [1:i-1, i+1:N];
        xS(i) = (b(i) - A(i,idx) * xS(idx)) / A(i,i);
    end
    sel = max(abs(xS - lama));
    fprintf('%4d | %9.5f | %9.5f | %9.5f | %9.5f | %9.5f | %.2e\n', it, xS, sel);
    iterGS = it;
    if sel < tol, break; end
end
if iterGS < maxit
    fprintf('Konvergen pada iterasi ke-%d\n', iterGS);
else
    fprintf('Mencapai batas maks iterasi (%d)\n', maxit);
end
