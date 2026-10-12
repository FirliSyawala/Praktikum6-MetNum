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

fprintf('\n===== METODE JACOBI =====\n');
fprintf('Iter |     p     |     q     |     r     |     s     |     t     | Selisih\n');
xJc = zeros(N,1);
iterJ = 0;
for it = 1:maxit
    xb = zeros(N,1);
    for i = 1:N
        idx = [1:i-1, i+1:N];
        xb(i) = (b(i) - A(i,idx) * xJc(idx)) / A(i,i);
    end
    sel = max(abs(xb - xJc));
    fprintf('%4d | %9.5f | %9.5f | %9.5f | %9.5f | %9.5f | %.2e\n', it, xb, sel);
    xJc = xb;
    iterJ = it;
    if sel < tol, break; end
end
if iterJ < maxit
    fprintf('Konvergen pada iterasi ke-%d\n', iterJ);
else
    fprintf('Mencapai batas maks iterasi (%d)\n', maxit);
end
