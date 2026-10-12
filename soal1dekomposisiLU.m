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

fprintf('\n===== DEKOMPOSISI LU =====\n');
L = eye(N); U = zeros(N);
for i = 1:N
    for j = i:N
        U(i,j) = A(i,j) - L(i,1:i-1) * U(1:i-1,j);
    end
    for j = i+1:N
        L(j,i) = (A(j,i) - L(j,1:i-1) * U(1:i-1,i)) / U(i,i);
    end
end
fprintf('Matriks L:\n'); disp(L)
fprintf('Matriks U:\n'); disp(U)
y = zeros(N,1);
for i = 1:N
    y(i) = b(i) - L(i,1:i-1) * y(1:i-1);
end
fprintf('Vektor y (Ly = b):\n'); disp(y)
xL = zeros(N,1);
for i = N:-1:1
    xL(i) = (y(i) - U(i,i+1:N) * xL(i+1:N)) / U(i,i);
end
fprintf('Solusi x (Ux = y):\n'); disp(xL)
