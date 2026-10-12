clc; clear;
A = [ 5 -1 -1  0  0;
     -1  5 -1 -1  0;
     -1 -1  4 -1 -1;
      0  0  1  4 -2;
      0  1 -1  1  4];
b = [-1; 2; 6; 2; -1];
N = 5; maxit = 20; tol = 0.001;
judul = {'Gauss','Gauss-Jordan','Dekomposisi LU','Jacobi','Gauss-Seidel'};

% --- Gauss ---
M = [A b];
for k = 1:N-1, M(k+1:N,:) = M(k+1:N,:) - M(k+1:N,k)/M(k,k)*M(k,:); end
x1 = triu(M(:,1:N)) \ M(:,end);

% --- Gauss-Jordan ---
J = [A b];
for k = 1:N
    J(k,:) = J(k,:)/J(k,k); idx = [1:k-1,k+1:N];
    J(idx,:) = J(idx,:) - J(idx,k)*J(k,:);
end
x2 = J(:,end);

% --- Dekomposisi LU ---
L = eye(N); U = A;
for k = 1:N-1
    L(k+1:N,k) = U(k+1:N,k)/U(k,k);
    U(k+1:N,:) = U(k+1:N,:) - L(k+1:N,k)*U(k,:);
end
x3 = U \ (L \ b);

% --- Jacobi ---
x4 = zeros(N,1); d = diag(A); R = A - diag(d);
for it = 1:maxit
    xb = (b - R*x4)./d; s = max(abs(xb-x4)); x4 = xb;
    if s < tol, break, end
end

% --- Gauss-Seidel ---
x5 = zeros(N,1);
for it = 1:maxit
    xb = tril(A) \ (b - triu(A,1)*x5); s = max(abs(xb-x5)); x5 = xb;
    if s < tol, break, end
end

% ===== SOAL 2: GALAT =====
xe = inv(A) * b;                    % nilai eksak (analitik)
hasil = [x1 x2 x3 x4 x5];
galat = abs(xe - hasil);            % tiap kolom = satu metode

fprintf('\nSOAL 2: Galat = |eksak - hampiran|\n');
fprintf('%-16s %10s %10s %10s %10s %10s\n','Metode','p','q','r','s','t');
for k = 1:5
    fprintf('%-16s %10.2e %10.2e %10.2e %10.2e %10.2e\n', judul{k}, galat(:,k));
end
