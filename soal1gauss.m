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

fprintf('\n===== METODE GAUSS =====\n');
M = [A b];
for k = 1:N-1
    [~, m] = max(abs(M(k:N,k))); m = m + k - 1;
    if m ~= k, M([k m],:) = M([m k],:); end
    for i = k+1:N
        f = M(i,k) / M(k,k);
        M(i,:) = M(i,:) - f * M(k,:);
    end
end
fprintf('Matriks augmented setelah eliminasi (segitiga atas):\n');
disp(M)
xG = zeros(N,1);
for i = N:-1:1
    xG(i) = (M(i,end) - M(i,i+1:N) * xG(i+1:N)) / M(i,i);
end
fprintf('Hasil substitusi mundur:\n');
for i = N:-1:1
    fprintf('  %s = %.6f\n', nama{i}, xG(i));
end
