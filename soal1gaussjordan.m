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

fprintf('\n===== METODE GAUSS-JORDAN =====\n');
J = [A b];
for k = 1:N
    [~, m] = max(abs(J(k:N,k))); m = m + k - 1;
    if m ~= k, J([k m],:) = J([m k],:); end
    J(k,:) = J(k,:) / J(k,k);
    for i = [1:k-1, k+1:N]
        J(i,:) = J(i,:) - J(i,k) * J(k,:);
    end
end
fprintf('Matriks akhir [I | x]:\n');
disp(J)
xJ = J(:,end);
fprintf('Himpunan penyelesaian: {%s}\n', ...
    strjoin(arrayfun(@(i) sprintf('%s=%.6f', nama{i}, xJ(i)), 1:N, 'UniformOutput', false), ', '));
