function soal3runningtime()
A=[5 -1 -1 0 0;-1 5 -1 -1 0;-1 -1 4 -1 -1;0 0 1 4 -2;0 1 -1 1 4]; b=[-1;2;6;2;-1];
N=5; mx=20; tol=0.001; u=500; w=zeros(1,5);
m={'Gauss','Gauss Jordan','Dekomposisi LU','Jacobi','Gauss Seidel'};

tic; for r=1:u  % Gauss
 M=[A b]; for k=1:N-1, M(k+1:N,:)=M(k+1:N,:)-M(k+1:N,k)/M(k,k)*M(k,:); end
 x=zeros(N,1); for i=N:-1:1, x(i)=(M(i,N+1)-M(i,i+1:N)*x(i+1:N))/M(i,i); end
end; w(1)=toc/u;

tic; for r=1:u  % Gauss Jordan
 M=[A b]; for k=1:N, M(k,:)=M(k,:)/M(k,k); i=[1:k-1,k+1:N]; M(i,:)=M(i,:)-M(i,k)*M(k,:); end
end; w(2)=toc/u;

tic; for r=1:u  % LU
 L=eye(N); U=A; for k=1:N-1, L(k+1:N,k)=U(k+1:N,k)/U(k,k); U(k+1:N,:)=U(k+1:N,:)-L(k+1:N,k)*U(k,:); end
 y=zeros(N,1); for i=1:N, y(i)=b(i)-L(i,1:i-1)*y(1:i-1); end
 x=zeros(N,1); for i=N:-1:1, x(i)=(y(i)-U(i,i+1:N)*x(i+1:N))/U(i,i); end
end; w(3)=toc/u;

tic; for r=1:u  % Jacobi
 x=zeros(N,1); for it=1:mx, xb=x; for i=1:N, j=[1:i-1,i+1:N]; xb(i)=(b(i)-A(i,j)*x(j))/A(i,i); end
 s=max(abs(xb-x)); x=xb; if s<tol, break, end, end
end; w(4)=toc/u;

tic; for r=1:u  % Gauss Seidel
 x=zeros(N,1); for it=1:mx, o=x; for i=1:N, j=[1:i-1,i+1:N]; x(i)=(b(i)-A(i,j)*x(j))/A(i,i); end
 if max(abs(x-o))<tol, break, end, end
end; w(5)=toc/u;

fprintf('Komparasi Running Time\n\n%-16s | Running Time (s)\n%s\n','Metode',repmat('-',1,36));
for i=1:5, fprintf('%-16s | %.8f\n',m{i},w(i)); end
[a,p]=min(w); [c,q]=max(w);
fprintf('\nTercepat : %s (%.8f detik)\nTerlambat: %s (%.8f detik)\n',m{p},a,m{q},c);
end
