[nc,nd]=size(cD1);
[CC,DD]=gKernel(nc,nd,cD1);
A=pca_energy((circsim1+CC)/2,0.7);
B=pca_energy((dissim1+DD)/2,0.9);
M_Omega=cD1;
Omega_linear=find(cD1==1);
lambda=1;
max_iter=1000;
r_a=size(A,2);   % matrix A column number
r_b=size(B,2);
L=1;	% parameter s
gamma=2;
Z0=zeros(r_a,r_b);
Z=Z0;
alpha0=1;
alpha=1;
i=0;
convergence=zeros(max_iter,1);

M_Omega_linear=full(M_Omega(Omega_linear))';     
[n,m]=size(M_Omega);
[row,column]=index2spa(Omega_linear,n);

svdt3=A'*M_Omega*B;
AZ0BOmega=xumm(A*Z0,B',row,column);
AZBOmega=AZ0BOmega;

while i<max_iter	% APG
    i=i+1;
    Y=Z+alpha*(1/alpha0-1)*(Z-Z0);
    Z0=Z;
    AYBOmega=(1+alpha*(1/alpha0-1))*AZBOmega-(alpha*(1/alpha0-1))*AZ0BOmega;   
    AZ0BOmega=AZBOmega;
	
    svdt2=A'*(sparse(row,column,AYBOmega,n,m))*B;	%
	
    Z=sidesvd2Threshold(Y,svdt2,svdt3,L,lambda);     
	
    AZBOmega=xumm(A*Z,B',row,column);   
	
    qlpl1=norm(AYBOmega-M_Omega_linear,'fro')^2/2;
    qlpl2=svdt2-svdt3;
    DiffL2=norm(AZBOmega-M_Omega_linear,'fro')^2/2;	% approximation error
	
    while DiffL2>Qlpl(Z,Y,L,qlpl1,qlpl2)	
        L=L*gamma;
        Z=sidesvd2Threshold(Y,svdt2,svdt3,L,lambda);	% SVT to sovle approximation of Z
        AZBOmega=xumm(A*Z,B',row,column);               
        DiffL2=norm(AZBOmega-M_Omega_linear,'fro')^2/2;
    end
	
    alpha0=alpha;
    alpha=(sqrt(alpha^4+4*alpha^2)-alpha^2)/2;    
    convergence(i,1)=sideobjectCalc(Z,lambda,DiffL2); 
      if i>1
          if abs(convergence(i,1)-convergence(i-1,1))<(1e-5)*convergence(i,1)
                break;
          end
      end
end
M_recover=A*Z*B';
