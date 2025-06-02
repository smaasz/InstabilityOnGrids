% Eady setup, rectangular C-grid, baroclinic and symmetric 
% instability
  a=12500.0;               % -- side of quad
  
  theta=0;              % -- l=sqrt(3)Ksin(theta)/2, k=Kcos(theta)
  thetaw=0;                % -- mean flow direction
  
  f0=-0.0001;
  g=100000000;               % To effectively impose the rigid lid
  N=0.001; N2=N*N;
  Ri=100;                  % -- Richardson number
  M2=abs(N*f0)/sqrt(Ri);        % -- M^2, i. e. the horizontal stratification
  Nz=64;                   % -- the number of vertical layers
  H0=4000;                 % -- fluid depth
  dz=H0/Nz;
  U=-(M2*dz/f0)*cos(thetaw)*[-(1:Nz)+(Nz+1)/2];
  V=-(M2*dz/f0)*sin(thetaw)*[-(1:Nz)+(Nz+1)/2];
  bx=-M2*sin(thetaw);
  by=M2*cos(thetaw);
  % ==================
  % Auxiliary matrices:
  % ==================
  Az=zeros([Nz,Nz+1]);    % -- Vertical averaging of w
  for nz=1:Nz,
      Az(nz,nz)=1/2; Az(nz,nz+1)=1/2;
  end
  Wz=zeros([Nz+1,Nz]); 
  for nz=1:Nz,
      Wz(nz,nz:Nz)=1; 
  end
  W=dz*Az*Wz;            % -- Returns vertically av. w if applied to -gx*u-gy*v
  Bz=zeros([Nz,Nz]);
  for nz=1:Nz-1,
      Bz(nz,nz:nz+1)=0.5;
  end
  Pz=zeros([Nz,Nz]);
  for nz=1:Nz-1,
      Pz(nz,nz:nz+1)=[1,-1];
  end
  Pz(Nz,:)=1;
  P=dz*inv(Pz)*Bz;      % -- Returns pressure when multiplied with buoyancy
  
 % ================
 % Cycle
 % ================

 Kvec=0.025:0.05:pi;
 NN=length(Kvec);
 w=zeros([1 NN]);      % array for frequencies 
 n=0;
 for K=Kvec,   % for K=0 sometimes there are problems with
     k=K*cos(theta);      % ordering; k here == ka and l==lh
     l=K*sin(theta);
     n=n+1
     
     %%%% Horizontal gradients
     gx=(2*i/a)*sin(k/2);
     gy=(2*i/a)*sin(l/2);
     gx=i*k/a; gy=i*l/a;     %%%
     %%%% Horizontal averaging
     ax=cos(k/2); ay=cos(l/2);
     %ax=1; ay=1;
     %%%% Coefficients in the derivatives of b in b equation
     
     %%%% Matrices acting on u and v in the u and v equations
     Au=diag(ax*gx*U+ay*gy*V)-(M2/f0)*cos(thetaw)*ax*(-gx)*W;
     %Au=diag(i*k*U/a+i*l*V/a)-(M2/f0)*cos(thetaw)*ax*(-gx)*W;
     
     Av=-ax*ay*f0*diag(ones([1,Nz]))-(M2/f0)*cos(thetaw)*ax*(-gy)*W;
     % +gx*P b +g*gx eta
     Bv=diag(ax*gx*U+ay*gy*V)-(M2/f0)*sin(thetaw)*ay*(-gy)*W;
     %Bv=diag(i*k*U/a+i*l*V/a)-(M2/f0)*sin(thetaw)*ay*(-gy)*W;
     
     Bu=ax*ay*f0*diag(ones([1,Nz]))-(M2/f0)*sin(thetaw)*ay*(-gx)*W;
     % +gy*P b +g*gy eta
     Cb=diag(ax*gx*U+ay*gy*V);
     %Cb=diag(i*k*U/a+i*l*V/a);
     Cu=N2*W*(-gx)+(bx*ax)*diag(ones([1,Nz]));
     Cv=N2*W*(-gy)+(by*ay)*diag(ones([1,Nz]));
     %%% sea surface height matrix
     Eu=g*gx*ones([Nz,1]);
     Ev=g*gy*ones([Nz,1]);
     Du=dz*gx*ones([1,Nz]);
     Dv=dz*gy*ones([1,Nz]);
     %%% Construct system matrix:
     S=-[Au,Av,gx*P,Eu; Bu,Bv,gy*P,Ev;Cu, Cv,Cb,zeros([Nz,1]);Du,Dv,zeros([1,Nz+1])];
    %S=-[Au,Av,i*k*P/a,Eu; Bu,Bv,i*l*P/a,Ev;Cu, Cv,Cb,zeros([Nz,1]);Du,Dv,zeros([1,Nz+1])];
    
     %%% Find maximum unstable eigenvalue
     E=eig(S);
     w(n)=max(real(E));
 end

 df.f0 = f0;
 df.a = a/1000;
 df.type = "standard";
 df.grid = "quad-C";
 df.N = N;
 df.Ri = Ri;
 df.theta = theta;
 df.thetaU = thetaw;
 df.date = datestr(now, 'yy-mm-dd-HH:MM');
 df.Nz = Nz;
 df.ks = Kvec/a;
 df.vs = w;

 data = jsonencode(df);

 f = fopen('/Users/stmaas001/Projects/OceanFlows/SemiAnalyticInstabilityAnalysis/data/data.jsonl','a');
 fprintf(f, '%s\n', data);
 fclose(f);
 
plot(Kvec/a, w*N/abs(M2))
