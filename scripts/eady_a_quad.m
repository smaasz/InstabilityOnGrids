% Eady setup, triangular A-grid, baroclinic and symmetric 
% instability
  a=12500.0;                % -- side of triangle
  h=a*sqrt(3)/2;           % -- the height of triangle 
  
  theta=0;                 % -- l=sqrt(3)Ksin(theta)/2, k=Kcos(theta)
  thetaw=0;                % -- mean flow direction
  
  f0=-0.0001;
  g=100000000;             % To effectively impose the rigid lid
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

 Kvec=0.025:0.05:pi*2/sqrt(3);
 NN=length(Kvec);
 w=zeros([1 NN]);      % array for frequencies 
 n=0;
 for K=Kvec,   % for K=0 sometimes there are problems with
     k=K*cos(theta);      % ordering; k here == ka and l==lh
     l=K*sqrt(3)*sin(theta)/2;
     n=n+1
     
     %%%% Horizontal gradients
     gx=(2*i/3/a)*(sin(k)+0.5*sin(k/2+l)-0.5*sin(-k/2+l));
     gy=(i/2/h)*(sin(k/2+l)+sin(-k/2+l));
     %gx=i*k/a; gy=i*l/a;     %%%
     
     %%%% Coefficients in the derivatives of b in b equation
     axu=(2/3)*(0.5*(1+cos(k))+0.125*(1+cos(k/2+l))+0.125*(1+cos(-k/2+l)));
     axv=(sqrt(3)/12)*(cos(k/2+l)-cos(-k/2+l));
     ayu=axv;
     ayv=0.25*(2+cos(k/2+l)+cos(-k/2+l));
     %axu=1; ayv=1; axv=0.;ayu=0.;    %%%
     
     %%%% Matrices acting on u and v in the u and v equations
     Au=diag(gx*U+gy*V)-(M2/f0)*cos(thetaw)*(-gx)*W;
     Av=-f0*diag(ones([1,Nz]))-(M2/f0)*cos(thetaw)*(-gy)*W;
     % +gx*P b +g*gx eta
     Bv=diag(gx*U+gy*V)-(M2/f0)*sin(thetaw)*(-gy)*W;
     Bu=f0*diag(ones([1,Nz]))-(M2/f0)*sin(thetaw)*(-gx)*W;
     % +gy*P b +g*gy eta
     Cb=diag(gx*U+gy*V);
     Cu=N2*W*(-gx)+(bx*axu+by*ayu)*diag(ones([1,Nz]));
     Cv=N2*W*(-gy)+(bx*axv+by*ayv)*diag(ones([1,Nz]));
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
 df.a = a * 1e-3;
 df.type = "standard";
 df.grid = "quad-A";
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
