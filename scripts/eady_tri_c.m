% Eady setup, triangular mimetic C-grid, baroclinic and symmetric 
% instability
  a=12500.0;               % -- side of triangle
  h=a*sqrt(3)/2;           % -- the height of triangle 
  
  theta=0;              % -- l=sqrt(3)Ksin(theta)/2, k=Kcos(theta)
  thetaw=0;                % -- mean flow direction
  
  f0=-0.0001;
  g=100000000;               % To effectively impose the rigid lid
  %g=10;
  N=0.001; N2=N*N;
  Ri=100;                  % -- Richardson number
  M2=abs(N*f0)/sqrt(Ri);        % -- M^2, i. e. the horizontal stratification
  Nz=64;                   % -- the number of vertical layers
  H0=4000;                 % -- fluid depth
  dz=H0/Nz;
  U=-(M2*dz/f0)*cos(thetaw)*[-(1:Nz)+(Nz+1)/2];
  V=-(M2*dz/f0)*sin(thetaw)*[-(1:Nz)+(Nz+1)/2];
  %U=-(M2*dz/f0)*cos(thetaw)*[-(1:Nz)+1/2];   % non-symmetric
  %V=-(M2*dz/f0)*sin(thetaw)*[-(1:Nz)+1/2];   % non-symmetric
  
  U0=-(M2*dz/f0)*cos(thetaw); V0=-(M2*dz/f0)*sin(thetaw);
  UV=[-(1:Nz)+(Nz+1)/2];   % Profile 
  %UV=[-(1:Nz)+1/2];       % non-symmetric
  
  bx=-M2*sin(thetaw);
  by=M2*cos(thetaw);
  bua=bx*a/4+by*h/6; bub=-bx*a/4+by*h/6; buc=-by*h/3;  % Buoyancy increments
                                                       % on u-triangle
  Ua=U*sqrt(3)/2+V/2; Ub=-U*sqrt(3)/2+V/2; Uc=-V;      % Three normal velocities
  M2a=-(M2/f0)*(cos(thetaw)*sqrt(3)/2+sin(thetaw)/2);  % Their contrib. in M2-term
  M2b=-(M2/f0)*(-cos(thetaw)*sqrt(3)/2+sin(thetaw)/2); % in momentum
  M2c=-(M2/f0)*(-sin(thetaw));
     
  
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
 for K=Kvec,   
     k=K*cos(theta);   % k here == ka and l==lh
     l=K*sqrt(3)*sin(theta)/2;
     n=n+1
     %%%% Phase shifts
     aa=exp(i*k/4+i*l/6); bb=exp(-i*k/4+i*l/6); cc=exp(-i*l/3);
     %%%% Horizontal gradients
     ga=(3/2/h)*[-conj(aa),aa]; gb=(3/2/h)*[-conj(bb),bb]; gc=(3/2/h)*[-conj(cc),cc];
     %%%% Perot reconstruction
     Pe=(2/h)*[aa*a/4, -bb*a/4, 0; aa*h/6, bb*h/6, -cc*h/3];
     PP=[Pe;conj(Pe)];
     %%%% Transpose Perot
     PT=(3/2/h)*[conj(aa)*a/4, conj(aa)*h/6; -conj(bb)*a/4,conj(bb)*h/6; ...
                 0, -conj(cc)*h/3];
     PT=[PT,conj(PT)]; Me=PT*PP;    % This is a P^TP matrix, a kind
                                    % of elementary mass matrix.
     %%%% Fourier symbol of horizontal divergence
     Dh=(2/h)*[aa,bb,cc;-conj(aa),-conj(bb),-conj(cc)];
     %%%% Mean divergence at a,b,c edges:
     Ae=(1/2)*[conj(aa),aa; conj(bb),bb; conj(cc),cc]; 
     De=Ae*Dh*Me;   %%% This is 3 by 3 matrix that acts on velocities to produce 
                    %%% the averaged divergence at velocity locations.  
     De=diag([M2a,M2b,M2c])*De; 
                    %%% Weighted with M2 term amplitudes (still 3 by 3)
     % Vector-invariant advection 
     Om=i*(4/3/a)*[-sin(-k/4+l/2), sin(k/4+l/2), -sin(k/2)];
     % Multiply with (U,V) rotated to (-V,U) and project on normal direction     
     Ut=[-cos(k/4-l/2)*(-U0/2+V0*sqrt(3)/2); cos(k/4+l/2)*(U0/2+V0*sqrt(3)/2); ...
         -cos(k/2)*(U0)];   % This is a second operator of reciprocal Perot reconstruction 
                            % (from vertices) applied to mean velocity.   
     
     K=[U0,V0,0,0;0,0,U0,V0]; 
     G=[ga;gb;gc];
     AU=Ut*Om+G*K*PP;       % This is a Fourier symbol of advection. 
                            % It is 3 by 3 matrix, and profile has to be
                            % added in vertical.
     % Coriolis 1: reconstruct to vertices, and then back to edges   
    A1=exp(i*(k/4-l/2));
    B1=exp(i*(-k/4-l/2));
    C1=exp(i*k/2);
   
    PPx=[-(A1+conj(A1))*a/4, -(B1+conj(B1))*a/4, (C1+conj(C1))*a/2]*2/3/a;
    PPy=[(A1+conj(A1))*h/2, -(B1+conj(B1))*h/2, 0 ]*2/3/a;                       
     % Coriolis 2: average to mid-edges and project on normals 
    PPc=(C1+conj(C1))*[ 0 -1]/2;
    PPa=(A1+conj(A1))*[a/4 h/6]*3/2/h;
    PPb=(B1+conj(B1))*[-a/4 h/6]*3/2/h;
    
    Cor=f0*[PPa;PPb;PPc]*[PPx;PPy];    
                            
     %%%% Coefficients of the derivatives of b in b equation
     % Velocities are obtained through PT*PP, i.e. Me on the 2D
     % level, and the result is De*diag([bua,bub,buc])*Me
     %%%% The M2 term in the buoyancy equation
     BM=[Dh(1,:);-Dh(2,:)]*diag([bua,bub,buc])*Me;
     %BM=[[bx,by]*Pe; [bx,by]*conj(Pe)];

     %%%% Advection of buoyancy
     % (1/2)a/(ah/2)[Ua*aa^2+Ub*bb^2+Uc*cc^2] u-triangle
     % - conj    d-triangle
     Ab=(Ua*aa^2+Ub*bb^2+Uc*cc^2)/h;
     %Ab=i*k*U/a+i*l*V/h;
     %%% sea surface height matrix
     Ea=g*ones([Nz,1])*ga;
     Eb=g*ones([Nz,1])*gb;
     Ec=g*ones([Nz,1])*gc;
     
     D2=Dh*Me; % Horizontal divergence 
     Du=dz*[D2(1,1)*ones([1,Nz]),D2(1,2)*ones([1,Nz]),D2(1,3)*ones([1,Nz])];
     Dd=dz*[D2(2,1)*ones([1,Nz]),D2(2,2)*ones([1,Nz]),D2(2,3)*ones([1,Nz])];
     
     %%%% Vertical divergence on u and d locations weighted with N2:
     D2=D2*N2;
    

     % Viscosity operators -curl curl + grad div
     Curl2=(1/a)*[conj(A1)-A1; B1-conj(B1);C1-conj(C1)];
     nu=0.005*a^3;
     Lapl=G*Dh-Curl2*Om;
     Visc=nu*Lapl*Lapl;
     Cor=Cor+Visc;
     % Diffusion operator
     LaplS=Dh*G;
     kappa=0.001*a^3;
     Dif=kappa*LaplS*LaplS;

     %%% Construct system matrix:
     S=-[diag(AU(1,1)*UV+Cor(1,1))-De(1,1)*W, diag(AU(1,2)*UV+Cor(1,2))-De(1,2)*W, ...
         diag(AU(1,3)*UV+Cor(1,3))-De(1,3)*W, ga(1)*P, ga(2)*P, Ea; ...
         diag(AU(2,1)*UV+Cor(2,1))-De(2,1)*W, diag(AU(2,2)*UV+Cor(2,2))-De(2,2)*W, ...
         diag(AU(2,3)*UV+Cor(2,3))-De(2,3)*W, gb(1)*P, gb(2)*P, Eb; ...
         diag(AU(3,1)*UV+Cor(3,1))-De(3,1)*W, diag(AU(3,2)*UV+Cor(3,2))-De(3,2)*W, ...
         diag(AU(3,3)*UV+Cor(3,3))-De(3,3)*W, gc(1)*P, gc(2)*P, Ec; ...
         -D2(1,1)*W+BM(1,1)*diag(ones([1,Nz])),-D2(1,2)*W+BM(1,2)*diag(ones([1,Nz])), ...
         -D2(1,3)*W+BM(1,3)*diag(ones([1,Nz])),...
         Dif(1,1)*diag(ones([1,Nz])),diag(Ab+Dif(1,2)),zeros([Nz,2]); ...
         -D2(2,1)*W+BM(2,1)*diag(ones([1,Nz])),-D2(2,2)*W+BM(2,2)*diag(ones([1,Nz])), ...
         -D2(2,3)*W+BM(2,3)*diag(ones([1,Nz])),...
         diag(-conj(Ab)+Dif(2,1)),Dif(2,2)*diag(ones([1,Nz])),zeros([Nz,2]); ...
         Du,zeros([1,2*Nz+2]); ...
         Dd,zeros([1,2*Nz+2])];
    
     %%% Find maximum unstable eigenvalue
     E=eig(S);
     w(n)=max(real(E));
 end

 df.f0 = f0;
 df.a = a * 1e-3;
 df.type = "ICON-mimetic";
 df.u0 = norm([U*ones(Nz,1)/Nz, V*ones(Nz,1)/Nz]);
 df.grid = "tri-C";
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

% Summary: (i) There is a spurious instability on the baroclinic axis. (ii)
% There is a spurious instability on the symmetric axis (Hollingsworth ?).
% (i) can be damped by viscosity and diffusion, the latter is important,
% without it viscosity should be rather strong. The behavior on the symmetric 
% axis is not as good as for the quadrilateral C grid. Furthermore, if the
% velocity profile starts from zero at the bottom, the spurious 
% instability gets stronger. (No Galilean invariance on the discrete level.)
% Triangular B grid has similar problems.





