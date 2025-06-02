% Eady setup, triangular hexagonal C-grid, baroclinic and symmetric 
% instability
  a=12500.0;               % -- side of triangle
  h=a*sqrt(3)/2;           % -- the height of triangle 
  
  theta=0;              % -- l=sqrt(3)Ksin(theta)/2, k=Kcos(theta)
  thetaw=0;             % -- mean flow direction
  
  f0=-0.0001;
  g=100000000;             % To effectively impose the rigid lid
  %g=10;
  N=0.001; N2=N*N;
  Ri=100;                  % -- Richardson number
  M2=abs(N*f0)/sqrt(Ri);        % -- M^2, i. e. the horizontal stratification
  Nz=64;                   % -- the number of vertical layers
  H0=4000;                 % -- fluid depth
  dz=H0/Nz;
  u0=0.0
  U=-(M2*dz/f0)*cos(thetaw)*[-(1:Nz)+(Nz+1)/2];
  V=-(M2*dz/f0)*sin(thetaw)*[-(1:Nz)+(Nz+1)/2];
  %U=-(M2*dz/f0)*cos(thetaw)*[-(1:Nz)+1/2];   % non-symmetric
  %V=-(M2*dz/f0)*sin(thetaw)*[-(1:Nz)+1/2];   % non-symmetric
  %U=-(M2*dz/f0)*cos(thetaw)*[-(1:Nz)+(Nz+1)/2]+u0*cos(thetaw); % shifted
  %V=-(M2*dz/f0)*sin(thetaw)*[-(1:Nz)+(Nz+1)/2]+u0*sin(thetaw); % shifted
  
  U0=-(M2*dz/f0)*cos(thetaw); V0=-(M2*dz/f0)*sin(thetaw);
  UV=[-(1:Nz)+(Nz+1)/2];   % Profile 
  %UV=[-(1:Nz)+1/2];       % non-symmetric
    
  bx=-M2*sin(thetaw);
  by=M2*cos(thetaw);

  Ua=U0; Ub=U0/2+sqrt(3)*V0/2;Uc=-U0/2+sqrt(3)*V0/2;    % Three normal velocities
                                                       % (*UV)
  M2a=-(M2/f0)*(cos(thetaw));                           % Their contrib. in M2-term
  M2b=-(M2/f0)*(cos(thetaw)/2+sin(thetaw)*sqrt(3)/2);  % in momentum
  M2c=-(M2/f0)*(-cos(thetaw)/2+sin(thetaw)*sqrt(3)/2);
     
  
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
     %%% phase multipliers:
     aa=sin(k/2); bb=sin(k/4+l/2); cc=sin(-k/4+l/2);
     %%% gradient
     Ga=2*i*aa/a; Gb=2*i*bb/a; Gc=2*i*cc/a;
     G=[Ga;Gb;Gc];
     %%% horizontal divergence
     Dh=(4*i/3/a)*[aa,bb,cc]; 

     %%% Momentum advection:
     % definitions for the grad K part, see in notes;
     Ki=(2/3)*[Ua*cos(k/2),Ub*cos(k/4+l/2),Uc*cos(-k/4+l/2)];
     AK=G*Ki; 

     % Energy-conserving xi\times u part: for historical reasons
     % the c-direction is treated as opposite here.
     % Uc=-Uc; However, change in the direction also enforce the change 
     % of sign with Uc, so seemingly no change is needed. 
     
     ah=exp(-i*l/3); bh=exp(-i*k/4+i*l/6); ch=exp(i*k/4+i*l/6);
     xiu=(2/h)*[ah,-bh,ch]; xid=-conj(xiu);     % relative vorticity   %% -ch for -Uc
     xia=(xiu*conj(ah)+xid*ah)/2;               % edge values
     xib=(xiu*conj(bh)+xid*bh)/2;
     xic=(xiu*conj(ch)+xid*ch)/2;
     
     % a-edge:
     al=exp(-i*k/4+i*l/2); be=exp(-i*k*3/4+i*l/2);
     ga=exp(-i*k*3/4-i*l/2); de=exp(-i*k/4-i*l/2); 
     %Sa=2*Ub*xib*al-Uc*xic*be+Ub*xib*ga-2*Uc*xic*de;
     %Sa=Sa-2*Uc*xic*conj(de)+Ub*xib*conj(ga)-Uc*xic*conj(be)+2*Ub*xib*conj(al);
     
     Sa=2*Ub*xib*al+Uc*xic*be+Ub*xib*ga+2*Uc*xic*de;
     Sa=Sa+2*Uc*xic*conj(de)+Ub*xib*conj(ga)+Uc*xic*conj(be)+2*Ub*xib*conj(al);
     
     
     % b-edge:
     al=exp(-i*k/2); be=exp(-i*k*3/4-i*l/2);
     ga=exp(-i*l); de=exp(i*k/4-i*l/2); 
     %Sb=-2*Uc*xic*al-Ua*xia*be-Uc*xic*ga-2*Ua*xia*de;
     %Sb=Sb-2*Ua*xia*conj(de)-Uc*xic*conj(ga)-Ua*xia*conj(be)-2*Uc*xic*conj(al);
     
     Sb=2*Uc*xic*al-Ua*xia*be+Uc*xic*ga-2*Ua*xia*de;
     Sb=Sb-2*Ua*xia*conj(de)+Uc*xic*conj(ga)-Ua*xia*conj(be)+2*Uc*xic*conj(al);
     
     % c-edge:
     al=exp(i*k/4+i*l/2); be=exp(i*l);
     ga=exp(-i*k*3/4+i*l/2); de=exp(-i*k/2); 
     Sc=2*Ua*xia*al+Ub*xib*be+Ua*xia*ga+2*Ub*xib*de;
     Sc=Sc+2*Ub*xib*conj(de)+Ua*xia*conj(ga)+Ub*xib*conj(be)+2*Ua*xia*conj(al);
     
     Axi=-[Sa;Sb;-Sc]/6/sqrt(3); %Uc=-Uc;  
     % This is energy-conserving variant
     
     % PV-conserving:
     xiu=(2/h)*[ah,-bh,ch]; xid=-conj(xiu);     % relative vorticity   %%
     xia=(xiu*conj(ah)+xid*ah)/2;                % edge values
     xib=(xiu*conj(bh)+xid*bh)/2;
     xic=(xiu*conj(ch)+xid*ch)/2;

     Axipv=[-V0*xia;(-V0/2+sqrt(3)*U0/2)*xib;(V0/2+sqrt(3)*U0/2)*xic];
     AU=AK+ Axipv;          % Advection of momentum (PV conserving)
     AU=AK+(Axi+Axipv)/2;   % Advection of momentum (energy conserving). Test
     %%%% I forgot that in the energy conserving form I need (xi_e+xi_e')/2
     %%%% for symmetry! Because of linearity of AU, this is achieved in
     %%%% this way.

     %%% Coriolis 
     fff=f0*2*h/3/a;
     CB=fff*(2*cos(0.25*k+l/2)+cos(0.75*k-l/2))/3;  
     CC=fff*(2*cos(0.25*k-l/2)+cos(0.75*k+l/2))/3;  
     CA=fff*(2*cos(0.5*k)+cos(l))/3;

     Cor=-[0,CC,CB;-CC,0,CA;-CB, -CA,0];
     %Test Coriolis - done, - is needed.
     
     %%%% Mean divergence at a,b,c edges:
     Ae=[cos(k/2); cos(k/4+l/2); cos(-k/4+l/2)]; 
     De=Ae*Dh;       %%% This is 3 by 3 matrix that acts on velocities to produce 
                    %%% the averaged divergence at velocity locations.  
     De=diag([M2a;M2b;M2c])*De; 
                    %%% Weighted with M2 term amplitudes (still 3 by 3)
     %%%% The M2 term in the buoyancy equation
     BM=(2/3)*[bx*cos(k/2),(bx/2+by*sqrt(3)/2)*cos(k/4+l/2), ...
               (-bx/2+by*sqrt(3)/2)*cos(-k/4+l/2)];
     %%%% Advection of buoyancy
     
     Ab=(2*i/3/a)*(Ua*sin(k)+Ub*sin(k/2+l)+Uc*sin(-k/2+l));
     ll1=-(2/3)*(sin(k/2))^2;
     ll2=-(2/3)*(sin(k/4+l/2))^2;
     ll3=-(2/3)*(sin(-k/4+l/2))^2;
     
     Ab=(2*i/3/a)*(Ua*sin(k)*(1-ll1)+Ub*sin(k/2+l)*(1-ll2)+Uc*sin(-k/2+l)*(1-ll3));
     %Ab=i*k*U0/a+i*l*V0/h;

     %%% sea surface height matrix
     Ea=g*ones([Nz,1])*Ga;
     Eb=g*ones([Nz,1])*Gb;
     Ec=g*ones([Nz,1])*Gc;
     
     % Horizontal divergence 
     Du=dz*[Dh(1)*ones([1,Nz]),Dh(2)*ones([1,Nz]),Dh(3)*ones([1,Nz])];
     
     %%%% Vertical divergence on u and d locations weighted with N2:
     D2=Dh*N2;
    

     % Viscosity operators -curl curl + grad div
     Curl=(2/h)*[ah,-bh,ch;-conj(ah),-conj(-bh),-conj(ch)]; 
     Curl2=-(3/2/h)*[conj(ah), -ah; conj(bh), -bh; conj(ch), -ch];
     %Curl2=(1/a)*[conj(A1)-A1; B1-conj(B1);C1-conj(C1)];
     nu=0.002*a^3;
     Lapl=G*Dh-Curl2*Curl;
     Visc=nu*Lapl*Lapl;
     Cor=Cor+Visc;
     % Diffusion operator
     LaplS=Dh*G;
     kappa=0.001*a^3;
     Dif=kappa*LaplS*LaplS;
     %%% Construct system matrix:
     S=-[diag(AU(1,1)*UV+Cor(1,1))-De(1,1)*W, diag(AU(1,2)*UV+Cor(1,2))-De(1,2)*W, ...
         diag(AU(1,3)*UV+Cor(1,3))-De(1,3)*W, Ga*P, Ea; ...
         diag(AU(2,1)*UV+Cor(2,1))-De(2,1)*W, diag(AU(2,2)*UV+Cor(2,2))-De(2,2)*W, ...
         diag(AU(2,3)*UV+Cor(2,3))-De(2,3)*W, Gb*P, Eb; ...
         diag(AU(3,1)*UV+Cor(3,1))-De(3,1)*W, diag(AU(3,2)*UV+Cor(3,2))-De(3,2)*W, ...
         diag(AU(3,3)*UV+Cor(3,3))-De(3,3)*W, Gc*P, Ec; ...
         -D2(1)*W+BM(1)*diag(ones([1,Nz])),-D2(2)*W+BM(2)*diag(ones([1,Nz])), ...
         -D2(3)*W+BM(3)*diag(ones([1,Nz])), diag(Ab*UV+Dif),zeros([Nz,1]); ...
          Du,zeros([1,Nz+1])];
    
     %%% Find maximum unstable eigenvalue
     E=eig(S);
     w(n)=max(real(E));
 end

 df.f0 = f0;
 df.a = a * 1e-3;
 df.type = "MPAS";
 df.grid = "hex-C";
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
