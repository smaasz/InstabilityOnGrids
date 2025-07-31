% Eady setup, triangular B-grid, baroclinic and symmetric 
% instability
a=3.125e3;                % -- side of triangle
  h=a*sqrt(3)/2;           % -- the height of triangle 
  
  theta= 0;              % -- l=sqrt(3)Ksin(theta)/2, k=Kcos(theta)
  thetaw= 0;                % -- mean flow direction
  
  mom_a=1;

  f0=-0.0001;
  g=100000000;               % To effectively impose the rigid lid
  %g=10;
  N=0.001; N2=N*N;
  Ri=100;                  % -- Richardson number
  M2=abs(N*f0)/sqrt(Ri);        % -- M^2, i. e. the horizontal stratification
  Nz=8;                   % -- the number of vertical layers
  H0=4000;                 % -- fluid depth
  dz=H0/Nz;
  U=-(M2*dz/f0)*cos(thetaw)*[-(1:Nz)+(Nz+1)/2];
  V=-(M2*dz/f0)*sin(thetaw)*[-(1:Nz)+(Nz+1)/2];
  %U=-(M2*dz/f0)*cos(thetaw)*[-(0:Nz-1)-1/2+Nz];
  %V=-(M2*dz/f0)*sin(thetaw)*[-(0:Nz-1)-1/2+Nz];
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
  
  %P=zeros([Nz,Nz]);
  %for nz=1:Nz
  %P(nz,nz)=-0.5*dz; P(nz,1:nz-1)=-dz;
  %end
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
     
     % Phase shifts
     aa=exp(i*k/2+i*l/3); bb=exp(-i*k/2+i*l/3); cc=exp(-2*i*l/3);
     
     %%%% Horizontal gradients of velocities on u-triangles
     gxu=(1/a)*(aa-bb); gyu=(1/h)*(0.5*(aa+bb)-cc);
     gxd=-conj(gxu); gyd=-conj(gyu);
     
     %%%% Scalar gradients on u-triangles; they enter pressure grad and div
     gx=-conj(gxu); gy=-conj(gyu);
     %%%% d triangles: -conj(..)
     %%%% divergence (1/2)[-conj(gx),-conj(gy), gx, gy];
     
     %%%% Gradients from A-grid; they enter advection of buyoancy
     gxA=(2*i/3/a)*(sin(k)+0.5*sin(k/2+l)-0.5*sin(-k/2+l));
     gyA=(i/2/h)*(sin(k/2+l)+sin(-k/2+l));
     %gx=i*k/a; gy=i*l/a;     %%%
     
     %%%% Coefficients in the derivatives of b in b equation
     axu=(aa+bb+0.25*(aa+cc)+0.25*(bb+cc))/6;
     axv=(sqrt(3)/24)*(aa-bb);
     ayu=axv;
     ayv=0.125*(aa+2*cc+bb);
     %axu=0.5; ayv=0.5; axv=0.;ayu=0.;    %%%
     %gx=i*k/a; gy=i*l/a; 
     %gxA=i*k/a; gyA=i*l/a;
     M=conj(aa+bb+cc)/3;        % averaging to u triangle
     %M=1/3;
     %%%% Matrices acting on u and v in the u and v equations
     Au=diag(gxu*U+gyu*V);                  % on u d vel
     %Au=diag(i*k*U/a+i*l*V/a);
     Auu=-(M2/f0/2)*cos(thetaw)*(-gx)*W;    % on u d vel
     Cor=-f0*diag(ones([1,Nz]));            % on v u  
     Auv=-(M2/f0/2)*cos(thetaw)*(-gy)*W;    % on v d
     % +gx*P b +g*gx eta
     
     %Bv=diag(gx*U+gy*V)
     Bvv=-(M2/f0/2)*sin(thetaw)*(-gy)*W;     % on v d 
     Buv=-(M2/f0/2)*sin(thetaw)*(-gx)*W;     % on u d
     % Coriolis and advection can be taken from the A-part
     % +gy*P b +g*gy eta
     
     %Cb=diag(gxA*U+gyA*V);                   % on buoyancy
     % Lapl. correction
     ll=(1/9)*(2*cos(k)+2*cos(k/2+l)+2*cos(-k/2+l)-6);   % ll=Lapl*a^2/6;
     Cb=diag(gxA*U+gyA*V)*(1-ll);    % Correction to 4th order
                                     
     %Cb=diag(i*k*U/a+i*l*V/a);      % exact
     Cu=0.5*N2*W*(-gx);              % on u d
     Cuu=(bx*axu+by*ayu)*diag(ones([1,Nz])); % on u u
     Cv=0.5*N2*W*(-gy);              % on v d
     Cvv=(bx*axv+by*ayv)*diag(ones([1,Nz])); % on v u
     
     %%% sea surface height matrix
     Eu=g*gx*ones([Nz,1]);
     Ev=g*gy*ones([Nz,1]);
     Div=0.5*dz*[-conj(gx)*ones([1,Nz]),-conj(gy)*ones([1,Nz]), gx*ones([1,Nz]),gy*ones([1,Nz])];
     
     
     % Averaged pressure grad:
     gx1=0.5*(gx-conj(M*gx)); gy1=0.5*(gy-conj(M*gy));
     gx2=0.5*(-conj(gx)+M*gx); gy2=0.5*(-conj(gy)+M*gy);
     
     
     lll=(gxu*gxd+gyu*gyd);   % Laplacian acting on velocities
     Lapl=ll*6/a^2;           % Laplacian acting on scalars
     kappa=0.0007*a^3*lll^2;  % +kappad=0.0007 sufficient for mom_a=1 
     
     kappad=0.0005*a^3*Lapl^2;
     % if [-1,M^*;M,-1] is used, 12/a^2 gives a Laplacian
     %
     kappas=0.00005*(12)^2/a;  % corresponds to biharmonic;0.0005 for mom_adv=1
     kappa=0.;
     kappas=0.;
     kappad=0.;
     Visc=kappa*diag(ones([1,Nz]));
     Dif=kappad*diag(ones([1,Nz]));
     VV=diag(ones([1,Nz])); 
     v1=kappas*(1+M*conj(M)); v2=-2*kappas*conj(M);  v3=-2*kappas*M;
     
     %%%%%%%%%%%%% This is the inverted easy backscatter only for (mom_a=1)
     uu=kappas*(1-0.5*M*conj(M)+0.5*conj(M)^3); 
     uv=kappas*(-0.5*M^2-conj(M)+0.5*M*conj(M)^2);
     vv=kappas*(1-0.5*M*conj(M)+0.5*M^3); 
     vu=kappas*(-0.5*conj(M)^2-M+0.5*conj(M)*M^2);
     %%%%%%%%%%%%% This is our biharmonic viscosity
     uu=v1; vv=v1; uv=v2; vu=v3;
     %%%%%%%%%%%%%
     
     
     if mom_a==0,
     %%% Construct system matrix:
     S=-[-conj(Auu)*M+Visc+v1*VV, Cor-conj(Auv)*M, Au+Auu*M+v2*VV, Auv*M, gx*P, Eu; ... 
         -conj(Buv)*M-Cor,-conj(Bvv)*M+Visc+v1*VV, Buv*M, Bvv*M+Au+v2*VV, gy*P, Ev; ...
         -conj(Au)-conj(Auu*M)+v3*VV, -conj(Auv*M), Auu*conj(M)+Visc+v1*VV, Auv*conj(M)+Cor, -conj(gx)*P, -conj(Eu); ...
         -conj(Buv*M), -conj(Au)-conj(Bvv*M)+v3*VV, Buv*conj(M)-Cor, Bvv*conj(M)+Visc+v1*VV, -conj(gy)*P, -conj(Ev); ...
         Cuu-conj(Cu), Cvv-conj(Cv),Cu+conj(Cuu), Cv+conj(Cvv),Cb+Dif,zeros([Nz,1]); ....
         Div,zeros([1,Nz]),Cb(1,1)];
     end;

     if mom_a==1,
        F=0.5*(U*gxu+V*gyu); %*(1-ll*0.75); %+(ll*0.75)^2); 
        Au=diag(F);
     %%% Construct system matrix:
     S=-[-conj(Auu)*M+Au*M+Visc+uu*VV, Cor-conj(Auv)*M, -conj(Au)*M+Auu*M+uv*VV, Auv*M, gx*P, Eu; ... 
         -conj(Buv)*M-Cor,-conj(Bvv)*M+Au*M+Visc+uu*VV, Buv*M, Bvv*M-conj(Au)*M+uv*VV, gy*P, Ev; ...
         Au*conj(M)+vu*VV-conj(Auu*M), -conj(Auv*M), -conj(Au*M)+Auu*conj(M)+Visc+vv*VV, Auv*conj(M)+Cor, -conj(gx)*P, -conj(Eu); ...
         -conj(Buv*M), Au*conj(M)-conj(Bvv*M)+vu*VV, Buv*conj(M)-Cor, -conj(Au*M)+Bvv*conj(M)+Visc+vv*VV, -conj(gy)*P, -conj(Ev); ...
         Cuu-conj(Cu), Cvv-conj(Cv),Cu+conj(Cuu), Cv+conj(Cvv),Cb+Dif,zeros([Nz,1]); ....
         Div,zeros([1,Nz+1])];
     end;

     if mom_a==2,
     % U,V change sign with depth!!! Upwind is incorrect!!!    
     b1=(U*sqrt(3)/2+V/2)'*[0,(sqrt(3)*gxu/2+gyu/2)*h/3];   
     b2=(-U*sqrt(3)/2+V/2)'*[-conj(bb)+(sqrt(3)*gxd/2-gyd/2)*h/3, 1]*bb;
     b3=(-V)'*[-conj(cc)+(gyd)*h/3, 1]*cc;
     
     
     
     
     %b1=0.5*(U*sqrt(3)/2+V/2)'*[1-aa*(sqrt(3)*gxd/2+gyd/2)*h/3,aa+(sqrt(3)*gxu/2+gyu/2)*h/3];   
     %b2=0.5*(-U*sqrt(3)/2+V/2)'*[1+bb*(sqrt(3)*gxd/2-gyd/2)*h/3,bb-(sqrt(3)*gxu/2-gyu/2)*h/3];
     %b3=0.5*(-V)'*[1+cc*(gyd)*h/3, cc-(gyu)*h/3];
     
     %b1=(U*sqrt(3)/2+V/2)'*[1,aa]*0.5;   
     %b2=(-U*sqrt(3)/2+V/2)'*[1, bb]*0.5;
     %b3=(-V)'*[1,cc]*0.5;
     

     A1=(2/h)*diag(b1(:,1)+b2(:,1)+b3(:,1)); 
     A2=(2/h)*diag(b1(:,2)+b2(:,2)+b3(:,2));
     
     b1=-(U*sqrt(3)/2+V/2)'*[1,(sqrt(3)*gxu/2+gyu/2)*h/3-aa]*conj(aa);   
     b2=-(-U*sqrt(3)/2+V/2)'*[(sqrt(3)*gxd/2-gyd/2)*h/3, 0];
     b3=-(-V)'*[(gyd)*h/3, 0];
     
     %b1=-0.5*(U*sqrt(3)/2+V/2)'*[conj(aa)-(sqrt(3)*gxd/2+gyd/2)*h/3,1+conj(aa)*(sqrt(3)*gxu/2+gyu/2)*h/3];   
     %b2=-0.5*(-U*sqrt(3)/2+V/2)'*[conj(bb)+(sqrt(3)*gxd/2-gyd/2)*h/3, 1-conj(bb)*(sqrt(3)*gxu/2-gyu/2)*h/3];
     %b3=-0.5*(-V)'*[conj(cc)+(gyd)*h/3, 1-conj(cc)*(gyu)*h/3];
     %b1=-(U*sqrt(3)/2+V/2)'*[conj(aa),1]*0.5;   
     %b2=-(-U*sqrt(3)/2+V/2)'*[conj(bb), 1]*0.5;
     %b3=-(-V)'*[conj(cc), 1]*0.5;
     

     A3=(2/h)*diag(b1(:,1)+b2(:,1)+b3(:,1)); 
     A4=(2/h)*diag(b1(:,2)+b2(:,2)+b3(:,2));
     % This is the linear reconstruction upwind scheme. 
      
     %%% Construct system matrix:
     S=-[-conj(Auu)*M+A1+Visc+v1*VV, Cor-conj(Auv)*M, A2+Auu*M+v2*VV, Auv*M, gx*P, Eu; ... 
         -conj(Buv)*M-Cor,-conj(Bvv)*M+A1+Visc+v1*VV, Buv*M, Bvv*M+A2+v2*VV, gy*P, Ev; ...
         A3-conj(Auu*M)+v3*VV, -conj(Auv*M), A4+Auu*conj(M)+Visc+v1*VV, Auv*conj(M)+Cor, -conj(gx)*P, -conj(Eu); ...
         -conj(Buv*M), A3-conj(Bvv*M)+v3^VV, Buv*conj(M)-Cor, A4+Bvv*conj(M)+Visc+v1*VV, -conj(gy)*P, -conj(Ev); ...
         Cuu-conj(Cu), Cvv-conj(Cv),Cu+conj(Cuu), Cv+conj(Cvv),Cb+Dif,zeros([Nz,1]); ....
         Div,zeros([1,Nz+1])];
     end;
     %%% Find maximum unstable eigenvalue
     E=eig(S);
     w(n)=max(real(E));
 end

 df.f0 = f0;
 df.a = a * 1e-3;
 df.type = "flux";
 df.u0 = norm([U*ones(Nz,1)/Nz, V*ones(Nz,1)/Nz]);
 df.grid = "tri-B";
 df.N = N;
 df.Ri = Ri;
 df.theta = theta;
 df.thetaU = thetaw;
 df.date = datestr(now, 'yy-mm-dd-HH:MM');
 df.Nz = Nz;
 df.ks = Kvec/a;
 df.vs = w;

 data = jsonencode(df);

 % f = fopen('/Users/stmaas001/Projects/OceanFlows/InstabilityOnGrids/data/data.jsonl','a');
 % fprintf(f, '%s\n', data);
 % fclose(f);
 
plot(Kvec/a, w*N/abs(M2))
