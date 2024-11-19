% Script file run_acoustics2d_free
%
% Solve free space with CRBC on four sides - use potential version of CRBC 
%

ndg=9;
alpha=1;
nx=20;
ny=20;
c=1;
ttot=20;
nsteps=10000;

x=linspace(-1,1,nx+1);
hx=2/nx;
y=linspace(-1,1,ny+1);
hy=2/ny;
dt=ttot/nsteps;

[DD,Lift,s]=DGsetup(ndg);
Dx=DD/hx;
Liftx=Lift/hx;
Dy=DD/hy;
Lifty=Lift/hy;

p=zeros(ndg+1,ndg+1,nx,ny);
u=zeros(ndg+1,ndg+1,nx,ny);
v=zeros(ndg+1,ndg+1,nx,ny);

for k=1:ny
  for j=1:nx
    for kk=1:ndg+1
      for jj=1:ndg+1
         yloc=.5*(y(k)+y(k+1))+(.5*hy)*s(kk);
         xloc=.5*(x(j)+x(j+1))+(.5*hx)*s(jj);
         p(jj,kk,j,k)=10*exp(-40*(xloc^2+yloc^2));
      end
    end
  end
end

% Set up crbc - we just make nbc=7 and pretend a point source at the
% origin created the solution

nbc=7;
delta=1;
eta=delta/(c*ttot);
[a emax]=optimal_cosinesP(eta,nbc-1);
sig=(1-a.^2)./(ttot*a);

disp(sprintf('emax=%3.2e',emax))

% Factor the corner matrix

[cmL cmU pc]=acoustics2d_CRBC_cornermat(c,a,nbc);

pE=zeros(ndg+1,nbc+1,ny);
uE=zeros(ndg+1,nbc+1,ny);
vE=zeros(ndg+1,nbc+1,ny);
pW=zeros(ndg+1,nbc+1,ny);
uW=zeros(ndg+1,nbc+1,ny);
vW=zeros(ndg+1,nbc+1,ny);
pN=zeros(ndg+1,nbc+1,nx);
uN=zeros(ndg+1,nbc+1,nx);
vN=zeros(ndg+1,nbc+1,nx);
pS=zeros(ndg+1,nbc+1,nx);
uS=zeros(ndg+1,nbc+1,nx);
vS=zeros(ndg+1,nbc+1,nx);

pEdat=zeros(ndg+1,ny);
uEdat=zeros(ndg+1,ny);
pWdat=zeros(ndg+1,ny);
uWdat=zeros(ndg+1,ny);
pNdat=zeros(ndg+1,nx);
vNdat=zeros(ndg+1,nx);
pSdat=zeros(ndg+1,nx);
vSdat=zeros(ndg+1,nx);

pNE=zeros(nbc+1,nbc+1);
uNE=zeros(nbc+1,nbc+1);
vNE=zeros(nbc+1,nbc+1);
pNW=zeros(nbc+1,nbc+1);
uNW=zeros(nbc+1,nbc+1);
vNW=zeros(nbc+1,nbc+1);
pSE=zeros(nbc+1,nbc+1);
uSE=zeros(nbc+1,nbc+1);
vSE=zeros(nbc+1,nbc+1);
pSW=zeros(nbc+1,nbc+1);
uSW=zeros(nbc+1,nbc+1);
vSW=zeros(nbc+1,nbc+1);

pNEdatN=zeros(nbc+1,1);
uNEdatN=zeros(nbc+1,1);
pNEdatE=zeros(nbc+1,1);
vNEdatE=zeros(nbc+1,1);
pSEdatS=zeros(nbc+1,1);
uSEdatS=zeros(nbc+1,1);
pSEdatE=zeros(nbc+1,1);
vSEdatE=zeros(nbc+1,1);
pNWdatN=zeros(nbc+1,1);
uNWdatN=zeros(nbc+1,1);
pNWdatW=zeros(nbc+1,1);
vNWdatW=zeros(nbc+1,1);
pSWdatS=zeros(nbc+1,1);
uSWdatS=zeros(nbc+1,1);
pSWdatW=zeros(nbc+1,1);
vSWdatW=zeros(nbc+1,1);

% March with old reliable RK4

crk=dt*[1/2; 1/2; 1];
brk=dt*[1/6; 1/3; 1/3; 1/6];

% Necessary arrays  

pn=zeros(ndg+1,1);
un=zeros(ndg+1,1);
vn=zeros(ndg+1,1);
ppx=zeros(2,ndg+1);
ppy=zeros(2,ndg+1);
upx=zeros(2,ndg+1);
vpy=zeros(2,ndg+1);
ppcrbc=zeros(2,nbc+1);
upcrbc=zeros(2,nbc+1);
vpcrbc=zeros(2,nbc+1);
dpN=zeros(nbc+1);
duN=zeros(nbc+1);
dpE=zeros(nbc+1);
dvE=zeros(nbc+1);

  pt=zeros(size(p));
  ut=zeros(size(u));
  vt=zeros(size(v));
  pEt=zeros(size(pE));
  uEt=zeros(size(uE));
  vEt=zeros(size(vE));
  pWt=zeros(size(pW));
  uWt=zeros(size(uW));
  vWt=zeros(size(vW));
  pNt=zeros(size(pN));
  uNt=zeros(size(uN));
  vNt=zeros(size(vN));
  pSt=zeros(size(pS));
  uSt=zeros(size(uS));
  vSt=zeros(size(vS));
  pNEt=zeros(size(pNE));
  uNEt=zeros(size(uNE));
  vNEt=zeros(size(vNE));
  pSEt=zeros(size(pSE));
  uSEt=zeros(size(uSE));
  vSEt=zeros(size(vSE));
  pNWt=zeros(size(pNW));
  uNWt=zeros(size(uNW));
  vNWt=zeros(size(vNW));
  pSWt=zeros(size(pSW));
  uSWt=zeros(size(uSW));
  vSWt=zeros(size(vSW));

for it=1:nsteps    
  pnew=p;
  unew=u;
  vnew=v;
  pEnew=pE;
  uEnew=uE;
  vEnew=vE;
  pWnew=pW;
  uWnew=uW;
  vWnew=vW;
  pNnew=pN;
  uNnew=uN;
  vNnew=vN;
  pSnew=pS;
  uSnew=uS;
  vSnew=vS;
  pNEnew=pNE;
  uNEnew=uNE;
  vNEnew=vNE;
  pNWnew=pNW;
  uNWnew=uNW;
  vNWnew=vNW;
  pSEnew=pSE;
  uSEnew=uSE;
  vSEnew=vSE;
  pSWnew=pSW;
  uSWnew=uSW;
  vSWnew=vSW;
  pst=p;
  ust=u;
  vst=v;
  pEst=pE;
  uEst=uE;
  vEst=vE;
  pWst=pW;
  uWst=uW;
  vWst=vW;
  pNst=pN;
  uNst=uN;
  vNst=vN;
  pSst=pS;
  uSst=uS;
  vSst=vS;
  pNEst=pNE;
  uNEst=uNE;
  vNEst=vNE;
  pNWst=pNW;
  uNWst=uNW;
  vNWst=vNW;
  pSEst=pSE;
  uSEst=uSE;
  vSEst=vSE;
  pSWst=pSW;
  uSWst=uSW;
  vSWst=vSW;
%
  for ist=1:4
%
% compute time derivatives - first corners then edges then volume 
%
    dpN=pNst(ndg+1,:,nx);
    duN=uNst(ndg+1,:,nx);
    dpE=pEst(ndg+1,:,ny);
    dvE=vEst(ndg+1,:,ny);
    [pNEt uNEt vNEt pNEdatN uNEdatN pNEdatE vNEdatE]=acoustics2d_CRBC_corner(dpN,duN,dpE,dvE,pNEst,uNEst,vNEst,c,sig,nbc,cmL,cmU,pc); 
%
    dpN=pSst(ndg+1,:,nx);
    duN=uSst(ndg+1,:,nx);
    dpE=pEst(1,:,1);
    dvE=-vEst(1,:,1);
    [pSEt uSEt vSEt pSEdatS uSEdatS pSEdatE vSEdatE]=acoustics2d_CRBC_corner(dpN,duN,dpE,dvE,pSEst,uSEst,vSEst,c,sig,nbc,cmL,cmU,pc); 
%
    dpN=pNst(1,:,1);
    duN=-uNst(1,:,1);
    dpE=pWst(ndg+1,:,ny);
    dvE=vWst(ndg+1,:,ny);
    [pNWt uNWt vNWt pNWdatN uNWdatN pNWdatW vNWdatW]=acoustics2d_CRBC_corner(dpN,duN,dpE,dvE,pNWst,uNWst,vNWst,c,sig,nbc,cmL,cmU,pc); 
%
    dpN=pSst(1,:,1);
    duN=-uSst(1,:,1);
    dpE=pWst(1,:,1);
    dvE=-vWst(1,:,1);
    [pSWt uSWt vSWt pSWdatS uSWdatS pSWdatW vSWdatW]=acoustics2d_CRBC_corner(dpN,duN,dpE,dvE,pSWst,uSWst,vSWst,c,sig,nbc,cmL,cmU,pc); 
%
    for k=1:ny
      if (k==1)
        ppcrbc(1,1:nbc+1)=pSWdatW(1:nbc+1,1)';     
        vpcrbc(1,1:nbc+1)=-vSWdatW(1:nbc+1,1)';     
      else
        ppcrbc(1,1:nbc+1)=(squeeze(pWst(ndg+1,1:nbc+1,k-1)));     
        vpcrbc(1,1:nbc+1)=(squeeze(vWst(ndg+1,1:nbc+1,k-1)));     
      end
      if (k==ny)
        ppcrbc(2,1:nbc+1)=pNWdatW(1:nbc+1,1)';     
        vpcrbc(2,1:nbc+1)=vNWdatW(1:nbc+1,1)';     
      else
        ppcrbc(2,1:nbc+1)=(squeeze(pWst(1,1:nbc+1,k+1)));     
        vpcrbc(2,1:nbc+1)=(squeeze(vWst(1,1:nbc+1,k+1)));     
      end
      pn=pst(1,:,1,k);
      un=-ust(1,:,1,k);
      [pWt(:,:,k) uWt(:,:,k) vWt(:,:,k) pWdat(:,k) uWdat(:,k)]=acoustics2d_CRBC_face(pn,un,ppcrbc,vpcrbc, ...
							       pWst(:,:,k),uWst(:,:,k),vWst(:,:,k),c,Dy,Lifty,alpha,ndg,a,sig,nbc);
    end
%
    for k=1:ny
      if (k==1)
        ppcrbc(1,1:nbc+1)=pSEdatE(1:nbc+1,1)';     
        vpcrbc(1,1:nbc+1)=-vSEdatE(1:nbc+1,1)';     
      else
        ppcrbc(1,1:nbc+1)=(squeeze(pEst(ndg+1,1:nbc+1,k-1)));     
        vpcrbc(1,1:nbc+1)=(squeeze(vEst(ndg+1,1:nbc+1,k-1)));     
      end
      if (k==ny)
        ppcrbc(2,1:nbc+1)=pNEdatE(1:nbc+1,1)';     
        vpcrbc(2,1:nbc+1)=vNEdatE(1:nbc+1,1)';     
      else
        ppcrbc(2,1:nbc+1)=(squeeze(pEst(1,1:nbc+1,k+1)));     
        vpcrbc(2,1:nbc+1)=(squeeze(vEst(1,1:nbc+1,k+1)));     
      end
      pn=pst(ndg+1,:,nx,k);
      un=ust(ndg+1,:,nx,k);
      [pEt(:,:,k) uEt(:,:,k) vEt(:,:,k) pEdat(:,k) uEdat(:,k)]=acoustics2d_CRBC_face(pn,un,ppcrbc,vpcrbc, ...
							       pEst(:,:,k),uEst(:,:,k),vEst(:,:,k),c,Dy,Lifty,alpha,ndg,a,sig,nbc);
    end
%
    for j=1:nx 
      if (j==1)
        ppcrbc(1,1:nbc+1)=pSWdatS(1:nbc+1,1)';     
        upcrbc(1,1:nbc+1)=-uSWdatS(1:nbc+1,1)';     
      else
        ppcrbc(1,1:nbc+1)=(squeeze(pSst(ndg+1,1:nbc+1,j-1)));     
        upcrbc(1,1:nbc+1)=(squeeze(uSst(ndg+1,1:nbc+1,j-1)));     
      end
      if (j==nx)
        ppcrbc(2,1:nbc+1)=pSEdatS(1:nbc+1,1)';     
        upcrbc(2,1:nbc+1)=uSEdatS(1:nbc+1,1)';     
      else
        ppcrbc(2,1:nbc+1)=(squeeze(pSst(1,1:nbc+1,j+1)));     
        upcrbc(2,1:nbc+1)=(squeeze(uSst(1,1:nbc+1,j+1)));     
      end
      pn=pst(:,1,j,1)';
      vn=-vst(:,1,j,1)';
      [pSt(:,:,j) vSt(:,:,j) uSt(:,:,j) pSdat(:,j) vSdat(:,j)]=acoustics2d_CRBC_face(pn,vn,ppcrbc,upcrbc, ...
							       pSst(:,:,j),vSst(:,:,j),uSst(:,:,j),c,Dx,Liftx,alpha,ndg,a,sig,nbc);
    end
%
    for j=1:nx 
      if (j==1)
        ppcrbc(1,1:nbc+1)=pNWdatN(1:nbc+1,1)';     
        upcrbc(1,1:nbc+1)=-uNWdatN(1:nbc+1,1)';     
      else
        ppcrbc(1,1:nbc+1)=(squeeze(pNst(ndg+1,1:nbc+1,j-1)));     
        upcrbc(1,1:nbc+1)=(squeeze(uNst(ndg+1,1:nbc+1,j-1)));     
      end
      if (j==nx)
        ppcrbc(2,1:nbc+1)=pNEdatN(1:nbc+1,1)';     
        upcrbc(2,1:nbc+1)=uNEdatN(1:nbc+1,1)';     
      else
        ppcrbc(2,1:nbc+1)=(squeeze(pNst(1,1:nbc+1,j+1)));     
        upcrbc(2,1:nbc+1)=(squeeze(uNst(1,1:nbc+1,j+1)));     
      end
      pn=pst(:,ndg+1,j,ny)';
      vn=vst(:,ndg+1,j,ny)';
      [pNt(:,:,j) vNt(:,:,j) uNt(:,:,j) pNdat(:,j) vNdat(:,j)]=acoustics2d_CRBC_face(pn,vn,ppcrbc,upcrbc, ...
							       pNst(:,:,j),vNst(:,:,j),uNst(:,:,j),c,Dx,Liftx,alpha,ndg,a,sig,nbc);
    end
%
    for k=1:ny
    for j=1:nx
      if (k==1)
        ppy(1,1:ndg+1)=pSdat(1:ndg+1,j)';     
        vpy(1,1:ndg+1)=-vSdat(1:ndg+1,j)';     
      else
        ppy(1,1:ndg+1)=(squeeze(pst(1:ndg+1,ndg+1,j,k-1)))';     
        vpy(1,1:ndg+1)=(squeeze(vst(1:ndg+1,ndg+1,j,k-1)))';     
      end
      if (k==ny)
        ppy(2,1:ndg+1)=pNdat(1:ndg+1,j)';     
        vpy(2,1:ndg+1)=vNdat(1:ndg+1,j)';     
      else
        ppy(2,1:ndg+1)=(squeeze(pst(1:ndg+1,1,j,k+1)))';     
        vpy(2,1:ndg+1)=(squeeze(vst(1:ndg+1,1,j,k+1)))';     
      end
      if (j==1)
        ppx(1,1:ndg+1)=pWdat(1:ndg+1,k)'; 
        upx(1,1:ndg+1)=-uWdat(1:ndg+1,k)';     
      else
        ppx(1,1:ndg+1)=(squeeze(pst(ndg+1,1:ndg+1,j-1,k)));     
        upx(1,1:ndg+1)=(squeeze(ust(ndg+1,1:ndg+1,j-1,k)));     
      end
      if (j==nx)
        ppx(2,1:ndg+1)=pEdat(1:ndg+1,k)'; 
        upx(2,1:ndg+1)=uEdat(1:ndg+1,k)';     
      else
        ppx(2,1:ndg+1)=(squeeze(pst(1,1:ndg+1,j+1,k)));     
        upx(2,1:ndg+1)=(squeeze(ust(1,1:ndg+1,j+1,k)));     
      end
      [pt(:,:,j,k),ut(:,:,j,k),vt(:,:,j,k)]=acoustics2d(pst(:,:,j,k),ust(:,:,j,k),vst(:,:,j,k), ...
                                                        ppx,upx,ppy,vpy,c,Dx,Liftx,Dy,Lifty,alpha,ndg);
    end
    end
%
    if (ist < 4)
      pnew=pnew+brk(ist)*pt;
      unew=unew+brk(ist)*ut;
      vnew=vnew+brk(ist)*vt;
      pWnew=pWnew+brk(ist)*pWt;
      uWnew=uWnew+brk(ist)*uWt;
      vWnew=vWnew+brk(ist)*vWt;
      pEnew=pEnew+brk(ist)*pEt;
      uEnew=uEnew+brk(ist)*uEt;
      vEnew=vEnew+brk(ist)*vEt;
      pNnew=pNnew+brk(ist)*pNt;
      uNnew=uNnew+brk(ist)*uNt;
      vNnew=vNnew+brk(ist)*vNt;
      pSnew=pSnew+brk(ist)*pSt;
      uSnew=uSnew+brk(ist)*uSt;
      vSnew=vSnew+brk(ist)*vSt;
      pNEnew=pNEnew+brk(ist)*pNEt;
      uNEnew=uNEnew+brk(ist)*uNEt;
      vNEnew=vNEnew+brk(ist)*vNEt;
      pNWnew=pNWnew+brk(ist)*pNWt;
      uNWnew=uNWnew+brk(ist)*uNWt;
      vNWnew=vNWnew+brk(ist)*vNWt;
      pSEnew=pSEnew+brk(ist)*pSEt;
      uSEnew=uSEnew+brk(ist)*uSEt;
      vSEnew=vSEnew+brk(ist)*vSEt;
      pSWnew=pSWnew+brk(ist)*pSWt;
      uSWnew=uSWnew+brk(ist)*uSWt;
      vSWnew=vSWnew+brk(ist)*vSWt;
      pst=p+crk(ist)*pt;
      ust=u+crk(ist)*ut;
      vst=v+crk(ist)*vt;
      pWst=pW+crk(ist)*pWt;
      uWst=uW+crk(ist)*uWt;
      vWst=vW+crk(ist)*vWt;
      pEst=pE+crk(ist)*pEt;
      uEst=uE+crk(ist)*uEt;
      vEst=vE+crk(ist)*vEt;
      pNst=pN+crk(ist)*pNt;
      uNst=uN+crk(ist)*uNt;
      vNst=vN+crk(ist)*vNt;
      pSst=pS+crk(ist)*pSt;
      uSst=uS+crk(ist)*uSt;
      vSst=vS+crk(ist)*vSt;
      pNEst=pNE+crk(ist)*pNEt;
      uNEst=uNE+crk(ist)*uNEt;
      vNEst=vNE+crk(ist)*vNEt;
      pNWst=pNW+crk(ist)*pNWt;
      uNWst=uNW+crk(ist)*uNWt;
      vNWst=vNW+crk(ist)*vNWt;
      pSEst=pSE+crk(ist)*pSEt;
      uSEst=uSE+crk(ist)*uSEt;
      vSEst=vSE+crk(ist)*vSEt;
      pSWst=pSW+crk(ist)*pSWt;
      uSWst=uSW+crk(ist)*uSWt;
      vSWst=vSW+crk(ist)*vSWt;
    else
      p=pnew+brk(ist)*pt;
      u=unew+brk(ist)*ut;
      v=vnew+brk(ist)*vt;
      pW=pWnew+brk(ist)*pWt;
      uW=uWnew+brk(ist)*uWt;
      vW=vWnew+brk(ist)*vWt;
      pE=pEnew+brk(ist)*pEt;
      uE=uEnew+brk(ist)*uEt;
      vE=vEnew+brk(ist)*vEt;
      pN=pNnew+brk(ist)*pNt;
      uN=uNnew+brk(ist)*uNt;
      vN=vNnew+brk(ist)*vNt;
      pS=pSnew+brk(ist)*pSt;
      uS=uSnew+brk(ist)*uSt;
      vS=vSnew+brk(ist)*vSt;
      pNE=pNEnew+brk(ist)*pNEt;
      uNE=uNEnew+brk(ist)*uNEt;
      vNE=vNEnew+brk(ist)*vNEt;
      pNW=pNWnew+brk(ist)*pNWt;
      uNW=uNWnew+brk(ist)*uNWt;
      vNW=vNWnew+brk(ist)*vNWt;
      pSE=pSEnew+brk(ist)*pSEt;
      uSE=uSEnew+brk(ist)*uSEt;
      vSE=vSEnew+brk(ist)*vSEt;
      pSW=pSWnew+brk(ist)*pSWt;
      uSW=uSWnew+brk(ist)*uSWt;
      vSW=vSWnew+brk(ist)*vSWt;
    end
  end
%
  if (mod(it,50)==0)
    jwr=100+it/50;
    cj=num2str(jwr,'%i');
    fname=strcat('freeup_nbc7.',cj);
    fid=fopen(fname,'w');
    energy=0;
    vfact=hx*hy/((ndg+1)^2);
    t=dt*it;
    for k=1:ny
    for j=1:nx 
    for kk=1:ndg+1
    for jj=1:ndg+1
      energy=energy+vfact*(p(jj,kk,j,k)^2+u(jj,kk,j,k)^2+v(jj,kk,j,k)^2);
      fprintf(fid,'%16.14f %16.14f %16.14f \n',p(jj,kk,j,k),u(jj,kk,j,k),v(jj,kk,j,k));
    end
    end
    end
    end
    fclose(fid);
    disp(sprintf('t=%4.2f energy=%4.2f',t,sqrt(energy)))
  end
%
end



