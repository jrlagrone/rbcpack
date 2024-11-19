% Script file run_acoustics2d
%
% Solve in a duct with CRBC at both ends - use potential version of CRBC 
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
eta=1/(c*ttot);
[a emax]=optimal_cosinesP(eta,nbc-1);
sig=(1-a.^2)./(ttot*a);

disp(sprintf('emax=%3.2e',emax))

pr=zeros(ndg+1,nbc+1,ny);
ur=zeros(ndg+1,nbc+1,ny);
vr=zeros(ndg+1,nbc+1,ny);
pl=zeros(ndg+1,nbc+1,ny);
ul=zeros(ndg+1,nbc+1,ny);
vl=zeros(ndg+1,nbc+1,ny);
prdat=zeros(ndg+1,ny);
urdat=zeros(ndg+1,ny);
pldat=zeros(ndg+1,ny);
uldat=zeros(ndg+1,ny); 


% March with old reliable RK4

crk=dt*[1/2; 1/2; 1];
brk=dt*[1/6; 1/3; 1/3; 1/6];

% Necessary arrays

pn=zeros(ndg+1,1);
un=zeros(ndg+1,1);
ppx=zeros(2,ndg+1);
ppy=zeros(2,ndg+1);
upx=zeros(2,ndg+1);
vpy=zeros(2,ndg+1);
ppcrbc=zeros(2,nbc+1);
vpcrbc=zeros(2,nbc+1);

for it=1:nsteps 
  pnew=p;
  unew=u;
  vnew=v;
  prnew=pr;
  urnew=ur;
  vrnew=vr;
  plnew=pl;
  ulnew=ul;
  vlnew=vl;
  pst=p;
  ust=u;
  vst=v;
  prst=pr;
  urst=ur;
  vrst=vr;
  plst=pl;
  ulst=ul;
  vlst=vl;
  pt=zeros(size(p));
  ut=zeros(size(u));
  vt=zeros(size(v));
  prt=zeros(size(pr));
  urt=zeros(size(ur));
  vrt=zeros(size(vr));
  plt=zeros(size(pl));
  ult=zeros(size(ul));
  vlt=zeros(size(vl));
%
  for ist=1:4
%
% compute time derivatives - for the potential version we need to do CRBC first 
%
    for k=1:ny
      if (k==1)
        ppcrbc(1,1:nbc+1)=(squeeze(plst(1,1:nbc+1,k)));     
        vpcrbc(1,1:nbc+1)=-(squeeze(vlst(1,1:nbc+1,k)));     
      else
        ppcrbc(1,1:nbc+1)=(squeeze(plst(ndg+1,1:nbc+1,k-1)));     
        vpcrbc(1,1:nbc+1)=(squeeze(vlst(ndg+1,1:nbc+1,k-1)));     
      end
      if (k==ny)
        ppcrbc(2,1:nbc+1)=(squeeze(plst(ndg+1,1:nbc+1,k)));     
        vpcrbc(2,1:nbc+1)=-(squeeze(vlst(ndg+1,1:nbc+1,k)));     
      else
        ppcrbc(2,1:nbc+1)=(squeeze(plst(1,1:nbc+1,k+1)));     
        vpcrbc(2,1:nbc+1)=(squeeze(vlst(1,1:nbc+1,k+1)));     
      end
      pn=pst(1,:,1,k);
      un=-ust(1,:,1,k);
      [plt(:,:,k) ult(:,:,k) vlt(:,:,k) pldat(:,k) uldat(:,k)]=acoustics2d_CRBC_face(pn,un,ppcrbc,vpcrbc, ...
                                                              plst(:,:,k),ulst(:,:,k),vlst(:,:,k),c,Dy,Lifty,alpha,ndg,a,sig,nbc);
    end
%
    for k=1:ny
      if (k==1)
        ppcrbc(1,1:nbc+1)=(squeeze(prst(1,1:nbc+1,k)));     
        vpcrbc(1,1:nbc+1)=-(squeeze(vrst(1,1:nbc+1,k)));     
      else
        ppcrbc(1,1:nbc+1)=(squeeze(prst(ndg+1,1:nbc+1,k-1)));     
        vpcrbc(1,1:nbc+1)=(squeeze(vrst(ndg+1,1:nbc+1,k-1)));     
      end
      if (k==ny)
        ppcrbc(2,1:nbc+1)=(squeeze(prst(ndg+1,1:nbc+1,k)));     
        vpcrbc(2,1:nbc+1)=-(squeeze(vrst(ndg+1,1:nbc+1,k)));     
      else
        ppcrbc(2,1:nbc+1)=(squeeze(prst(1,1:nbc+1,k+1)));     
        vpcrbc(2,1:nbc+1)=(squeeze(vrst(1,1:nbc+1,k+1)));     
      end
      pn=pst(ndg+1,:,nx,k);
      un=ust(ndg+1,:,nx,k);
      [prt(:,:,k) urt(:,:,k) vrt(:,:,k) prdat(:,k) urdat(:,k)]=acoustics2d_CRBC_face(pn,un,ppcrbc,vpcrbc, ...
                                                               prst(:,:,k),urst(:,:,k),vrst(:,:,k),c,Dy,Lifty,alpha,ndg,a,sig,nbc);
    end
!
    for k=1:ny
    for j=1:nx
      if (k==1)
        ppy(1,1:ndg+1)=(squeeze(pst(1:ndg+1,1,j,k)))';     
        vpy(1,1:ndg+1)=-(squeeze(vst(1:ndg+1,1,j,k)))';     
      else
        ppy(1,1:ndg+1)=(squeeze(pst(1:ndg+1,ndg+1,j,k-1)))';     
        vpy(1,1:ndg+1)=(squeeze(vst(1:ndg+1,ndg+1,j,k-1)))';     
      end
      if (k==ny)
        ppy(2,1:ndg+1)=(squeeze(pst(1:ndg+1,ndg+1,j,k)))';     
        vpy(2,1:ndg+1)=-(squeeze(vst(1:ndg+1,ndg+1,j,k)))';     
      else
        ppy(2,1:ndg+1)=(squeeze(pst(1:ndg+1,1,j,k+1)))';     
        vpy(2,1:ndg+1)=(squeeze(vst(1:ndg+1,1,j,k+1)))';     
      end
      if (j==1)
        ppx(1,1:ndg+1)=pldat(1:ndg+1,k)';   % CRBC left     
        upx(1,1:ndg+1)=-uldat(1:ndg+1,k)';
      else
        ppx(1,1:ndg+1)=(squeeze(pst(ndg+1,1:ndg+1,j-1,k)));     
        upx(1,1:ndg+1)=(squeeze(ust(ndg+1,1:ndg+1,j-1,k)));     
      end
      if (j==nx)
        ppx(2,1:ndg+1)=prdat(1:ndg+1,k)';     % CRBC right 
        upx(2,1:ndg+1)=urdat(1:ndg+1,k)';     
      else
        ppx(2,1:ndg+1)=(squeeze(pst(1,1:ndg+1,j+1,k)));     
        upx(2,1:ndg+1)=(squeeze(ust(1,1:ndg+1,j+1,k)));     
      end
      [pt(:,:,j,k),ut(:,:,j,k),vt(:,:,j,k)]=acoustics2d(pst(:,:,j,k),ust(:,:,j,k),vst(:,:,j,k), ...
                                                        ppx,upx,ppy,vpy,c,Dx,Liftx,Dy,Lifty,alpha,ndg);
    end
    end
%
%
    if (ist < 4)
      pnew=pnew+brk(ist)*pt;
      unew=unew+brk(ist)*ut;
      vnew=vnew+brk(ist)*vt;
      plnew=plnew+brk(ist)*plt;
      ulnew=ulnew+brk(ist)*ult;
      vlnew=vlnew+brk(ist)*vlt;
      prnew=prnew+brk(ist)*prt;
      urnew=urnew+brk(ist)*urt;
      vrnew=vrnew+brk(ist)*vrt;
      pst=p+crk(ist)*pt;
      ust=u+crk(ist)*ut;
      vst=v+crk(ist)*vt;
      plst=pl+crk(ist)*plt;
      ulst=ul+crk(ist)*ult;
      vlst=vl+crk(ist)*vlt;
      prst=pr+crk(ist)*prt;
      urst=ur+crk(ist)*urt;
      vrst=vr+crk(ist)*vrt;
    else
      p=pnew+brk(ist)*pt;
      u=unew+brk(ist)*ut;
      v=vnew+brk(ist)*vt;
      pl=plnew+brk(ist)*plt;
      ul=ulnew+brk(ist)*ult;
      vl=vlnew+brk(ist)*vlt;
      pr=prnew+brk(ist)*prt;
      ur=urnew+brk(ist)*urt;
      vr=vrnew+brk(ist)*vrt;
    end
  end
%
  if (mod(it,50)==0)
    jwr=100+it/50;
    cj=num2str(jwr,'%i');
    fname=strcat('duct_nbc7.',cj);
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


