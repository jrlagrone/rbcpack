function [pCRBC_t uCRBC_t vCRBC_t pNdat uNdat pEdat vEdat]=acoustics2d_CRBC_cornerP(dpN,duN,dpE,dvE,pCRBC,uCRBC,vCRBC,c,sig,nbc,cmL,cmU,pc)

% Compute time derivatives of the auxiliary fields on a corner element. Also return normal derivative data
% which is used to compute fluxes. 
%
% We write as though this is the northeast corner and the first index corresponds to the east direction and
% the second to the north direction. The inputs are outgoing data along the faces and we assume the
% velocities are into the face. Thus dpN, duN are data from the north face with u normal to the east face
% and dpE, dvE are data from the east face with v normal to the north face. North will always be north or south
% and East will always be east or west. These variables all have dimension nbc+1. 
% pCRBC(nbc+1,nbc+1),uCRBC(nbc+1,nbc+1),vCRBC(nbc+1,nbc+1) are current values of the auxiliary variables
% c speed of sound (p_t+cu_x+cv_y=0, u_t+cp_x=0, v_t+cp_y=0 where x is normal and v tangential)
% sig(2*nbc) is a crbc parameter 
% cmL and cmU are precomputed LU factors of the corner matrix computed by acoustics2d_CRBC_cornermat.
% pc is the column permutation computed in acoustics2d_CRBC_cornermat.

pCRBC_t=zeros(nbc+1,nbc+1);
uCRBC_t=zeros(nbc+1,nbc+1);
vCRBC_t=zeros(nbc+1,nbc+1);
pNdat=zeros(nbc+1,1);
uNdat=zeros(nbc+1,1);
pEdat=zeros(nbc+1,1);
vEdat=zeros(nbc+1,1);
nn=4*(nbc+1)^2;
bb=zeros(nn,1);
xx=zeros(nn,1);
x=zeros(nn,1);

irow=1;

% Incoming data from the North

j=1;

for k=1:nbc+1
  bb(irow)=dpN(k)+duN(k);
  irow=irow+1;
end

% Incoming data from the East

k=1;

for j=1:nbc+1
  bb(irow)=dpE(j)+dvE(j);
  irow=irow+1;
end

% x-recursions c p_jx-c a_j (u_jx+v_jy) +sig_j p_j = -cp_(j+1)x-c abar_j (u_(j+1)x + v_(j+1)y) + sigbar_j p_(j+1)
%              c u_jx -ca_j p_jx + sig_j u_j = -cu_(j+1)x - c abar_j p_(j+1)x + sigbar_j u_(j+1) 

for k=1:nbc+1
for j=1:nbc
  bb(irow)=-sig(2*j-1)*pCRBC(j,k)+sig(2*j)*pCRBC(j+1,k);
  irow=irow+1;
  bb(irow)=-sig(2*j-1)*uCRBC(j,k)+sig(2*j)*uCRBC(j+1,k);
  irow=irow+1;
end
end

% y-recursions c p_ky-c a_k (u_kx+v_ky)+ sig_k p_k = -cp_(k+1)y-c abar_k (u_(k+1)x + v_(k+1)y)+ sigbar_k p_(k+1)
%              c v_ky -ca_k p_ky +sig_k v_k = -cv_(k+1)y - c abar_k p_(k+1)y + sigbar_k p_(k+1)

for k=1:nbc
for j=1:nbc+1
  bb(irow)=-sig(2*k-1)*pCRBC(j,k)+sig(2*k)*pCRBC(j,k+1);
  irow=irow+1;
  bb(irow)=-sig(2*k-1)*vCRBC(j,k)+sig(2*k)*vCRBC(j,k+1);
  irow=irow+1;
end
end

% No data at terminations 

xx=cmU\(cmL\bb);
x(pc)=xx;

% The pdes - indexing as in cornermat

for k=1:nbc+1
for j=1:nbc+1
  ipx=1+4*(j-1)+4*(nbc+1)*(k-1);
  pCRBC_t(j,k)=-c*x(ipx+2)-c*x(ipx+3);
  uCRBC_t(j,k)=-c*x(ipx);
  vCRBC_t(j,k)=-c*x(ipx+1);
end
end

for k=1:nbc+1
  ipx=1+4*(nbc+1)*(k-1);
  pNdat(k)=x(ipx);
  uNdat(k)=x(ipx+2);
end

for j=1:nbc+1
  ipx=1+4*(j-1);
  pEdat(j)=x(ipx+1);
  vEdat(j)=x(ipx+3);
end



