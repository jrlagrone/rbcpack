function [pCRBC_t uCRBC_t vCRBC_t pCRBC_dat uCRBC_dat]=acoustics2d_CRBC_face(pn,un,pCRBCp,vCRBCp,pCRBC,uCRBC,vCRBC,c,Dtan,Lift,alpha,ndg,a,sig,nbc)

% Compute time derivatives of the auxiliary fields on a face element and the normal derivatives which are used for
% volume flux data 
%
% Here p is pressure u is velocity outward normal to the face and v is tangential velocity in increasing index direction
% pn(ndg+1), un(ndg+1) are p and the outward normal velocity from the interior
% pCRBCp(2,nbc+1), vCRBCp(2,nbc+1) are outside states for the face element 
% pCRBC(ndg+1,nbc+1),uCRBC(ndg+1,nbc+1),vCRBC(ndg+1,nbc+1) are current values of the auxiliary variables
% c speed of sound (p_t+cu_x+cv_y=0, u_t+cp_x=0, v_t+cp_y=0 where x is normal and v tangential)
% Dtan(ndg+1,ndg+1) tangential derivative matrix
% Lift(ndg+1,2) lift matrix
% alpha upwinding parameter - see definitions of pflux and vflux 
% a(2*nbc),sig(2*nbc) - CRBC parameters
%

pCRBC_t=zeros(ndg+1,nbc+1);
uCRBC_t=zeros(ndg+1,nbc+1);
vCRBC_t=zeros(ndg+1,nbc+1);
pflux=zeros(2,nbc+1);
vflux=zeros(2,nbc+1); 

% tangential derivatives

ptan=Dtan*pCRBC;
vtan=Dtan*vCRBC;
pflux(1,1:nbc+1)=.5*(pCRBC(1,1:nbc+1)-pCRBCp(1,1:nbc+1))-.5*alpha*(vCRBCp(1,1:nbc+1)-vCRBC(1,1:nbc+1));
pflux(2,1:nbc+1)=-.5*(pCRBC(ndg+1,1:nbc+1)-pCRBCp(2,1:nbc+1))-.5*alpha*(vCRBCp(2,1:nbc+1)-vCRBC(ndg+1,1:nbc+1));
vflux(1,1:nbc+1)=.5*(vCRBC(1,1:nbc+1)-vCRBCp(1,1:nbc+1))-.5*alpha*(pCRBCp(1,1:nbc+1)-pCRBC(1,1:nbc+1));
vflux(2,1:nbc+1)=-.5*(vCRBC(ndg+1,1:nbc+1)-vCRBCp(2,1:nbc+1))-.5*alpha*(pCRBCp(2,1:nbc+1)-pCRBC(ndg+1,1:nbc+1));
ptan=ptan+Lift*pflux;
vtan=vtan+Lift*vflux;

% normal derivatives by recursion

rpn=zeros(ndg+1,nbc+1);
rmn=zeros(ndg+1,nbc+1);

rpn(:,1)=pn+un;  % Volume data 
for j=1:nbc
  rpn(:,j+1)=(c*(a(2*j-1)-1)*rpn(:,j)+c*a(2*j-1)*vtan(:,j)-sig(2*j-1)*(pCRBC(:,j)+uCRBC(:,j)) ...
              -c*a(2*j)*vtan(:,j+1)+sig(2*j)*(pCRBC(:,j+1)+uCRBC(:,j+1)))/(c*(1+a(2*j)));
end
rmn(:,nbc+1)=zeros(ndg+1,1);
for j=nbc:-1:1
  rmn(:,j)=(c*(a(2*j)-1)*rmn(:,j+1)-c*a(2*j)*vtan(:,j+1)+sig(2*j)*(pCRBC(:,j+1)-uCRBC(:,j+1)) ...
              +c*a(2*j-1)*vtan(:,j)-sig(2*j-1)*(pCRBC(:,j)-uCRBC(:,j)))/(c*(1+a(2*j-1)));
end

% Finally the equations

pCRBC_t=-(c/2)*(rpn-rmn)-c*vtan; 
uCRBC_t=-(c/2)*(rpn+rmn);
vCRBC_t=-c*ptan;
pCRBC_dat=(rpn(:,1)+rmn(:,1))/2;
uCRBC_dat=(rpn(:,1)-rmn(:,1))/2;





