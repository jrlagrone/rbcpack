function [p_t u_t v_t]=acoustics2d(p,u,v,ppx,upx,ppy,vpy,c,Dx,Liftx,Dy,Lifty,alpha,ndg)

% Compute time derivatives of the auxiliary fields on a face element 
%
% Here p is pressure u is x-velocity and v is y-velocity
% ppx(2,ndg+1),ppy(2,ndg+1),upx(2,ndg+1),vpy(2,ndg+1) are outside states
% c speed of sound (p_t+cu_x+cv_y=0, u_t+cp_x=0, v_t+cp_y=0 where x is normal and v tangential)
% Dx(ndg+1,ndg+1), Dy(ndg+1,ndg+1) are derivative matrices
% Liftx(ndg+1,2),Lifty(ndg+1,2) derivative matrices
% alpha upwinding parameter - see definitions of pflux, uflux, vflux 
%
p_t=zeros(ndg+1,ndg+1);
u_t=zeros(ndg+1,ndg+1);
v_t=zeros(ndg+1,ndg+1);
pflux=zeros(2,ndg+1);
uflux=zeros(2,ndg+1); 
vflux=zeros(2,ndg+1); 

% x derivatives

px=Dx*p;
ux=Dx*u;
pflux(1,1:ndg+1)=.5*(p(1,1:ndg+1)-ppx(1,1:ndg+1))-.5*alpha*(upx(1,1:ndg+1)-u(1,1:ndg+1));
pflux(2,1:ndg+1)=-.5*(p(ndg+1,1:ndg+1)-ppx(2,1:ndg+1))-.5*alpha*(upx(2,1:ndg+1)-u(ndg+1,1:ndg+1));
uflux(1,1:ndg+1)=.5*(u(1,1:ndg+1)-upx(1,1:ndg+1))-.5*alpha*(ppx(1,1:ndg+1)-p(1,1:ndg+1));
uflux(2,1:ndg+1)=-.5*(u(ndg+1,1:ndg+1)-upx(2,1:ndg+1))-.5*alpha*(ppx(2,1:ndg+1)-p(ndg+1,1:ndg+1));
px=px+Liftx*pflux;
ux=ux+Liftx*uflux;

% y derivatives

py=Dy*p';
vy=Dy*v';
pflux(1,1:ndg+1)=.5*(p(1:ndg+1,1)'-ppy(1,1:ndg+1))-.5*alpha*(vpy(1,1:ndg+1)-v(1:ndg+1,1)');
pflux(2,1:ndg+1)=-.5*(p(1:ndg+1,ndg+1)'-ppy(2,1:ndg+1))-.5*alpha*(vpy(2,1:ndg+1)-v(1:ndg+1,ndg+1)');
vflux(1,1:ndg+1)=.5*(v(1:ndg+1,1)'-vpy(1,1:ndg+1))-.5*alpha*(ppy(1,1:ndg+1)-p(1:ndg+1,1)');
vflux(2,1:ndg+1)=-.5*(v(1:ndg+1,ndg+1)'-vpy(2,1:ndg+1))-.5*alpha*(ppy(2,1:ndg+1)-p(1:ndg+1,ndg+1)');
py=py+Lifty*pflux;
vy=vy+Lifty*vflux;

p_t=-c*(ux+vy');
u_t=-c*px; 
v_t=-c*py';


