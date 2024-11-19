function [DD Lift s]=DGsetup(ndg)
  s=JacobiGL(0,0,ndg);
  V1D=Vandermonde1D(ndg,s); 
  DD=Dmatrix1D(ndg,s,V1D); 
  Emat=zeros(ndg+1,2);
  Emat(1,1)=1;
  Emat(ndg+1,2)=1;
  Lift=V1D*(V1D'*Emat);

