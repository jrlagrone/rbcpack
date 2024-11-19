function [cmL cmU pc]=acoustics2d_CRBC_cornermat(c,a,nbc)

% Form and factor the corner matrix - use colamd to reorder
%
% a(2*nbc) are the crbc parameters and
% c is the sound speed
% We write as though this is the northeast corner and the first index corresponds to the east direction and
% the second to the north direction. The inputs are outgoing normal derivatives along the faces and we assume the
% velocities are into the face.

nn=4*(nbc+1)^2;

% Unknowns are p_x(j,k) p_y(j,k) u_x(j,k) v_y(j,k)
%
% Indices: p_x 1+4(j-1)+4(nbc+1)(k-1)
%          p_y 2+4(j-1)+4(nbc+1)(k-1)
%          u_x 3+4(j-1)+4(nbc+1)(k-1)
%          v_y 4j+4(nbc+1)(k-1)
% For simplicity I will create a full matrix and then sparsify
%

B=zeros(nn,nn);


% Create indices
  p_x=zeros(nbc+1,nbc+1);
  p_y=zeros(nbc+1,nbc+1);
  u_x=zeros(nbc+1,nbc+1);
  v_y=zeros(nbc+1,nbc+1);
  for k=1:nbc+1
  for j=1:nbc+1
    p_x(j,k)=1+4*(j-1)+4*(nbc+1)*(k-1);
    p_y(j,k)=2+4*(j-1)+4*(nbc+1)*(k-1);
    u_x(j,k)=3+4*(j-1)+4*(nbc+1)*(k-1);
    v_y(j,k)=4*j+4*(nbc+1)*(k-1);
  end
  end
%
irow=1;

% Incoming data from the North

j=1;

for k=1:nbc+1
  B(irow,p_x(j,k))=1;
  B(irow,u_x(j,k))=1;
  irow=irow+1;
end

% Incoming data from the East

k=1;

for j=1:nbc+1
  B(irow,p_y(j,k))=1;
  B(irow,v_y(j,k))=1;
  irow=irow+1;
end

% x-recursions c p_jx-c a_j (u_jx+v_jy) = -cp_(j+1)x-c abar_j (u_(j+1)x + v_(j+1)y)
%              c u_jx -ca_j p_jx = -cu_(j+1)x - c abar_j p_(j+1)x

for k=1:nbc+1
for j=1:nbc
  B(irow,p_x(j,k))=c;
  B(irow,p_x(j+1,k))=c;
  B(irow,u_x(j,k))=-c*a(2*j-1);
  B(irow,u_x(j+1,k))=c*a(2*j);
  B(irow,v_y(j,k))=-c*a(2*j-1);
  B(irow,v_y(j+1,k))=c*a(2*j);
  irow=irow+1;
  B(irow,u_x(j,k))=c;
  B(irow,u_x(j+1,k))=c;
  B(irow,p_x(j,k))=-c*a(2*j-1);
  B(irow,p_x(j+1,k))=c*a(2*j);
  irow=irow+1;
end
end

% y-recursions c p_ky-c a_k (u_kx+v_ky) = -cp_(k+1)y-c abar_k (u_(k+1)x + v_(k+1)y)
%              c v_ky -ca_k p_ky = -cv_(k+1)y - c abar_k p_(k+1)y

for k=1:nbc
for j=1:nbc+1
  B(irow,p_y(j,k))=c;
  B(irow,p_y(j,k+1))=c;
  B(irow,u_x(j,k))=-c*a(2*k-1);
  B(irow,u_x(j,k+1))=c*a(2*k);
  B(irow,v_y(j,k))=-c*a(2*k-1);
  B(irow,v_y(j,k+1))=c*a(2*k);
  irow=irow+1;
  B(irow,v_y(j,k))=c;
  B(irow,v_y(j,k+1))=c;
  B(irow,p_y(j,k))=-c*a(2*k-1);
  B(irow,p_y(j,k+1))=c*a(2*k);
  irow=irow+1;
end
end

% Terminations in j recursion

j=nbc+1;

for k=1:nbc+1
  B(irow,p_x(j,k))=1;
  B(irow,u_x(j,k))=-1;
  irow=irow+1;
end

% Terminations in k recursion

k=nbc+1;

for j=1:nbc+1
  B(irow,p_y(j,k))=1;
  B(irow,v_y(j,k))=-1;
  irow=irow+1;
end

% Reality check

%disp(sprintf('%3.0f %3.0f',nn,irow-1))

A=sparse(B);

% Sparsity check

%disp(sprintf('%7.0f %7.0f',nn^2,nnz(A)))

pc=colamd(A);

[cmL cmU]=lu(A(:,pc)); 


