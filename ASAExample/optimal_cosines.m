function [a p errest]=optimal_cosines(eta,tol,pmax)

%
% Uses the Remez algorithm to find optimal
% cosines {a(j), j=1 ... 2p} which minimize
% the maximum norm on 0 < x < 1 of:
%
% e(x) = exp(-eta/x)*((1-x)/(1+x))*prod_{j=1}^{2p} (a(j)-x)/(a(j)+x)   
%
% Returns: {a(j) j=1, ... ,2p}, and the smallest p such that the tolerance
% is met - p>pmax not allowed 
%
% The relevance of these numbers to the construction of optimal local
% radiation boundary conditions for the wave equation and related
% first order systems is discussed in ``Complete radiation boundary
% conditions: minimizing the long-time error growth of local methods''
% by Tom Hagstrom and Tim Warburton. A preprint may be downloaded from
% this site.
%
% Written by T. Hagstrom, February 18, 2011.
%
% This code is freely distributed for research purposes. We provide no
% warranties of its performance or suitability for any purpose.
%

global aloc etaloc q

ntol=10^(-10); % Tolerance for fzero when computing local maxima
rtol=10^(-12); % Tolerance for Remez algorithm

% Check for reasonable inputs

if eta > 0.1
  warning(' eta=0.1 is the maximum allowed: results returned for eta=0.1 ')
  etaloc=1;
elseif eta < 10^(-7)   warning(' eta=10^(-7) is the minimum allowed: results returned for eta=10^(-7) ')
  etaloc=10^(-7);
else
  etaloc=eta;
end 

if tol > 0.1
  warning(' tol=.1 is the maximum allowed: results returned for tol=.1 ')
  ttol=0.1;
elseif tol < 10^(-8) 
  warning(' tol=10^(-8) is the minimum allowed: results returned for tol=10^(-8) ')
  ttol=10^(-8);
else
  ttol=tol;
end 


if pmax > 40
  warning(' pmax=40 is the maximum allowed: results returned for pmax=40 ')
  ptmax=40;
elseif pmax < 2 
  warning(' pmax=2 is the minimum allowed: results returned for pmax=2 ')
  ptmax=2;
else
  ptmax=ceil(pmax);
end 

% Start Remez - begin with p=2 and increment until done

p=1;
err=2*ttol;

while (err > ttol) && (p < ptmax)

p=p+1;
[a err]=optimal_cosinesP(eta,p-1);

end 

if (p==ptmax) && (err > ttol)
  warning(' Tolerance unmet - returning optimal cosines with p=pmax')
  errest=err;
end 

errest=err;

return;

function [a emax]=optimal_cosinesP(eta,p)

%
% Uses the Remez algorithm to find optimal
% cosines {a(j), j=1 ... 2p+2} which minimize
% the maximum norm on 0 < x < 1 of:
%
% e(x) = exp(-eta/x)*((1-x)/(1+x))*prod_{j=1}^{2p+2} (a(j)-x)/(a(j)+x)   
%
% Returns: {a(j) j=1, ... ,2p+2}, and max_{0<x<1} |e(x)|
%
% The relevance of these numbers to the construction of optimal local
% radiation boundary conditions for the wave equation and related
% first order systems is discussed in ``Complete radiation boundary
% conditions: minimizing the long-time error growth of local methods''
% by Tom Hagstrom and Tim Warburton. A preprint may be downloaded from
% this site.
%
% Written by T. Hagstrom, January 1, 2009.
%
% This code is freely distributed for research purposes. We provide no
% warranties of its performance or suitability for any purpose.
%

global aloc etaloc q

ntol=10^(-10); % Tolerance for fzero when computing local maxima
rtol=10^(-12); % Tolerance for Remez algorithm 

etaloc=eta; 

xmin=etaloc/log(1/rtol); % For x<xmin |e|<rtol 

q=2*p+2;
aloc=(etaloc).^([1:q]'/q);  % Initial guess - geometric distribution
zopt=optimset('TolX',ntol);
remmax=100; % Maximum number of Remez iterations
it=1;
dn=1;

while (it < remmax) && (dn > ntol)
    ap=[1; aloc; xmin];
    rm=zeros(2*p+3,1);
    sgn=1;
    sgns=rm;
    xm=rm;
    for k=1:2*p+3
        xm(k)=fzero(@dref,ap(k:k+1),zopt);
        sgns(k)=sgn;
        sgn=-sgn;
    end
    rm=ref(xm);
    refmax=max(abs(rm));
    itn=1;
    dan=1;
    a0=aloc;
    lam=min(abs(rm));
    while (itn < 20) && (dan > rtol)
        jac=zeros(2*p+3,2*p+3);
        dz=zeros(2*p+3,1);
        rm=ref(xm);
        res=lam*sgns-rm;
        for k=1:q
            jac(:,k)=-2*rm.*xm./(xm.*xm-aloc(k)*aloc(k));
        end
        jac(:,2*p+3)=-sgns;
        dz=jac\res;
        
        sm=1;
        aap=[1; aloc; xmin];
        ahp=.5*(aap(1:q+1)+aap(2:q+2));
        for k=1:q
            atry=aloc(k)+sm*dz(k);
            if (atry > ahp(k))
                sm=(ahp(k)-aloc(k))/dz(k);
            end
            if (atry < ahp(k+1))
                sm=(ahp(k+1)-aloc(k))/dz(k);
            end
        end
        if ((lam+sm*dz(q+1)) < lam/2)
            sm=-lam/(2*dz(q+1));
        end
        aloc=aloc+sm*dz(1:q);
        lam=lam+sm*dz(q+1);
        dan=sm*max(abs(dz));
        itn=itn+1;
    end
    aloc=sort(aloc,'descend');
    dn=max(abs(aloc-a0));
    it=it+1;
end
   
a=aloc;

ap=[1; aloc; xmin];
for k=1:2*p+3
    xm(k)=fzero(@dref,ap(k:k+1),zopt);
end
rm=ref(xm);
emax=max(abs(rm));

return;

function y=ref(x)

global aloc etaloc q

nx=length(x);
y=zeros(size(x));

for j=1:nx
    if (x(j)==0)
        y(j)=0;
    else
        y(j)=exp(-etaloc/x(j))*(1-x(j))/(1+x(j));
        for k=1:q
            y(j)=y(j)*(x(j)-aloc(k))/(x(j)+aloc(k));
        end
    end
end

return;

function y=dref(x)

global aloc etaloc q

nx=length(x);
y=zeros(size(x));

for j=1:nx
    [dmin,kmin]=min(abs(aloc-x(j)));
    if (x(j)==0)
        y(j)=-10^(-50);
    elseif (x(j)==1)
        y(j)=-exp(-etaloc)/2;
        for k=1:q
            y(j)=y(j)*(1-aloc(k))/(1+aloc(k));
        end
    elseif (dmin==0)
        y(j)=exp(-etaloc/x(j))*(1-x(j))/((1+x(j))*(x(j)+aloc(kmin)));
        for k=1:q
            if (k ~= kmin)
                y(j)=y(j)*(x(j)-aloc(k))/(x(j)+aloc(k));
            end
        end
    else
        w=exp(-etaloc/x(j))*(1-x(j))/(1+x(j));
        dw=w*(etaloc/(x(j)*x(j))-2/(1-x(j)*x(j)));
        for k=1:q
            dw=dw*(x(j)-aloc(k))/(x(j)+aloc(k));
            w=w*(x(j)-aloc(k))/(x(j)+aloc(k));
        end
        for k=1:q
            dw=dw+w*2*aloc(k)/(x(j)^2-aloc(k)^2);
        end
        y(j)=dw;
    end
end

return;

