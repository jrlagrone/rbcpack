% solve the Klein-Gordon equation with a double absorbing boundary layer as
% described in "The Double Absorbing Boundary method" by Hagstrom, Givoli,
% et. al.
%
% d^2 u / dt^2 - c^2 d^2 u / dx^2 + su = 0,
% s > 0, c > 0
% u(0) = 0, DAB at x = 1

% domain parameters
n = 500;           % number of grid points in domain [0,1]
m = 600;           % number of grid points to extent the domain by for a reference 
                   % solution using a larger simulation
s = 25;            % dispersion parameter
c = 1;             % wave speed
nsteps = 1500;     % number of time steps
cfl = 0.99;        % cfl ratio 

% compute grid spacing / time step
dx = 1.0 / (n - 1);
dt = cfl * dx / c;

% source paramters
tw = 25*dt;             % pulse width
t0 = 5*tw;              % pulse delay (roughly need t0 > 4*tw to be smooth)
amp = 1;                % pulse "amplitude"
sloc = 90;              % source location

% DAB parameters
p = 3;                  % DAB/CRBC order
ndab = 3;               % DAB width
a = ones(p,1);          % choose all the cosines to be one for simplicity
ab = ones(p,1);         % choose all the cosines to be one for simplicity
sig = zeros(p,1);       % choose all the sigmas to be zero for simplicity
sigb = zeros(p,1);      % choose all the sigmas to be zero for simplicity

% ... or use optimal cosines
T = nsteps * dt;
eta = (n - sloc)*dx / (c * T);
if (p>0)
    [at errest] = optimal_cosines(eta, p-1);
    a = at(1:2:2*p);
    ab = at(2:2:2*p);
    sig = 0.5*dt*(1 - a.*a) ./ (T*a);
    sigb = 0.5*dt*(1 - ab.*ab) ./ (T*ab);
end

% allocate storage
unew = zeros(n,1); % field values
ucur = zeros(n,1);
uold = zeros(n,1);

runew = zeros(n+m,1); % for larger reference simulation 
rucur = zeros(n+m,1);
ruold = zeros(n+m,1);

udabnew = zeros(ndab, p+1); % dab aux. variables
udabcur = zeros(ndab, p+1);
udabold = zeros(ndab, p+1);

% time step
for t=1:nsteps
    
    % internal updates --- eqn 54, in DAB paper
    unew(2:n-1) = 2*ucur(2:n-1) - uold(2:n-1) + ((c*dt)/dx)^2 * (ucur(3:n) ...
        - 2*ucur(2:n-1) + ucur(1:n-2)) - s*dt^2*ucur(2:n-1);
    
    % reference solution
    runew(2:m+n-1) = 2*rucur(2:m+n-1) - ruold(2:m+n-1) + ((c*dt)/dx)^2 * (rucur(3:m+n) ...
        - 2*rucur(2:m+n-1) + rucur(1:m+n-2)) - s*dt^2*rucur(2:m+n-1);
    
    % add source 
    unew(sloc) = unew(sloc) - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);
    runew(sloc) = runew(sloc) - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);
    
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper
    for q=1:p+1
       udabnew(2:ndab-1, q) = 2*udabcur(2:ndab-1, q) - udabold(2:ndab-1, q) + ...
           ((c*dt)/dx)^2 * (udabcur(3:ndab, q)- 2*udabcur(2:ndab-1, q) +...
           udabcur(1:ndab-2, q)) - s*dt^2*udabcur(2:ndab-1, q); 
    end
    
    % left boundary is correctly set to zero, copy data to DAB boundary for
    % the right boundary
    udabnew(1,1) = unew(n-1);
    
    % run the "forward" recursion --- from eqn. 60-61 (a=ab=1,sig=sigb=0)
    w = 1/dt + c/dx;

    % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        udabnew(1,q+1) = ...
            (ab(q) - c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(1,q+1) ...
            +(ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2,q+1) ...
            +(-a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2,q) ...
            +(-a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(1,q) ...
            +(-ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2,q+1) ...
            +(a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2,q) ...
            +(a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(1,q);   
    end
    
    % apply the termination conditon, sommerfeld --- from eqn 56-57
    udabnew(ndab, p+1) = ((udabcur(ndab-1, p+1) - udabnew(ndab-1, p+1) + udabcur(ndab, p+1)) / dt ...
        + c*(udabcur(ndab-1,p+1) - udabcur(ndab, p+1) + udabnew(ndab-1, p+1))/dx)/w;

    % run the "backward" recursions --- from eqn. 58-59, generalized\
    for q=p:-1:1
        udabnew(ndab,q) = ...
            (a(q) - c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(ndab,q) ...
            +(a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(ndab-1,q) ...
            +(-ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(ndab-1,q+1) ...
            +(-ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(ndab,q+1) ...
            +(-a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(ndab-1,q) ...
            +(ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(ndab-1,q+1) ...
            +(ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(ndab,q+1);   
    end
    
    % copy DAB value back into the field
    unew(n) = udabnew(2,1);
    
    % figures
    
    % field and comparison to larger simulation
    figure(1)
    subplot(1,2,1)
    plot(1:n, unew);
    title('field values')
    subplot(1,2,2)
    plot(1:n, unew - runew(1:n))
    title('Error compared to larger simulation')
    drawnow
    
    % field and auxiliary fields
%     figure(2)
%     subplot(1, p+4, 1:3)
%     plot(1:n, unew);
%     title('field values')
%     for i=1:p+1
%         subplot(1, p+4, i+3)
%         plot(1:ndab, udabnew(:,i));
%         title(sprintf('p = %i', i-1))
%     end
%     drawnow
    
    % swap old, new, and current values
    uold = ucur;
    ucur = unew;
    
    ruold = rucur;
    rucur = runew;
    
    udabold = udabcur;
    udabcur = udabnew;
    
end

