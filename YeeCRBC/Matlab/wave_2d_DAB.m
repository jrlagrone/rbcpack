% solve the 2d scalar wave equation with a double absorbing boundary layer 
% as described in "The Double Absorbing Boundary method" by Hagstrom, Givoli,
% et. al.

% domain parameters
n = 300;           % number of grid points in domain [0,0.1]
m = 300;           % number of grid points to extent the domain by for a reference 
                   % solution using a larger simulation
c = 1;             % wave speed
nsteps = 700;     % number of time steps
cfl = 0.99; % cfl ratio (1 is exact dispersion relation, but num. unstable)

% compute grid spacing / time step
dx = 0.1 / (n - 1); % we'll just use dx=dy
dt = cfl / (c*sqrt(2/dx^2));

% source paramters
tw = 25*dt;             % pulse width
t0 = 5*tw;              % pulse delay (roughly need t0 > 4*tw to be smooth)
amp = 1;                % pulse "amplitude"
sloc = [210, 175];      % source location

% DAB parameters
p = 5;                  % DAB/CRBC order
ndab = 3;               % DAB width
a = ones(p,1);          % choose all the cosines to be one for simplicity
ab = ones(p,1);         % choose all the cosines to be one for simplicity
sig = zeros(p,1);       % choose all the sigmas to be zero for simplicity
sigb = zeros(p,1);      % choose all the sigmas to be zero for simplicity

% or choose optimal cosines
T = nsteps * dt;
eta = (n - sloc(2))*dx / (c * T);
if (p>0)
    [at, errest] = optimal_cosines(eta, p-1);
    a = at(1:2:2*p);
    ab = at(2:2:2*p);
    sig = 0.5*dt*(1 - a.*a) ./ (T*a);
    sigb = 0.5*dt*(1 - ab.*ab) ./ (T*ab);
end

% allocate storage
unew = zeros(n); % field values
ucur = zeros(n);
uold = zeros(n);

runew = zeros(n, n+m); % for larger reference simulation 
rucur = zeros(n, n+m);
ruold = zeros(n, n+m);

udabnew = zeros(n, ndab, p+1); % dab aux. variables
udabcur = zeros(n, ndab, p+1);
udabold = zeros(n, ndab, p+1);

% time step
for t=1:nsteps
    
    % internal updates --- eqn 54, in DAB paper
    unew(2:n-1, 2:n-1) = 2*ucur(2:n-1, 2:n-1) - uold(2:n-1, 2:n-1) ...
        + ((c*dt)/dx)^2 * (ucur(3:n, 2:n-1) - 4*ucur(2:n-1,2:n-1) ...
        + ucur(1:n-2, 2:n-1) + ucur(2:n-1, 1:n-2) + ucur(2:n-1, 3:n));
    
    % reference solution (same thing on a bigger domain)
    runew(2:n-1, 2:n+m-1) = 2*rucur(2:n-1, 2:n+m-1) - ruold(2:n-1, 2:n+m-1) ...
        + ((c*dt)/dx)^2 * (rucur(3:n, 2:n+m-1) - 4*rucur(2:n-1,2:n+m-1) ...
        + rucur(1:n-2, 2:n+m-1) + rucur(2:n-1, 1:n+m-2) + rucur(2:n-1, 3:n+m));
    
    % add source
    unew(sloc(1), sloc(2)) = unew(sloc(1), sloc(2)) ...
        - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);
    runew(sloc(1), sloc(2)) = runew(sloc(1), sloc(2)) ...
        - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);
    
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper 
    udabnew(2:n-1, 2:ndab-1,:) = 2*udabcur(2:n-1, 2:ndab-1,:) - udabold(2:n-1, 2:ndab-1,:) ...
        + ((c*dt)/dx)^2 * (udabcur(3:n, 2:ndab-1,:) - 4*udabcur(2:n-1,2:ndab-1,:) ...
        + udabcur(1:n-2, 2:ndab-1,:) + udabcur(2:n-1, 1:ndab-2,:) + udabcur(2:n-1, 3:ndab,:));
  
    
    % copy data to DAB boundary for
    % the right boundary in the y-direction
    udabnew(:,1,1) = unew(:,n-1);
    
    % run the "forward" recursion --- from eqn. 60-61 (a=ab=1,sig=sigb=0)
    w = 1/dt + c/dx;

    % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        udabnew(2:n-1,1,q+1) = ...
            (ab(q) - c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2:n-1,1,q+1) ...
            +(ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2:n-1,2,q+1) ...
            +(-a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2:n-1,2,q) ...
            +(-a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2:n-1,1,q) ...
            +(-ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2:n-1,2,q+1) ...
            +(a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2:n-1,2,q) ...
            +(a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2:n-1,1,q);   
    end
    
    % apply the termination conditon, sommerfeld --- from eqn 56-57
    udabnew(2:n-1,ndab, p+1) = ((udabcur(2:n-1,ndab-1, p+1) - udabnew(2:n-1,ndab-1, p+1) + udabcur(2:n-1,ndab, p+1)) / dt ...
        + c*(udabcur(2:n-1,ndab-1,p+1) - udabcur(2:n-1,ndab, p+1) + udabnew(2:n-1,ndab-1, p+1))/dx)/w;
    
    % run the "backward" recursions --- from eqn. 58-59, generalized
    for q=p:-1:1
        udabnew(2:n-1,ndab,q) = ...
            (a(q) - c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(2:n-1,ndab,q) ...
            +(a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(2:n-1,ndab-1,q) ...
            +(-ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(2:n-1,ndab-1,q+1) ...
            +(-ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(2:n-1,ndab,q+1) ...
            +(-a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(2:n-1,ndab-1,q) ...
            +(ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(2:n-1,ndab-1,q+1) ...
            +(ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(2:n-1,ndab,q+1);   
    end
    
    % copy DAB value back into the field
    unew(:,n) = udabnew(:,2,1);
    
    % figures
    
    % field and comparison to larger simulation
    figure(1)
    subplot(1,2,1)
    surf(unew');
    view(2)
    shading interp
    colorbar;
    title('field values')
    subplot(1,2,2)
    surf(unew' - runew(1:n,1:n)')
    colorbar;
    view(2)
    shading interp
    title('Error compared to larger simulation')
    drawnow
    
    % field and auxiliary fields
%     figure(2)
%     subplot(1, p+4, 1:3)
%     surf(unew);
%     view(2)
%     shading interp
%     colorbar;
%     title('field values')
%     for i=1:p+1
%         subplot(1, p+4, i+3)
%         surf(udabnew(:,:,i));
%         shading interp
%         view(2)
%         colorbar;
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

