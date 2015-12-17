% solve the Klein-Gordon equation with a double absorbing boundary layer as
% described in "The Double Absorbing Boundary method" by Hagstrom, Givoli,
% et. al.
%
% d^2 u / dt^2 - c^2 (d^2 u / dx^2 + d^2 u / dy^2) + su = 0,
% s > 0, c > 0
% u(0,y) = 0, u(x,0) = 0, DAB at x = .1, y = .1

% domain parameters
n = 250;           % number of grid points in domain [0,0.1]
m = 250;           % number of grid points to extent the domain by for a reference
% solution using a larger simulation
c = 1;             % wave speed
s = 20;            % dispersion parameter
nsteps = 1000;     % number of time steps
cfl = 0.99; % cfl ratio (1 is exact dispersion relation, but num. unstable)

% plotting
figure(1)
plot_freq = 10;   % plot frequency (in time steps)

% compute grid spacing / time step
dx = 0.1 / (n - 1); % we'll just use dx=dy
dt = cfl / (c*sqrt(2/dx^2));

% source paramters
tw = 25*dt;             % pulse width
t0 = 5*tw;              % pulse delay (roughly need t0 > 4*tw to be smooth)
amp = 1;                % pulse "amplitude"
sloc = [60, 25];      % source location

% DAB parameters
p = 5;                  % DAB/CRBC order
ndab = 3;               % DAB width

% or choose optimal cosines
T = nsteps * dt;        % simulation length
delta = max(n-sloc)*dx; % Seperation of DAB from source
eta = delta / (c * T);
if (p>0)
    [at, errest] = optimal_cosines(eta, p-1);
    a = at(1:2:2*p);
    ab = at(2:2:2*p);
    sig = 0.5*dt*(1 - a.*a) ./ (T*a);
    sigb = 0.5*dt*(1 - ab.*ab) ./ (T*ab);
end

fprintf('The reflection coefficient is %e \n', errest);

% allocate storage
unew = zeros(n); % field values
ucur = zeros(n);
uold = zeros(n);

runew = zeros(n+m); % for larger reference simulation
rucur = zeros(n+m);
ruold = zeros(n+m);

udabnew_east = zeros(ndab, n, p+1); % east dab aux. variables
udabcur_east = zeros(ndab, n, p+1);
udabold_east = zeros(ndab, n, p+1);

udabnew_north = zeros(n, ndab, p+1); % north dab aux. variables
udabcur_north = zeros(n, ndab, p+1);
udabold_north = zeros(n, ndab, p+1);

udabnew_ne = zeros(ndab, ndab, p+1, p+1); % northeast corner dab aux. variables
udabcur_ne = zeros(ndab, ndab, p+1, p+1);
udabold_ne = zeros(ndab, ndab, p+1, p+1);

% time step
for t=1:nsteps
    
    % internal updates --- eqn 54, in DAB paper
    unew(2:n-1, 2:n-1) = 2*ucur(2:n-1, 2:n-1) - uold(2:n-1, 2:n-1) ...
        + ((c*dt)/dx)^2 * (ucur(3:n, 2:n-1) - 4*ucur(2:n-1,2:n-1) ...
        + ucur(1:n-2, 2:n-1) + ucur(2:n-1, 1:n-2) + ucur(2:n-1, 3:n)) ...
        - s*dt^2*ucur(2:n-1, 2:n-1);
    
    % reference solution (same thing on a bigger domain)
    runew(2:n+m-1, 2:n+m-1) = 2*rucur(2:n+m-1, 2:n+m-1) - ruold(2:n+m-1, 2:n+m-1) ...
        + ((c*dt)/dx)^2 * (rucur(3:n+m, 2:n+m-1) - 4*rucur(2:n+m-1,2:n+m-1) ...
        + rucur(1:n+m-2, 2:n+m-1) + rucur(2:n+m-1, 1:n+m-2) + rucur(2:n+m-1, 3:n+m)) ...
        - s*dt^2*rucur(2:n+m-1, 2:n+m-1);
    
    % add source
    unew(sloc(1), sloc(2)) = unew(sloc(1), sloc(2)) ...
        - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);
    runew(sloc(1), sloc(2)) = runew(sloc(1), sloc(2)) ...
        - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);
    
    %% North DAB
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper
    udabnew_north(2:n-1, 2:ndab-1,:) = 2*udabcur_north(2:n-1, 2:ndab-1,:) - udabold_north(2:n-1, 2:ndab-1,:) ...
        + ((c*dt)/dx)^2 * (udabcur_north(3:n, 2:ndab-1,:) - 4*udabcur_north(2:n-1,2:ndab-1,:) ...
        + udabcur_north(1:n-2, 2:ndab-1,:) + udabcur_north(2:n-1, 1:ndab-2,:) + udabcur_north(2:n-1, 3:ndab,:)) ...
        - s*dt^2*udabcur_north(2:n-1, 2:ndab-1,:);
    
    
    % copy data to DAB boundary for
    % the top boundary in the y-direction
    udabnew_north(:,1,1) = unew(:,n-1);
    
    % run the "forward" recursion --- from eqn. 60-61 (a=ab=1,sig=sigb=0)
    w = 1/dt + c/dx;
    
    % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        udabnew_north(2:n-1,1,q+1) = ...
            (ab(q) - c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_north(2:n-1,1,q+1) ...
            +(ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_north(2:n-1,2,q+1) ...
            +(-a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_north(2:n-1,2,q) ...
            +(-a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_north(2:n-1,1,q) ...
            +(-ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_north(2:n-1,2,q+1) ...
            +(a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_north(2:n-1,2,q) ...
            +(a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_north(2:n-1,1,q);
    end
    
    % apply the termination conditon, sommerfeld --- from eqn 56-57
    udabnew_north(2:n-1,ndab, p+1) = ((udabcur_north(2:n-1,ndab-1, p+1) - udabnew_north(2:n-1,ndab-1, p+1)...
        + udabcur_north(2:n-1,ndab, p+1)) / dt + c*(udabcur_north(2:n-1,ndab-1,p+1) ...
        - udabcur_north(2:n-1,ndab, p+1) + udabnew_north(2:n-1,ndab-1, p+1))/dx)/w;
    
    % run the "backward" recursions --- from eqn. 58-59, generalized
    for q=p:-1:1
        udabnew_north(2:n-1,ndab,q) = ...
            (a(q) - c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_north(2:n-1,ndab,q) ...
            +(a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_north(2:n-1,ndab-1,q) ...
            +(-ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_north(2:n-1,ndab-1,q+1) ...
            +(-ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_north(2:n-1,ndab,q+1) ...
            +(-a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_north(2:n-1,ndab-1,q) ...
            +(ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_north(2:n-1,ndab-1,q+1) ...
            +(ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_north(2:n-1,ndab,q+1);
    end
    
    %% East DAB
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper
    udabnew_east(2:ndab-1, 2:n-1,:) = 2*udabcur_east(2:ndab-1, 2:n-1,:) - udabold_east(2:ndab-1, 2:n-1, :) ...
        + ((c*dt)/dx)^2 * (udabcur_east(3:ndab, 2:n-1, :) - 4*udabcur_east(2:ndab-1, 2:n-1,:) ...
        + udabcur_east(1:ndab-2, 2:n-1,:) + udabcur_east(2:ndab-1, 1:n-2,:) + udabcur_east(2:ndab-1, 3:n,:))...
        - s*dt^2*udabcur_east(2:ndab-1, 2:n-1,:);
    
    % copy data to DAB boundary for
    % the right boundary in the x-direction
    udabnew_east(1,:,1) = unew(n-1,:);
    
    % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        udabnew_east(1,2:n-1,q+1) = ...
            (ab(q) - c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_east(1,2:n-1,q+1) ...
            +(ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_east(2,2:n-1,q+1) ...
            +(-a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_east(2,2:n-1,q) ...
            +(-a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_east(1,2:n-1,q) ...
            +(-ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_east(2,2:n-1,q+1) ...
            +(a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_east(2,2:n-1,q) ...
            +(a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_east(1,2:n-1,q);
    end
    
    % apply the termination conditon, sommerfeld --- from eqn 56-57
    udabnew_east(ndab,2:n-1, p+1) = ((udabcur_east(ndab-1,2:n-1, p+1) - udabnew_east(ndab-1,2:n-1, p+1)...
        + udabcur_east(ndab,2:n-1, p+1)) / dt + c*(udabcur_east(ndab-1,2:n-1,p+1) ...
        - udabcur_east(ndab,2:n-1, p+1) + udabnew_east(ndab-1,2:n-1, p+1))/dx)/w;
    
    % run the "backward" recursions --- from eqn. 58-59, generalized
    for q=p:-1:1
        udabnew_east(ndab,2:n-1,q) = ...
            (a(q) - c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_east(ndab,2:n-1,q) ...
            +(a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_east(ndab-1,2:n-1,q) ...
            +(-ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_east(ndab-1,2:n-1,q+1) ...
            +(-ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_east(ndab,2:n-1,q+1) ...
            +(-a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_east(ndab-1,2:n-1,q) ...
            +(ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_east(ndab-1,2:n-1,q+1) ...
            +(ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_east(ndab,2:n-1,q+1);
    end
    
    %% Northeast Corner DAB
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper
    udabnew_ne(2:ndab-1, 2:ndab-1,:,:) = 2*udabcur_ne(2:ndab-1, 2:ndab-1,:,:) - udabold_ne(2:ndab-1, 2:ndab-1, :,:) ...
        + ((c*dt)/dx)^2 * (udabcur_ne(3:ndab, 2:ndab-1, :,:) - 4*udabcur_ne(2:ndab-1, 2:ndab-1,:,:) ...
        + udabcur_ne(1:ndab-2, 2:ndab-1,:,:) + udabcur_ne(2:ndab-1, 1:ndab-2,:,:) + udabcur_ne(2:ndab-1, 3:ndab,:,:))...
        - s*dt^2*udabcur_ne(2:ndab-1, 2:ndab-1,:,:);
    
    
    % copy data to DAB boundary from the North and East DAB layers
    udabnew_ne(2:ndab-1, 1,:,1) = udabnew_east(2:ndab-1, n-1, :); % first aux index corresponds to x-direction
    udabnew_ne(1, 2:ndab-1,1,:) = udabnew_north(n-1, 2:ndab-1, :); % second aux index corresponds to y-direction
    
    % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        
        % recursions in the x-direction
        udabnew_ne(1,2:ndab-1,q+1,:) = ...
            (ab(q) - c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_ne(1,2:ndab-1,q+1,:) ...
            +(ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_ne(2,2:ndab-1,q+1,:) ...
            +(-a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_ne(2,2:ndab-1,q,:) ...
            +(-a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_ne(1,2:ndab-1,q,:) ...
            +(-ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_ne(2,2:ndab-1,q+1,:) ...
            +(a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_ne(2,2:ndab-1,q,:) ...
            +(a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_ne(1,2:ndab-1,q,:);
        
        % recursions in the y-direction
        udabnew_ne(2:ndab-1,1,:,q+1) = ...
            (ab(q) - c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_ne(2:ndab-1,1,:,q+1) ...
            +(ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_ne(2:ndab-1,2,:,q+1) ...
            +(-a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_ne(2:ndab-1,2,:,q) ...
            +(-a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur_ne(2:ndab-1,1,:,q) ...
            +(-ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_ne(2:ndab-1,2,:,q+1) ...
            +(a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_ne(2:ndab-1,2,:,q) ...
            +(a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew_ne(2:ndab-1,1,:,q);
    end
    
    % apply the termination conditon, sommerfeld --- from eqn 56-57
    % x-direction
    udabnew_ne(ndab,2:ndab-1,p+1,:) = ((udabcur_ne(ndab-1,2:ndab-1,p+1,:) - udabnew_ne(ndab-1,2:ndab-1,p+1,:)...
        + udabcur_ne(ndab,2:ndab-1,p+1,:)) / dt + c*(udabcur_ne(ndab-1,2:ndab-1,p+1,:) ...
        - udabcur_ne(ndab,2:ndab-1,p+1,:) + udabnew_ne(ndab-1,2:ndab-1,p+1,:))/dx)/w;
    
    % y-direction
    udabnew_ne(2:ndab-1,ndab,:,p+1) = ((udabcur_ne(2:ndab-1,ndab-1,:,p+1) - udabnew_ne(2:ndab-1,ndab-1,:,p+1)...
        + udabcur_ne(2:ndab-1,ndab,:,p+1)) / dt + c*(udabcur_ne(2:ndab-1,ndab,:,p+1) ...
        - udabcur_ne(2:ndab-1,ndab,:,p+1) + udabnew_ne(2:ndab-1,ndab-1,:,p+1))/dx)/w;
    
    % run the "backward" recursions --- from eqn. 58-59, generalized
    for q=p:-1:1
        
        % recursions in the x-direction
        udabnew_ne(ndab,2:ndab-1,q,:) = ...
            (a(q) - c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_ne(ndab,2:ndab-1,q,:) ...
            +(a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_ne(ndab-1,2:ndab-1,q,:) ...
            +(-ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_ne(ndab-1,2:ndab-1,q+1,:) ...
            +(-ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_ne(ndab,2:ndab-1,q+1,:) ...
            +(-a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_ne(ndab-1,2:ndab-1,q,:) ...
            +(ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_ne(ndab-1,2:ndab-1,q+1,:) ...
            +(ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_ne(ndab,2:ndab-1,q+1,:);
        
        % recursions in the y-direction
        udabnew_ne(2:ndab-1,ndab,:,q) = ...
            (a(q) - c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_ne(2:ndab-1,ndab,:,q) ...
            +(a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_ne(2:ndab-1,ndab-1,:,q) ...
            +(-ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_ne(2:ndab-1,ndab-1,:,q+1) ...
            +(-ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur_ne(2:ndab-1,ndab,:,q+1) ...
            +(-a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_ne(2:ndab-1,ndab-1,:,q) ...
            +(ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_ne(2:ndab-1,ndab-1,:,q+1) ...
            +(ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew_ne(2:ndab-1,ndab,:,q+1);
    end
    
    % copy values from the Northest DAB corner to the North and East DAB
    % layers
    udabnew_east(2:ndab-1, n, :) = udabnew_ne(2:ndab-1, 2,:,1);
    udabnew_north(n, 2:ndab-1, :) = udabnew_ne(2, 2:ndab-1,1,:);
    
    % copy values from the DAB layers into the internal solver
    unew(n,:) = udabnew_east(2,:,1);
    unew(:,n) = udabnew_north(:,2,1);
    
    % figures
    
    % field and comparison to larger simulation
    if (mod(t,plot_freq) == 0)
        subplot(1,2,1)
        surf(unew'); % transpose because matlab plots x and y "backwards"
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
    end
    
    % swap old, new, and current values
    uold = ucur;
    ucur = unew;
    
    ruold = rucur;
    rucur = runew;
    
    udabold_east = udabcur_east;
    udabcur_east = udabnew_east;
    
    udabold_north = udabcur_north;
    udabcur_north = udabnew_north;
    
    udabold_ne = udabcur_ne;
    udabcur_ne = udabnew_ne;
    
end