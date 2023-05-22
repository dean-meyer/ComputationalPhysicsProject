% Initialization of stuff and things. Clear everything out.
clear all;

% How many objects do we want.
n_obj=399;
randomvar = 500;
plot_mode=1;

% All objects are assigned a random mass from 0 to 1.
for i=(1:n_obj)
    L = rand(1);
    mass(i) = L;
end

% Create the radii.
rad = sqrt(mass);

% Period and radius.
p_jup=200;
r_jup=15; 

% Time steps.
dt = 0.5;
irt=1;

% It's time to do some preallocation.
x=zeros(1,n_obj);
y=zeros(1,n_obj);
vx=zeros(1,n_obj);
vy=zeros(1,n_obj);
px = zeros(1,n_obj);
py = zeros(1,n_obj);

% Randomly create planetoid periods as a fraction of the above period
% as well as a rotational angle.

irat1 = 3.0*rand(1,1);
irat=3.0*rand(1,n_obj);
rangl=2*pi*rand(1,n_obj);

% The "Kuldge" factor helps us speed the simulation, it is included in the 
% initial velocity calculations. This can result in increased error, but we
% can somewhat offset by modifying our dt value.s
kludj=0.8;

% Create the initial positions, velocities, momentums. Also set prior
% values the same as current values for the purposes of Euler-Cromer.
x(1:n_obj)=-cos(rangl).*(p_jup*irat).^(2/3);
y(1:n_obj)=-sin(rangl).*(p_jup*irat).^(2/3);
vx(1:n_obj)=kludj*2*pi*sin(rangl).*((p_jup*irat).^(-1/3));
vy(1:n_obj)=-kludj*2*pi*cos(rangl).*((p_jup.*irat).^(-1/3));
px(1:n_obj) = vx .* mass;
py(1:n_obj) = vy .* mass;
ox=x;
oy=y;
ovx=vx;
ovy=vy;
opx = px;
opy = py;

% Time steps.
tt=0;
step=1;

% Create the figure and set a figure handle.
figure(1);
g=figure(1);

str='';
slen=length(str);

% Loop to cover the orbits one at a time.
while(randomvar > 15)
    
    % Acceleration due to gravitational force.
    fx=[];
    fy=[];
    for jj=1:length(mass)
        
        fx=0;
        fy=0;
        for ii=1:length(mass)
            if(jj~=ii)
                c_r3=mass(ii)./((ox(ii)-ox(jj)).^2 + ...
                                (oy(ii)-oy(jj)).^2).^(1.5);
                fx=fx-(ox(jj)-ox(ii)).*c_r3;
                fy=fy-(oy(jj)-oy(ii)).*c_r3;
            end
        end
        
        % Update the velocities, positions, and momentums.
        vx(jj)=ovx(jj)+fx*dt;
        vy(jj)=ovy(jj)+fy*dt;
        x(jj)=ox(jj)+vx(jj)*dt;
        y(jj)=oy(jj)+vy(jj)*dt;
        px(jj)=vx(jj)*mass(jj);
        py(jj)=vy(jj)*mass(jj);
        
    end
    
    % Let's plot.
    if(~mod(step,1) & plot_mode)
        figure(1);
        spinor(x(:),y(:),mass(:),zeros(1,numel(mass))); 
        hold off;
        xlim([-100 100]);
        ylim([-100 100]);
        drawnow;
    end
    
    % Advance everything by one step.
    tt=tt+dt;
    step=step+1;
    
    % Pass everything to the Collision Detection funcation for processing
    % and then send it all back out to be used next time through the loop.
    [n_obj,mass,x,y,vx,vy,px,py,rad] = CollisionDetection(n_obj,mass,x,y,vx,vy,px,py,rad);
    ox=x;
    oy=y;
    ovx=vx;
    ovy=vy;
    opx = px;
    opy = py;

    % We have to update our stop condition here so we can track it to end
    % the code execution.
    randomvar = length(mass);
end