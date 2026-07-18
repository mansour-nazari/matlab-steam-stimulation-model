% % General assumptions
% % In order to illustrate the conditions under which this new model
% % could be used, some basic assumptions are made below:
% % (a) the reservoir is homogeneous, and fluid flows radially into the oil layer;
% % (b) the thickness of the oil layer is large enough, so that the heat loss to surrounding layers can be neglected;
% % (c) the rock achieves the fluid temperature instantaneously;
% % (d) the thermal energy is fully absorbed by the oil layer and energy conservation law is functional throughout the CSS process;
% % (e) the reservoir temperature is constant in steam chamber, decreasing in hot liquid zone, and equal to the initial reservoir temperature in cold oil zone;
% % (f) steam heats the formation by means of heat convection and heat conduction;
% % (g) the fluctuations of fluid properties and flow parameters along the horizontal section like temperature, pressure and steam quality are neglected
clc
clear
close all
% %% input data
Tr=50; % initial reservoir temperature(C)
Pr=9.25*10^6; % reservoir pressure(pa)
K=1.4*10^-12; % absolute permeability(m^2)
Sw=0.35;
So=0.65;
PHI=30;
HCw=4.2*10^6; %heat capacity of water(J/(m^3*C))
HCo=2.1*10^6; %heat capacity of oil(J/(m^3*C))
HCr=2.34*10^6; %heat capacity of reservoir(J/(m^3*C))
HEATC=2.2; %heat conductivity of reservoir(W/(m*C))
RHOo=0.95*10^3; % OIL DENSITY(KG/m^3)
Lt=300; % length of horizontal section(m)
Rw=0.1; % radius of wellbore(m)
Ps=10.5*10^6; % steam pressure(pa)
Ts=300; % steam temprature(C)
X=0.7; % steam quality
LATENTs=1.4*10^6; % latanet heat of steam(J/KG)
ENTHALs=1.34*10^6; % enthalpy of steam(J/KG)
INJs=16*1000/3600; % steam injection rate(kg/s)
INJp=5*86400; % injection period(s)
SOAK=5*86400; % soaking period(s)
prodt=120*86400; % production period(s)
Hws=2.34*10^6;
Hwr=0.4*10^6;
%%
HCl=HCw*Sw+HCo*So;
Rs=sqrt(INJs*X*LATENTs*INJp/(3.14*HCr*(Ts-Tr)*Lt)+Rw^2); % steam chamber radius
Qw=INJs/998; % condensate water rate(rho water = 998 kg/m^3)
tD=HEATC*INJp/(HCr*Rw^2);
M=INJs*HCw/(2*3.14*HEATC*998*Lt); alpha=HEATC/HCr;
Dr=0.7;
tolerance=10;
while tolerance>=10
    Rh=Rs+Dr;
    Dr=Dr+0.05;
    %%
    %Diffusion equation - 1D - Explicit Method - Cylindrical / Polar 
    % Coordinates
    % Dirichlet BC conditions - Constant temperatures at boundaries
    % Inputs
    Rs=2;
    M=3.3;
    alpha=1.38*10^-6;
    C = -M; % Thermal diffusivity
%     tD=70;
    % boundary conditions
    Tr_in= 1; % BC1, Temperature at r_in, deg.C
    Tr_out=0; % BC2, Temperature at r_out, deg.C
    Ti = 0; % IC, Initial Temperature at r_in <r<r_out, deg.C 
    
    % grid in r direction
 
    rin = Rs/Rw; % Inside Radius , r_in, say m
    rout = Rh/Rw; % Outside Radius , r_out, say m
    mr = 12; % no. of sections divided between r_in and r_out eg 5, 10
    delta_r= (rout - rin)/mr; % section length, m
    nr=mr+1; % total no. of radial points
    % phi linspace(0, 2*pi, 500);
    % total angle = 2*pi
    
    % grid in theta direction
    ntheta = 72; % no. of angle steps
    dtheta=360/ntheta;
    delta_phi = 2*pi/ntheta; % angle step, rad 
    
    % time step discretizing
    time=tD; % total dimesionless time, s eg 200, 9000
    nt=400; % no. of time steps eg 2, 300
    delta_t = time/nt; % timestep, s


    % Solution
    %{
    delta_phi phi(1) = 0; for j = 1:nphi+1
    phi(j) = phi(1)+(-1)*delta_phi;
    end phi
    %}
     
    d = delta_t/delta_r^2; % diffusion number
    d1 = C*delta_t/(2*delta_r);
    
    if d <0.5
        fprintf('solution stable\nd = %10.7f\n', d)
    else
        fprintf('solution unstable\nd = %10.7f\n', d)
    end

    for i = 1:nr
        r(i)=rin +(i-1)* delta_r;
    end
    
    theta_degree=0:dtheta:360+dtheta; 
    theta=theta_degree.*pi/180;

    % Creating initial and boundary conditions T= zeros(nr,nphi+2,nt);
    for k = 1:nt+1
        for j = 1:ntheta+2 
            for i = 1:nr
                if (i == 1) 
                    T(i,j,k) = Tr_in;
                elseif (i == nr)
                    T(i,j,k) = Tr_out;
                else
                    T(i,j,k) = Ti;
                end
            end
        end
    end
    T;
    % Creating T matrix using explicit method
    for k = 1:nt
        for j = 1:ntheta+2
            for i = 2:nr-1
                T(i,j,k+1)= T(i,j,k) + d*(T(i+1,j,k) - 2* T(i,j,k) + T(i-1,j,k)) + (d1/r(i)) * (T(i+1,j,k) - T(i-1,j,k));
            end
        end
    end
    
    
R=r'.*ones(mr+1,ntheta+2);
T1=R.*T(:,:,end);
R=R(:,1);
T1=T1(:,1);
Int=0;
% trapozied rule for integration
    for i=1:mr+1
        if i==1 && i==mr+1
           Int=T1(i)+Int;
        else
            Int=2*T1(i)+Int;
        end
    end

Int=0.5*delta_r*Int;
tolerance=abs((INJs*tD*(Hws-Hwr)/alpha)/(2*3.14*Lt*HCr*(Ts-Tr))-Int);


end
% Convert Polar/Cylindrical Coordinates to Cartesian Coordinates
% Initial Temperature Plot
T1 = T(:,:,1); 
subplot(2,2,1)
r = r;
[Phi,R] = meshgrid(theta,r);
Z= T1; % T1* Phi./Phi 
[X,Y,Z]=pol2cart(Phi,R,Z);
% figure(1)
[C,h] = contourf(X, Y, Z, 300);
set(h,'LineColor', 'none')
grid on
axis equal
colorbar
title({'Initial Temperature Profile - 1D - Cylindrical Coordinates'; 'tD=0'})
xlabel('rD ')

% Final Temperature Plot 
T2 = T(:,:,nt+1); 
subplot(2,2,3)
r=r;
[Phi,R] = meshgrid(theta,r);
Z=T2.* Phi./Phi;
[X,Y,Z]=pol2cart(Phi,R,Z);
% figure(1)
[C,h] = contourf(X, Y, Z, 300);
set(h,'LineColor', 'none')
grid on
axis equal
colorbar
title({'Final Temperature Profile - 1D - Cylindrical Coordinates '; 'tD=40 '})
xlabel('rD ')

%last time step
subplot(2,2,4)
r;
theta;
Z = Z';
[R,PHI] = meshgrid(r,theta+1);
h = surf(R.*cos(PHI), R.*sin(PHI), Z);
    shading interp
    colormap('jet')
    daspect([60 60 1])
    colorbar;
    title('at final time step tD=40')
    xlabel('rD')
    ylabel('rD')
    zlabel('TD')
    xlim([-60 60])
    ylim([-60 60])
% Surface Plot - Animated
subplot(2,2,2)

r;
theta;
[R,PHI] = meshgrid(r,theta+1);

Z=T(:,:,1);
for k = 1:nt+1
    Z = Z';
    h = surf(R.*cos(PHI), R.*sin(PHI), Z);
    shading interp
    colormap('jet')
    daspect([60 60 1])
    colorbar;
    title({['Transient Heat Conduction']; ['tD (\itt) = ',num2str((k-1)*delta_t)]})
    drawnow; 
    pause(0.0001); 
    refreshdata(h);
    if  k~=nt+1
        Z = T(:,:,k+1);
    else
        break;
    end
end

figure
plt=T(:,:,end);
plt=plt(:,end);
plt=(Ts-Tr).*plt+Tr;
pltx=Rw.*R;
plot(pltx,plt,'marker','P')
grid on
title('temprature versus radius at final time step')
xlabel('r(m)')
ylabel('T(C)')