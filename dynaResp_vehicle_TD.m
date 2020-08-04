function [Do,Ao,Aoz_car] = dynaResp_vehicle_TD(Bridge,Wind,varargin)
% function [Do] = dynaResp_vehicle_TD(Bridge,Wind,nCoeff,varargin) Computes the PSD
% of the bridge displacement to wind and or traffic loading. 
%
% INPUT
%
% Bridge: structure variables that contains the information about the bridge
% Wind : structure variables that contains the information about the wind turbulence
%  varargin:
%   - rho: air density (default value is 1.25)
%   - g: gravitational constant (9.81 m∕s^2)
%   - k: term introducing an aerodynamic damping to the torsional deck motion.
%     default value is 0 (0.25 for flat plate aerodynamic)
%   - 'vehicle': structure including the vehicle aprameters to study
%   traffic-induced vibrations (cf. example1.mlx)
%   - 'comments': "on" or "off" to have some comments on possible  issues
% 
% 
% Output
%
% Do: [3xMxN] matrix: Bridge displacement response, where M is the number
%   of nodes and N is the number of time steps. "3" is the number of degrees
%   of freedom considered: lateral, vertial and torsional.
% Ao: [3xMxN] matrix: Bridge acceleration response
% Aoz_car: [PxN] matrix: vehicle acceleration response considering P
%   vehicles
%
% Author: E Cheynet - UiB - last modified: 03/08/2020
%

%% Input parser
p = inputParser();
p.CaseSensitive = false;
p.addOptional('rho',1.25);
p.addOptional('k',1/4);
p.addOptional('g',9.81);
p.addOptional('vehicle',[]);
p.addOptional('comments','on');
p.parse(varargin{:});
%% shorthen the variables name
rho  = p.Results.rho ;
k  = p.Results.k ;
g  = p.Results.g ;
vehicle = p.Results.vehicle ;
comments = p.Results.comments ;
% shortened variable
D = Bridge.D;
B = Bridge.B;
Cd = Bridge.Cd;
Cl = Bridge.Cl;
Cm = Bridge.Cm;
dCd = Bridge.dCd;
dCl = Bridge.dCl;
dCm =Bridge.dCm;
phi = Bridge.phi;
wn = Bridge.wn;
zetaStruct = Bridge.zetaStruct;
u = Wind.u;
w = Wind.w;
t = Wind.t;
%% Case of a zero mean wind speed
if strcmpi(comments,'on')
    if abs(nanmean(u(1,:)))<0.1
        warning('u should have a non-zero mean value !')
    end
end
%% Fundamental aprameters definition
dt = median(diff(t));
N = numel(t);
% Definition of Nyy, and Nmodes:
[~,Nmodes,Nyy]= size(Bridge.phi);
% Check bridge normalization
if max(Bridge.x)== 1
    y = Bridge.x.*Bridge.L;
else
    warning('Bridge.x is not normalized');
end
%%  Check the vehicle parameters
Ncar = numel(vehicle);
if ~isempty(vehicle)
    if strcmpi(comments,'on'), fprintf('movingLoad model used \n');end
end

%% Bridge modal mass, stifness and damping calculation
%  we assume the vehicle mass to be negligible compared to the mass of the
%  bridge
Mtot = diag([Bridge.m+2*Bridge.mc,Bridge.m+2*Bridge.mc,Bridge.m_theta]);
phi0 = reshape(phi,[],Nyy);
phi0_N = bsxfun(@times,phi0,1./max(abs(phi0),[],2));
Nm = size(phi0,1);

Mtot = repmat(Mtot,round(Nm/3),round(Nm/3));
M = zeros(Nm+Ncar);
K = zeros(Nm+Ncar);
C = zeros(Nm+Ncar);
for pp=1:Nm
    for qq=1:Nm/3
        if (pp+3*qq)<=Nm
            Mtot(pp,pp+3*qq) = 0;
            Mtot(pp+3*qq,pp) = 0;
        end
    end
    for qq=1:Nm
        M(pp,qq)  = trapz(y,phi0_N(pp,:).*phi0_N(qq,:).*Mtot(pp,qq));
    end
end
K(1:Nm,1:Nm) = diag(wn(:)).^2*M(1:Nm,1:Nm);
C(1:Nm,1:Nm) = 2.*diag(wn(:))*M(1:Nm,1:Nm)*diag(zetaStruct(:));

%% aerodynamic properties
CoeffAero = @(C,alpha) C(1)+alpha.*C(2); % linear
C1=[Cd,dCd];
C2=[Cl,dCl];
C3=[Cm,dCm];
%% INITIALISATION
Do = zeros(3,Nyy,N);
Ao = zeros(3,Nyy,N);
Aoz_car = zeros(Ncar,N);
rt=zeros(Nyy,1);
dry=zeros(Nyy,1);
drz=zeros(Nyy,1);
drt=zeros(Nyy,1);

% Case of no wind
if mean(nanstd(u,0,2))<0.01 && mean(nanstd(w,0,2))<0.01 && nanmean(nanmean(u,2))<0.1
    noWind=1;
else
    noWind=0;
end

%% Main loop
for idt=1:N
    % get  modal forces for wind
    if noWind==0 % if wind turbulence is acounted for
        [Fmodal_wind(1:Nm,1:Nm)]= getFmodal(CoeffAero,C1,C2,C3,y,...
            u(:,idt),w(:,idt),dry,drz,drt,D,B,rt,k,rho,Nmodes,phi0_N);
        if Ncar>0
            Fmodal_wind(Nm+1:end,Nm+1:end)=0; % ensure that there is no wind load on the car
        end
    else
        Fmodal_wind = zeros(Nm+Ncar);
    end
    
    if Ncar>0
        [Fmodal_car] = getVehicleLoad(t,vehicle,Bridge,Nm,y,phi0_N,idt);
        Fmodal = Fmodal_wind+ Fmodal_car; % load on the bridge
    else
        Fmodal = Fmodal_wind;
    end
        
    if idt ==1 % initial acceleration
        DoM = zeros(Nm);
        VoM = zeros(Nm);
        AoM = M(1:Nm,1:Nm)\(Fmodal(1:Nm,1:Nm)-C(1:Nm,1:Nm).*VoM-K(1:Nm,1:Nm).*DoM);
        [DoM,VoM,AoM,Do(:,:,idt),Vo,Ao(:,:,idt),Aoz_car(:,idt),~] =...
            Newmark(dt,DoM,VoM,AoM,Fmodal(1:Nm,1:Nm),M(1:Nm,1:Nm),K(1:Nm,1:Nm),C(1:Nm,1:Nm),phi0,Ncar,Nyy,Nm);
    else
        [DoM,VoM,AoM,Do(:,:,idt),Vo,Ao(:,:,idt),Aoz_car(:,idt),~] =...
            Newmark(dt,DoM,VoM,AoM,Fmodal(1:Nm,1:Nm),M(1:Nm,1:Nm),K(1:Nm,1:Nm),C(1:Nm,1:Nm),phi0,Ncar,Nyy,Nm);
    end
    rt = Do(3,:,idt)';
    dry = Vo(1,:)';
    drz = Vo(2,:)';
    drt = Vo(3,:)';
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Fmodal]= getFmodal(CoeffAero,C1,C2,C3,y,u,w,dry,drz,drt,D,B,rt,k,rho,Nmodes,phi0_N)
        % [Fmodal]= getFmodal(y,u,w,dry,drz,drt,D,B,rt) computes the modal
        %  forces acting on the main-span of a suspension bridge bridge deck
        %%
        Vrel = (u-dry);
        W = (w-drz-k*B.*drt);
        Vrel = sqrt(Vrel.^2+W.^2); % is [Nyy x 1]
        beta = atan(W./Vrel); % is [Nyy x 1]
        alpha = rt+beta; % is [Nyy x N]
        COEFF = 1/2*rho*B.*Vrel.^2;     
        Ftot(:,1)=COEFF.*(D/B.*CoeffAero(C1,alpha).*cos(beta) - CoeffAero(C2,alpha).*sin(beta));
        Ftot(:,2)=COEFF.*(D/B.*CoeffAero(C1,alpha).*sin(beta) +  CoeffAero(C2,alpha).*cos(beta));
        Ftot(:,3)=COEFF.*(B.*CoeffAero(C3,alpha));
        F = repmat(Ftot,1,Nmodes)';
        Fmodal = trapz(y,F.*phi0_N,2);
        Fmodal = diag(Fmodal(:));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [x1,dx1,ddx1,Do,Vo,Ao,Aoz_car,Doz_car] =...
            Newmark(dt,x0,dx0,ddx0,F,M,K,C,phi0,Ncar,Nyy,Nm,varargin)
        
        
        % options: default values
        inp = inputParser();
        inp.CaseSensitive = false;
        inp.addOptional('alpha',1/4);
        inp.addOptional('beta',1/2);
        inp.parse(varargin{:});
        % shorthen the variables name
        alphaCoeff = inp.Results.alpha ;
        beta = inp.Results.beta;
        
        aDT2 = (alphaCoeff.*dt.^2);
        aDT = (alphaCoeff.*dt);
        
        A = (1./aDT2.*M+beta/aDT*C+K);
        B1 = F+M.*(1./aDT2*x0+1./aDT*dx0+(1/(2*alphaCoeff)-1)*ddx0);
        B2 = C.*(beta/aDT*x0+(beta/alphaCoeff-1).*dx0);
        B3 = (beta/alphaCoeff-2)*dt/2*ddx0;
        
        x1 = A\(B1+B2+B3);
        ddx1= 1/aDT2.*(x1-x0)-1/aDT.*dx0-(1/(2*alphaCoeff)-1).*ddx0;
        dx1= dx0+(1-beta).*dt*ddx0+beta.*dt*ddx1;
        
        x2=diag(x1);
        dx2=diag(dx1);
        ddx2=diag(ddx1);
        Vo = zeros(3,Nyy);
        Do = zeros(3,Nyy);
        Ao = zeros(3,Nyy);
        %         nu = zeros(3,Nmodes);
        
        for oo = 1:3
            Ao(oo,:) = phi0(oo:3:Nm,:)'*ddx2(oo:3:Nm);
            Do(oo,:) = phi0(oo:3:Nm,:)'*x2(oo:3:Nm);
            Vo(oo,:) = phi0(oo:3:Nm,:)'*dx2(oo:3:Nm);
        end
        Aoz_car = zeros(Ncar,1);
        Doz_car = zeros(Ncar,1); 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Fmodal_car] = getVehicleLoad(t,vehicle,Bridge,Nm,y,phi0_N,idt)
        Ny = numel(Bridge.x);
        Nv = numel(vehicle);
        F1 = zeros(Ny,1);
        for rr=1:Nv
            F0 = movingLoad(t(idt),vehicle(rr).tStart,Bridge.L,vehicle(rr).speed,Bridge.x,vehicle(rr).mass,vehicle(rr).direction);
            F1 = F1 + F0(:); % add concentrated load for each car in a given position on the bridge
        end
        F1 = repmat(F1,1,Nm);
%       no force applied on the car since it is simply a moving mass
        Fmodal_car = zeros(Nm+Nv,1); % preallocation 
        Fmodal_car(1:Nm) = trapz(y,F1'.*phi0_N,2); % transform nodal forces to modal forces
        Fmodal_car(1:3:Nm)=0; % no forces on lateral modes
        Fmodal_car(3:3:Nm)=0; % no forces on torsional modes
        Fmodal_car = diag(Fmodal_car); % no coupling for the car
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
