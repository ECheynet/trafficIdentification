function [wn,phi,phi_cables] = eigenBridge(Bridge,Nmodes)
% function [wn,phi,phi_cables] = eigenBridge(Bridge) computes
% the mode shape and eigen-frequency of a single-span suspension bridge.
% It is designed to be fast and straitghforward, although it might slightly
% less accurate than a Finite element model.
%
% INPUT:
%   Bridge: type: structure (see file studyCase.m)
%   Nyy: type: float [1 x 1] : number of "nodes" consituting the deck
%   Nmodes: type: float [1 x 1] : number of "modes" to be computed
%
% OUTPUT
% wn: type: 2D matrix [3 x Nmodes] : eigen-modes of the suspension bridge
%   wn(1,:) --> all the eigen-frequencies for the lateral displacement (x axis)
%   wn(2,:) --> all the eigen-frequencies for the vertical displacement (z axis)
%   wn(3,:) --> all the eigen-frequencies for the torsional angle (around y axis)
%
% phi: type: 3D matrix [3 x Nmodes x Nyy] : modes shapes of the bridge girder ( =deck)
%   phi(1,i,:) --> mode shape for the i-th eigen frequency for the lateral
%                  bridge displacement ( x axis)
%   phi(2,i,:) --> mode shape for the i-th eigen frequency for the vertical
%                  bridge displacement  ( z axis)
%   phi(3,i,:) --> mode shape for the i-th eigen frequency for the torsional
%                  bridge angle (around y axis)
%
% phi_cables: type: 2D matrix [Nmodes x Nyy]:
% modes shapes of the bridge main cables along x axis. No mode
% shapes for vertical and torsional motion of the cable are calculated.
%
% Author info:
% E. Cheynet , University of Stavanger. last modified: 31/12/2016
%
% see also LysefjordBridge hardangerBridge
%

%%
% preallocation
if isfield(Bridge,'Nyy'),
    Nyy = Bridge.Nyy;
else
    error(' ''Nyy'' is not a field of the structure ''Bridge'' ')
end

wn = zeros(3,Nmodes); % eigen frequency matrix
phi = zeros(3,Nmodes,Nyy); % mode shapes of bridge girder
phi_cables = zeros(Nmodes,Nyy); % mode shapes of cables

% discretisation of bridge span
x = linspace(0,Bridge.L,Nyy); % vector bridge span
x = x./Bridge.L ;% reduced length


%% LATERAL MOTION
% preallocation
alpha = zeros(2*Nmodes,2*Nmodes); % reduced variable
beta = zeros(2*Nmodes,2*Nmodes); % reduced variable
gamma = zeros(2*Nmodes,2*Nmodes); % reduced variable

m_tilt= Bridge.m;
mc_tilt = 2.*Bridge.mc;

% calculation of alpha, beta and gamma
% alpha and beta
% cf E.N. Str�mmen "STRUCTURAL DYNAMICS" for explanations
for nn=1:2*Nmodes
    alpha(nn,nn) = (Bridge.E.*Bridge.Iz).*(nn*pi./Bridge.L).^4;
    beta(nn,nn) = 2*Bridge.H_cable.*(nn*pi./Bridge.L).^2;
end

% gamma
% cf E.N. Str�mmen "STRUCTURAL DYNAMICS" for explanations
for pp=1:2*Nmodes
    for nn=1:2*Nmodes,
        if and(rem(pp,2)==1,rem(nn,2)==0)||and(rem(pp,2)==0,rem(nn,2)==1)
            gamma(pp,nn)=0;
        else
            gamma(pp,nn)= 2*Bridge.m*Bridge.g/(Bridge.L*Bridge.ec).*trapz(x.*Bridge.L,sin(pp.*pi.*x).*sin(nn.*pi.*x)./...
                (1+Bridge.hm/Bridge.ec-4.*(x).*(1-x)));
        end
    end
end

% matrix mass
M = diag(repmat([m_tilt;mc_tilt],[Nmodes,1]));

% Matrix stiwness K
for p=1:2:2*Nmodes
    for nn=1:2:2*Nmodes
        if p==nn
            clear Omega
            Omega(1,1) = alpha((p+1)/2,(nn+1)/2)+gamma((p+1)/2,(nn+1)/2);
            Omega(1,2) = -gamma((p+1)/2,(nn+1)/2);
            Omega(2,1) = -gamma((p+1)/2,(nn+1)/2);
            Omega(2,2) = beta((p+1)/2,(nn+1)/2)+gamma((p+1)/2,(nn+1)/2);
            K(p:p+1,nn:nn+1) = Omega;
        else
            clear V
            V = gamma((p+1)/2,(nn+1)/2).*[1,-1;-1,1];
            K(p:p+1,nn:nn+1) = V;
        end
    end
end

% eigen-value problem solved for non-trivial solutions
[vector,lambda]=eig(K,M,'chol');

wn(1,:) = sqrt(diag(lambda(1:Nmodes,1:Nmodes))); % filling the matrix wn

% Normalization
for ii=1:Nmodes
    % deck mode shape construction using series expansion
    phi(1,ii,:) = vector(1:2:end,ii)'*sin([1:1:Nmodes]'.*pi*x);
    % cables mode shape construction using series expansion
    phi_cables(ii,:) = vector(2:2:end,ii)'*sin([1:1:Nmodes]'.*pi*x); 
    
    modeMax = max([max(abs(phi_cables(ii,:))),max(abs(phi(1,ii,:)))]);
    phi(1,ii,:) = phi(1,ii,:)./modeMax; % normalisation
    phi_cables(ii,:) = phi_cables(ii,:)./modeMax; % normalisation
end


%% VERTICAL MOTION

clear K M
% INITIALISATION
kappa = zeros(Nmodes,Nmodes); %reduced variable
lambda = zeros(Nmodes,Nmodes); %reduced variable
mu = zeros(Nmodes,Nmodes); %reduced variable

le = Bridge.L*(1+8*(Bridge.ec/Bridge.L)^2); % effective length
% cf page 124  "STRUCTURAL DYNAMICS" of E.N. Str�mmen

for nn=1:Nmodes,
    kappa(nn,nn) = Bridge.E*Bridge.Iy.*(nn*pi./Bridge.L).^4;
    lambda(nn,nn) = 2*Bridge.H_cable.*(nn*pi./Bridge.L).^2;
end

for p=1:Nmodes,
    for nn=1:Nmodes,
        if and(rem(p,2)==1,rem(nn,2)==1) % sont impaires
            mu(p,nn)=(32*Bridge.ec/(pi*Bridge.L)).^2*(Bridge.Ec*Bridge.Ac)/(Bridge.L*le)/(p*nn);
        else
            mu(p,nn)= 0;
        end
    end
end

M  = (2*Bridge.mc+Bridge.m).*eye(Nmodes); % mass matrix
K= kappa+lambda +mu; % stiwness matrix


% eigen-value problem solved for non-trivial solutions
[vector,lambda]=eig(K,M,'chol');

wn(2,:) = sqrt(diag(lambda(1:Nmodes,1:Nmodes))); % filling of matrix wn

for ii=1:Nmodes,
    phi(2,ii,:) = vector(1:Nmodes,ii)'*sin([1:Nmodes]'.*pi*x); % mode shape construction using series expansion
    phi(2,ii,:) = phi(2,ii,:)./max(abs(phi(2,ii,:))); % normalisation
end

%% TORSIONAL MOTION

clear K M
omega = zeros(Nmodes,Nmodes);
v = zeros(Nmodes,Nmodes);
V = zeros(Nmodes,Nmodes);
xi = zeros(Nmodes,Nmodes);
m_tilt_theta_tot = zeros(Nmodes,Nmodes);
m_tilt = 2*Bridge.mc+Bridge.m;


m_theta_tot = Bridge.m_theta + Bridge.mc*(Bridge.bc^2/2);

for nn=1:Nmodes,
    omega(nn,nn) = (nn*pi./Bridge.L).^2.*(Bridge.GIt+(nn*pi/Bridge.L)^2*(Bridge.E*Bridge.Iw));
    V(nn,nn) = Bridge.H_cable * (Bridge.bc^2/2)*(nn*pi./Bridge.L).^2;
    v(nn,nn) = (m_tilt)*Bridge.g*Bridge.hr;
    m_tilt_theta_tot(nn,nn) = m_theta_tot;
end

for p=1:Nmodes,
    for nn=1:Nmodes,
        if and(rem(p,2)==1,rem(nn,2)==1) % sont impaires
            xi(p,nn)=(16*Bridge.ec*Bridge.bc/(pi*Bridge.L)).^2*Bridge.Ec*Bridge.Ac/(Bridge.L*le)*1/(p*nn);
        else
            xi(p,nn)= 0;
        end
    end
end


K= omega+V+v+xi; % stiwness matrix
M =  diag(diag(m_tilt_theta_tot)); % mass matrix

% eigen-value problem solved for non-trivial solutions
[vector,lambda]=eig(K,M,'chol');
wn(3,:) = sqrt(diag(lambda(1:Nmodes,1:Nmodes))); % filling the wn matrix

for ii=1:Nmodes,
    phi(3,ii,:) = vector(1:Nmodes,ii)'*sin([1:Nmodes]'.*pi*x); % mode shape construction using series expansion
    phi(3,ii,:) = phi(3,ii,:)./max(abs(phi(3,ii,:))); % normalisation
end


% positiv max value
for ii=1:3,
    for jj=1:Nmodes,
        if abs(min(squeeze(phi(ii,jj,:))))> abs(max(squeeze(phi(ii,jj,:)))),
            phi(ii,jj,:)=-phi(ii,jj,:);
        end
    end
end




end

