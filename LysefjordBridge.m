function [Bridge] = LysefjordBridge(Nyy)
% [Bridge] = LysefjordBridge(Nyy,Nmodes) computes the structural
% properties of the Lysefjord bridge, for a given number of mode Nmodes,
% and a number of nodes Nnodes for the discretized span.
%
% % Input
% * Nyy: scalar: Number of nodes for the discretization of the span.
%
% Output:
% * bridge: structures with the fields:
%         - B : deck width
%         - D : Deck height
%         - L : length of main span (m)
%         - x : discretisation of bridge axis into normalized coordinates
%         - E : young modulus steel (Pa)
%         - Ec : young modulus steel (Pa)
%         - Ac : cross section main cable (m^2)
%         - g : acceleration of gravity
%         - m : lineic mass of girder (kg/m)
%         - mc : lineic mass of cable (kg/m)
%         - ec : sag (m)*
%         - hm : hanger length at mid span (m)
%         - hr : distance between shear center and hanger attachment
%         - bc : distance betweem main cable (m)
%         - H_cable : horizontal tension force in the main cables
%         - Cd : drag coefficient
%         - dCd : first derivative of drag coefficient
%         - Cl : lift coefficient
%         - dCl : first derivative of lift coefficient
%         - Cm : pitching moment coefficient
%         - dCm : first derivative of pitching moment coefficient
%         - Iz : Moment of inertia with respect to bending about y axis (used for lateral bridge analysis)
%         - Iy : Moment of inertia with respect to bending about z axis (used for vertical bridge analysis)
%         - m_theta :
%         - Iw : WARPING RESISTANCE
%         - GIt : TORSIONAL STIFFNESS
%         - k : bridge constant . cf aerodynamic of streamlined bridge
%
% Author information
%
% Author: Etienne Cheynet -- University of Stavanger, Norway.
%  last modified: 21/05/2016
%
% See also eigenBridge
%
%

% GENERAL INPUTS
Bridge.B = 12.3; % deck width
Bridge.D = 2.76; % Deck height
Bridge.L = 446 ; % length of main span (m) *
Bridge.Nyy=Nyy;% Discretisation of bridge main span in Nyy points
Bridge.x = linspace(0,1,Bridge.Nyy); % discretisation of bridge axis into normalized coordinates


Bridge.E = 210000e6; % young modulus steel (Pa) *
Bridge.Ec = 180000e6; % young modulus steel (Pa) *

Bridge.Ac = 0.038 ;% cross section main cable (m^2) *
Bridge.g = 9.81;
Bridge.m =5350 ; % lineic mass of girder (kg/m)*
Bridge.mc =408 ; % lineic mass of cable (kg/m)*
Bridge.ec= 45; % sag (m)*
Bridge.hm = 3 ; % hanger length at mid span (m)*
Bridge.hr =0.400; % distance between shear center and hanger attachment
Bridge.bc = 10.2500; % distance betweem main cable (m)
Bridge.H_cable = Bridge.m*Bridge.g*Bridge.L^2/(16*Bridge.ec)*(1+2*Bridge.mc/Bridge.m*(1+4/3*(Bridge.ec/Bridge.L)^2));
% aerodynamic coefficient (quasi steady terms)
Bridge.Cd = 1;% drag coefficient
Bridge.dCd = 0;% first derivative of drag coefficient
Bridge.Cl = 0.1;% lift coefficient
Bridge.dCl = 3;% first derivative of lift coefficient
Bridge.Cm = 0.02;% pitching moment coefficient
Bridge.dCm = 1.12;% first derivative of pitching moment coefficient
% ADDITIONAL INPUTS FOR LATERAL MODES
Bridge.Iz = 4.952; % Moment of inertia with respect to bending about y axis (used for lateral bridge analysis)
% ADDITIONAL INPUTS FOR VERTICAL MODES
Bridge.Iy = 0.429; % Moment of inertia with respect to bending about z axis (used for vertical bridge analysis)
% ADDITIONAL INPUTS FOR TORSIONAL MODES
% Bridge.m_theta = 82430; %kg.m^2/m*
Bridge.m_theta = 58730; %kg.m^2/m*

Bridge.Iw = 4.7619; % WARPING RESISTANCE
Bridge.GIt = 0.75e11; % TORSIONAL STIFFNESS
Bridge.k = 1/4; % bridge constant . cf aerodynamic of streamlined bridge
end