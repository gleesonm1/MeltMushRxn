%% melt-mush reaction. Each step of the model starts with a pre-defined mush composition and the melt phase determined from the previous step of the model
clear all
close all

%% load in appropriate packages and functions
addpath('MELTS')
addpath('function')
warning('off','all')

%% set initial mineral compositions
Ol_As = zeros(19, 1); % array for Si, Ti, Al, Fe3+, Cr, Fe, Mn, Mg, Ni, Co, Ca, Na, K, P, H, C, S, Cl, F oxides respectively
Ol_As([1 6 8 11]) = [39.40 16.87 42.96 0.05]; % SiO2, FeOt, MgO, CaO
Ol_As=Ol_As/sum(Ol_As)*100; % make sure the composition (in g) is normalised to 100.

Pl_As = zeros(19, 1); % array for Si, Ti, Al, Fe3+, Cr, Fe, Mn, Mg, Ni, Co, Ca, Na, K, P, H, C, S, Cl, F oxides respectively
Pl_As([1 3 6 11:13]) = [50.34 31.31 0.14 14.49 3.66 0.04]; % SiO2, Al2O3, FeOt, CaO, Na2O, K2O
Pl_As=Pl_As/sum(Pl_As)*100; % make sure the composition (in g) is normalised to 100.

Cp_As = zeros(19, 1); % array for Si, Ti, Al, Fe3+, Cr, Fe, Mn, Mg, Ni, Co, Ca, Na, K, P, H, C, S, Cl, F oxides respectively
Cp_As([1:3 5 6 8 11 12]) = [51.67 0.64 3.34 0 4.94 16.52 21.64 0.37]; % SiO2, TiO2, Al2O3, Cr2O3, FeOt, MgO, CaO, Na2O
Cp_As=Cp_As/sum(Cp_As)*100; % make sure the composition (in g) is normalised to 100.

%% set initial melt composition
bulk = zeros(19, 1);
bulk([1:6 8 11:15]) = [48.25 0.88 17.77 0.84 0 6.79 9.47 11.74 2.79 0.05 0.17 0.2]; % SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MgO, CaO, Na2O, K2O, P2O5, H2O
liq=bulk/sum(bulk)*100; % make sure the composition (in g) is normalised to 100.

%% define initial conditions and carry out FC calculation
pressure=1000; % pressure (bars)
T_start = 1230; % temperature to start FC calculation (oC)
T_end = 1050; % temperature to end FC calculation (oc)

phi=0.20; % initial porosity of the system

R = 0.5; % define recharge rate

M=0.30; % mass of solid reactant relative to mass of melt

% set mineralogy of the mush (must sum to 1)
X_Mush=struct('olivine1', 0.05, 'plagioclase1', 0.55, 'clinopyroxene1', 0.40);

% set composition of the mush components
C_Mush=struct('olivine1', Ol_As, 'plagioclase1', Pl_As, 'clinopyroxene1', Cp_As);

% define the modal mineralogy of the reacted assemblage
X_React=struct('olivine1', 0.05, 'plagioclase1', 0.25, 'clinopyroxene1', 0.70);
X_React_new=struct('olivine1', 0.05, 'plagioclase1', 0.70, 'clinopyroxene1', 0.25);

T_mush=1210; % temperature of the mush at the start of reaction (oC)
deltaT=00; % offset in temperature between melt and mush (oC)

%% carry out fractional crystallisation modelling
Cond_start=[pressure T_start]; % conditions at start of model
Cond_end=[pressure T_end]; % conditions at end of model
dc=[0 1]; % pressure and temperature increments
Frac=true; % fractionate solids
liquidus=false; % don't find the liquidus prior to calculation of LLD
Fluids=false; % don't fractionate fluids
fO2="None"; % no fO2 buffer set

% carry out FC modelling
[ptpath,CompositionFC,MassFC,ConditionsFC]=Crystalise(liq,Cond_start,Cond_end,dc,Frac,liquidus,fO2,Fluids); 

% save the enthalpy of the melt phase at each model step
%ConditionsFC.H=100.*ptpath.getListProperty('h','liquid1')'./ptpath.getListProperty('mass','liquid1')'; 

%% Melt-mush reaction calculations
% set composition and enthalpy of new liquid phase for reaction calculations
liq=CompositionFC.liquid1{ConditionsFC.T==T_mush+deltaT,:};

% run Melt-mush reaction calculations
[melts,MeltMushTables]=MeltMushRxn(20,[pressure T_mush],phi,M,R,X_Mush,C_Mush,X_React,X_React_new,liq, deltaT);

%% plot results
% normalise reacted assemblage
MeltMushTables.React_Mass.liquid1=MeltMushTables.Liq_Mass.liquid1;
React_Mass_norm=MeltMushTables.React_Mass;
React_Mass_norm{:,:}=MeltMushTables.React_Mass{:,:}./sum(MeltMushTables.React_Mass{:,:},2);

% get bulk assemblage
Bulk_System=MeltMushTables.Unreact_Mass;
Bulk_System{:,:}=MeltMushTables.Unreact_Mass{:,:}+MeltMushTables.React_Mass{:,:};

% area plots to show how the modal proportions change with time
x=linspace(0,20,20+1);
Y=[React_Mass_norm.olivine1,React_Mass_norm.plagioclase1,React_Mass_norm.clinopyroxene1,React_Mass_norm.liquid1];

figure('rend','painters','pos',[20 10 850 450])
subaxis(1,2,2,'SpacingVert',0.04,'SpacingHoriz',0.04)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'YAxisLocation','right')
hold on
box on
a = area(x,Y);
a(1).FaceColor = [0 0.56 0];
a(2).FaceColor = [0 0 0.56];
a(3).FaceColor = [0 0.2 0];
a(4).FaceColor = [1 0.6 .6];
ylim([0 1])
ylabel('Modal Proportions','FontSize',16)
xlabel('Reaction no.','FontSize',16)
title('Reacted Assemblage','FontSize',16)

Y=[Bulk_System.olivine1,Bulk_System.plagioclase1,Bulk_System.clinopyroxene1,Bulk_System.liquid1];

subaxis(1,2,1,'SpacingVert',0.04,'SpacingHoriz',0.04)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
a = area(x,Y);
a(1).FaceColor = [0 0.56 0];
a(2).FaceColor = [0 0 0.56];
a(3).FaceColor = [0 0.2 0];
a(4).FaceColor = [1 0.6 .6];
ylim([0 1])
ylabel('Modal Proportions','FontSize',16)
xlabel('Reaction no.','FontSize',16)
title('Bulk Assemblage','FontSize',16)

