%% melt-mush reaction. Each step of the model starts with a pre-defined mush composition and the melt phase determined from the previous step of the model
clear all
close all

%% load in appropriate packages and functions
addpath('MELTS')
addpath('function')
warning('off','all')

%% load SWIR data
SWIR = readtable('Results/SWIR_glass.xlsx');

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
T_end = 1100; % temperature to end FC calculation (oc)

phi=0.20; % initial porosity of the system

M=0.30; % mass of solid reactant relative to mass of melt

% set mineralogy of the mush (must sum to 1)
X_Mush=struct('olivine1', 0.05, 'plagioclase1', 0.55, 'clinopyroxene1', 0.40);

% set composition of the mush components
C_Mush=struct('olivine1', Ol_As, 'plagioclase1', Pl_As, 'clinopyroxene1', Cp_As);

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

%% set up element plots
figure('rend','painters','pos',[20 10 850 800])
subaxis(2,2,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
plot(SWIR.MgO, SWIR.SiO2, 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerSize', 5)
plot(CompositionFC.liquid1.MgO, CompositionFC.liquid1.SiO2, '-k', 'LineWidth', 2)
xlabel('MgO (wt%)', 'FontSize', 16)
ylabel('SiO_{2} (wt%)', 'FontSize', 16)

subaxis(2,2,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12, 'YAxisLocation','right')
hold on
box on
plot(SWIR.MgO, SWIR.Al2O3, 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerSize', 5)
plot(CompositionFC.liquid1.MgO, CompositionFC.liquid1.Al2O3, '-k', 'LineWidth', 2)
xlabel('MgO (wt%)', 'FontSize', 16)
ylabel('Al_{2}O_{3} (wt%)', 'FontSize', 16)

subaxis(2,2,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
ylabel('H (J/kg)','FontSize',16)
xlabel('Melt Mass Ratio','FontSize',16)

for i = 1:3
    % define the modal mineralogy of the reacted assemblage
    if i == 1
        X_React=struct('olivine1', 0.05, 'plagioclase1', 0.75, 'clinopyroxene1', 0.20);
    elseif i ==2
        X_React=struct('olivine1', 0.05, 'plagioclase1', 0.55, 'clinopyroxene1', 0.40);
    else
        X_React=struct('olivine1', 0.05, 'plagioclase1', 0.25, 'clinopyroxene1', 0.70);
    end
        

    T_mush=1180; % temperature of the mush at the start of reaction (oC)
    deltaT=0; % offset in temperature between melt and mush (oC)

    

    % save the enthalpy of the melt phase at each model step
    %ConditionsFC.H=100.*ptpath.getListProperty('h','liquid1')'./ptpath.getListProperty('mass','liquid1')'; 

    %% Melt-mush reaction calculations
    % set composition and enthalpy of new liquid phase for reaction calculations
    liq=CompositionFC.liquid1{ConditionsFC.T==T_mush+deltaT,:};

    % run melt-mush reaction calculations
    [melts,MeltMushTables]=MeltMushRxn_MeltPath(6,[pressure T_mush],phi,M,X_Mush,C_Mush,X_React,liq,deltaT);

    % determine % porosity change
    MeltMassRatio=MeltMushTables.Liq_Mass.liquid1(2:end)./MeltMushTables.Liq_Mass.liquid1(1);
    
%     X_Mu=struct2table(X_React);
%     H_react=H_Mush*X_Mu{:,:}';
    
    figure(1)
    subaxis(2,2,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
    if i == 1
        plot(MeltMushTables.React_Composition.liquid1.MgO, MeltMushTables.React_Composition.liquid1.SiO2, '-r', 'LineWidth', 2)
    elseif i == 2
        plot(MeltMushTables.React_Composition.liquid1.MgO, MeltMushTables.React_Composition.liquid1.SiO2, '-b', 'LineWidth', 2)
    else
        plot(MeltMushTables.React_Composition.liquid1.MgO, MeltMushTables.React_Composition.liquid1.SiO2, '-y', 'LineWidth', 2)
    end 
    
    subaxis(2,2,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
    if i == 1
        plot(MeltMushTables.React_Composition.liquid1.MgO, MeltMushTables.React_Composition.liquid1.Al2O3, '-r', 'LineWidth', 2)
    elseif i == 2
        plot(MeltMushTables.React_Composition.liquid1.MgO, MeltMushTables.React_Composition.liquid1.Al2O3, '-b', 'LineWidth', 2)
    else
        plot(MeltMushTables.React_Composition.liquid1.MgO, MeltMushTables.React_Composition.liquid1.Al2O3, '-y', 'LineWidth', 2)
    end 
    
    % determine initial enthalpy of assimilated component
    React=fieldnames(X_React);
    assimilate = MELTSdynamic(1);
    assimilate.engine.pressure=MeltMushTables.Conditions.P(1);
    assimilate.engine.temperature=MeltMushTables.Conditions.T(1);
    
    C_React=struct();
    H_As=zeros(1,length(React));
    for j=1:length(React)  
        C_React.(React{j})=zeros(19,1);
        C_React.(React{j})=MeltMushTables.Unreact_Composition.(React{j}){1,:}';      
        assimilate.engine.calcPhaseProperties(React{j},MeltMushTables.Unreact_Composition.(React{j}){1,:});
        H_As(j)=assimilate.engine.getProperty('h',React{j});
    end

    % convert struct to table to extract data
    C_As=struct2table(C_React);
    X_As=struct2table(X_React);

    % determine bulk composition of reacted component
    C_Assim=C_As{:,:}*X_As{:,:}';

    % determine enthalpy of reacted component
    H_Assim=H_As*X_As{:,:}';

    % ensure enthalpy is relative to 100g
    H_react=100.*H_Assim./sum(C_Assim);

    subaxis(2,2,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
    if i == 1
        plot(MeltMassRatio, MeltMushTables.Conditions.H(2:end), '-r', 'LineWidth', 2)
        plot([0.7 1.2], [H_react H_react], ':r')
    elseif i == 2
        plot(MeltMassRatio, MeltMushTables.Conditions.H(2:end), '-b', 'LineWidth', 2)
        plot([0.7 1.2], [H_react H_react], ':b')
    else
        plot(MeltMassRatio, MeltMushTables.Conditions.H(2:end), '-y', 'LineWidth', 2)
        plot([0.7 1.2], [H_react H_react], ':y')
    end 
end

plot([0.7 1.2], [MeltMushTables.Conditions.H(1) MeltMushTables.Conditions.H(1)], '-k', 'LineWidth', 2)


