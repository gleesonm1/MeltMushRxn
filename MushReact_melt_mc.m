%% code to simulate reactive porous flow at a single horizon in a mush system
clear all
close all

%% load in MELTS modules
addpath('MELTS')
addpath('function')
warning('off','all')

%% load data
SWIR_Mineral_major=readtable('Results/AB_mineral_major.xlsx');
SWIR_glass=readtable('Results/SWIR_glass.xlsx');

%% define initial mush composition
% set assimilant comp
Lith="Troctolite";
Mineral="Ol";
depth=500;

Ol_As = zeros(19, 1); % array for Si, Ti, Al, Fe3+, Cr, Fe, Mn, Mg, Ni, Co, Ca, Na, K, P, H, C, S, Cl, F oxides respectively
Ol_As([1 6 8 11]) = [mean(SWIR_Mineral_major.SiO2((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.FeO((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.MgO((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.CaO((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth)))];% [39 9 51 0.0];

Mineral="Pl";
Pl_As = zeros(19, 1); % array for Si, Ti, Al, Fe3+, Cr, Fe, Mn, Mg, Ni, Co, Ca, Na, K, P, H, C, S, Cl, F oxides respectively
Pl_As([1 3 6 11:13]) = [mean(SWIR_Mineral_major.SiO2((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.Al2O3((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.FeO((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.CaO((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.Na2O((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.K2O((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth)))];%[47.06 34.34 0 18.11 1.4 0.04];%[52.52 30.35 0.0 13.66 3.84 0.03]; %0.19

Mineral="Cpx";
Cp_As = zeros(19, 1); % array for Si, Ti, Al, Fe3+, Cr, Fe, Mn, Mg, Ni, Co, Ca, Na, K, P, H, C, S, Cl, F oxides respectively
Cp_As([1:3 5 6 8 11 12]) = [mean(SWIR_Mineral_major.SiO2((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    0.67.*mean(SWIR_Mineral_major.TiO2((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.Al2O3((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    0.*mean(SWIR_Mineral_major.Cr2O3((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.FeO((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.MgO((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.CaO((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth))) ...
    mean(SWIR_Mineral_major.Na2O((SWIR_Mineral_major.Mineral==Mineral) & (SWIR_Mineral_major.Lithology==Lith) & (SWIR_Mineral_major.GrainArea=="core") & (SWIR_Mineral_major.Depth_mbsf_>depth)))];%[52.49 0.36 3.51 0 3.74 16.80 21.94 0.38];

% set initial melt composition
Num=[1 2 4];
Fe3=0.10*80/72*mean(SWIR_glass.FeOt(Num)-1);
Fe2=(1-0.10)*mean(SWIR_glass.FeOt(Num)-1);
bulk = zeros(19, 1);

% Primitive SWIR basalt composition
bulk([1:6 8 11:15]) = [mean(SWIR_glass.SiO2(Num)) mean(SWIR_glass.TiO2(Num)) ...
    mean(SWIR_glass.Al2O3(Num)) Fe3 0 Fe2 mean(SWIR_glass.MgO(Num)) mean(SWIR_glass.CaO(Num)) ...
    mean(SWIR_glass.Na2O(Num)) mean(SWIR_glass.K2O(Num)) mean(SWIR_glass.P2O5(Num)) 0.20]; %0.03 Ti 1.1

% Composition from Dick et al. 2000
% bulk([1:6 8 11:15]) = [51.0 1.57 15.7 0.1*9.07*(159.69/2)/71.844 0 0.9*9.07 ...
%     8.09 11 3.01 0.08 0.16 0.2];
%% define initial conditions and carry out FC calculation
As_Cond(1)=1000;

liq=bulk/sum(bulk)*100;

Cond_start=[As_Cond(1) 1230];
Cond_end=[As_Cond(1) 1050];
dc=[0 1];
Frac=true;
liquidus=false;
Fluids=false;
fO2="None";

[ptpath,CompositionFC,MassFC,ConditionsFC]=Crystalise(liq,Cond_start,Cond_end,dc,Frac,liquidus,fO2,Fluids);

ConditionsFC.H=100.*ptpath.getListProperty('h','liquid1')'./ptpath.getListProperty('mass','liquid1')';

%% define figures
Input = load('input_scenario1.mat');
%% extract compositional information
Mineral="Cpx";

%% define key parameters for models
len=2;

Results = [];

for model=1:75
    % define mass fraction of melt in the initial system
    phi = Input.Results(model,2);
    
    % define M^{Ass} parameter
    M = Input.Results(model,3);
    
    % define the modal mineralogy of the bulk mush
    X_Mush=struct('olivine1', 0.05, 'plagioclase1', 0.55, 'clinopyroxene1', 0.40); % 0.05, 0.55, 0.4
    C_Mush=struct('olivine1', Ol_As, 'plagioclase1', Pl_As, 'clinopyroxene1', Cp_As);
    
    % set up mush assimilants
    O = Input.Results(model,4);
    P = Input.Results(model,5);
    C = Input.Results(model,6);
    
    % define the modal mineralogy of the reacted assemblage
    X_React=struct('olivine1', O, 'plagioclase1', P, 'clinopyroxene1', C);

    % set initial conditions
    As_Cond(2)=round(Input.Results(model,1));
    deltaT=0;

    % set composition and enthalpy of new liquid phase for reaction
    % calculations
    liq=CompositionFC.liquid1{ConditionsFC.T==As_Cond(2)+deltaT,:};

    %% Melt-mush reaction calculations
    [melts,MeltMushTables]=MeltMushRxn_MeltPath(len,As_Cond,phi,M,X_Mush,C_Mush,X_React,liq,deltaT);

    %% save results
    % determine % porosity change
    MeltMass=MeltMushTables.Liq_Mass.liquid1(2:end)./MeltMushTables.Liq_Mass.liquid1(1);
    ReactMass=zeros(length(MeltMass),1);

    for i=1:length(ReactMass)
        if i==1
            Mliq=1;
            RM(i)=M*Mliq;
        else
            Mliq=Mliq*(MeltMass(i));
            RM(i)=M*Mliq;
        end        
    end
    RM=cumsum(RM);

    R_Min=struct();
    React=fieldnames(X_React);
    for i=1:length(React)
        R_Min.(React{i})=MeltMushTables.React_Mass.(React{i})(2:end)./MeltMushTables.React_Mass_In.(React{i})(1);    
    end
    
    Results((model-1)*len+1:(model-1)*len+len,1) = MeltMass';
    Results((model-1)*len+1:(model-1)*len+len,2) = O+zeros(len,1);
    Results((model-1)*len+1:(model-1)*len+len,3) = P+zeros(len,1);
    Results((model-1)*len+1:(model-1)*len+len,4) = C+zeros(len,1);
    Results((model-1)*len+1:(model-1)*len+len,5) = phi+zeros(len,1);
    Results((model-1)*len+1:(model-1)*len+len,6) = M+zeros(len,1);
    Results((model-1)*len+1:(model-1)*len+len,7) = As_Cond(2)+zeros(len,1);
    off = MeltMushTables.React_Mass.olivine1./sum(MeltMushTables.React_Mass{:,:},2);
    Results((model-1)*len+1:(model-1)*len+len,8) = off(2:end);
    off = MeltMushTables.React_Mass.plagioclase1./sum(MeltMushTables.React_Mass{:,:},2);
    Results((model-1)*len+1:(model-1)*len+len,9) = off(2:end);
    off = MeltMushTables.React_Mass.clinopyroxene1./sum(MeltMushTables.React_Mass{:,:},2);
    Results((model-1)*len+1:(model-1)*len+len,10) = off(2:end);
    Results((model-1)*len+1:(model-1)*len+len,11) = MeltMushTables.Conditions.H(1)-MeltMushTables.Conditions.H(2:end);
end

%% plot results
figure('rend','painters','pos',[20 10 1250 320])
subaxis(1,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(Results(:,1), 20, 'Normalization', 'probability', 'FaceColor', 'r')
xlabel('Melt Mass Ratio', 'FontSize', 16)
ylabel('Density', 'FontSize', 16)

subaxis(1,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(Results(:,1), Results(:,3), 40, Results(:,4), 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
ylabel('X_{Plagioclase}^{Reacted assemblage}', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

%colormap viridis

subaxis(1,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(Results(:,1), -Results(:,11), 40, Results(:,4), 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
ylabel('\Delta Plag', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Olivine}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

%%
figure('rend','painters','pos',[20 10 450 420])
% sample correlation matrix
r = corrcoef(Results)

% labels
labels = ["MMR", "X_{Olivine}", "X_{Plagioclase}", "X_{Clinopyroxene}", "phi", "M", "T", "X_{OlCryst}", "X_{PlagCryst}", "X_{ClinoCryst}",  "delta H"];
% scatter plot
n = size(r, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
scatter(x(:), y(:), 300.*abs(r(:)), r(:), 'filled', 'MarkerFaceAlpha', 0.6)
% enclose markers in a grid
xl = [1:n+1;repmat(n+1, 1, n+1)];
xl = [xl(:, 1), xl(:, 1:end-1)];
yl = repmat(n+1:-1:1, 2, 1);
line(xl, yl, 'color', 'k') % horizontal lines
line(yl, xl, 'color', 'k') % vertical lines
% show labels
text(1:n, (n:-1:1) + 0.5, labels, 'HorizontalAlignment', 'right')
text((1:n) + 0.5, repmat(n + 1, n, 1), labels, ...
    'HorizontalAlignment', 'right', 'Rotation', 270)
h = gca;
colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.9;
axis(h, 'equal')
%colormap('viridis')
