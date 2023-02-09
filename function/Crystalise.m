function [ptpath,Composition,Mass,Conditions]=Crystalise(bulk,Cond_start,Cond_end,dc,Frac,liquidus,fO2,Fluids)
%% Crystallisation function to simulate the differentation of magma.
%   Required Inputs: 
%           'bulk' (19-by-1 vector with initial liquid composition)
%           'Cond_start' (2-by-1 vector containing pressure (1) and
%           temperature (2) at the start of the model)
%           'Cond_end' (2-by-1 vector containing pressure (1) and
%           temperature (2) at the end of the model)
%           'dc' (2-by-1 vector of pressure (1) and temperature (2) with each step)
%           'Frac' (true or false - determines whether frac or equil crystallisation)
%   Optional Inputs:
%           'liquidus' (true or false; finds liquidus before calculation performed. If absent will assume false)
%           'fO2' (string outlining how fO2 varies from QFM, can set to 'None')
%           'Fluids' (true or false; choice to frac fluids)
%
%   Outputs:
%           'ptpath' (MELTSdynamic for the crystallisation path)
%           'Composition' (Stacked table containing the composition of each
%           phase that is present in the system; for example,
%           Composition.liquid1.MgO would return the MgO content of the
%           liquid at each step)
%           'Mass' (Table containing the mass (per 100g) of each phase at
%           each step)
%           'Conditions' (Table containing the pressure and temperature at
%           each model step)
%       Note: If only assimilate and melts are specified as outputs the
%       code will terminate after the main for loop.

    ptpath=MELTSdynamic(1);
    ptpath.engine.pressure=Cond_start(1);
    ptpath.engine.temperature=Cond_start(2);
    
    % set bulk composition and normalise to 100g
    ptpath.engine.set('BulkComposition',100*bulk/sum(bulk));
    
    if nargin>6 && fO2~='None'
        ptpath.engine.setSystemProperties("Log fO2 Path",fO2)
    end
    
    if nargin>5 && liquidus==true
        ptpath.engine.findLiquidus;
    end
        
    if Frac==true
        ptpath.engine.setSystemProperties("Mode", "Fractionate Solids");
    end
    
    if nargin>7 && Fluids == true
        ptpath.engine.setSystemProperties("Mode", "Fractionate Fluids");
    end
        
    ptpath.engine.calcEquilibriumState(1,1);
    
    while ptpath.engine.temperature>=Cond_end(2) && ptpath.engine.pressure>=Cond_end(1)
        % add new node
        ptpath=ptpath.addNodeAfter;
        
        % Change P-T conditions
        ptpath.engine.pressure = ptpath.engine.pressure - dc(1);
        ptpath.engine.temperature = ptpath.engine.temperature - dc(2);
        
        % calculate the equilibrium
        ptpath.engine.calcEquilibriumState(1,1);
        
        liq=ptpath.engine.getProperty('dispComposition','liquid1');
        
        if liq(19)~=0
            break
        end
        
        %ptpath.engine.set('BulkComposition',liq);

    end
    
    if nargout>1
        len=length(ptpath.getListProperty('pressure'))-1;
        
        % extract conditions of model
        N={'P','T'};
        Conditions=table(zeros(len+1,1),zeros(len+1,1));
        Conditions.Properties.VariableNames=N;
        Conditions.P=ptpath.getListProperty('pressure')';
        Conditions.T=ptpath.getListProperty('temperature')';
        
        % Define names of all possible phases in the MELTS models (this can be
        % cut down to save time
        Names={'liquid1','fluid1','olivine1','olivine2','clinopyroxene1',...
            'clinopyroxene2','plagioclase1','plagioclase2',...
            'spinel1','spinel2','orthopyroxene1','orthopyroxene2',...
            'kfeldspar1','kfeldspar2','apatite1','rhmoxide1',...
            'quartz1','biotite1','whitlockite1'};
        Mass=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
        Mass.Properties.VariableNames=Names;

        SolidPhases=zeros(len+1,19,width(Mass));
        Elements={'SiO2','TiO2','Al2O3','Fe2O3','Cr2O3','FeO','MnO','MgO',...
            'NiO','CoO','CaO','Na2O','K2O','P2O5','H2O','CO2','SO2','Cl','F'};
        Phase=table(zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),...
            zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),...
            zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),...
            zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),...
            zeros(len+1,1),zeros(len+1,1),zeros(len+1,1),zeros(len+1,1));
        Phase.Properties.VariableNames=Elements;
        Composition=table(Phase,Phase,Phase,Phase,...
            Phase,Phase,Phase,Phase,Phase,...
            Phase,Phase,Phase,Phase,Phase,...
            Phase,Phase,Phase,Phase,Phase);
        Composition.Properties.VariableNames=Names;

        for j=1:width(Mass)
            Mass{1:len+1,j}=ptpath.getListProperty('mass',Names(j))';
            SolidPhases(:,:,j)=ptpath.getListProperty('dispComposition',Names(j))';
            Phase{:,:}=SolidPhases(:,:,j);
            Composition(:,j)=table(Phase);
            if Names(j)=="rhmoxide1"
                Mass{1:len+1,j}=ptpath.getListProperty('mass','rhm-oxide1')';
                SolidPhases(:,:,j)=ptpath.getListProperty('dispComposition','rhm-oxide1')';
                Phase{:,:}=SolidPhases(:,:,j);
                Composition(:,j)=table(Phase);
            end 
            if Names(j)=="kfeldspar1"
                Mass{1:len+1,j}=ptpath.getListProperty('mass','k-feldspar1')';
                SolidPhases(:,:,j)=ptpath.getListProperty('dispComposition','k-feldspar1')';
                Phase{:,:}=SolidPhases(:,:,j);
                Composition(:,j)=table(Phase);
            end  
            if Names(j)=="kfeldspar2"
                Mass{1:len+1,j}=ptpath.getListProperty('mass','k-feldspar2')';
                SolidPhases(:,:,j)=ptpath.getListProperty('dispComposition','k-feldspar1')';
                Phase{:,:}=SolidPhases(:,:,j);
                Composition(:,j)=table(Phase);
            end  
        end
    end 
end