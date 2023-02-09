function [melts,MeltMushTables]=MeltMushRxn_MeltPath(len,As_Cond,phi,M,X_Mush,C_Mush,X_React,liq,deltaT)
%% Function to track the evolution of a melt phase through several steps of melt-mush reaction.
%   Required Inputs:
%       'len' (integer stating the number of melt-mush reactions to
%       simulate)
%       'As_Cond' (2-by-1 vector containing the pressure (1) and initial
%       temperature (2) of melt-mush reaction.
%       'phi' (initial porosity of the mush)
%       'M' (relative mass of solid material reacted compared to the mass
%       of melt in the chemical system)
%       'X_Mush' (Struct containing the relative mineral proportions and
%       names in the mush system)
%       'C_Mush (Struct containing the composition of each phase specified
%       in X_Mush as 19-by-1 vectors)
%       'X_React' (Struct containing the mineral proportions in the reacted
%       solid assemblage)
%       'liq' (19-by-1 vector containing the initial liquid composition
%       used in the reaction calculations)
%       'deltaT' (offset between the mush temperature and the melt
%       temperature)
%
%   Outputs:
%       'melts' (MELTSdynamic for the reaction path)
%       'MeltMushTables' (Struct containing results tables for the reaction
%       calculations)


    [MeltMushTables]=CreateTables(len);

    MeltMushTables.Liq_Mass.liquid1(1)=phi;

    MeltMushTables.Conditions.P(1)=As_Cond(1);
    MeltMushTables.Conditions.T(1)=As_Cond(2)+deltaT;

    % Mass of minerals in initial unreacted component
    Minerals=fieldnames(X_Mush);
    React=fieldnames(X_React);
    for i=1:length(Minerals)
        MeltMushTables.Unreact_Mass.(Minerals{i})(1)=X_Mush.(Minerals{i})*(1-phi);
        MeltMushTables.Unreact_Composition.(Minerals{i}){1,:}=C_Mush.(Minerals{i})';
        
        MeltMushTables.React_Mass_In.(React{i})(1)=phi*M*X_React.(React{i});
    end

    % Set up MELTS dynamics for the AFC calculations
    assimilate = MELTSdynamic(1);

    % set dynamic for the bulk mush
    mush = MELTSdynamic(1);
    mush.engine.pressure=MeltMushTables.Conditions.P(1);
    mush.engine.temperature=MeltMushTables.Conditions.T(1)-deltaT;

    % set dynamic for the melts calculations
    melts=MELTSdynamic(1);
    melts.engine.pressure=MeltMushTables.Conditions.P(1);
    melts.engine.temperature=MeltMushTables.Conditions.T(1);

    %% calculate initial enthalpy of the unreacted assemblage (at P and T of interest)
    H_Mush=zeros(1,length(Minerals));
    for i=1:length(Minerals)
        mush.engine.calcPhaseProperties(Minerals{i},MeltMushTables.Unreact_Composition.(Minerals{i}){1,:})
        H_Mush(i)=mush.engine.getProperty('h',Minerals{i});
    end
    X_Mu=struct2table(X_Mush);
    H_bulkMush=H_Mush*X_Mu{:,:}';
    H_bulkMush_In=H_bulkMush;

    %% determine composition and enthalpy of the initial melt phase
    % set liquid composition in the melts dynamic
    % normalise bulk composition to 100g
    liq=liq/sum(liq)*100;

    % set bulk composition as initial liquid
    melts.engine.setBulkComposition(liq);

    % calc equilibrium
    melts.engine.calcEquilibriumState(1,0);

    % calc enthaly of system
    Melt_H=melts.engine.h('bulk');
    
    MeltMushTables.React_Composition.liquid1{1,:} = liq;
    MeltMushTables.Conditions.H(1) = Melt_H;
    
    %% enthalpy of assimilated components
    % calculate composition and enthalpy of the initial reacted assemblage
    assimilate.engine.pressure=MeltMushTables.Conditions.P(1);
    assimilate.engine.temperature=MeltMushTables.Conditions.T(1)-deltaT;

    % set composition of the reacted mineral phases
    C_React=struct();
    H_As=zeros(1,length(React));
    for i=1:length(React)
        C_React.(React{i})=zeros(19,1);
        C_React.(React{i})=MeltMushTables.Unreact_Composition.(React{i}){1,:}';
        
        assimilate.engine.calcPhaseProperties(React{i},MeltMushTables.Unreact_Composition.(React{i}){1,:});
        H_As(i)=assimilate.engine.getProperty('h',React{i});
    end

    % convert struct to table to extract data
    C_As=struct2table(C_React);
    X_As=struct2table(X_React);

    % determine bulk composition of reacted component
    C_Assim=C_As{:,:}*X_As{:,:}';

    % determine enthalpy of reacted component
    H_Assim=H_As*X_As{:,:}';

    % ensure enthalpy is relative to 100g
    H_Assim=100.*H_Assim./sum(C_Assim);

    %% for loop time
    for j=1:len
        mush.engine.temperature=As_Cond(2);
        
        % determine bulk enthalpy
        H=phi*(Melt_H)+M*phi*H_Assim+(1-phi-M*phi)*H_bulkMush_In;

        % determine the chemical composition of the local system
        C=(phi/(phi+M*phi))*(liq)+((M*phi)/(phi+M*phi))*C_Assim';

        % carry out AFC calculation
        melts=melts.addNodeAfter;

        % set new bulk 
        melts.engine.setBulkComposition(100*C/sum(C));

        % calc equilibrium state
        melts.engine.calcEquilibriumState(1,0);

        % Indentifying common messages that indicate the code is about to crash!
        tf = isequal(melts.engine.status.message,'Quadratic iterations exceeded.');
        if tf ==1
            break
        end
        tf = isequal(melts.engine.status.message,'MELTS call calcEquilibriumState failed. Cleaning up...');
        if tf ==1
            break
        end
        if melts.engine.getProperty('X',"liquid1","SiO2")<0
            break
        end

        % determine enthalpy of local chemical system
        H1=melts.engine.h("bulk");

        H2=H1*(phi+M*phi)+H_bulkMush_In*(1-phi-M*phi);
        
        % initial adjust to make the system at a single temperature
        mush.engine.temperature=melts.engine.temperature;

        melts.engine.calcEquilibriumState(1,0);
        H1=melts.engine.h("bulk");

        H_Mush=zeros(1,length(Minerals));
        X_Unreact=zeros(1,length(Minerals));
        for i=1:length(Minerals)
            mush.engine.calcPhaseProperties(Minerals{i},MeltMushTables.Unreact_Composition.(Minerals{i}){1,:})
            H_Mush(i)=mush.engine.getProperty('h',Minerals{i});
            X_Unreact(i)=MeltMushTables.Unreact_Mass.(Minerals{i})(1)-M*phi*X_React.(Minerals{i});
        end
        H_bulkMush=sum(H_Mush.*X_Unreact./sum(X_Unreact));

        H2=H1*(phi+M*phi)+H_bulkMush*(1-phi-M*phi);
        
        
        if H2<H
            while H2<H
                melts.engine.temperature=melts.engine.temperature+0.1;
                mush.engine.temperature=mush.engine.temperature+0.1;

                melts.engine.calcEquilibriumState(1,0);
                H1=melts.engine.h("bulk");

                H_Mush=zeros(1,length(Minerals));
                X_Unreact=zeros(1,length(Minerals));
                for i=1:length(Minerals)
                    mush.engine.calcPhaseProperties(Minerals{i},MeltMushTables.Unreact_Composition.(Minerals{i}){1,:})
                    H_Mush(i)=mush.engine.getProperty('h',Minerals{i});
                    X_Unreact(i)=MeltMushTables.Unreact_Mass.(Minerals{i})(1)-M*phi*X_React.(Minerals{i});
                end
                H_bulkMush=sum(H_Mush.*X_Unreact./sum(X_Unreact));

                H2=H1*(phi+M*phi)+H_bulkMush*(1-phi-M*phi);
            end
        else
            while H2>H
                melts.engine.temperature=melts.engine.temperature-0.1;
                mush.engine.temperature=mush.engine.temperature-0.1;
                
                melts.engine.calcEquilibriumState(1,0);
                H1=melts.engine.h("bulk");
                
                H_Mush=zeros(1,length(Minerals));
                X_Unreact=zeros(1,length(Minerals));
                for i=1:length(Minerals)
                    mush.engine.calcPhaseProperties(Minerals{i},MeltMushTables.Unreact_Composition.(Minerals{i}){1,:})
                    H_Mush(i)=mush.engine.getProperty('h',Minerals{i});
                    X_Unreact(i)=MeltMushTables.Unreact_Mass.(Minerals{i})(1)-M*phi*X_React.(Minerals{i});
                end
                H_bulkMush=sum(H_Mush.*X_Unreact./sum(X_Unreact));

                H2=H1*(phi+M*phi)+H_bulkMush*(1-phi-M*phi);
            end
        end

        % calculate new composition and enthalpy of the melt phase
        liq=melts.engine.getProperty('dispComposition','liquid1');
        liq=liq'/(sum(liq))*100;

        % determine entropy of the liquid phase and set to enthalpy per 100g
        Melt_H=melts.engine.h("liquid1")/melts.engine.getProperty('mass','liquid1')*100;

        % assign results to relevant tables
        % set new bulk proportion of liquid
        MeltMushTables.Liq_Mass.liquid1(j+1)=(phi+M*phi)*(melts.engine.getProperty('mass','liquid1')/melts.engine.getProperty('mass','bulk'));

        % set new reacted mass and composition
        for i=1:length(React)
            A=melts.engine.getProperty('mass',React{i});
            A(isnan(A))=0;
            B=melts.engine.getProperty('dispComposition',React{i});
            B(isnan(B))=0;
            MeltMushTables.React_Composition.(React{i}){j+1,:}=B';
            MeltMushTables.React_Mass.(React{i})(j+1)=(phi+M*phi)*A/melts.engine.getProperty('mass','bulk');
        end
        
        MeltMushTables.React_Composition.liquid1{j+1,:} = liq;

        % set conditions
        MeltMushTables.Conditions.T(j+1)=melts.engine.getProperty('temperature');
        MeltMushTables.Conditions.P(j+1)=melts.engine.getProperty('pressure');
        MeltMushTables.Conditions.H(j+1)=Melt_H;
    end
end