classdef Bcell
% Instance of a GC B cell
% Methods:
%   Bcell - Initialization
%   AgCaptureBcell - Antigen capture by B cell
%   AgCaptureBcellAnalytic - Antigen capture by B cell using analytical
%                            formula. This function is called by AgCaptureBcell when needed.
%   Mutate - Mutate B cell
% --------------------------------------------
    %%
    properties
        %Constants 
        Constants;
        
        %Properties of B Cell
        Lineage;            %Clone lineage, between 1 and 100
        TargetEp;           %Epitope the B cell is targeting, between 1 and 4
        E_AgAb;             %Affintiy between B cell and target Ag: Scalar(Strain-specific) or Vector(BnAb)
        AgCaptured;         %1 by N_Ag array: number of antigen i captured
        NumMut;             %Number of total mutations
        NumBenMut;          %Number of affinity-increasing mutations
        NumAffMut;          %Number of affinity-affecting mutations
        MutHistory;         %History of mutation
    end   
    
    %%
    methods
        %% Initialization
        function obj = Bcell(Lineage, TargetEp)
            if(nargin>0)
                obj.Constants = UseSharedConstants.Constants;
                obj.Lineage = Lineage;
                obj.TargetEp = TargetEp;
                if ismember(obj.TargetEp, obj.Constants.Ep_specific)
                    obj.E_AgAb = obj.Constants.E_initial;
                else
                    obj.E_AgAb = [obj.Constants.E_initial,obj.Constants.E_initial,obj.Constants.E_initial];
                end
                obj.NumMut = 0;
                obj.NumAffMut = 0;
                obj.NumBenMut = 0;
                obj.MutHistory = zeros(1,10);
            else %default initialization
                
            end
        end
        
        %% Ag Capture
        function obj = AgCaptureBcell(obj, W)
            if ~exist('W', 'var')
                W = 0;
            end
            if W<=0 %Explicit Ag capture
                obj.AgCaptured = zeros(1,obj.Constants.N_Ag); %Number of each captured antigen
                if ismember(obj.TargetEp, obj.Constants.Ep_specific) % Strain-specific B cell
                    if obj.Constants.N_Ag == 1 %Chimeric
                        obj.AgCaptured = E2Ag(4, obj.E_AgAb, 0, 0);
                    elseif obj.Constants.N_Ag == 3 %Cocktail
                        Energy = [0, 0, 0];
                        Energy(obj.TargetEp) = obj.E_AgAb;
                       obj.AgCaptured = E2Ag(1, Energy(1), Energy(2), Energy(3)); 
                    end
                else % BnAb precursor
                    if obj.Constants.N_Ag == 1 %Chimeric
                        obj.AgCaptured = E2Ag(2, obj.E_AgAb(1), obj.E_AgAb(2), obj.E_AgAb(3));
                    elseif obj.Constants.N_Ag == 3 %Cocktail
                        obj.AgCaptured = E2Ag(3, obj.E_AgAb(1), obj.E_AgAb(2), obj.E_AgAb(3));                    
                    end
                end
            else % Analytical Ag Capture
                obj = AgCaptureBcellAnalytic(obj, W);
            end
        end
        
        %% Analytical Ag Capture
        function obj = AgCaptureBcellAnalytic(obj, W)
            obj.AgCaptured = zeros(1,obj.Constants.N_Ag); %Number of each captured antigen
            if ismember(obj.TargetEp, obj.Constants.Ep_specific) %If specific B cell
                if obj.E_AgAb <= obj.Constants.E_initial %Capture Ag only if above this affinity
                    obj.AgCaptured(obj.TargetEp) = exp(-W*(obj.E_AgAb-obj.Constants.E_initial));
                end
            else %If bnAb B cell
                for j=1:obj.Constants.N_Ag %for each epitope
                   if obj.E_AgAb(j) <= obj.Constants.E_initial
                       obj.AgCaptured(j) = exp(-W*(obj.E_AgAb(j)-obj.Constants.E_initial));
                   end
                end
            end
        end
        
        %% Mutate
        function obj = Mutate(obj)
            if rand*0.7 < 0.5 %silent mutation
                obj.NumMut = obj.NumMut + 1; 
            else %nonsynonymous mutation
                if ismember(obj.TargetEp, obj.Constants.Ep_specific) %if specific B cell
                    obj.NumMut = obj.NumMut + 1; % Increas number of total mutations
                    obj.NumAffMut = obj.NumAffMut + 1; %Affinity-affecting mutations
                    dE = getAffinityChange('specific');
                    if dE<0 %If mutation is affinity-increasing
                        obj.NumBenMut = obj.NumBenMut + 1;
                        obj.Constants.NumBenMut = obj.Constants.NumBenMut+1;
                        obj.Constants.BenMuts{obj.Constants.NumBenMut} = dE;
                        obj.MutHistory(sum(obj.NumBenMut)) = obj.Constants.NumBenMut;
                    end                    
                    obj.E_AgAb = obj.E_AgAb + dE;
                else %if bnAb B cell
                    obj.NumMut = obj.NumMut+1;
                    dE = getAffinityChange('bnab');
                    obj.E_AgAb = obj.E_AgAb + dE;
                    obj.NumAffMut = obj.NumAffMut + 1;
                    if any(dE(1,:)<0)
                        obj.NumBenMut = obj.NumBenMut + 1;
                        obj.Constants.NumBenMut = obj.Constants.NumBenMut+1;
                        obj.Constants.BenMuts{obj.Constants.NumBenMut} = dE;
                        obj.MutHistory(sum(obj.NumBenMut)) = obj.Constants.NumBenMut;
                    end
                end    
            end     
        end
    end
end