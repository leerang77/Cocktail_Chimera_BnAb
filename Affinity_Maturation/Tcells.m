classdef Tcells < handle
%% Header
% Class Tcells provides methods for selection by T cells
% methods summary:
%   Tcells - Initialization
%   Selection - T cell selection of B cells
%   Internalized
%   HelpFunction
% ------------------------------
    %%
    properties
        Constants       % Constants
        TcellsFraction  % 1 by N_Tp array specifying the fraction of T cells
                            % that targets specific epitope. Sum to 1. 
    end
    
    %%
    methods
        %% Initialization
        function obj = Tcells(TcellsFraction)
            if(nargin>0)
                obj.Constants = UseSharedConstants.Constants;
                obj.TcellsFraction = TcellsFraction;
            end
        end
        
        %% T cell selection
        function [p_Sel,SelectionIndex] = Selection(obj, BcellsArray, NonLinearity) %Coarse Grained
            CurrentBcellNum = 0;
            for i=1:length(BcellsArray)
                if ~isempty(BcellsArray(i).Lineage)
                    CurrentBcellNum = CurrentBcellNum+1;
                else
                    break
                end
            end
            TcellHelp = zeros(CurrentBcellNum,obj.Constants.N_Tp);
            p_Sel = zeros(1,CurrentBcellNum);
            SelectionIndex = zeros(1,CurrentBcellNum);
            I_Tp = Internalized(obj, CurrentBcellNum, BcellsArray); %Number of pMHC of each type presented by each B cell 
            I_Tp_Avg = zeros(1,obj.Constants.N_Tp);
            BcellAgCapIdx = cell(1,obj.Constants.N_Tp);
            
            %Get the average amount of pMHC presented
            for i=1:obj.Constants.N_Tp
                n = sum(I_Tp(:,i)~=0);
                total = sum(I_Tp(:,i));
                if(n>0)
                    I_Tp_Avg(1,i) = total/n;
                else
                    I_Tp_Avg(1,i) = 0;
                end
                BcellAgCapIdx{i} = find(I_Tp(:,i)); 
            end
            
            %Get the availability of T cells 
            TcellNum = obj.TcellsFraction*obj.Constants.TotalTcell;
            for i=1:obj.Constants.N_Tp
                if TcellNum(i)>0 && ~isempty(BcellAgCapIdx{i})
                    TcellAvailability(i) = TcellNum(i)/length(BcellAgCapIdx{i});
                else
                    TcellAvailability(i) = 0;
                end
            end
            
            %Get the amount of T cell help and selection probability
            for i=1:CurrentBcellNum
                for j=1:obj.Constants.N_Tp
                    if(I_Tp(i,j)>0)
                        TcellHelp(i,j) = TcellAvailability(j)*HelpFunction(I_Tp(i,j)/I_Tp_Avg(1,j), NonLinearity);
                    else
                        TcellHelp(i,j) = 0;
                    end
                end
                p_Sel(i) = sum(TcellHelp(i,:))/(1+sum(TcellHelp(i,:)));
            end
            p_Sel = obj.Constants.pmax*p_Sel;
                      
            %Selection
            idx = 0;
            for i=1:CurrentBcellNum
                if rand<p_Sel(1,i)
                    idx = idx+1;
                    SelectionIndex(idx) = i;
                end
            end
            SelectionIndex = SelectionIndex(SelectionIndex>0);
            
            %%
            function I_Tp = Internalized(obj, CurrentBcellNum, BcellsArray)
                % Determines the Langmuir constant and Mean contact time
                I_Ag = zeros(CurrentBcellNum, obj.Constants.N_Ag); %Number of Internalized Antigens
                for j=1:CurrentBcellNum
                    AgCaptured = BcellsArray(j).AgCaptured;
                    I_Ag(j,:) = AgCaptured;
                end
                I_Tp = I_Ag*obj.Constants.AgtoTp;
            end

            %%    
            function Help = HelpFunction(AgCapturedRatio, NonLinearity)
                Help = AgCapturedRatio^NonLinearity;
            end
        end
    end
end
