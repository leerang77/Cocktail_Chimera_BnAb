classdef SharedConstants < handle
    properties
        %constants related to GC dynamics
        N_Seed;             %Number of initial seeding B cells
        N_Max;              %Max GC capacity
        GC_Length;          %Number of GC cycles before termination
        PrecursorFreq;      %1 by #of specific epitope; Fraction of germinal B cells targeting epitope i
        
        %constants related to affinities
        E_initial;          %Initial affinity of B cells       
        Sigma;              %Covariance matrix
        MutationPDF;        %PDF function for mutation
        
        %constants related to antigen
        AgtoEp;             %N_Ag by N_Ep array: AgtoEp(i,k) = number of epitopes k
                            %present on antigen i
        N_Ag;               %Number of different antigen types
        N_Ep;               %Number of different epitope types
        AgType;             %string; 'chimeric' or 'cocktail'
        Ep_specific;        %Indices of strain-specific epitopes
        Ep_bnab;            %Indices of bnab epitopes
        
        %constants related to T cells
        N_Tp;               % Number of different types of T epitopes
        AgtoTp;             % N_Ag by N_Tp array: Number of each epitope in each antigen
        TotalTcell;         % Total availability of T cells 
        pmax;               % Maximum positive selection probability
        pdeath;             % Death probability for unselected B cells
        pdiff;              % Differentiation probability for selected B cells        
        
        %variables
        BenMuts;            %Cell array that stores all beneficial mutations that occur
        NumBenMut;          %Number of all beneficial mutations that occur
    end
end
