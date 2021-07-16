function InitializeSharedConstants(AgType,rho,bnab,h)
%% Header
%  This function is used to initialize parameters of the simulation. 
%  See Class SharedConstants for description.

%constants related to GC dynamics
h.N_Seed = 100;
h.N_Max = 4000;
h.GC_Length = 28;
if bnab==0
    h.PrecursorFreq = [1,0,0,0];
else
    h.PrecursorFreq = [0.33, 0.33, 0.33, 0.01];
end

%constants related to affinities
h.E_initial = -13.8;
h.Sigma = [1, rho, rho; rho, 1, rho; rho, rho, 1];
h.MutationPDF = [3.1, 1.4, 2.3]; %lognormal - mu, sigma, translation

%Constants related to antigen
h.AgType = AgType;
h.Ep_specific = [1,2,3];        %Indices of strain-specific epitopes
h.Ep_bnab = 4;            %Indices of bnab epitopes
h.N_Ep = 4;
if strcmp(AgType, 'cocktail')
    h.N_Ag = 3;             
    h.AgtoEp = [3, 0, 0, 3; 
    0, 3, 0, 3;
    0, 0, 3, 3;];
    h.AgtoTp = [3, 0, 0, 3; 
        0, 3, 0, 3;
        0, 0, 3, 3;];
elseif strcmp(AgType, 'chimeric')
    h.N_Ag = 1;
    h.AgtoEp = [1, 1, 1, 3];
    h.AgtoTp = [1, 1, 1, 3];
end

%constants related to T cells
h.N_Tp = 4;
h.TotalTcell = 1500;
h.pmax = 0.6;
h.pdeath = 1;
h.pdiff = 0.1;

%variables
h.BenMuts = cell(1,2000);
h.NumBenMut = 0;
end