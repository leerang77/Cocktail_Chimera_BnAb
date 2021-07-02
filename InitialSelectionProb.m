X = linspace(0.4,1.5,50);
InitialProb = zeros(2,50);
for AgTypeIdx = 1:2
    InitialProb(AgTypeIdx,:) = getInitialSelectionProb(AgTypeIdx,0.7,0,X,0,1);
end
dlmwrite('InitialBnAbSelectionProb.csv', [X;InitialProb]);

function p = getInitialSelectionProb(AgTypeIdx,rho,SharedTcellFraction,X,W,bnab)
    %% Run a single GC simulation
    % Input
    %X - Float, Stringency of T cell selection

    % Output
    %GCStat - Struct with following fields:
    %   BenMuts: All affinity-increasing mutations that occur in GC
    %   BnAbMutations: BnAb mutations that become prevalent in the population
    %   BenMutsByTime: Number of affinity-increasing mutations by time
    %   SpecificMutations: Specific mutations that become prevalent in the population
    %   BcellsByEp: Number of B cells targeting each epitope by time
    %   SurvivingLineages: Number of surviving lineages by time
    %   DomOccupancy: Occupancy of the most dominant lineage by time 
    %   SpecificAff: Mean affinity and number of positively selected strain-specific B cells 
    %   BnAbAff: Mean affinity and number of positively selected strain-specific B cells
    %   BnAbAffAll: Array of affinities of bnAb precursors
    %   SpecificNumMut: Average number of mutations (total, affinity affecting) on SS B cells
    %   BnAbNumMut: Average number of mutations (total, affinity affecting) on BnAb precursors

    %--------------------------------
    Constants = UseSharedConstants.Constants; %get handle for the shared constants
    if AgTypeIdx==1, AgType='chimeric';
    elseif AgTypeIdx==2, AgType='cocktail';
    end
    
    if ~exist('W','var') %If W is not passed, then explicit Ag capture
        W=0; bnab=1;
    end
    InitializeSharedConstants(AgType,rho,bnab,Constants); %Initialize the shared constants
    
    %Define the similarity of the T cell epitopes
    SpecificTcellFraction = (1-SharedTcellFraction);
    TcellsFraction = SpecificTcellFraction*[1/3,1/3,1/3,0]+[0,0,0,SharedTcellFraction];
    TfhC = Tcells(TcellsFraction);
    
    %Initializations
    Constants.BenMuts = cell(1,2000); % for storing all affinity-increasing mutations
    Constants.NumBenMut = 0;
    totalBcells(Constants.GC_Length,Constants.N_Max) = Bcell(); % for B cells at each time
    selectedBcells(Constants.GC_Length,Constants.N_Max/2) = Bcell(); % for selected B cells at each time
    PlasmaArray(1,4000) = Bcell(); % for differentiated plasma cells
    BcellsArray(1,Constants.N_Max) = Bcell(); % array of B cells currently in the GC 
    CurrentBcellNum = 0; % number of non-empty instances in BcellsArray 
    CurrentPlasmaNum = 0;% number of non-empty instances in PlasmaArray 
    BenMutsNumByTime = zeros(1,Constants.GC_Length);
    BnAbSelectionProb = cell(1,Constants.GC_Length);
    %Determine the targeting epitopes of seeding B cells
    %specific: 1,2,3 and bnab: 4
    BnAbPrecursorNum = round(Constants.N_Seed * ...
                       Constants.PrecursorFreq(1,Constants.Ep_bnab));
    SpecificPrecursorNum = zeros(1,Constants.N_Ep-1);
    r = rand(1,Constants.N_Seed-BnAbPrecursorNum)...
        *(1-Constants.PrecursorFreq(1,Constants.Ep_bnab));
    cmf = cumsum(Constants.PrecursorFreq);
    for i=Constants.Ep_specific
        idx = r<cmf(i);
        SpecificPrecursorNum(1,i) = sum(idx);
        r(idx) = 1;
    end
    TargetEpArray = zeros(1,Constants.N_Seed);
    TargetEpArray(1:BnAbPrecursorNum) = Constants.Ep_bnab;
    for i=Constants.Ep_specific
        idx = find(TargetEpArray==0,1);
        TargetEpArray(idx:idx+SpecificPrecursorNum(i)-1) = i;
    end
    clear BnAbPrecursorNum SpecificPrecursorNum r cmf idx

    % Seeding the GC
    for i=1:Constants.N_Seed
        BcellsArray(1,i) = Bcell(i,TargetEpArray(i));
    end
    clear TargetEpArray
    CurrentBcellNum = Constants.N_Seed;
    for RepNum=1:4
        BcellsArray(CurrentBcellNum+1:CurrentBcellNum*2) = BcellsArray(1:CurrentBcellNum);
        CurrentBcellNum = CurrentBcellNum*2;
    end

% Competitive Phase
    for j=1:CurrentBcellNum
        BcellsArray(j) = AgCaptureBcell(BcellsArray(j),W);
    end
    %Help and Selection
    p = zeros(1,length(X));
    for i=1:length(X)
        [pSel, SelectionIdx] = Selection(TfhC, BcellsArray, X(i));
        p(i) = pSel(1);
    end
end