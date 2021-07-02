function runMultipleGCCombined_edit(GCNum,AgTypeIdx,rho,SharedTcellFraction,X,W,bnab)
%% HEADER
 %This code can run affinity maturation simulations with both the analytical antigen capture 
 %and cocktail/chimeric antigen capture.
 %To run analytical antigen capture, use AgTypeIdx=2 and SharedTcellFraction=1.
 %To run chimeric antigen, use AgTypeIdx=1, SharedTcellFraction=1, and do
 %not pass W and bnab.
 %To run cocktail antigen, use AgTypeIdx=2, SharedTcellFraction=0, and do
 %not pass W and bnab.
 %The code runs number of GCs specified by GCNum, and derive statistics.
 
 % Inputs
 %GCNum - Int, Total number of GCs (i.e. number of independent repeats)
 %AgTypeIdx - If 1, chimeric antigen
 %            If 2, cocktail antigen
 %rho       - Float between 0 and 1, Antigenic distance between the different strains
 %            Specifically, the mutation is derived from a 3-D Gaussian
 %            with rho as the correlation and then converted to shifted
 %            log-normal
 %SharedTcellFraction - Float between 0 and 1, The fraction of T cell epitopes shared among the antigens
 %E0 - Float, The initial affinity of the B cells
 %X  - Float, Stringency of T cell selection in the GC
 %W (only for analytical Ag capture) - Float, Stringency of Ag capture in the GC
 %bnab (only for analytical Ag capture) - 0 or 1, Whether there is 0 or 1 bnAb
 %                                        precursor among the founders

 %  Outputs
 %Saves the following statistics in a *.mat file.
 %SingleGCStat - Statistics from a single GC run. See the function
 %               runSingleGC for more description.
 %TotalBcellsByEp - GC_Length x N_Ep array, Number of B cells targeting
 %                  each epitope after each cycle
 %MeanSurvivingLineages - GC Length x 1 array,
 %MeanDomOccupancy - GC Length x 1 array,
 %MeanSpecificAff - GC Length x 1 array,
 %MeanBnAbAff - GC Length x 4 array,
 %MeanSpecificNumMut - GC Length x 2 array,
 %MeanBnAbNumMut - GC Length x 3 array,
 %%  -------- Main Script --------
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
    
    %Run the simulations 
    SingleGCStat = cell(1,GCNum);
    for gcnumber=1:GCNum
        disp(gcnumber);
        rng(gcnumber);
        SingleGCStat{gcnumber} = runSingleGC(X);
    end

    %Analyze to derive mean statistics
    analyzeMultipleGC(W)
    
%% Subfunctions
    function analyzeMultipleGC(W)
     %% Use the results of individual GCs to calculate overall statistics
     % Input
     %W - Ag capture stringency. If 0, explicit Ag capture; If >0,
     %    analytical Ag capture
     
     % Output
     %None, saves the result into *mat file
     %----------------------------------------
        TotalBcellsByEp = zeros(Constants.GC_Length, Constants.N_Ep);
        MeanSurvivingLineages = zeros(Constants.GC_Length,1);
        MeanDomOccupancy = zeros(Constants.GC_Length,1);
        MeanSpecificAff = zeros(Constants.GC_Length,1);
        MeanBnAbAff = zeros(Constants.GC_Length,3);
        MeanSpecificNumMut = zeros(Constants.GC_Length,2);
        MeanBnAbNumMut = zeros(Constants.GC_Length,2);
        totalspecific = zeros(Constants.GC_Length,1);
        totalbnab = zeros(Constants.GC_Length,1);
        Detectable = zeros(Constants.GC_Length,2);
        DetectableSelected = zeros(Constants.GC_Length,2);
        Day3BnAbSelectionProb = cell(1,3);
        Day5BnAbSelectionProb = cell(1,3);
        Day10BnAbSelectionProb = cell(1,3);
        Day12BnAbSelectionProb = cell(1,3);
        % Calculation of mean statistics
        % For affinity and number of mutations, B cells are pooled from all GCs then mean is calculated
        
        for i=1:GCNum
            for j=1:3
                Day3BnAbSelectionProb{j}{i} = SingleGCStat{i}.BnAbSelectionProb{6}{j};
                Day5BnAbSelectionProb{j}{i} = SingleGCStat{i}.BnAbSelectionProb{10}{j};
                Day10BnAbSelectionProb{j}{i} = SingleGCStat{i}.BnAbSelectionProb{20}{j};               
                Day12BnAbSelectionProb{j}{i} = SingleGCStat{i}.BnAbSelectionProb{24}{j};
            end
            TotalBcellsByEp = TotalBcellsByEp + SingleGCStat{i}.BcellsByEp;
            MeanSurvivingLineages = MeanSurvivingLineages + SingleGCStat{i}.SurvivingLineages/GCNum;
            MeanDomOccupancy = MeanDomOccupancy + SingleGCStat{i}.DomOccupancy/GCNum;
            totalSpecificAff = MeanSpecificAff + ...
                               SingleGCStat{i}.SpecificAff(:,1).*SingleGCStat{i}.SpecificAff(:,2);
            totalspecific = totalspecific + SingleGCStat{i}.SpecificAff(:,2);
            totalBnAbAff = MeanBnAbAff + ...
                           min(SingleGCStat{i}.BnAbAff(:,1:3),-11).*SingleGCStat{i}.BnAbAff(:,4);
            totalbnab = totalbnab + SingleGCStat{i}.BnAbAff(:,4);
            totalSpecificNumMut = MeanSpecificNumMut + ...
                                  SingleGCStat{i}.SpecificNumMut(:,1:2).* ...
                                  SingleGCStat{i}.SpecificNumMut(:,3);
            totalBnAbNumMut = MeanBnAbNumMut + ...
                              SingleGCStat{i}.BnAbNumMut(:,1:2).*SingleGCStat{i}.BnAbNumMut(:,3);
            Detectable = Detectable + SingleGCStat{i}.Detectable;
            DetectableSelected = DetectableSelected + SingleGCStat{i}.DetectableSelected;
        end
        for j=1:3
            Day3BnAbSelectionProb{j} = [Day3BnAbSelectionProb{j}{:}];
            Day5BnAbSelectionProb{j} = [Day5BnAbSelectionProb{j}{:}];
            Day10BnAbSelectionProb{j} = [Day10BnAbSelectionProb{j}{:}];
            Day12BnAbSelectionProb{j} = [Day12BnAbSelectionProb{j}{:}];
        end
        MeanSpecificAff = totalSpecificAff./totalspecific;       
        MeanBnAbAff = totalBnAbAff./totalbnab;       
        MeanSpecificNumMut = totalSpecificNumMut./totalspecific; 
        MeanBnAbNumMut = totalBnAbNumMut./totalbnab;
        MeanSpecificAff(isnan(MeanSpecificAff))=0; %change nan to 0
        MeanBnAbAff(isnan(MeanBnAbAff))=0;
        MeanSpecificNumMut(isnan(MeanSpecificNumMut))=0;
        MeanBnAbNumMut(isnan(MeanBnAbNumMut))=0;
        
        % Define filename for output
        if bnab==1, prefix = 'BnAb';
        elseif bnab==0, prefix = 'NoBnAb';
        end
        if W>0 % Analytical Ag capture
            filename = strcat(prefix,sprintf('Analytical_rho_%.1f_E0_%.1f_X_%.2f_W_%.2f_f_%.1f.mat'...
                ,rho,Constants.E_initial,X,W,SharedTcellFraction));
        else % Explicit Ag capture
            filename = strcat(AgType,sprintf('Explicit_rho_%.1f_E0_%.1f_X_%.2f_f_%.1f.mat',...
                rho,Constants.E_initial,X,SharedTcellFraction));
        end
        [~,fn,~] = fileparts(filename);
        for j=1:3
            dlmwrite(strcat(fn,'_Day3BnAbSelectionProb.csv'), Day3BnAbSelectionProb{j}, '-append')
            dlmwrite(strcat(fn,'_Day5BnAbSelectionProb.csv'), Day5BnAbSelectionProb{j}, '-append')
            dlmwrite(strcat(fn,'_Day10BnAbSelectionProb.csv'), Day10BnAbSelectionProb{j}, '-append')
            dlmwrite(strcat(fn,'_Day12BnAbSelectionProb.csv'), Day12BnAbSelectionProb{j}, '-append')
        end
        save(filename,'SingleGCStat','TotalBcellsByEp','MeanSurvivingLineages',...
        'MeanDomOccupancy','MeanSpecificAff','MeanBnAbAff','MeanSpecificNumMut','MeanBnAbNumMut',...
        'Detectable','DetectableSelected')
    end

    function GCStat = runSingleGC(X)
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
        for i=1:Constants.GC_Length
            %Antigen Capture
            for j=1:CurrentBcellNum
                BcellsArray(j) = AgCaptureBcell(BcellsArray(j),W);
            end
            %Help and Selection
            [pSel, SelectionIdx] = Selection(TfhC, BcellsArray, X);
            BnAbSelectionProb{i} = getBnAbSelectionProb(pSel, BcellsArray);
            selectedBcells(i,1:length(SelectionIdx)) = BcellsArray(SelectionIdx);
            survivalidx = setdiff(1:CurrentBcellNum, SelectionIdx);
            survivalidx = survivalidx(rand(1,length(survivalidx))<(1-Constants.pdeath));
            %Export 10% of selected cells as plasma cells
            PlasmaIdx = zeros(1, length(SelectionIdx));
            for j=1:length(SelectionIdx)
                if rand<Constants.pdiff
                    PlasmaIdx(j) = SelectionIdx(j);
                end
            end
            PlasmaIdx = PlasmaIdx(PlasmaIdx~=0);
            PlasmaArray(CurrentPlasmaNum+1:CurrentPlasmaNum+length(PlasmaIdx)) = BcellsArray(PlasmaIdx);
            CurrentPlasmaNum = CurrentPlasmaNum + length(PlasmaIdx);
            %Apoptosis of unselected cells
            NextGenBcells(1,Constants.N_Max) = Bcell();
            NextGenBcells(1:length(survivalidx)) = BcellsArray(survivalidx);
            NextGenBcellNum = length(survivalidx);
            %Selected cells proliferate and mutate
            %first division
            ParentIdx = setdiff(SelectionIdx, PlasmaIdx);
            DaughterIdx = ParentIdx(rand(1,length(ParentIdx))<0.7); % Mutation leads to death 30% of time
            n1 = length(ParentIdx); n2 = length(DaughterIdx);
            NextGenBcells(NextGenBcellNum+1:NextGenBcellNum+n1) = BcellsArray(ParentIdx);
            NextGenBcells(NextGenBcellNum+n1+1:NextGenBcellNum+n1+n2) = BcellsArray(DaughterIdx);
            for j=NextGenBcellNum+n1+1:NextGenBcellNum+n1+n2 % One of the two cells mutate after division
                NextGenBcells(j) = Mutate(NextGenBcells(j));
            end
            clear ParentIdx DaughterIdx
            %second division
            ParentIdx = NextGenBcellNum+1:NextGenBcellNum+n1+n2;
            DaughterIdx = ParentIdx(rand(1,length(ParentIdx))<0.7); % Mutation leads to death 30% of time
            NextGenBcellNum = NextGenBcellNum+n1+n2;
            n1 = length(ParentIdx); n2 = length(DaughterIdx);
            NextGenBcells(NextGenBcellNum+1:NextGenBcellNum+n1) = NextGenBcells(ParentIdx);
            NextGenBcells(NextGenBcellNum+n1+1:NextGenBcellNum+n1+n2) = NextGenBcells(DaughterIdx);
            for j=NextGenBcellNum+n1+1:NextGenBcellNum+n1+n2 % One of the two cells mutate after division
                NextGenBcells(j) = Mutate(NextGenBcells(j));
            end
            BcellsArray = NextGenBcells;
            clear ParentIdx DaughterIdx n1 n2 NextGenBcells
            CurrentBcellNum = 0;
            for j=1:length(BcellsArray) % Get the number of GC B cells
                if ~isempty(BcellsArray(j).Lineage)
                    CurrentBcellNum = CurrentBcellNum+1;
                else
                    break
                end
            end
            BenMutsNumByTime(i) = Constants.NumBenMut;
            totalBcells(i,:) = BcellsArray;
        end
        for i=1:length(PlasmaArray)
            if isempty(PlasmaArray(i).TargetEp)
                break
            else
                PlasmaArray(i).Constants=[];
            end
        end
        i=i-1;
        PlasmaArray = PlasmaArray(1:i);
        % Analysis
        [BcellsByEp, Detectable, SurvivingLineages, DomOccupancy, BnAbAffAll] = getBcellsByEp(totalBcells);
        [SpecificAff,BnAbAff,BnAbAffSelected,DetectableSelected,SpecificNumMut,BnAbNumMut,BnAbMutations,SpecificMutations]...
            = getBcellsAff(selectedBcells);
        GCStat.BnAbSelectionProb = BnAbSelectionProb;
        GCStat.BenMuts = Constants.BenMuts;
        GCStat.BnAbMutations = BnAbMutations;
        GCStat.BenMutsByTime = BenMutsNumByTime;
        GCStat.SpecificMutations = SpecificMutations;
        GCStat.BcellsByEp = BcellsByEp;
        GCStat.SurvivingLineages = SurvivingLineages;
        GCStat.DomOccupancy = DomOccupancy;
        GCStat.SpecificAff = SpecificAff;
        GCStat.BnAbAff = BnAbAff;
        GCStat.BnAbAffAll = BnAbAffAll;
        GCStat.BnAbAffSelected = BnAbAffSelected;
        GCStat.SpecificNumMut = SpecificNumMut;
        GCStat.BnAbNumMut = BnAbNumMut;
        GCStat.Detectable = Detectable;
        GCStat.DetectableSelected = DetectableSelected;
        GCStat.PlasmaArray = PlasmaArray;
    end

    function bnAbSelectionProb = getBnAbSelectionProb(pSel, BcellsArray)
        bnAbSelectionProb = cell(1,3);
        for i=1:3
            bnAbSelectionProb{i} = zeros(1,length(BcellsArray));
        end
        for i=1:length(pSel)
            if BcellsArray(i).TargetEp==4
                breadth = sum(BcellsArray(i).E_AgAb<-13.8);
                if breadth>0
                    bnAbSelectionProb{breadth}(i) = pSel(i);
                end
            end
        end
        for i=1:3
            bnAbSelectionProb{i} = bnAbSelectionProb{i}(bnAbSelectionProb{i}~=0);
            bnAbSelectionProb{i} = unique(bnAbSelectionProb{i});
        end
    end

    function [BcellsByEp, Detectable, SurvivingLineages, DomOccupancy, BnAbAffAll] = getBcellsByEp(BcellsArray)
        %% Analyze the epitope specificity and clonal lineages of GC B cells
        % Input
        %BcellsArray: GC Length x N_max array. Contains GC B cells by time
        
        % Outputs
        %BcellsByEp: GC Length x N_Ep array. Number of B cells targeting each epitope
        %SurvivingLineages: GC Length x 1 array. Number of surviving lineages
        %DomOccpancy: GC Length x 1 array. Occupancy of the dominant clone
        %-------------------------

      
        [gclength,number] = size(BcellsArray);
        BnAbAffAll = cell(1,gclength);  
        BcellsByEp = zeros(gclength,Constants.N_Ep);
        Detectable = zeros(gclength,2);
        SurvivingLineages = zeros(gclength,1);
        DomOccupancy = zeros(gclength,1);
        for dim1 = 1:gclength
            BnAbAffAll{dim1} = zeros(number,3);
            Lineages = zeros(1,Constants.N_Seed);
            for dim2 = 1:number
                if isempty(BcellsArray(dim1,dim2).TargetEp)
                    BnAbAffAll{dim1} = BnAbAffAll{dim1}(sum(BnAbAffAll{dim1},2)~=0,:);
                    break
                end
                Ep = BcellsArray(dim1,dim2).TargetEp;
                if any(BcellsArray(dim1,dim2).E_AgAb<=-13.8)
                    Detectable(dim1,1) = Detectable(dim1,1)+1;
                end
                if any(BcellsArray(dim1,dim2).E_AgAb<-13.8)
                    Detectable(dim1,2) = Detectable(dim1,2)+1;
                end
                BcellsByEp(dim1,Ep) = BcellsByEp(dim1,Ep)+1;
                if ~ismember(Ep,Constants.Ep_specific)
                    E = BcellsArray(dim1,dim2).E_AgAb;
                    E_var = sort(E(1,:),'descend');
                    BnAbAffAll{dim1}(dim2,1:3) = E_var;
                end
                Lineages(BcellsArray(dim1,dim2).Lineage) = Lineages(BcellsArray(dim1,dim2).Lineage)+1;
            end
            SurvivingLineages(dim1,1) = sum(Lineages>0);
            DomOccupancy(dim1,1) = max(Lineages)/sum(Lineages);
        end
    end

    function [SpecificAff,BnAbAff,BnAbAffAll,Detectable,SpecificNumMut,BnAbNumMut,BnAbMutations,SpecificMutations] = getBcellsAff(BcellsArray)
        %% Analyze the affinities and mutations of selected B cells
        %SpecificAff: GC Length x 2 array. First column: Mean affinity,
        %    Second column: Cells Number
        %BnAbAff: GC Length x 4 array. Columns 1-3: Mean affinity towards epitopes 1-3, sorted in
        %    decreasing order. Column 4: Cells Number
        %BnAbAffAll: 1xGC Length cell array. Each cell is N_max x 4 array.
        %    Each row of this array contains the affinity of a bnAb
        %    precursor.
        %SpecificNumMut: GCLength x 3 array. First column: Mean total number
        %    of mutations. Second column: Mean total number of
        %    affinity-affecting mutations. column 3: Cells Number
        %BnAbNumMut: GCLength x 3 array. First column: Mean total number of
        %    mutations. Second column: Mean total number of affinity-affecting
        %    mutations(conserved region). Second column: Mean total number of
        %    affinity-affecting mutations(variable region). Column 3: Number of B Cells
        %BnAbMutations: 2 x GCLength cell. BnAb mutations that are prevalent.
        %    Each row corresponds to a different definition of "prevalence". 
        %SpecificMutations: 2 x GCLength cell. Strain-specific mutations that are prevalent.
        [gclength,number] = size(BcellsArray);
        SpecificAff = zeros(gclength,2);
        BnAbAff = zeros(gclength,4);
        SpecificNumMut = zeros(gclength,3);
        BnAbNumMut = zeros(gclength,3);
        BnAbMutations = cell(2,gclength);
        SpecificMutations = cell(2,gclength);
        BnAbAffAll = cell(1,gclength);
        Detectable = zeros(gclength,2);
        for dim1 = 1:gclength
            BnAbAffAll{dim1} = zeros(number,3);
            BnAbMutations{1,dim1} = zeros(number,10);
            BnAbMutations{2,dim1} = zeros(number,10);
            SpecificMutations{1,dim1} = zeros(number,10);
            SpecificMutations{2,dim1} = zeros(number,10);
            bnabcnt = 0;
            sscnt = 0;
            for dim2 = 1:number
                if isempty(BcellsArray(dim1,dim2).TargetEp)
                    BnAbAffAll{dim1} = BnAbAffAll{dim1}(sum(BnAbAffAll{dim1},2)~=0,:);
                    break
                end
                if any(BcellsArray(dim1,dim2).E_AgAb<=-13.8)
                    Detectable(dim1,1) = Detectable(dim1,1)+1;
                end
                if any(BcellsArray(dim1,dim2).E_AgAb<-13.8)
                    Detectable(dim1,2) = Detectable(dim1,2)+1;
                end                
                Ep = BcellsArray(dim1,dim2).TargetEp;
                if ismember(Ep,Constants.Ep_specific)
                    sscnt = sscnt + 1;
                    SpecificAff(dim1,1) = SpecificAff(dim1,1) + BcellsArray(dim1,dim2).E_AgAb;
                    SpecificAff(dim1,2) = SpecificAff(dim1,2) + 1;
                    SpecificNumMut(dim1,1) = SpecificNumMut(dim1,1) + BcellsArray(dim1,dim2).NumMut;
                    SpecificNumMut(dim1,2) = SpecificNumMut(dim1,2) + sum(BcellsArray(dim1,dim2).NumAffMut);
                    SpecificNumMut(dim1,3) = SpecificNumMut(dim1,3) + 1;
                    SpecificMutations{1,dim1}(dim2,:) = BcellsArray(dim1,dim2).MutHistory;                    
                else
                    bnabcnt = bnabcnt + 1;
                    E = BcellsArray(dim1,dim2).E_AgAb;
                    E_var = sort(E(1,:),'descend');
                    BnAbAff(dim1,1:3) = BnAbAff(dim1,1:3) + E_var;
                    BnAbAff(dim1,4) = BnAbAff(dim1,4) + 1;
                    BnAbAffAll{dim1}(dim2,1:3) = E_var;
                    BnAbMutations{1,dim1}(dim2,:) = BcellsArray(dim1,dim2).MutHistory;
                    BnAbNumMut(dim1,1) = BnAbNumMut(dim1,1) + BcellsArray(dim1,dim2).NumMut;
                    BnAbNumMut(dim1,2) = BnAbNumMut(dim1,2) + BcellsArray(dim1,dim2).NumAffMut;
                    BnAbNumMut(dim1,3) = BnAbNumMut(dim1,3) + 1;
                end
            end
            BnAbAffAll{dim1} = BnAbAffAll{dim1}(sum(BnAbAffAll{dim1},2)~=0,:);
            SpecificAff(dim1,1) = SpecificAff(dim1,1)/SpecificAff(dim1,2);
            SpecificAff(isnan(SpecificAff))=0;
            BnAbAff(dim1,1:3) = BnAbAff(dim1,1:3)/BnAbAff(dim1,4);
            BnAbAff(isnan(BnAbAff))=0;
            SpecificNumMut(dim1,1:2) = SpecificNumMut(dim1,1:2)/SpecificNumMut(dim1,3);
            SpecificNumMut(isnan(SpecificNumMut))=0;
            BnAbNumMut(dim1,1:2) = BnAbNumMut(dim1,1:2)/BnAbNumMut(dim1,3);
            BnAbNumMut(isnan(BnAbNumMut))=0;
            [cnt_unique, unique_mut] = hist(BnAbMutations{1,dim1}(BnAbMutations{1,dim1}>0)...
                ,unique(BnAbMutations{1,dim1}(BnAbMutations{1,dim1}>0)));
            prevalent = cnt_unique/bnabcnt > 0.5;          
            BnAbMutations{1,dim1} = unique_mut(prevalent);
            prevalent = cnt_unique > 100;
            BnAbMutations{2,dim1} = unique_mut(prevalent);
            [cnt_unique, unique_mut] = hist(SpecificMutations{1,dim1}(SpecificMutations{1,dim1}>0)...
                ,unique(SpecificMutations{1,dim1}(SpecificMutations{1,dim1}>0)));
            prevalent = cnt_unique/sscnt > 0.5;          
            SpecificMutations{1,dim1} = unique_mut(prevalent);
            prevalent = cnt_unique > 100;
            SpecificMutations{2,dim1} = unique_mut(prevalent);
        end
    end
end