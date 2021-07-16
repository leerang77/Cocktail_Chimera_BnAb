function runMultipleGC(beginning,GCNum,AgTypeIdx,rho,SharedTcellFraction,X,W,bnab)
%% HEADER
 %This code can run affinity maturation simulations with both the analytical antigen capture 
 %and cocktail/chimeric antigen capture.
 %To run chimeric antigen, use AgTypeIdx=1, SharedTcellFraction=0
 %To run cocktail antigen, use AgTypeIdx=2, SharedTcellFraction=0
 %Pass W=0 for antigen capture from explicit simulation
 %The code runs number of GCs indexed between beginning - beginning+GCNum
 % e.g. 1 - 50; this way, many simulations can be run in parallel. 
 
 % Inputs
 %beginning - Int, index of the first GC
 %GCNum - Int, Number of GCs to run (i.e. number of independent repeats)
 %AgTypeIdx - If 1, chimeric antigen
 %            If 2, cocktail antigen
 %rho       - Float between 0 and 1, Antigenic distance between the different strains
 %            Specifically, the mutation is derived from a 3-D Gaussian
 %            with rho as the correlation and then converted to shifted
 %            log-normal
 %SharedTcellFraction - Float between 0 and 1, The fraction of T cell epitopes shared among the antigens
 %W - If W=0, antigen capture based on explicit simulation is used. If W>0,
 %      simple analytical formula is used.
 %X  - Float, Stringency of T cell selection in the GC
 %bnab- 0 or 1, Whether there is 0 or 1 bnAb precursor among the founders

 %  Outputs
 %Saves the following statistics in a *.mat file.
 %SingleGCStat - Statistics from a single GC run. See the function
 %               runSingleGC for more description.
 %TotalBcellsByEp - GC_Length x N_Ep array, Number of B cells targeting
 %                  each epitope after each cycle
 %MeanSurvivingLineages - GC Length x 1 array
 %MeanDomOccupancy - GC Length x 1 array
 %%  -------- Main Script --------
    Constants = UseSharedConstants.Constants; %get handle for the shared constants
    if AgTypeIdx==1, AgType='chimeric';
    elseif AgTypeIdx==2, AgType='cocktail';
    end
    InitializeSharedConstants(AgType,rho,bnab,Constants); %Initialize the shared constants
    
    %Define the similarity of the T cell epitopes
    SpecificTcellFraction = (1-SharedTcellFraction);
    TcellsFraction = SpecificTcellFraction*[1/3,1/3,1/3,0]+[0,0,0,SharedTcellFraction];
    TfhC = Tcells(TcellsFraction);
    
    %Run the simulations 
    SingleGCStat = cell(1,GCNum);
    for gcnumber=beginning:(beginning+GCNum-1)
        disp(gcnumber);
        rng(gcnumber+1234); %Seed random numbers for reproducibility
        SingleGCStat{gcnumber-beginning+1} = runSingleGC(X);
    end

    %Analyze to derive mean statistics
    analyzeMultipleGC(W,beginning,gcnumber)
    
%% Subfunctions
    function analyzeMultipleGC(W,beginning,gcnumber)
     %% Use the results of individual GCs to calculate overall statistics
     % Input
     %W - Ag capture stringency. If 0, explicit Ag capture; If >0,
     %    analytical Ag capture
     %beginning - Int, index of the first GC
     %GCNum - Int, Number of GCs to run (i.e. number of independent repeats)
     % Output
     %None, saves the result into *mat file
     %----------------------------------------
        TotalBcellsByEp = zeros(Constants.GC_Length, Constants.N_Ep);
        Detectable = zeros(Constants.GC_Length,2); % First row - Number of B cells with affinity <= E0, 
                                                   % Second row - Number of B cells with affinity < E0
        DetectableSelected = zeros(Constants.GC_Length,2); % Number of selected B cells
        Day3BnAbSelectionProb = cell(1,3); % Selection probabilities of bnAb lineage B cells on day 3
        Day5BnAbSelectionProb = cell(1,3);
        Day10BnAbSelectionProb = cell(1,3);
        Day12BnAbSelectionProb = cell(1,3);

        for i=1:GCNum
            for j=1:3
                Day3BnAbSelectionProb{j}{i} = SingleGCStat{i}.BnAbSelectionProb{6}{j};
                Day5BnAbSelectionProb{j}{i} = SingleGCStat{i}.BnAbSelectionProb{10}{j};
                Day10BnAbSelectionProb{j}{i} = SingleGCStat{i}.BnAbSelectionProb{20}{j};               
                Day12BnAbSelectionProb{j}{i} = SingleGCStat{i}.BnAbSelectionProb{24}{j};
            end
            TotalBcellsByEp = TotalBcellsByEp + SingleGCStat{i}.BcellsByEp;
            Detectable = Detectable + SingleGCStat{i}.Detectable;
            DetectableSelected = DetectableSelected + SingleGCStat{i}.DetectableSelected;
        end
        for j=1:3
            Day3BnAbSelectionProb{j} = [Day3BnAbSelectionProb{j}{:}];
            Day5BnAbSelectionProb{j} = [Day5BnAbSelectionProb{j}{:}];
            Day10BnAbSelectionProb{j} = [Day10BnAbSelectionProb{j}{:}];
            Day12BnAbSelectionProb{j} = [Day12BnAbSelectionProb{j}{:}];
        end
        
        % Define filename for output
        if bnab==1, prefix = 'BnAb';
        elseif bnab==0, prefix = 'NoBnAb';
        end
        if W>0 % Analytical Ag capture
            filename = strcat(prefix,sprintf('Analytical_rho_%.1f_E0_%.1f_X_%.2f_W_%.2f_f_%.1f_gcnum_%d_to_%d.mat'...
                ,rho,Constants.E_initial,X,W,SharedTcellFraction,beginning,gcnumber));
        else % Explicit Ag capture
            filename = strcat(AgType,sprintf('Explicit_rho_%.1f_E0_%.1f_X_%.2f_f_%.1f_gcnum_%d_to_%d.mat',...
                rho,Constants.E_initial,X,SharedTcellFraction,beginning,gcnumber));
        end
        [~,fn,~] = fileparts(filename);
        writeBnAbSelectionProb = 0; % Change to 1 if want the csv files containing bnAb selection probabilities
        if writeBnAbSelectionProb == 1
            for j=1:3
                dlmwrite(strcat(fn,'_Day3BnAbSelectionProb.csv'), Day3BnAbSelectionProb{j}, '-append')
                dlmwrite(strcat(fn,'_Day5BnAbSelectionProb.csv'), Day5BnAbSelectionProb{j}, '-append')
                dlmwrite(strcat(fn,'_Day10BnAbSelectionProb.csv'), Day10BnAbSelectionProb{j}, '-append')
                dlmwrite(strcat(fn,'_Day12BnAbSelectionProb.csv'), Day12BnAbSelectionProb{j}, '-append')
            end
        end
        save(filename,'SingleGCStat','TotalBcellsByEp','Detectable','DetectableSelected') 
    end

    function GCStat = runSingleGC(X)
        %% Run a single GC simulation - includes seeding GC, competitive phase, and analysis
        % Input
        %X - Float, Stringency of T cell selection
        
        % Output
        %GCStat - Struct with following fields:
        %   BenMuts: Cell array; All affinity-increasing mutations that occur in GC,
        %            either off-target of bnAb
        %   BnAbMutations: 2 x GC length Cell array; BnAb mutations that
        %                  are prevalent in the population at given time, based on two
        %                  different definitions of "prevalence" (each row). 
        %   BenMutsByTime: 1 x GC length;
        %                  Number of total affinity-increasing mutations that occur in the GC by time
        %   SpecificMutations: 2 x GC length Cell array; Specific mutations that are prevalent in the population
        %   BcellsByEp: GC length x 4; Number of B cells targeting each
        %               epitope by time (three strain-specific epitopes and the bnAb epitope)
        %   SurvivingLineages: GC length x 1; Number of surviving lineages by time
        %   DomOccupancy: GC length x 1; Occupancy of the most dominant lineage by time 
        %   SpecificAff: GC length x 2 array; 
        %                At each time, highest and lowest affinity of selected strain-specific B cells 
        %   BnAbAff: GC length x 4 array; Mean affinity and number of positively selected bnAb B cells
        %   BnAbAffAll: GC length x 1 cell; Each cell is an array of
        %               affinities of all bnAb B cells at the time
        %   SpecificNumMut: Average number of mutations (total, affinity affecting) on SS B cells
        %   BnAbNumMut: Average number of mutations (total, affinity affecting) on BnAb precursors
        %   Detectable: GC length x 2 array; Number of B cells that have affinity above threshold
                                             % First row - Number of B cells with affinity <= E0, 
                                             % Second row - Number of B cells with affinity < E0
        %   DetectableSelected: GC length x 2 array; Same as Detectable,
        %                                            but only counting positively selected B cells
        %   PlasmaArray: Cell Array; Contains all plasma cells
        
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
        end % End of competitive phase
        
        % To save storage space, clean up PlasmaArray cell array
        for i=1:length(PlasmaArray)
            if isempty(PlasmaArray(i).TargetEp)
                break
            else
                PlasmaArray(i).Constants=[]; 
            end
        end
        i=i-1;
        PlasmaArray = PlasmaArray(1:i);
        
        % Analysis of the GC
         %Use helper functions to get desired statistics
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
    %% Determine selection probabilities of all bnAb B cells
    %Input
    % pSel - Array containing selection probabilities of all B cells
    % BcellsArray - Array of B cells
    %Output
    % bnAbSelectionProb - Array containing selection probabilities of only
    %                     bnAb B cells
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
        %Detectable: Detectable B cells 
        %SurvivingLineages: GC Length x 1 array. Number of surviving lineages
        %DomOccpancy: GC Length x 1 array. Occupancy of the dominant clone
        %BnAbAffAll: Cell array containing affinities of all bnAb B cells
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
        SpecificAff(:,2) = SpecificAff(:,2)-100;
        BnAbAff = zeros(gclength,4);
        SpecificNumMut = zeros(gclength,3);
        BnAbNumMut = zeros(gclength,3);
        BnAbMutations = cell(2,gclength);
        SpecificMutations = cell(2,gclength);
        BnAbAffAll = cell(1,gclength);
        Detectable = zeros(gclength,2);
        for dim1 = 1:gclength
            Aff = zeros(number,1);
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
                    SpecificNumMut(dim1,1) = SpecificNumMut(dim1,1) + BcellsArray(dim1,dim2).NumMut;
                    SpecificNumMut(dim1,2) = SpecificNumMut(dim1,2) + sum(BcellsArray(dim1,dim2).NumAffMut);
                    SpecificNumMut(dim1,3) = SpecificNumMut(dim1,3) + 1;
                    SpecificMutations{1,dim1}(dim2,:) = BcellsArray(dim1,dim2).MutHistory;  
                    Aff(dim2) = BcellsArray(dim1,dim2).E_AgAb;
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
            Aff = Aff(Aff~=0);
            SpecificAff(dim1,1) = prctile(Aff,10);
            SpecificAff(dim1,2) = prctile(Aff,90);           
            BnAbAffAll{dim1} = BnAbAffAll{dim1}(sum(BnAbAffAll{dim1},2)~=0,:);
%             SpecificAff(dim1,1) = SpecificAff(dim1,1)/SpecificAff(dim1,2);
%             SpecificAff(isnan(SpecificAff))=0;
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