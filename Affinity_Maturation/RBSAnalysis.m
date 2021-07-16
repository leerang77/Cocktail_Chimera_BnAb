%% Analyze RBS-directed B cell statistics
cutoff = -13.8; index = 2; % This makes the analysis focus on B cells that have affinities at least greater than E0
X = [0.4, 0.6, 1.0, 1.5];
chimericbindN = cell(1,length(X));
cocktailbindN = cell(1,length(X));
chimericbindNSelected = cell(1,length(X));
cocktailbindNSelected = cell(1,length(X));
detectablech = cell(1,length(X));
detectableco = cell(1,length(X));
totalBcellch = cell(1,length(X));
totalBcellco = cell(1,length(X));
for rho=[0.7]
    for i=1:length(X)
        [chimericbindN{i},detectablech{i}, totalBcellch{i}] = getBindNCombined(cutoff, 1, rho, 0, -13.8, X(i));
        [cocktailbindN{i},detectableco{i}, totalBcellco{i}] = getBindNCombined(cutoff, 2, rho, 0, -13.8, X(i));
        if index==1
            prefix = 'precursors_included_';
        else
            prefix = '';
        end
        %Fraction of GC B cells that are RBS-directed
        dlmwrite(strcat(prefix,sprintf('chimeric_bnAb_fraction_rho%.1f.csv',rho)), sum(chimericbindN{i}(1:3,:))./detectablech{i}(index,:), '-append')
        dlmwrite(strcat(prefix,sprintf('cocktail_bnAb_fraction_rho%.1f.csv',rho)), sum(cocktailbindN{i}(1:3,:))./detectableco{i}(index,:), '-append')
        %Fraction of RBS-directed B cells that are cross-reactive
        dlmwrite(strcat(prefix,sprintf('chimeric_bnAb_crossreactive_rho%.1f.csv',rho)), sum(chimericbindN{i}(2:3,:))./sum(chimericbindN{i}(1:3,:)), '-append')
        dlmwrite(strcat(prefix,sprintf('cocktail_bnAb_crossreactive_rho%.1f.csv',rho)), sum(cocktailbindN{i}(2:3,:))./sum(cocktailbindN{i}(1:3,:)), '-append')
    end
end

function [bindN,total,detectableBcell,totalBcell] = getBindNCombined(cutoff,AgTypeIdx,rho,SharedTcellFraction,E0,X,W,bnab)
first = 1;
gcnum = 50;
last = 1000;
for index=first:gcnum:(last-gcnum+1)
    if index==1
        [bindN,total,detectableBcell,totalBcell] = getBindN(cutoff,index,gcnum,AgTypeIdx,rho,SharedTcellFraction,E0,X,W,bnab);
    else
        [a,b,c,d] = getBindN(cutoff,index,gcnum,AgTypeIdx,rho,SharedTcellFraction,E0,X,W,bnab);
        bindN = bindN + a;
        total = total + b;
        detectableBcell = detectableBcell + c;
        totalBcell = totalBcell + d;
    end
end
end

function [bindN,detectableBcell,totalBcell] = getBindN(cutoff,beginning,gcnum,AgTypeIdx,rho,SharedTcellFraction,E0,X,W,bnab)
%First load the data
    Constants = UseSharedConstants.Constants; %get handle for the shared constants
    switch AgTypeIdx
        case 1
            AgType = 'chimeric';
        case 2
            AgType = 'cocktail';
    end
    if ~exist('W','var')
        W=0;
        bnab=1;
    end
    InitializeSharedConstants(AgType,rho,E0,Constants); %initialize the shared constants
    if bnab==1
        prefix = 'BnAb';
    elseif bnab==0
        prefix = 'NoBnAb';
    end    
    if W>0
        filename = strcat(prefix,sprintf('Analytical_rho_%.1f_E0_%.1f_X_%.2f_W_%.2f_f_%.1f.mat',rho,E0,X,W,SharedTcellFraction));
        data = load(filename);
    else
        filename = strcat(AgType,sprintf('Explicit_rho_%.1f_E0_%.1f_X_%.2f_f_%.1f_gcnum_%d_to_%d.mat',rho,E0,X,SharedTcellFraction,beginning,beginning+gcnum-1));
        data = load(filename);
    end

%Derive the desired statistics
    bindN = zeros(Constants.GC_Length,3);
    for i=1:length(data.SingleGCStat)
       BnAbAffAll = data.SingleGCStat{i}.BnAbAffAll;
       for j=1:Constants.GC_Length
          epnum = sum(BnAbAffAll{j}<cutoff, 2); %Number of variants that the B cell recognizes
          for idx = epnum'
              if idx>0
                  bindN(j,idx) = bindN(j,idx)+1;
              end
          end
       end
    end
    bindN = bindN';
    detectableBcell = data.Detectable';
    totalBcell = data.TotalBcellsByEp';
end