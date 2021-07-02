% cutoff = -13.8; index = 2;
cutoff = -13.79; index = 1;
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
        [chimericbindN{i},~,detectablech{i}, totalBcellch{i}] = getBindN(cutoff, 1, rho, 0, -13.8, X(i));
        [chimericbindNSelected{i},~,detectableselected] = getBindNSelected(cutoff, 1, rho, 0, -13.8, X(i));
        [cocktailbindN{i},~,detectableco{i}, totalBcellco{i}] = getBindN(cutoff, 2, rho, 0, -13.8, X(i));
        [cocktailbindNSelected{i},~,detectableselected] = getBindNSelected(cutoff, 2, rho, 0, -13.8, X(i));
        if index==1
            prefix = 'precursors_included_';
        else
            prefix = '';
        end
        dlmwrite(strcat(prefix,sprintf('chimeric_bnAb_fraction_rho%.1f.csv',rho)), sum(chimericbindN{i}(1:3,:))./detectablech{i}(index,:), '-append')
        dlmwrite(strcat(prefix,sprintf('cocktail_bnAb_fraction_rho%.1f.csv',rho)), sum(cocktailbindN{i}(1:3,:))./detectableco{i}(index,:), '-append')
        dlmwrite(strcat(prefix,sprintf('chimeric_bnAb_crossreactive_rho%.1f.csv',rho)), sum(chimericbindN{i}(2:3,:))./sum(chimericbindN{i}(1:3,:)), '-append')
        dlmwrite(strcat(prefix,sprintf('cocktail_bnAb_crossreactive_rho%.1f.csv',rho)), sum(cocktailbindN{i}(2:3,:))./sum(cocktailbindN{i}(1:3,:)), '-append')
    end
end

% %% Compare selected vs. all (detectable) cross-reactive fraction of bnAbs 
% for i=1:2
%    figure
%    plot(1:28, sum(chimericbindN{i}(2:3,:))./(sum(chimericbindN{i}(1:3,:))), 'DisplayName', sprintf('All-X=%.1f',X(i)))
%    hold on
%    plot(1:28, sum(chimericbindNSelected{i}(2:3,:))./(sum(chimericbindNSelected{i}(1:3,:))), 'DisplayName', sprintf('Selected-X=%.1f',X(i)))
%     title('Cross-reactive fraction of bnAbs (detectable precursors only)')
%     legend('show')
% end
% 
% 
% %% Compare bnAb fraction in GC - all vs. detectable
% for i=1:2
%    figure
%    plot(1:28, sum(chimericbindN{i}(1:3,:))./detectablech{i}(1,:), 'DisplayName', sprintf('Detectable-X=%.1f',X(i)))
%    hold on
%    plot(1:28, totalBcellch{i}(4,:)./(sum(totalBcellch{i},1)), 'DisplayName', sprintf('All-X=%.1f',X(i)))
%    title('BnAb fraction of GC')
%    legend('show')
% end
% 
% %%
% for i
% 
% %% Cross-reactive bnAb fraction in GC - all vs. detectable
% for i=1:2
%    figure
%    plot(1:28, sum(chimericbindN{i}(2:3,:))./detectablech{i}(1,:), 'DisplayName', sprintf('Detectable-X=%.1f',X(i)))
%    hold on
%    plot(1:28, totalBcellch{i}(4,:)./(sum(totalBcellch{i},1)).*sum(chimericbindNSelected{i}(2:3,:))./(sum(chimericbindNSelected{i}(1:3,:))), 'DisplayName', sprintf('Current-X=%.1f',X(i)))
%    title('Cross-reactive BnAb fraction of GC')
%    legend('show')
% end

% figure
% for i=1:4
%    plot(1:28, sum(chimericbindN{i}(2:3,:))./detectablech{i}(index,:), 'DisplayName', sprintf('Chimeric-X=%.1f',X(i)))
%    hold on
%    plot(1:28, sum(cocktailbindN{i}(2:3,:))./detectableco{i}(index,:), 'DisplayName', sprintf('Cocktail-X=%.1f',X(i)))
%    sprintf('Chimeric-X=%.1f',X(i))
%    a = sum(chimericbindN{i}(2:3,:))./detectablech{i}(index,:);
%    disp([a(5), a(19)])
%    sprintf('Cocktail-X=%.1f',X(i))
%    a = sum(cocktailbindN{i}(2:3,:))./detectableco{i}(index,:);
%    disp([a(5), a(19)])
%    title('Cross-reactive BnAb fraction of GC')
%    legend('show')
% end
% 
% figure
% for i=1:4
%    plot(1:28, sum(chimericbindN{i}(2:3,:))./(sum(chimericbindN{i}(1:3,:))), 'DisplayName', sprintf('Chimeric-X=%.1f',X(i)))
%    hold on
%    plot(1:28, sum(cocktailbindN{i}(2:3,:))./(sum(cocktailbindN{i}(1:3,:))), 'DisplayName', sprintf('Cocktail-X=%.1f',X(i)))  
%    title('Cross-reactive fraction of bnAbs (detectable precursors only)')
%    legend('show')
% end


function [bindN,total,detectableBcell,totalBcell] = getBindN(cutoff,AgTypeIdx,rho,SharedTcellFraction,E0,X,W,bnab)
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
        filename = strcat(AgType,sprintf('Explicit_rho_%.1f_E0_%.1f_X_%.2f_f_%.1f.mat',rho,E0,X,SharedTcellFraction));
        data = load(filename);
    end

    bindN = zeros(Constants.GC_Length,3);
    total = zeros(Constants.GC_Length,1);
    for i=1:length(data.SingleGCStat)
       BnAbAffAll = data.SingleGCStat{i}.BnAbAffAll;
       for j=1:Constants.GC_Length
          total(j) = total(j) + length(BnAbAffAll{j});
          epnum = sum(BnAbAffAll{j}<cutoff, 2);
          for idx = epnum'
              if idx>0
                  bindN(j,idx) = bindN(j,idx)+1;
              end
          end
       end
    end
    bindN = bindN';
    total = total';
    detectableBcell = data.Detectable';
    totalBcell = data.TotalBcellsByEp';
end

function [bindN,total,totalBcell] = getBindNSelected(cutoff,AgTypeIdx,rho,SharedTcellFraction,E0,X,W,bnab)
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
        filename = strcat(AgType,sprintf('Explicit_rho_%.1f_E0_%.1f_X_%.2f_f_%.1f.mat',rho,E0,X,SharedTcellFraction));
        data = load(filename);
    end

    bindN = zeros(Constants.GC_Length,3);
    total = zeros(Constants.GC_Length,1);
    for i=1:length(data.SingleGCStat)
       BnAbAffSelected = data.SingleGCStat{i}.BnAbAffSelected;
       for j=1:Constants.GC_Length
          total(j) = total(j) + length(BnAbAffSelected{j});
          epnum = sum(BnAbAffSelected{j}<cutoff, 2);
          for idx = epnum'
              if idx>0
                  bindN(j,idx) = bindN(j,idx)+1;
              end
          end
       end
    end
    bindN = bindN';
    total = total';
    totalBcell = data.DetectableSelected';
end