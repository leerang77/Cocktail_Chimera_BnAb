%% Main
chimeric = zeros(28,4);
cocktail = zeros(28,4);
X = [0.4,0.6,1,1.5];
for i=1:4
    chimeric(:,i) = getGCsContainingBnAb(1000, 1, 0.7, 0, -13.8, X(i))/1000;
    cocktail(:,i) = getGCsContainingBnAb(1000, 2, 0.7, 0, -13.8, X(i))/1000;
end
dlmwrite('chimericGCsContainingBnAb.csv', chimeric');
dlmwrite('cocktailGCsContainingBnAb.csv', cocktail');


function GCs = getGCsContainingBnAb(GCNum,AgTypeIdx,rho,SharedTcellFraction,E0,X,W,bnab)
    Constants = UseSharedConstants.Constants; %get handle for the shared constants
    if AgTypeIdx==1, AgType='chimeric';
    elseif AgTypeIdx==2, AgType='cocktail';
    end
    if ~exist('W','var') %If W is not passed, then explicit Ag capture
        W=0; bnab=1;
    end
    InitializeSharedConstants(AgType,rho,E0,Constants); %initialize the shared constants
    if bnab==1
        prefix = 'BnAb';
    elseif bnab==0
        prefix = 'NoBnAb';
    end    
    if W>0
        filename = strcat(prefix,sprintf('Analytical_rho_%.1f_E0_%.1f_X_%.2f_W_%.2f_f_%.1f.mat'...
            ,rho,E0,X,W,SharedTcellFraction));
        data = load(fullfile('AM_data',filename));
    else
        filename = strcat(AgType,sprintf('Explicit_rho_%.1f_E0_%.1f_X_%.2f_f_%.1f.mat'...
            ,rho,E0,X,SharedTcellFraction));
        data = load(fullfile('AM_data',filename));
    end
    
    GCs = zeros(1,28);
    for i=1:length(data.SingleGCStat)
        for j=1:length(GCs)
            if data.SingleGCStat{i}.BcellsByEp(j,4)>0
                GCs(j) = GCs(j) + 1;
            end
        end
    end
end