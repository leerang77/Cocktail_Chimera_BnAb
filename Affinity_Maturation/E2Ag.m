function Ag = E2Ag(antigentype, E1, E2, E3, w)
% From the antigen type and binding energies, determine the amount of
% antigen captured. Can handle both analytical and explicit antigen
% capture, and both strain-specific and bnAb B cells

% Inputs
%antigentype: 1 - Cocktail, strain-specific
%             2 - Chimeric, bnAb precursor
%             3 - Cocktail, bnAb precursor
%             4 - Chimeric, strain-specific
%E1, E2, E3: Binding energy towards the three antigen variants
%w (optional): Antigen capture stringency. Only required for analytical
%              capture. If w is either 0 or not given, explicit Ag capture 

% Outputs
%Ag: If antigentype is 1 or 4, a scalar; If 2 or 3, a 1x3 vector
%----------------------------------------------------------
    E0 = -13.8;
    if ~exist('w','var')
        w = 0;
    end
    %% Explicit Ag capture
    if w == 0
        [Energy, I] = sort([E1, E2, E3]); % sort affinities and remember order
        Energy = max(Energy, -20.8); % can't exceed maximum
        E1 = Energy(1); E2 = Energy(2); E3 = Energy(3);
        persistent data
        if isempty(data)
            data = load('E2AgTable.mat'); % load the standard interpolation functions
        end
        if antigentype==1 %cocktail, strain-specific
            if E1>E0 %non-binding
                Ag = [0, 0, 0];
            else %binding
                ag = data.F{1}(E1);
                Ag(I) = [ag, 0, 0];
            end
        elseif antigentype==2 %chimeric, bnAb precursor
            if all([E1,E2,E3]>E0) %non-binding
                Ag = 0;
            elseif all([E2,E3]>E0) %bind 1
                Ag = data.F{2}{3}(E1);
            elseif E3>E0 %bind 2
                Ag = data.F{2}{2}(E1,E2);
            else %bind 3
                Ag = data.F{2}{1}(E1, E2, E3);
            end   
        elseif antigentype==3 %cocktail, bnAb precursor
            if all([E1,E2,E3]>E0) %non-binding
                Ag = [0, 0, 0];
            elseif all([E2,E3]>E0) %bind 1
                ag1 = data.F{3}{3,1}(E1);
                Ag(I) = [ag1, 0, 0];
            elseif E3>E0 %bind 2
                ag1 = data.F{3}{2,1}(E1,E2);
                ag2 = data.F{3}{2,2}(E2,E2);
                Ag(I) = [ag1, ag2, 0];
            else %bind 3
                ag1 = data.F{3}{1,1}(E1, E2, E3);
                ag2 = data.F{3}{1,2}(E1, E2, E3);
                ag3 = data.F{3}{1,3}(E1, E2, E3);
                Ag(I) = [ag1, ag2, ag3];
            end
        elseif antigentype==4 %chimeric, strain-specific
            if E1>E0 %non-binding
                Ag = 0;
            else %binding
                Ag = data.F{4}(E1);
            end    
        end
    else
    %% Analytical Ag capture  
        [Energy, I] = sort([E1, E2, E3]);
        E1 = Energy(1); E2 = Energy(2); E3 = Energy(3);
        if ismember(antigentype, [1,4]) %strain-specific
            Ag(I) = [exp(-w*(E1-E0)), 0, 0];
        else %bnab precursor
            Ag(I) = [exp(-w*(E1-E0)), exp(-w*(E2-E0)), exp(-w*(E3-E0))];
        end
    end
end