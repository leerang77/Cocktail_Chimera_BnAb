function extracted = analyzeAgCaptured(antigentype, E1, E2, E3, R0, L0, Emem, rt, totalrunnum)
%% Header
% Determines the amounts of antigen extracted from the simulations
% output:
%   extracted - If strain-specific or chimeric+bnAb, then runnum x 1 vector
%               If cocktail+bnAb, then runnum x 3 vector

%Initialization
if ismember(antigentype, [1,2,4])
    extracted = zeros(totalrunnum,1);
else
    extracted = zeros(totalrunnum,3); %this is bnAb + cocktail
end
Z = 0; 

%Iterate over repeats
for runnum=1:totalrunnum
    dir = 'Data';
    if ismember(antigentype, [1,4])
        fnm = sprintf('Antigen%d_E1%.1f_R0%d_L0%d_Extracted_Emem%.1f_rt%.2f_run%d.mat',antigentype,E1,R0,L0,Emem,rt,runnum);
    else
        fnm = sprintf('Antigen%d_E1%.1f_E2%.1f_E3%.1f_R0%d_L0%d_Extracted_Emem%.1f_rt%.2f_run%d.mat',antigentype,E1,E2,E3,R0,L0,Emem,rt,runnum);
        dir = sprintf('Antigen%d_E1%.1f_E2%.1f_E3%.1f_R0%d_L0%d',antigentype,E1,E2,E3,R0,L0);
    end
    if exist(fullfile('Data',dir,fnm), 'file') == 2 %check if the *.mat file exists, then try loading
        try
            load(fullfile('Data',dir,fnm))
            Z = Z+1;
        catch
            disp(strcat('could not load the file',fnm))
            continue
        end
    else %if file does not exist
        disp(strcat(fnm, ' does not exist'))
        fid = fopen('missing_parameters.txt', 'a+'); %record the parameters for which the file doesn't exist
        fprintf(fid, '%d %.2f %.2f %.2f %d %d %d\n', [antigentype, E1, E2, E3, R0, L0, runnum]);
        fclose(fid);
        continue
    end
    if ismember(antigentype, [1,2,4]) 
        extracted(runnum) = length(Ag.Position) - length(Ag.Live);
    else %for cocktail + bnab, need to consider what types of Ag are captured 
        extracted(runnum,1) = length(setdiff(find(Ag.Type==1),Ag.Live));
        extracted(runnum,2) = length(setdiff(find(Ag.Type==2),Ag.Live));
        extracted(runnum,3) = length(setdiff(find(Ag.Type==3),Ag.Live));
    end
end
    if Z>0
        extracted = sum(extracted,1)/Z;
    end
end