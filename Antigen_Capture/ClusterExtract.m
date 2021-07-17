function ClusterExtract(antigen, E1, E2, E3, R0, L0, runnum)
%% Header
% This function simulates explicit antigen capture through clustering and
% extraction. It calls two subfunctions Cluster and Extract in order.

% Inputs:
%   antigentype: 1 - Cocktail, strain-specific
%                2 - Chimeric, bnAb precursor
%                3 - Cocktail, bnAb precursor
%                4 - Chimeric, strain-specific
%   E1, E2, E3 - Binding energy towards the three antigen variants
%   R0 - Number of BCR receptors in the B cell-APC interface
%   L0 - Number of antigen in the B cell-APC interface
%   runnum - index for repeat runs

% Outputs:
%   Saves the simulation result to a .mat file

% Summary of subfunctions:
%   Initialize - Randomly initialize the positions of R0 BCRs and L0
%                antigens on circular lattice
%   Cluster    - Simulate diffusion, bond forming, and bond breaking of
%                BCRs and antigens
%   Extract    - Simulate force-based extraction of Ags
%   Diffusion  - Take a step of free diffusion of Ag and BCR clusters
%   analyzeClusters  - Analyzes the number of BCRs and Ags contained in the clusters 
%   Internalization  - During the extraction phase, internalize any Ags
%                      that are bound to BCRs and not to APC; also, retract BCRs that are not
%                      bound to any Ags
%   BondBreakExtraction  - Breaking of bonds during the extraction phase
%                          according to the off-rates
%   BondBreakClustering  - Breaking of bonds during the clustering phase
%                          according to the off-rates
%   BondForm  - Formation of bonds during the clustering/extraction phase 
%   Visualize  - Visualization of the B cell-APC interface
%   Clustering  - Analysis of the cluster structures based on the Ag and
%               BCR bonds

% Key variables:
% Ag - Struct with the following fields:
%   Epitopes - 1 x 6 vector of epitopes that the antigen contains. For
%              example, [1 1 1 4 4 4] refers to homotrimer of variant 1,
%              and [1 2 3 4 5 6] refers to a heterotrimer. 
%              Exception: 3 x 6 vector for the cocktail Ag + bnAb
%              combination
%   Lattice - 2601 x 2 double.
%             First column - occupying molecule type:
%               0-None, 2-Ag, 3-Blocked
%             Second column - occupying molecule index
%   Position - L0 x 1 double, The position of Ag molecules on the lattice
%   Binding - L0 x 6 double, each column contains the index of the BCR that
%             the corresponding epitope of the Ag is bind to, 
%             or 0 if the given epitope is not bound.
%   Type - L0 x 1 double, Identity of the Ag molecules
%   Live - n x 1, where n is the number of Ag molecules still in the interface.
%          Contains the indices of live Ags (not yet internalized).
%   Mem - L0 x 1 double, 1 if the Ag is bound to the membrane, 0 otherwise.
%   Force - Force on the Ag when the Ag-Mem bond breaks
% BCR - Struct with the following fields:
%   Energy - Binding affinity towards the six epitopes (three
%            strain-specific, and three bnAB epitope). 0 if the BCR doesn't
%            bind the epitope
%   koff - Off-rate towards the six epitopes
%   Lattice - 2601 x 2 double. 
%             First column - occupying molecule type:
%               0-None, 1-BCR, 3-Blocked
%             Second column - occupying molecule index
%   Position - R0 x 1. The position of BCR molecules on the lattice
%   Binding - R0 x 4. The first two columns are the indices of Ags that are
%             bound to the BCR. The third and fourth columns are the identities of
%             the epitopes. For example, if a bnAb B cell it can be [10,
%             25, 4, 6].
%   Live - n x 1, where n is the number of BCR molecules still in the interface.
%          Contains the indices of live BCRs (not yet internalized).
%   Refractory - R0 x 1. 1 if the BCR is in refractory phase, 0 if otherwise.
%   Rebinding - R0 x 1. Total number of times that BCR rebinds. 
% Bond - L0 x R0 array. An entry is 1 if there is a bond between BCR and Ag.
% clusters - Struct with the following fields:
%   NumElem - Number of clusters x 2 array. First column - number of BCRs,
%             second column - number of Ags in the clusters.
%   Elem - Number of clusters x 2 array. First column - Indices of BCRs, 
%          Second column - Indices of Ags contained in the clusters.

%------- Main Script -------
fprintf('antigen=%d, E1=%.1f, E2=%.1f, E3=%.1f, R0=%d, L0=%d, runnum=%d\n',...
    antigen, E1, E2, E3, R0, L0, runnum);
%Constants
rt = 0.1;
Emem = -19;
%Define parameters according to the given scenario
switch antigen
    case 1 %Single heterotrimer
        fnm = sprintf('Antigen%d_E1%.1f_R0%d_L0%d',antigen,E1,R0,L0);
        Ag.Epitopes = [1, 1, 1, 4, 4, 4;
               2, 2, 2, 5, 5, 5;
               3, 3, 3, 6, 6, 6;];
        BCR.Energy = [E1, 0, 0, 0, 0, 0];
    case 2 %Chimeric
        fnm = sprintf('Antigen%d_E1%.1f_E2%.1f_E3%.1f_R0%d_L0%d',antigen,E1,E2,E3,R0,L0);
        Ag.Epitopes = [1, 2, 3, 4, 5, 6];
        BCR.Energy = [0, 0, 0, E1, E2, E3];      
    case 3 %Cocktail
        fnm = sprintf('Antigen%d_E1%.1f_E2%.1f_E3%.1f_R0%d_L0%d',antigen,E1,E2,E3,R0,L0);
        Ag.Epitopes = [1, 1, 1, 4, 4, 4;
               2, 2, 2, 5, 5, 5;
               3, 3, 3, 6, 6, 6;];
        BCR.Energy = [0, 0, 0, E1, E2, E3];
    case 4 %Single homotrimer
        fnm = sprintf('Antigen%d_E1%.1f_R0%d_L0%d',antigen,E1,R0,L0);
        Ag.Epitopes = [1,2,3,4,5,6];
        BCR.Energy = [E1, 0, 0, 0, 0, 0];
end
BCR.koff = 10^6./exp(-BCR.Energy);

if ~exist('Data', 'dir')
    mkdir('Data')
end
if ~exist(fullfile('Data',fnm),'dir')
    mkdir(fullfile('Data',fnm))
end
rng(runnum+1000)
if ~exist(fullfile('Data',fnm,strcat(fnm,sprintf('_Clustered_run%d.mat',runnum))), 'file') %Run Cluster if file doesn't exist
    Cluster(antigen, fnm, BCR, Ag, R0, L0, runnum);
end
rng(runnum+2000)
if ~exist(fullfile('Data',fnm,strcat(fnm,sprintf(... %Run Extract if file doesn't already exist
        '_Extracted_Emem%.1f_rt%.2f_run%d.mat',Emem,rt,runnum))), 'file')
    load(fullfile('Data',fnm,strcat(fnm,sprintf('_Clustered_run%d.mat',runnum)))) %load cluster file
    Extract(fnm, Bond, BCR, Ag, N, R0, L0, clustersbytime, index, runnum); %Run Extract
end
clear all
end
%% Subfunctions
%%
function [BCR, Ag, Bond] = Initialize(BCR, Ag, N, antigen, Dim, l, R0, L0)
% Randomly initialize the positions of R0 BCRs and L0
% antigens on circular lattice
% Inputs:
%   BCR, Ag -> structs
%   N -> Number of lattice sites along the diameter of the interface
%   antigen -> Antigen type
%   Dim -> Diameter of the interface
%   l -> lattice dimension
%   R0, L0 -> Number of BCRs and Ags
% Outputs:
%   BCR, Ag -> structs
%   Bond -> Adjacency matrix for BCRs and Ags

%First, default initialization
Ag.Lattice = zeros(N^2,2);
BCR.Lattice = zeros(N^2,2);
BCR.Position = zeros(R0,1);
BCR.Binding = zeros(R0,4);
Ag.Position = zeros(L0,1);
Ag.Binding = zeros(L0,6);
Ag.Type = zeros(L0,1);
Bond = zeros(R0,L0);
open = zeros(1,N^2);
BCR.Live = (1:R0)';
Ag.Live = (1:L0)';
Ag.Mem = ones(L0,1);
BCR.Refractory = zeros(R0,1);
% Define the circular region of the Ag - BCR interface
for i=1:N^2 
[loc1,loc2] = ind2sub([N,N], i);
d = sqrt((((N+1)/2-loc1)^2+((N+1)/2-loc2)^2)*l^2)/1000; %distance from center, in micrometers
if d > Dim/2 % block if further than the radius
  Ag.Lattice(i,1) = 3;
  BCR.Lattice(i,1) = 3;
else % store the index of open positions
  open(i)=1;
end
end
open = find(open==1);
% Randomly place the Ags and BCRs on open positions of the lattice
I = datasample(open,R0+L0,'Replace',false); 
for i=1:R0 % Place BCRs on the lattice
BCR.Position(i) = I(i);
BCR.Lattice(I(i),1) = 1;
BCR.Lattice(I(i),2) = i;
end
for i=R0+1:R0+L0 % Place Ags on the lattice
Ag.Position(i-R0) = I(i);
Ag.Lattice(I(i),1) = 2;
Ag.Lattice(I(i),2) = i-R0;
end
%Assign Ag types
if ismember(antigen,[1,2,4]) %If cocktail+specific or chimeric, there is only one Ag type
Ag.Type(:) = 1; 
elseif antigen==3 %If cocktail+bnAb, there are three Ag types. 
              %Three variants are equally represented.
Ag.Type(1:L0/3,1) = 1;
Ag.Type(L0/3+1:L0*2/3,1) = 2;
Ag.Type(L0*2/3+1:L0,1) = 3;
end
end

%%
function Cluster(antigen, fnm, BCR, Ag, R0, L0, runnum)
% Simulate diffusion, bond forming, and bond breaking of
% BCRs and antigens
% Inputs:
%   antigen -> Ag type
%   fnm -> file name for saving
%   BCR, Ag -> structs
%   R0, L0 -> BCR and Ag number in interface
% Outputs:
%   Saves the workspace to a *.mat file

% constants
l = 10; %nm
D0 = 5*10^4; %nm^2 s^-1
dt = 5*10^-4; %s
q0 = 10; %s^-1
Dim = 0.5; %micrometer
N = round(Dim*1000/l)+1; %Number of lattices within the diameter

% Initialization
disp('Initialization of molecules begins')
[BCR, Ag, Bond] = Initialize(BCR, Ag, N, antigen, Dim, l, R0, L0);
disp('Initialization of molecules completed')
% Tracking of clusters
Clustertypes_str = {'0-1','1-0','1-1','1-2','2-1','2-2','2-3','3-1','3-2','3-3','higher order'};
                %cluster description - number of BCR - number of Ag
clustersbytime = zeros(520,11); %tracks the number of clusters in the categories at 0.5 s interval
clusters = Clustering(Bond, BCR, Ag, R0, L0); %Define the clusters from the adjacency matrix
index = 1;
clustersbytime(index,:) = analyzeClusters(clusters);
index = index+1;
% Clustering phase
disp('Clustering Begins')
tic
for i=1:2*10^4 % 10 seconds
%Each step contains diffusion, bondbreaking, and bondformation
[BCR, Ag] = Diffusion(clusters, BCR, Ag, N, D0, dt, l, 'clustering');
[BCR,Ag,Bond] = BondBreakClustering(BCR, Ag, Bond, dt);
[BCR,Ag,Bond] = BondForm(BCR, Ag, Bond, N, q0, dt);
clusters = Clustering(Bond, BCR, Ag, R0, L0);
if rem(i,10*10^3)==0
    clustersbytime(index,:) = analyzeClusters(clusters);
    index = index+1;
    fprintf('SimTime=%.1f seconds, Elapsed=%.1f seconds\n', i*dt, toc)
%     Visualize(N,BCR,Ag,Bond);
end    
end
fprintf('Clustering Completed. Elapsed Time = %.1f seconds\n', toc)
save(fullfile('Data',fnm,strcat(fnm,sprintf('_Clustered_run%d.mat',runnum)))) %Save the result
end

%%
function Extract(fnm, Bond, BCR, Ag, N, R0, L0, clustersbytime, index, runnum)
%   Simulate force-based extraction of Ags
% Inputs:
%   fnm - file name for saving
%   Bond - adjacency matrix of BCR and Ag
%   BCR, Ag - structs
%   N - number of lattice sites along the diameter 
%   R0, L0 - number of BCRs and Ags in the interface
%   runnum - repeat number
% Outputs:
%   Saves workspace into a *.mat file

%constants
l = 10; %nm
D0 = 5*10^4; %nm^2 s^-1
dt = 5*10^-4; %s
rt = 0.1; %s - refractory time 
q0 = 10; %s^-1
xb = 1; %nm
F = 8; %pN
kbT = 4.1; %pN nm
Emem = -19; %kbT
stability = 2*D0*dt/(l^2); %should be less than 1
rng(runnum)

clusters = Clustering(Bond, BCR, Ag, R0, L0);
BCR.Rebinding = zeros(R0,2);
Ag.Force = zeros(L0,2);
% Extraction Phase
disp('Extraction Begins')
tic
for i=1:50*10^5 % 2500 seconds max
    [BCR, Ag] = Diffusion(clusters, BCR, Ag, N, D0, dt, l, 'extraction');
    [BCR, Ag, Bond] = BondBreakExtraction(BCR, Ag, Bond, dt, rt, Emem, xb, F, kbT);
    [BCR, Ag, Bond] = BondForm(BCR, Ag, Bond, N, q0, dt);
    clusters = Clustering(Bond, BCR, Ag, R0, L0);
    [BCR, Ag, Bond] = Internalization(Bond, clusters, BCR, Ag);
    if isempty(BCR.Live)
        fprintf('Extraction has been competed. SimTime=%.1f seconds\n', i*dt)
        break
    end
    if rem(i,10*10^3)==0
        fprintf('SimTime=%.1f seconds, Elapsed=%.1f seconds, LiveBCR=%d\n', i*dt, toc, length(BCR.Live))
        if index <= length(clustersbytime)
            clustersbytime(index,:) = analyzeClusters(clusters);
            index = index+1;
        end
    %     Visualize(N,BCR,Ag,Bond);
    end
end
%Save result
save(fullfile('Data',fnm,strcat(fnm,sprintf('_Extracted_Emem%.1f_rt%.2f_run%d.mat',Emem,rt,runnum))))
end

%%
function clusternum = analyzeClusters(clusters)
% Analyzes the number of BCRs and Ags contained in the clusters 
% Input:
%   clusters - struct containing the clusters
% Output:
%   clusternum - 1 x 11 array. Number of clusters in the 11 categories
clusternum = zeros(1,11);
%Clustertypes_str = {'0-1','1-0','1-1','1-2','2-1','2-2','2-3','3-1','3-2','3-3','higher order'};
Clustertypes_num = {[0,1],[1,0],[1,1],[1,2],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]};
for i=1:size(clusters.NumElem,1)
    for j=1:10
       if isequal(clusters.NumElem(i,:), Clustertypes_num{j})
          clusternum(1,j) = clusternum(1,j)+1;
       end
    end
end
clusternum(1,11) = size(clusters.NumElem,1) - sum(clusternum(1,1:10));
end

%%
function clusters = Clustering(Bond, BCR, Ag, R0, L0)
% From the adjacency matrix, define the clusters formed by BCRs and Ags
% Outputs:
%   clusters - Struct with the following fields:
%       NumElem - Number of clusters x 2 array. First column - number of BCRs,
%               second column - number of Ags in the clusters.
%       Elem - Number of clusters x 2 array. First column - Indices of BCRs, 
%               Second column - Indices of Ags contained in the clusters.

%First, default initialization
numcluster = 0;
newclusteridx = 0;
BCRcluster = zeros(R0,1);
Agcluster = zeros(L0,1);
clusters.NumElem = zeros(numcluster,2);
clusters.Elem = cell(numcluster,2);

%Define clusters by iterating over BCRs
for i=1:R0 %for each BCR
if ~ismember(i, BCR.Live) %Ensure the BCR is live
    continue
end
bonds = find(Bond(i,:)>0); %find the indices of all Ag bound to it
joining = Agcluster(bonds');%find all the clusters to be connected
joining = unique(joining(joining>0));
if isempty(joining) % if none of the antigens already belong to a cluster
    numcluster = numcluster+1; %create a new cluster
    newclusteridx = newclusteridx+1;
    BCRcluster(i) = newclusteridx;
    Agcluster(bonds') = newclusteridx;
else % if at least one antigen already belongs to a cluster
    s = min(joining); % get the lowest index of the clusters to be joined
    BCRcluster(ismember(BCRcluster,joining)) = s; % coalesce all clusters
    Agcluster(ismember(Agcluster,joining)) = s;
    Agcluster(bonds') = s; % ensure the Ags not appearing in other clusters are added
    BCRcluster(i) = s; % ensure the current BCR is added
    numcluster = numcluster - (length(joining)-1); % decrease the number of clusters accordingly
end
end

%Add free Ags
singleAgs = find(Agcluster==0); %Remaining Ags that don't belong to clusters yet
for i=1:length(singleAgs)
if ~ismember(singleAgs(i), Ag.Live) %Ensure the Ag is live
    continue
end
numcluster = numcluster + 1; %Each Ag becomes a new cluster
newclusteridx = newclusteridx + 1;
Agcluster(singleAgs(i)) = newclusteridx;
end

%Summarize the clusters
idx = 0;
for i=1:newclusteridx
if any([BCRcluster==i;Agcluster==i])
   idx = idx+1;
   clusters.NumElem(idx,1) = sum(BCRcluster==i);% Number of BCRs in the cluster
   clusters.NumElem(idx,2) = sum(Agcluster==i);% Number of Ags in the cluster
   clusters.Elem{idx,1} = find(BCRcluster==i); % Indices BCRs
   clusters.Elem{idx,2} = find(Agcluster==i); % Indices of Ags
end
end
end

%%
function [BCR, Ag, Bond] = Internalization(Bond, clusters, BCR, Ag)
%   During the extraction phase, internalize any clusters that connected to APC. That is, Ags
%   that are bound to BCRs but not to APCs, as well as free BCRs that are not
%   bound to any Ags
% Inputs:
%   Bond - BCR and Ag adjacency matrix
%   clusters - struct containing the composition of each cluster
%   BCR, Ag - structs
% Output:
%   BCR, Ag - structs, updated after internalization
%   Bond - updated adjacency matrix

intclusters = []; %clusters to be internalized
for i=1:size(clusters.Elem,1) % loop over all clusters
if clusters.NumElem(i,2) == 0 %If no Ag is in the cluster, i.e. Free BCR
    BCR.Live = setdiff(BCR.Live,clusters.Elem{i,1}); %Remove BCR from live index
    BCR.Lattice(BCR.Position(clusters.Elem{i,1}),[1,2]) = [0, 0]; %Set lattice occupancy to 0
    BCR.Position(clusters.Elem{i,1}) = 0; %Set BCR positions to 0
    intclusters = [intclusters;i]; 
elseif all(Ag.Mem(clusters.Elem{i,2})==0) %Or if all Ags are detached from the APC
    BCR.Live = setdiff(BCR.Live,clusters.Elem{i,1}); %Remove BCRs and Ags from live index
    Ag.Live = setdiff(Ag.Live,clusters.Elem{i,2});
    BCR.Lattice(BCR.Position(clusters.Elem{i,1}),[1,2]) = 0; %Set lattice occupancy to 0
    Ag.Lattice(Ag.Position(clusters.Elem{i,2}),[1,2]) = 0;
    BCR.Position(clusters.Elem{i,1}) = 0; %Set positions to 0
    Ag.Position(clusters.Elem{i,2}) = 0;
    Bond(clusters.Elem{i,1},:) = 0; %Remove bonds of the internalized molecules from adjacency matrix
    intclusters = [intclusters;i];
end
end
if ~isempty(intclusters) %Remove the internalized clusters
clusters.NumElem(intclusters,:) = [];
clusters.Elem(intclusters,:) = [];
end
end

%%
function [BCR, Ag, Bond] = BondForm(BCR, Ag, Bond, N, q0, dt)
% Probabilistically forms bonds between Ags and BCRs when they are at the same or adjacent
% lattice sites.
% Inputs:
%   BCR, Ag - structs
%   Bond - adjacency matrix
%   N - Number of lattice sites along diameter
%   q0 - base on-rate for one BCR arm and one Ag epitope
%   dt - time step
% Outputs:
%   BCR, Ag, Bond - updated

%Iterate over live BCRs
for i=1:length(BCR.Live)
% check if the BCR is in refractory period
if BCR.Refractory(BCR.Live(i))>0
    continue
end
% check if any of the BCR arms is empty
freebcr = find(BCR.Binding(BCR.Live(i),1:2)==0);
if isempty(freebcr)
    continue
end
% check if any Ag molecule in possible binding position
neighbors = [BCR.Position(BCR.Live(i));BCR.Position(BCR.Live(i))-1;...
    BCR.Position(BCR.Live(i))+1;BCR.Position(BCR.Live(i))-N;BCR.Position(BCR.Live(i))+N];
% exclude the site that contains already bound Ag
if length(freebcr)==1
   partner = BCR.Binding(BCR.Live(i),1);
   neighbors = setdiff(neighbors, Ag.Position(partner));
end
% define the targets for binding
targets.num = zeros(5,1); %number of target epitopes available at each neighbor site
targets.agnum = zeros(5,1); %Index of the antigens at each neighbor site (if exists)
targets.epitopes = []; %all target epitopes 
for j=1:length(neighbors) %for each neighbor
    if neighbors(j)>0 && neighbors(j)<N^2+1 %check if not out of bound
        if Ag.Lattice(neighbors(j),1)==2 %check if there is Ag
            agnum = Ag.Lattice(neighbors(j),2);
            agtype = Ag.Type(agnum);
            epitopes = Ag.Epitopes(agtype,:);
            epitopes = epitopes(Ag.Binding(agnum,:)==0); %free epitopes
            targetep = find(BCR.Energy<0); %get the BCR target epitopes
            epitopes = epitopes(ismember(epitopes,targetep)); %Ag epitopes that are targeted by the BCR
            targets.num(j) = length(epitopes);
            targets.agnum(j) = agnum;
            targets.epitopes = [targets.epitopes, epitopes];
        end
    end
end
% attempt to bind
q = q0*length(freebcr)*length(targets.epitopes); %total on-rate
Pon = 1-exp(-q*dt); %probability of binding
if rand < Pon %binding happens
    k = randi(length(targets.epitopes)); %sample random epitope
    epitope = targets.epitopes(k);
    s = find(cumsum(targets.num)>k-1,1); %find the Ag that the epitope belongs to
    agnum = targets.agnum(s);
    BCR.Binding(BCR.Live(i),[freebcr(1),freebcr(1)+2]) = [agnum, epitope]; %update binding
    %now also need to assign which of the Ag epitopes are bound to the BCR
    candidates = find(Ag.Epitopes(Ag.Type(agnum),:)==epitope); % find which of the 6 epitopes are matching
    for site = candidates 
        if Ag.Binding(agnum,site)==0 % if there is an empty epitope, assign BCR to it
            Ag.Binding(agnum,site) = BCR.Live(i);
            Bond(BCR.Live(i),agnum) = 1;
            break
        end
    end
end
end
end

%%
function [BCR, Ag, Bond] = BondBreakExtraction(BCR, Ag, Bond, dt, rt, Emem, xb, F, kbT)
%   Breaking of bonds during the extraction phase.
%   Off-rate depends on the force applied

%Iterate over all Live BCRs
for i=1:length(BCR.Live)
BCR.Refractory(BCR.Live(i)) = max(0,BCR.Refractory(BCR.Live(i))-dt); %Decrease refractory time
if any(Bond(BCR.Live(i),:)) %If bound to at least 1 Ag
    epitopes = BCR.Binding(BCR.Live(i),3:4); %Identify the epitopes bound to
    v = sum(epitopes>0); %Number of epitopes bound to - used to calculate force
    for j=1:2 %for each arm
        if epitopes(j)>0 %if the arm is bound to Ag
            koff = BCR.koff(epitopes(j))*exp(xb*F/(v*kbT)); 
            if rand < 1-exp(-koff*dt) %bond breaks
               agnum = BCR.Binding(BCR.Live(i),j); %get the Ag index
               BCR.Binding(BCR.Live(i),[j,j+2])=[0 0]; %BCR binding to 0
               Ag.Binding(agnum,(Ag.Binding(agnum,:)==BCR.Live(i))) = 0; %Ag binding to 0
               Bond(BCR.Live(i),agnum) = 0; %Update adjacency
               BCR.Refractory(BCR.Live(i)) = rt; %Update refractory time
               BCR.Rebinding(BCR.Live(i),1) = BCR.Rebinding(BCR.Live(i),1)+1; 
            end
        end
    end
    if BCR.Binding(BCR.Live(i),1)==0 && BCR.Binding(BCR.Live(i),2)>0
        BCR.Binding(BCR.Live(i),[1 3]) = BCR.Binding(BCR.Live(i),[2 4]);
        BCR.Binding(BCR.Live(i),[2 4]) = [0 0];
    end
end
end
%Iterate over all live Ags
for i=1:length(Ag.Live)
if any(Bond(:,Ag.Live(i))) %If bound to at least 1 BCR
    boundBCRs = Ag.Binding(Ag.Live(i),Ag.Binding(Ag.Live(i),:)>0); %Get the indices of bound BCRs
    %calculate the total force
    Ftot = 0;
    for j=1:length(boundBCRs)
        v = sum(BCR.Binding(boundBCRs(j),1:2)>0);
        Ftot = Ftot + F/v;
    end
    %calculate the Ag-Mem off-rate
    koffmem = 10^6/exp(-Emem)*exp(xb*Ftot/kbT);
    if rand < 1 - exp(-koffmem*dt) %Ag-mem bond breaks
        Ag.Force(Ag.Live(i)) = Ftot; %record the force at breaking
        Ag.Mem(Ag.Live(i)) = 0; %remove from live index
    end
end
end
end

%%
function [BCR, Ag, Bond] = BondBreakClustering(BCR, Ag, Bond, dt)
%   Breaking of bonds during the clustering phase.
%   Only Ag-BCR bonds (but not Ag-Mem bonds) break during this phase.

%Iterate over each BCR
for i=1:length(BCR.Live)
    if any(Bond(BCR.Live(i),:)) %If bound to at least 1 Ag
        epitopes = BCR.Binding(BCR.Live(i),3:4);
        for j=1:2
            if epitopes(j)>0
                koff = BCR.koff(epitopes(j));
                if rand < 1-exp(-koff*dt) %bond breaks
                   agnum = BCR.Binding(BCR.Live(i),j);
                   BCR.Binding(BCR.Live(i),[j,j+2])=[0 0];
                   Ag.Binding(agnum,(Ag.Binding(agnum,:)==BCR.Live(i))) = 0;
                   Bond(BCR.Live(i),agnum) = 0;
                end
            end
        end
        %If the second arm is bound but the first arm is not, swap
        if BCR.Binding(BCR.Live(i),1)==0 && BCR.Binding(BCR.Live(i),2)>0
            BCR.Binding(BCR.Live(i),[1 3]) = BCR.Binding(BCR.Live(i),[2 4]);
            BCR.Binding(BCR.Live(i),[2 4]) = [0 0];
        end
    end
end
end

%%
function [BCR, Ag] = Diffusion(clusters, BCR, Ag, N, D0, dt, l, phase)
%   Diffusion of molecules during both clustering and extraction phases
if strcmp(phase,'clustering')
    Mobile = find(sum(clusters.NumElem,2)<=3); %during clustering, clusters up to size 3 can diffuse
elseif strcmp(phase,'extraction')
    Mobile = find(clusters.NumElem(:,1)==0); %during extraction, BCRs do not diffuse
end
%   Iterate over all mobile clusters
for i=1:length(Mobile)
    % first identify the types of clusters
    Diffusivity = D0/sum(clusters.NumElem(Mobile(i),:)); %simplified assumption - 
                                                         %radius proportional to number of molecules
    p = 4*Diffusivity*dt/l^2; %probability of successful diffusion; 1 for a single molecule
    if rand > p
        continue
    end
    % get the positions of the BCRs and Ags
    BCRpos = BCR.Position(clusters.Elem{Mobile(i),1}); 
    Agpos = Ag.Position(clusters.Elem{Mobile(i),2});

    % randomly choose which direction to move
    switch randi(4)
        case 1 %up
            BCRposNew = BCRpos - 1;
            AgposNew = Agpos - 1;
        case 2 %down
            BCRposNew = BCRpos + 1;
            AgposNew = Agpos + 1;                
        case 3 %left
            BCRposNew = BCRpos - N;
            AgposNew = Agpos - N;                        
        case 4 %right
            BCRposNew = BCRpos + N;
            AgposNew = Agpos + N;
    end

    % move molecules
    if all([BCRposNew;AgposNew]>0) && all([BCRposNew;AgposNew]<N^2+1) % check if not out of bound
        if all(BCR.Lattice(setdiff(BCRposNew,BCRpos),1)==0) && ... 
                all(Ag.Lattice(setdiff(AgposNew,Agpos),1)==0) %check if the new position is empty 
                                                              % or will be emptied during move
            BCR.Position(clusters.Elem{Mobile(i),1}) = BCRposNew;
            Ag.Position(clusters.Elem{Mobile(i),2}) = AgposNew;
            BCR.Lattice(BCRpos,1:2) = repmat([0,0],length(BCRpos),1);
            Ag.Lattice(Agpos,1:2) = repmat([0,0],length(Agpos),1);
            BCR.Lattice(BCRposNew,1) = 1;
            BCR.Lattice(BCRposNew,2) = clusters.Elem{Mobile(i),1};
            Ag.Lattice(AgposNew,1) = 2;
            Ag.Lattice(AgposNew,2) = clusters.Elem{Mobile(i),2};
        end
    end
end
end

%%
function h = Visualize(N, BCR, Ag, Bond)
% Visualize the cluster interface by representing BCRs as blue and Ags as Red  
BCRlivepos = BCR.Position(BCR.Live);
Aglivepos = Ag.Position(Ag.Live);
[i1, j1] = ind2sub([N N], BCRlivepos);
[i2, j2] = ind2sub([N N], Aglivepos);
blocked = BCR.Lattice==3;
[i3,j3] = ind2sub([N,N],find(blocked));
h = figure;
axes('NextPlot','add','XLim',[1 N],'YLim',[1 N]);
axis square;
hold on
plot(i1, j1, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 1.5);
plot(i2, j2, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 1.5);   
plot(i3, j3, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 1.5);
end

function Ag = AbFeedback(Ag)
% A scenario in which Ab feedback randomly blocks some epitopes of the Ags
for i=1:length(Ag.Binding)
   for j=1:length(Ag.Binding(i,:))
       if rand < 1/3
           Ag.Binding(i,j) = -1;
       end
   end
end
end