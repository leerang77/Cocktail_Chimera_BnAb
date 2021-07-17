%% Summary
%  Performs standard interpolation to determine Ag captured amount as a function of affinity
%   based on the matlab function griddedInterpolant

%  antigentype - 1: Cocktail, strain-specific
%                2: Chimeric, bnAb precursor
%                3: Cocktail, bnAb precursor
%                4: Chimeric, strain-specific
%  outputs:
%   AgCaptured - multi-dimensional cell array
%               1st dimension: type of antigen + B cell type
%               2nd dimension: (for bnabs) bind 3, 2, 1 variants
%               3rd dimension: (for bnab+cocktail) 1st, 2nd, 3rd variant
%               Each cell contains an array that stores corresponding
%               amount of antigen captured.
%   F - gridded data interpolant object. Can be called to evaluate at
%   specific points. https://www.mathworks.com/help/matlab/ref/griddedinterpolant.html
% ---------------------------------------------------

E = flip(-13.8:-0.5:-20.3);
N = length(E);
rt = 0.1;
Emem = -19;
R0 = 120;
L0 = 120;
repeats = 30;
AgCaptured{1} = zeros(N,1);
AgCaptured{4} = zeros(N,1);
tic
for j = 1:N
  tic
  AgCaptured{1}(j) = max(analyzeAgCaptured(3, E(j), 0, 0, R0, L0, Emem, rt, repeats)); 
  AgCaptured{4}(j) = analyzeAgCaptured(2, E(j), 0, 0, R0, L0, Emem, rt, repeats); 
  toc
end
toc
AgCaptured{2}{1} = zeros(N,N,N);
AgCaptured{2}{2} = zeros(N,N);
AgCaptured{2}{3} = AgCaptured{4};
AgCaptured{3}{1,1} = zeros(N,N,N); AgCaptured{3}{1,2} = zeros(N,N,N); AgCaptured{3}{1,3} = zeros(N,N,N);
AgCaptured{3}{2,1} = zeros(N,N); AgCaptured{3}{2,2} = zeros(N,N);
AgCaptured{3}{3,1} = AgCaptured{1};
for j1 = 1:N
    for j2 = j1:length(E)
        for j3 = j2:length(E)
           captured = analyzeAgCaptured(2,E(j1),E(j2),E(j3),R0, L0,Emem,rt,repeats);
           ind = perms([j1,j2,j3]);
           AgCaptured{2}{1}(sub2ind(size(AgCaptured{2}{1}),ind(:,1),ind(:,2),ind(:,3))) = captured;
        end
         captured = analyzeAgCaptured(2,E(j1),E(j2),0,R0,L0,Emem,rt,repeats);
         ind = perms([j1,j2]);
         AgCaptured{2}{2}(sub2ind(size(AgCaptured{2}{2}),ind(:,1),ind(:,2))) = captured;
    end
end
for j1 = 1:N
    for j2 = j1:length(E)
        for j3 = j2:length(E)
           ag = analyzeAgCaptured(3,E(j1),E(j2),E(j3),R0,L0,Emem,rt,repeats);
           ind = perms([j1, j2, j3]);
           ag = perms(ag);
           for k=1:3
               AgCaptured{3}{1,k}(sub2ind(size(AgCaptured{3}{1,k}),ind(:,1),ind(:,2),ind(:,3))) = ag(:,k);
           end
        end
        ag = analyzeAgCaptured(3,E(j1),E(j2),0,R0,L0,Emem,rt,repeats);
        ind = perms([j1, j2]);
        ag = perms(ag(1:2));
        for k=1:2
            AgCaptured{3}{2,k}(sub2ind(size(AgCaptured{3}{2,k}),ind(:,1),ind(:,2))) = ag(:,k);
        end
    end
end
[x1, x2] = ndgrid(E, E);
[y1, y2, y3] = ndgrid(E, E, E);
F{1} = griddedInterpolant(E, AgCaptured{1});
F{4} = griddedInterpolant(E, AgCaptured{4});
F{2}{1} = griddedInterpolant(y1, y2, y3, AgCaptured{2}{1});
F{2}{2} = griddedInterpolant(x1, x2, AgCaptured{2}{2});
F{2}{3} = griddedInterpolant(E, AgCaptured{2}{3});
for j=1:3
F{3}{1,j} = griddedInterpolant(y1, y2, y3, AgCaptured{3}{1,j});
end
for j=1:2
F{3}{2,j} = griddedInterpolant(x1, x2, AgCaptured{3}{2,j});
end
F{3}{3,1} = griddedInterpolant(E, AgCaptured{3}{3,1});
save('E2AgTable.mat','AgCaptured', 'F')