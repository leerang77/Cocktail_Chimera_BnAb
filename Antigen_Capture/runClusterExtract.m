function runClusterExtract(antigen, E1, E2, E3, R0, L0, first, last)
for runnum=first:last
    ClusterExtract(antigen, E1, E2, E3, R0, L0, runnum)
end
end