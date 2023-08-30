%--------------------------------------------------------------------------
% This function takes the groups resulted from spectral clutsering and the
% ground truth to compute the misclassification rate.
% groups: [grp1,grp2,grp3] for three different forms of Spectral Clustering
% s: ground truth vector
% Missrate: 3x1 vector with misclassification rates of three forms of
% spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------


function [Missrate, index] = Misclassification(groups,s)

n = max(s);
index = zeros(size(groups, 2), n);
for i = 1:size(groups,2)
    [miss, index(i,:)] = missclassGroups( groups(:,i),s,n );
    Missrate(i,1) = miss ./ length(s); 
end