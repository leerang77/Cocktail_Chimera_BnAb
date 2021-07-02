function dE = getAffinityChange(Bcelltype)
% Returns energy change from a mutation
% Input
%Bcelltype: 'specific' or 'bnab'
% Output
%dE: If Bcelltype is 'specific', a scalar
%    If Bcelltype is 'bnab', a 1x3 vector
%----------------------------------------------
Constants = UseSharedConstants.Constants;
if ~isempty(Constants.MutationPDF) %sample dE from lognormal distribution 
    if strcmp(Bcelltype, 'specific')
        X = normrnd(Constants.MutationPDF(1),Constants.MutationPDF(2));
        dE = exp(X) - Constants.MutationPDF(3);
    elseif strcmp(Bcelltype, 'bnab')
        X = Constants.MutationPDF(1) + mvnrnd([0,0,0],Constants.MutationPDF(2)^2*Constants.Sigma);
        dE = exp(X) - Constants.MutationPDF(3);
    end
else %sample dE from normal distribution
    if strcmp(Bcelltype, 'specific')
        dE = normrnd(0,Constants.MutMagnitude);
    elseif strcmp(Bcelltype, 'bnab')
        dE = mvnrnd([0,0,0],Constants.MutMagnitude^2*Constants.Sigma);
    end
end
end
