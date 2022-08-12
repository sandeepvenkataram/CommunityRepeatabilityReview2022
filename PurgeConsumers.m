function [newC1, newC2] = PurgeConsumers(C1, C2, V)

% PURGECONSUMERS removes consumers that are represented by less than 1
% individual
%
% [newC1, newC2] = PurgeConsumers(C1, C2, V) checks the absolute abundances
% of consumer 1 and all consumer 2 strains (consumer 1 has only one
% strain). If consumer 1 is represented by less than 1 individual, its
% abundance is set to zero. Those consumer 2 strains that are represented
% by less than 1 individual are removed.
%
% INPUT
%
% C1 is a structure that holds all necessary information about current
% consumer state of consumer 1 with the following fields: g is an nC by nR
% matrix of per capita resource consumption rates, i.e., g(i,j) is the rate
% with which an individual of species i converts resource j into its own
% biomass; d is an nC by 1 vector of death rates (all d values have to be
% greater or equal to 1 to account for diluation); alpha is an nC by nR
% matrix of per capita resource production rates, i.e., alpha(i,j) is the
% rate of production of resource j by an individual of type i; N is an nC
% by 1 vector of current consumer abundances. C2 is the same for consumer
% 2. V is the volume of the chemostat.
%
% OUTPUT
%
% newC1 and newC2 arethe updated consumer structures of consumers 1 and 2
% (same as inputs C1 and C2) with low-abundance strains purged

if C1.N * V < 1
    % fprintf('Consumer 1 has gone extinct');
    newC1 = C1;
    newC1.N = 0;
else
    newC1 = C1;
end

TF = C2.N * V >= 1;

if all(TF)
    newC2 = C2;
    return;
else
    newC2.N = C2.N(TF );
    newC2.g = C2.g(TF,:);
    newC2.d = C2.d(TF,:);
    newC2.alpha = C2.alpha(TF,:);
end