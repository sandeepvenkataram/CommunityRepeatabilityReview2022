function newC = GenNewMutants(C, MutModel, V, U, dt)

% GENNEWMUTANTS generate new mutations
%
% newC = GenNewMutants(C, MutModel, V, U, dt) updates the vector of current
% consumer abundances by adding new mutants. Mutants are produced as a
% Poisson process with rate N * V * U where N is the current consumer
% density, V is the volume and U is the per capita per unit time mutation
% rate. All N, V and U are assumed static.
% 
% INPUT
%
% C is a structure that holds all necessary information about current
% consumer state with the following fields: g is an nC by nR matrix of per
% capita resource consumption rates, i.e., g(i,j) is the rate with which an
% individual of species i converts resource j into its own biomass; d is an
% nC by 1 vector of death rates (all d values have to be greater or equal
% to 1 to account for diluation); alpha is an nC by nR matrix of per capita
% resource production rates, i.e., alpha(i,j) is the rate of production of
% resource j by an individual of type i; N is an nC by 1 vector of
% current consumer abundances. MutModel is the mutation model generated by
% GenMutModel. V is the volume of the chemostat. U is the mutation rate per
% unit time. dt is the amout of time during which mutants are generated.
%
%
% OUTPUT
%
% newC is the updated consumer structure (same as input C) with mutants
% added


% Expected number of mutants in each strain:
ExpMut = C.N * V * U * dt;

% Actual number of mutants in each strain:
nMutVec = poissrnd(ExpMut);

% Total number of new mutants:
nMut = sum(nMutVec); 

if nMut == 0
    newC = C;
    return;
end

% Reduce the number of parent who got converted into mutants and add
% mutants with so far blank parameters
newC.N = [C.N - nMutVec/V; ones( nMut, 1)/V ];
newC.g = [C.g ; zeros(nMut, 2)];
newC.d = [C.d ; ones(nMut, 1)];
newC.alpha = [C.alpha ; zeros(nMut, 2)];

% Indices of strains that generated mutants:
parentix = find(nMutVec > 0);

% Index of the current mutant:
mutix = size(C.g,1)+1;

% Go through all parents that spawned mutants and populate mutants'
% traits
for ix = parentix'
    
    % Number of mutants spawned by the current parent:
    nMutCurr = nMutVec(ix);
    
    % These are the traits of the new spawed mutants:
    newC.g(mutix:mutix+nMutCurr-1,:) = repmat(C.g(ix,:), nMutCurr,1) + MutModel( nMutCurr );
    
    % Update the index of the current mutant:
    mutix = mutix+nMutCurr;
end


