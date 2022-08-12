function M = GenMutModel(Name, Params)

% MUTMODEL mutation-generating model

% M = MutModel(Name, Params) creates a generative model of mutations in a
% n-dimensional trait space.
% 
% INPUT
%
% Name (string) is the type of mutational model: 'Modular', 'Discrete' or
% 'Gaussian'. Params specifies the parameters of the model.

% For 'Modular', Params is an 2 by n matrix. Params(1,:) is a vector of
% probabilities that the mutation affects one of the species traits, where
% n is the total number of traits (sum(Params(1,:)) = 1). In other words,
% if a mutation happens, it affects trait i (and no other traits) with
% probability Params(1,i). The second row of Params specifies the increment
% by which each trait changes if it gets mutated;
%
% For 'Discrete', Params is a n+1 by k matrix. First row is the row of
% probabilities, such that sum(Params(1,:)) = 1. Each column starting from
% the second row specifies the coordinates of the corresponding mutational
% vector.
%
% For 'Gaussian', Params is an n by n variance-covariance matrix. In other
% words, each mutation affects all traits with the probabilities given by
% the Gaussian distribution with mean vector zero and variance-covariance
% matrix Params.
%
% OUTPUT
%
% M(N) is a function that generates and N by n matrix of mutants, where N
% is the number of mutants and each mutant is described by n traits whose
% values are drawn randomly from the specified model.


if strcmp(Name, 'Modular')
    M = @modular;
elseif strcmp(Name, 'Gaussian')
    M = @gaussian;
elseif strcmp(Name, 'Discrete')
    M = @discrete;
end



function r = modular(nmut)
    if Params(1,1) == 1
        r = [Params(2,1)*randn(nmut,1), zeros(nmut,1)];
    elseif Params(1,2) == 1
        r = [zeros(nmut,1), Params(2,2)*randn(nmut,1)];
    else
        r = zeros(nmut, 2);
        TF = rand(nmut,1) < Params(1,1);
        n1 = nnz(TF);
        n2 = nnz(~TF);
        r(TF,:) = [Params(2,1)*randn(n1,1), zeros(n1,1)];
        r(~TF,:) = [zeros(n2,1), Params(2,2)*randn(n2,1)];
    end
end


function r = discrete(nmut)
    kvec = mnrnd(nmut, Params(1,:));
    
    ixvec = find(kvec > 0);
    
    r = nan(nmut, size(Params,2));
    currix = 1;
    for ix = ixvec
        CurrNmut = kvec(ix);
        r(currix:currix+CurrNmut-1,:) = repmat(Params(2:end,ix)', CurrNmut, 1);
        currix = currix+CurrNmut;
    end
end


function r = gaussian(nmut)
    r = mvnrnd([0, 0], Params, nmut);
end

end