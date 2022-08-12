function Model = GenCRModel(Rin, g, d, alpha, type, gamma)

% GENCRMODEL consumer-resource model with substitutable resources at
% steady supply

% Model = GenCRModel(Rin, g, d, alpha, type, gamma) creates the model to
% calculate the derivatives of the number of consumer individuals for nC
% consumer species and the derivatives of the resource densities for nR
% resources in a chemostat-type model. Time is measured in the units of
% whole volume replacements.
% 
% INPUT
%
% Rin is a nR by 1 vector of influx resource densities; g is an nC by nR
% matrix of per capita resource consumption rates, i.e., g(i,j) is the rate
% with which an individual of species i converts resource j into its own
% biomass; d is an nC by 1 vector of death rates (all d values have to be
% greater or equal to 1 to account for diluation); alpha is an nC by nR
% matrix of per capita resource production rates, i.e., alpha(i,j) is the
% rate of production of resource j by an individual of type i; type is an
% nC by 1 vector of initial consumer numbers.
%
% Model = CR_subst_res(..., gamma), gamma is an nC by nR matrix of
% stoichiometry coefficients, i.e., gamma(i,j) is the amount of resource
% j that needs to be consumed by an species of type i to create one
% individual. By default this is a matrix of ones.
%
% OUTPUT
%
% Model is a structure with the following fields: nC is number of consumer
% types; nP is number of resource types; pcdNdt(R) is a function that
% calculates the PER CAPITA rate of change of the number of consumer
% individuals given the current resource abundances R; dNdt(N,R) is a
% function of N (vector of consumer species numbers) and R (vector of resource
% concentrations) that gives the rate of change in the number of consumers;
% dPdt(N,R) is a function that calculates the derivative of the resource
% abundances given the current numbers N and resource abundences R; J(N,R)
% is a function that returns the jacobian matrix; params holds all model
% parameters; types holds the current number of consumer individuals.


Model.nC = size(g,1);
Model.nR = size(g,2);

Model.params.Rin = Rin;
Model.params.g = g;
Model.params.d = d;
Model.params.alpha = alpha;
Model.type = type;

if nargin < 6
    gamma = ones( Model.nC, Model.nR );
end

Model.params.gamma = gamma;

Model.pcdNdt = @pcdNdt;
Model.dNdt = @dNdt;
Model.dRdt = @dRdt;
Model.J = @jacobian;


    function res = pcdNdt(R)
        res = g * R - d;
    end

    function res = dNdt(N,R)
        res = N .* pcdNdt(R);
    end

    function res = dRdt(N,R)
        res = Rin + alpha' * N - R .* (1 + (gamma .* g)' * N );
    end

    function res = jacobian(N,R)        
        NN = diag( pcdNdt(R) );
        NR = repmat(N,1,Model.nR) .* g;
        RR = diag( - (1 + (gamma .* g)' * N )  );
        RN = alpha' - g' .* gamma' .* repmat(R,1,Model.nC);
        res = [NN NR ; RN RR];
    end

end

