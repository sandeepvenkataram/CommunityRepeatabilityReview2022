function EcolEquil = GetEcolEquil(gMat, Rin)

% GETECOLEQUIL calculate the co-existence equilibrium with two consumers
%
% EcolEquil = GetEcolEquil(gMat, Rin) computes the co-existence equilibrium
% for the special case of two consumers, no resource production by
% consumers, no death in excess of dilution and 1-to-1 stoichiometry
% 
% INPUT
%
% gMat is an 2 by 2 matrix of per capita resource consumption rates,
% i.e., g(i,j) is the rate with which an individual of species i converts
% resource j into its own biomass; Rin is an 2 by 1 vector of input
% resource concentrations
%
% OUTPUT
%
% EcolEquil is a structure with two fields: R is a 2 by 1 vector of
% equilibrium resource concentrations; C is the 2 by 1 vector of consumer
% abundances


EcolEquil.R = gMat \ [1;1] ;
EcolEquil.C = gMat' \ (Rin./EcolEquil.R - 1);
