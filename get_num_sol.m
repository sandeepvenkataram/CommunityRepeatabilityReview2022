function [Nvec, Rvec, tvec] = get_num_sol(Model, N0, R0, Tf)

% GET_NUM_SOL numerically solves the system of consumer-resource model
% equations
%
% [Nvec, Rvec, tvec] = get_num_sol(Model, N0, R0, Tf) numerically solves
% the system of consumer-resource model equations.
% 
% INPUT
%
% Model is the consumer-resource model specified by function GenCRModel; N0
% is an nC by 1 vector of initial consumer abundances; R0 is an nR by 1
% vector of initial resource concentrations; Tf is the time horizon over
% which to solve the equations.
%
% OUTPUT
%
% Nvec is an nT by nC matrix of consumer abundances over time; Rvec is the nT by nR
% matrix of resource concentrations over time; tvec is an nT by 1 vector of
% time points at which the solution was evaluated.

nC = Model.nC;
% nR = Model.nR;
iR = nC + 1;
tspan = [0 Tf];
y0 = [N0; R0];

opts = odeset('Jacobian', @J, 'NonNegative', 1:length(y0));

[tvec, yvec] = ode15s(@(t,y) [Model.dNdt( y(1:nC), y(iR:end)); Model.dRdt( y(1:nC), y(iR:end))], tspan, y0, opts);

Nvec = yvec(:,1:nC);
Rvec = yvec(:,iR:end);

function r = J(t,y)
    N = y(1:nC);
    R = y(iR:end);
    
    r = Model.J(N,R);
end

end