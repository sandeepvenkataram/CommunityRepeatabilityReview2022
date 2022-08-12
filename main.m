%%% SIMULATION OF ECO-EVOLUTIONARY DYNAMICS IN A TWO-SPECIES CONSUMER
%%% RESOURCE MODEL WITH TWO SUBSTITUTABLE RESOURCES
%%%
%%% Venkataram S, Kryazhimskiy S (2022). Phil Trans B (SUBMITTED)
%%%
%%% Notes:
%%% (*) Only consumer 2 can evolve
%%% (*) When a new parameter set is chosen (as opposed to loaded from
%%% previous runs), a new run_id is automatically assigned. If the saving
%%% option is activated, all data and figures will be saved in the folder
%%% named by the run_id.
%%% (*) The ecological dynamics (w/o evolution) are plotted only the first
%%% time. If data is loaded, only the eco-evolutionary dynamics are plotted
%%%
%%% INPUT:
%%% Specify the variable parameters in the structures p and init
%%%
%%% OUTPUT:
%%% The main output variable that is generated is structure Output, which
%%% has the following fields.
%%%
%%% Output.tvec is a p.nT by 1 vector of sampling time points
%%% Output.Rvec is a p.nT by 2 by p.Nrep array of resource concentrations
%%% Output.Nvec is a p.nT by 2 by p.Nrep array of consumer densities
%%% Output.meanGve is a p.nT by 2 by p.Nrep array of mean trait values for consumer 2
%%% Output.pvec is a p.nT by 2 by p.Nrep array of relative abundances of both consumer species
%%% Output.s is a p.nT by 1 vector of state repeatabilities
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maindir = '.';

IFLOAD = true; % load a previous simulation (for plotting)?
IFSAVE = false; % save new data and plots?

if IFLOAD
    
    %%% UNCOMMENT TO USE EXISTING DATA:
    % run_id = "HR"; % high repeatability
    % run_id = "LR"; % low repeatability
   
    filename = sprintf('%s/figures/run_%s/data.mat', maindir, run_id);
    load(filename);
    close all;
else
    
    % Generate run id
    run_id = regexprep(string(datetime),':','');       
    
%% ===== Specify variable parameters =====
    p.U = 2e-4; % Rate of mutations (per unit time)
    p.Nrep = 5; % Number of replicate runs
    p.V = 1000;  % Volume of the chemostat (necessary to calculate abolute consumer numbers)
    p.Tecol = 30; % time to reach ecological equilibrium
    p.Tmax = 2e3; % duration of the experiment
    p.Rin = [0.5 ; 1]; % influx rate of resources
    p.tvec = (0:50:p.Tmax)'; % vector of sampling time points
    p.nT = length(p.tvec); % Number of time points
    p.nSample = 5; % Number of samples to compute repeatability
    
    %%% Consumer 1 traits and density at the initial time point (this consumer does NOT evolve, but its population size can change)
    init.C1.g = [3 1];     % C1 mostly limited by R2
    init.C1.d = 1;         % no death other than by dilution
    init.C1.alpha = [0 0]; % no resource production
    init.C1.N = 10/p.V;    % initial number of type 1 consumer species
    
    %%% Consumer 2 traits and density at the initial time point (this consumer can evolve)
    init.C2.g = [1 2];     % C2 consumes mostly R2
    init.C2.d = 1;         % no death other than by dilution
    init.C2.alpha = [0 0]; % no resource production
    init.C2.N = 10/p.V;    % initial number of type 2 consumer species
    
    %%% Initial resource concentrations
    init.R = [0; 0];
    
    
    %%% Mutation model (choose one)    
    %%%% Low repeatability:
    %     m1 = [0.08 0.0]';    % First mutation type (lambda increments)
    %     m2 = [-0.05 0.10]';  % Second mutation type (lambda increments)
    %     prob = [0.26 0.74];  % Probabilities of two mutation types
    %     p.MutParams = [prob; [m1 m2]];
    %     GenMut = GenMutModel('Discrete', p.MutParams);
    %     clear m1 m2 prob;
    
    
    %%%% High repeatability
    m1 = [0.08 0.0]';    % First mutation type (lambda increments)     
    m2 = [-0.05 0.10]';  % Second mutation type (lambda increments)
    prob = [0.8 0.2];    % Probabilities of two mutation types
    p.MutParams = [prob; [m1 m2]];
    GenMut = GenMutModel('Discrete', p.MutParams);
    clear m1 m2 prob;
%%% ===== END Specify variable parameters =====


%% ===== Define the ecological model =====
    Model = GenCRModel(p.Rin, [init.C1.g; init.C2.g], [init.C1.d; init.C2.d], ...
        [init.C1.alpha; init.C2.alpha], [1; 2*ones(size(init.C2.N))]);
%%% ==== END Define the ecological model =====    

    
%% ===== Relax the community towards its ecological equilibrium =====
%%% This is done so that evolution begins at or close to an ecological
%%% equilibrium.
%%%
%%% *** NOTE that, once the initial conditions are updated, the data
%%% generated in this section are not saved even if the IFSAVE option is
%%% TRUE.

    %%% Find ecological equilibrium analytically (for testing purposes)
    EcolEquil = GetEcolEquil([init.C1.g; init.C2.g], p.Rin); % Ecological equilibrium for two species with d = 1;
    
    %%% Propagating model up to p.Tecol
    [Nvec, Rvec, tvec] = get_num_sol(Model, [init.C1.N; init.C2.N], init.R, p.Tecol);
    
    %%% Update initial conditions:
    init.C1.N = Nvec(end,1);
    init.C2.N = Nvec(end,2:end)';
    init.R = Rvec(end,:)';
%%% ===== END Relax the community towards its ecological equilibrium =====
    
    
        
%% =====  Plotting ecological dynamics =====
    close all;
    
    %%% Specify the following dimensions:
    fdim.spwa = 6; % subplotwidth in cm
    fdim.spha = 5; % subplotheight in cm
    
    fdim.nx = 1; % number of panels along the horizontal dimension
    fdim.ny = 2; % number of panels along the vertical dimension
    
    fdim.xma = [1.5 0.5]; % left right horizontal margin in cm
    fdim.yma = [1.3 0.5]; % bottom top vertical margin cm
    
    fdim.dxa = 0.3; % horizontal distance between panels in cm
    fdim.dya = 0.4; % vertical distance between panels in cm
    
    fdim.tickfs = 8;
    fdim.labelfs = 10;
    
    %%% These will be computed automatically:
    fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
    fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);
    
    fdim.spwr = fdim.spwa / fdim.fw;
    fdim.sphr = fdim.spha / fdim.fh;
    fdim.xmr = fdim.xma / fdim.fw;
    fdim.ymr = fdim.yma / fdim.fh;
    fdim.dxr = fdim.dxa / fdim.fw;
    fdim.dyr = fdim.dya / fdim.fh;
    
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [0 0 fdim.fw fdim.fh]);
    
    cc = [
        0, 114, 178;    % blue
        213, 94, 0;     % vermillion
        86, 180, 233;   % sky blue
        230 159, 0;     % orange
        204, 121, 167;   % raddish purple
        0, 158, 115;    % bluish green
        240, 228, 66   % yellow
        ]./256;
    
    fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
    fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );
        
    %%% Plotting consumers:
    subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr]),
    hold on, box on;
    set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    
    for ic = 1:2
        plot( tvec, Nvec(:,ic)*p.V, '-', 'LineWidth', 3, 'Color', cc(ic,:) );
        plot( tvec, EcolEquil.C(ic)*ones(size(tvec))*p.V, '--', 'LineWidth', 1, 'Color', cc(ic+2,:));
        text(tvec(end), EcolEquil.C(ic)*p.V, sprintf('C_%d', ic),...
            'FontSize', fdim.labelfs, 'FontName', 'Helvetica', 'Color', cc(ic,:),...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
    clear ic;
        
    set(gca, 'XTickLabel', {}, 'YLim', [-30 630], 'XLim', [-2 32]);
    
    ylabel('Consumer population size', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
    
    
    %%% Plotting resources:
    subplot('Position', [fdim.spxvec(1) fdim.spyvec(2) fdim.spwr fdim.sphr]),
    hold on, box on;
    set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    
    for ic = 1:2
        plot( tvec, Rvec(:,ic), '-', 'LineWidth', 3, 'Color', cc(ic,:) );
        plot( tvec, EcolEquil.R(ic)*ones(size(tvec)), '--', 'LineWidth', 1, 'Color', cc(ic+2,:));
        text(tvec(end), EcolEquil.R(ic), sprintf('R_%d', ic),...
            'FontSize', fdim.labelfs, 'FontName', 'Helvetica', 'Color', cc(ic,:),...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
    clear ic
    
    set(gca, 'XLim', [-2 32], 'YLim', [-0.02 1.02]);
    
    ylabel('Resource cooncentration', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
    xlabel('Time', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');%, 'Position', [55, -0.4]);
    
    clear fdim cc;
    clear Nvec Rvec tvec;
    
    if IFSAVE
        mkdir(sprintf('%s/figures/run_%s', maindir, run_id));
        filename = sprintf('%s/figures/run_%s/EcolEquil.eps', maindir, run_id);
        % filename = sprintf('%s/figures/run_1_EcolEquil.eps', maindir);
        saveas(gcf, filename, 'epsc');
    end
%%% ===== END Plotting ecological dynamics =====
    



    
%% ===== Eco-evolutionary dynamics =====  
    Output = struct('tvec', p.tvec, 'Rvec', nan(p.nT, 2, p.Nrep),...
        'Nvec', nan(p.nT, 2, p.Nrep), 'meanGvec', nan(p.nT, 2, p.Nrep),...
        'pvec', nan(p.nT, 2, p.Nrep), 's', nan(p.nT, 1));
    
    for irep = 1:p.Nrep
        
        fprintf('Rep %d\n', irep);
    
        % Initial condition
        Output.Nvec(1,1,irep) = init.C1.N;
        Output.Nvec(1,2,irep) = init.C2.N;
        Output.meanGvec(1,:,irep) = init.C2.g;
        Output.Rvec(1,:,irep) = init.R';
        Output.pvec(1,:,irep) = Output.Nvec(1,:,irep) / sum(Output.Nvec(1,:,irep));
        
        % Current resources
        R = init.R;
        
        % Current consumer 1 traits and density
        C1.g = init.C1.g;
        C1.d = init.C1.d;
        C1.alpha = init.C1.alpha;
        C1.N = init.C1.N;
        
        % Current consumer 2 traits and density
        C2.g = init.C2.g;
        C2.d = init.C2.d;
        C2.alpha = init.C2.alpha;
        C2.N = init.C2.N;
        
        t = 0;
        tix = 2; % index of the next sampling time point
        
        while t < p.Tmax
            
            % Check that at least one species is still alive
            if C1.N == 0 && sum(C2.N) == 0
                printf('Extinction of both species');
                break;
            else
                % Define the current model (with possibly multiple strains of consumer 2)
                Model = GenCRModel(p.Rin, [C1.g; C2.g], [C1.d; C2.d], [C1.alpha; C2.alpha], [1; 2*ones(size(C2.N))]);
            end
            
            % Expected number of mutations per unit time:
            TotExpMutRate = sum(C2.N) * p.U * p.V;
            
            % Simulate the ecological model for this number of time units
            % (so that we have ~1 new mutant during this run)
            dt = 1/TotExpMutRate;
            
            % Book keeping so that we sample at the right time and stop the
            % simulation at the right time
            if t+dt > p.Tmax
                dt = p.Tmax - t;
            elseif t+dt > p.tvec(tix)
                dt = p.tvec(tix) - t;
            end
            
            % Propagating consumers and resources to t+dt
            [Nvec, Rvec, tvec] = get_num_sol(Model, [C1.N; C2.N], R, dt);
            C1.N = Nvec(end,1);
            C2.N = Nvec(end,2:end)';
            R = Rvec(end,:)';
            t = t + dt;
            
            % New mutations:
            C2 = GenNewMutants(C2, GenMut, p.V, p.U, dt);
            
            % Purging consumers that are below 1 individual per volume:
            [C1, C2] = PurgeConsumers(C1, C2, p.V);
            
            % Recording to output:
            if t == p.tvec(tix)
                Output.Rvec(tix,:,irep) = R';
                Output.Nvec(tix,1,irep) = C1.N;
                Output.Nvec(tix,2,irep) = sum(C2.N);
                Output.meanGvec(tix,:,irep) = mean(C2.g,1);
                Output.pvec(tix,:,irep) = Output.Nvec(tix,:,irep) / sum(Output.Nvec(tix,:,irep));
                
                tix = tix+1;
            end
        end
    end
    clear irep R C1 C2 t dt TotExpMutRate Nvec Rvec tvec;
    
    % Calculate repeatability index
    Output.s = GetRepeatability(Output.pvec, p.nSample);
    
    if IFSAVE
        filename = sprintf('%s/figures/run_%s/data.mat', maindir, run_id);
        save(filename, 'p', 'init', 'GenMut', 'Model', 'Output');
    end
%% ===== END Eco-evolutionary dynamics =====
end



%% ===== Plotting eco-evolutionary dynamics =====
figure;

%%% Specify the following dimensions:
fdim.spwa = 6; % subplotwidth in cm
fdim.spha = 5; % subplotheight in cm

fdim.nx = 1; % number of panels along the horizontal dimension
fdim.ny = 3; % number of panels along the vertical dimension

fdim.xma = [1.5 0.5]; % left right horizontal margin in cm
fdim.yma = [1.3 0.5]; % bottom top vertical margin cm

fdim.dxa = 0.3; % horizontal distance between panels in cm
fdim.dya = 0.4; % vertical distance between panels in cm

fdim.tickfs = 8;
fdim.labelfs = 10;

%%% These will be computed automatically:
fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);

fdim.spwr = fdim.spwa / fdim.fw;
fdim.sphr = fdim.spha / fdim.fh;
fdim.xmr = fdim.xma / fdim.fw;
fdim.ymr = fdim.yma / fdim.fh;
fdim.dxr = fdim.dxa / fdim.fw;
fdim.dyr = fdim.dya / fdim.fh;

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [8 0 fdim.fw fdim.fh]);

cc = [
    0, 114, 178;    % blue
    213, 94, 0;     % vermillion
    86, 180, 233;   % sky blue
    230 159, 0;     % orange
    204, 121, 167;   % raddish purple
    0, 158, 115;    % bluish green
    240, 228, 66   % yellow
    ]./256;

fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );


%%% Plotting consumers over time
subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');

plot( Output.tvec, squeeze(Output.Nvec(:,1,:)*p.V), '-', 'LineWidth', 1, 'Color', cc(1,:) );
plot( Output.tvec, squeeze(Output.Nvec(:,2,:)*p.V), '-', 'LineWidth', 1, 'Color', cc(2,:) );

clear irep;

set(gca, 'XTickLabel', {});

ylabel('Consumer population size', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
                


%%% Plotting resources over time:
subplot('Position', [fdim.spxvec(1) fdim.spyvec(2) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    
plot( Output.tvec, squeeze(Output.Rvec(:,1,:)), '-', 'LineWidth', 1, 'Color', cc(1,:) );
plot( Output.tvec, squeeze(Output.Rvec(:,2,:)), '-', 'LineWidth', 1, 'Color', cc(2,:) );

clear irep ic;

set(gca, 'XTickLabel', {});

ylabel('Resources density', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');






%%% Plotting traits over time
subplot('Position', [fdim.spxvec(1) fdim.spyvec(3) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    
plot( Output.tvec, squeeze(Output.meanGvec(:,1,:)*p.V), '-', 'LineWidth', 1, 'Color', cc(1,:) );
plot( Output.tvec, squeeze(Output.meanGvec(:,2,:)*p.V), '-', 'LineWidth', 1, 'Color', cc(2,:) );

clear irep ic;

ylabel('Trait values', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
xlabel('Time', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');%, 'Position', [55, -0.4]);

if IFSAVE
    filename = sprintf('%s/figures/run_%s/EvolTime.eps', maindir, run_id);
    saveas(gcf, filename, 'epsc');
end



%%% Plotting consumer phase space:
figure;

%%% Specify the following dimensions:
fdim.spwa = 6; % subplotwidth in cm
fdim.spha = 5; % subplotheight in cm

fdim.nx = 1; % number of panels along the horizontal dimension
fdim.ny = 1; % number of panels along the vertical dimension

fdim.xma = [1.5 0.5]; % left right horizontal margin in cm
fdim.yma = [1.3 0.5]; % bottom top vertical margin cm

fdim.dxa = 0.3; % horizontal distance between panels in cm
fdim.dya = 0.4; % vertical distance between panels in cm

fdim.tickfs = 8;
fdim.labelfs = 10;

%%% These will be computed automatically:
fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);

fdim.spwr = fdim.spwa / fdim.fw;
fdim.sphr = fdim.spha / fdim.fh;
fdim.xmr = fdim.xma / fdim.fw;
fdim.ymr = fdim.yma / fdim.fh;
fdim.dxr = fdim.dxa / fdim.fw;
fdim.dyr = fdim.dya / fdim.fh;

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [16 0 fdim.fw fdim.fh]);

fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );

fdim


subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    

plot( Output.Nvec(1,1,1)*p.V, Output.Nvec(1,2,1)*p.V, 'o',...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', cc(2,:), 'MarkerSize', 8);

plot( squeeze(Output.Nvec(:,1,:)*p.V), squeeze(Output.Nvec(:,2,:)*p.V), '-', 'LineWidth', 1, 'Color', 0.4*[1 1 1]);
plot( squeeze(Output.Nvec(end,1,:)*p.V), squeeze(Output.Nvec(end,2,:)*p.V), 's',...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor',   0.4*[1 1 1], 'MarkerSize', 8);


set(gca, 'XLim', [-30 330], 'YLim', [550 1100], 'XGrid', 'on', 'YGrid', 'on');

xlabel('C_1 population size', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');%, 'Position', [55, -0.4]);
ylabel('C_2 population size', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');

clear fdim cc;

if IFSAVE
    filename = sprintf('%s/figures/run_%s/EvolPhaseSpace.eps', maindir, run_id);
    saveas(gcf, filename, 'epsc');
end



%%% Plotting relative frequency over time and repeatability over time

figure;

%%% Specify the following dimensions:
fdim.spwa = 6; % subplotwidth in cm
fdim.spha = 5; % subplotheight in cm

fdim.nx = 1; % number of panels along the horizontal dimension
fdim.ny = 3; % number of panels along the vertical dimension

fdim.xma = [1.5 0.5]; % left right horizontal margin in cm
fdim.yma = [1.3 0.5]; % bottom top vertical margin cm

fdim.dxa = 0.3; % horizontal distance between panels in cm
fdim.dya = 0.4; % vertical distance between panels in cm

fdim.tickfs = 8;
fdim.labelfs = 10;

%%% These will be computed automatically:
fdim.fw = fdim.spwa * fdim.nx + fdim.dxa * (fdim.nx - 1) + sum(fdim.xma);
fdim.fh = fdim.spha * fdim.ny + fdim.dya * (fdim.ny - 1) + sum(fdim.yma);

fdim.spwr = fdim.spwa / fdim.fw;
fdim.sphr = fdim.spha / fdim.fh;
fdim.xmr = fdim.xma / fdim.fw;
fdim.ymr = fdim.yma / fdim.fh;
fdim.dxr = fdim.dxa / fdim.fw;
fdim.dyr = fdim.dya / fdim.fh;

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [24 0 fdim.fw fdim.fh]);

fdim.spxvec = fdim.xmr(1) + fdim.spwr * ( 0:(fdim.nx-1) ) + fdim.dxr * ( 0:(fdim.nx-1) );
fdim.spyvec = fdim.ymr(1) + fdim.sphr * ( (fdim.ny-1):-1:0 ) + fdim.dyr * ( (fdim.ny-1):-1:0 );

fdim



%%% Visualizing the mutation model
subplot('Position', [fdim.spxvec(1) fdim.spyvec(1) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');

for imut = 1:size(p.MutParams,2)
    plot( [0, p.MutParams(2,imut)], [0, p.MutParams(3,imut)], '-',...
        'LineWidth', 10*p.MutParams(1,imut), 'Color', 'k' );
end
    
set(gca, 'XLim', [-0.12 0.12], 'YLim', [-0.12 0.12], 'XTick', -0.1:0.1:0.1, 'YTick', -0.1:0.1:0.1);
set(gca, 'XGrid', 'on', 'YGrid', 'on');

xlabel('\Delta \lambda_{21}', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');
ylabel('\Delta \lambda_{22}', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');




%%% Plotting the relative abundance over time 
subplot('Position', [fdim.spxvec(1) fdim.spyvec(2) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    
plot( Output.tvec, squeeze(Output.pvec(:,1,:)), '-', 'LineWidth', 1, 'Color', 0.4*[1 1 1] );

set(gca, 'XLim', [-50 2050], 'YLim', [0 0.4], 'XTick', 0:500:2000, 'YTick', 0:0.1:0.4);
set(gca, 'XTickLabel', {});

ylabel('Relative abundance of C_1', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');




%%% Plotting repeatability over time 
subplot('Position', [fdim.spxvec(1) fdim.spyvec(3) fdim.spwr fdim.sphr]),
hold on, box on;
set(gca, 'FontName', 'Helvetica', 'FontSize', fdim.tickfs, 'Layer', 'top');
    
plot( Output.tvec, squeeze(Output.s), '-', 'LineWidth', 2, 'Color', 'k' );

set(gca, 'XLim', [-50 2050], 'YLim', [0.9 1]);

xlabel('Time', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');%, 'Position', [55, -0.4]);
ylabel('Repeatability <s>', 'FontSize', fdim.labelfs, 'FontName', 'Helvetica');


clear fdim cc;

if IFSAVE
    filename = sprintf('%s/figures/run_%s/EvolRepeatability.eps', maindir, run_id);
    saveas(gcf, filename, 'epsc');
end

%%% ===== END Plotting eco-evolutionary dynamics =====



