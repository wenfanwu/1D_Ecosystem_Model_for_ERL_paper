% ==================================================================================
% This demo aims to study the response of ecosystem to different SPM conditions using a 1D ecosystem model.
% It is also part of Wu et al. (2023, Weak local upwelling may elevate the risks of harmful algal blooms and hypoxia 
% in shallow waters during the warm season, Environ. Res. Lett. 18(11), 114031).
%
% Author Info: Wenfan Wu, Postdoc, Virginia Institute of Marine Science, William & Mary | wwu@vims.edu
% ==================================================================================
%% Step-1: Model Config
clc;clearvars
Mobj.nml_file = 'C:\Users\15641\OneDrive - vims.edu\GitHub_Projects\1D_Ecosystem_Model_for_ERL_paper\cosine.nml';  % NEED TO BE CHANGED!
Mobj.latitude = 39.5;  % site latitude

Mobj.dt = 0.025;  % day
Mobj.nyears = 0.3;  % year
Mobj.dz = 0.25;   % meter
Mobj.depmax = 25;  % meter

Mobj.depth = (Mobj.dz:Mobj.dz:Mobj.depmax)';  % depth layers
Mobj.nLevs = Mobj.depmax/Mobj.dz+1;   % # of vertical levels

%% Step-2: Initial Conditions (Bohai Sea)
InitCnd.PO4 = 0.3; InitCnd.NO3 = 8; InitCnd.SiO4 = 8;  InitCnd.NH4 = 0.8;
InitCnd.S1 = 0.1; InitCnd.S2 = 0.1; InitCnd.Z1 = 0.1; InitCnd.Z2 = 0.1;
InitCnd.DN = 0.1; InitCnd.DSi = 0.1; InitCnd.DOX = 250;

% InitCnd.SPM = ones(Mobj.nLevs, 1)*10.5;  % high SPM case
InitCnd.SPM = ones(Mobj.nLevs, 1)*2;    % low SPM case

nTimes = 365;
InitCnd.temp_time = (1:nTimes)-2;
InitCnd.temp = 22*ones(Mobj.nLevs, nTimes);   % constant temperature

InitCnd.mld_time = (1:nTimes)-2;
InitCnd.mld = 5*ones(1, nTimes);   % mixed layer depth (5 m)

InitCnd.bbl_time = (1:nTimes)-2;
InitCnd.bbl = (Mobj.depmax-5)*ones(1, nTimes);  % bottom boundary layer thickness (3 m)

%% Step-3: Model Run
TA = CoSiNE_1D_main(Mobj, InitCnd);

%% Step-4: Visualize
nyears = Mobj.nyears;                     % # of years
dt = Mobj.dt;                                     % time step (days)
nsteps = round(365*nyears/dt);        % number of time steps
time = (1:nsteps)*dt;
depth = 0:Mobj.dz:Mobj.depmax;  % depth for visualize
load('cmap_data.mat');  % colormap

varList = {'NO3', 'SiO4', 'NH4', 'PO4', 'S1', 'S2', 'Z1', 'Z2', 'DN', 'DSi', 'DOX'}; sn = 0;
figure('Position', [447,705,994,317], 'Color', 'w')
for ind = [1 2 4 5 6 9 10 11]
    sn = sn +1;
    varName = varList{ind};
    varTmp = squeeze(TA(:, ind,:));
    switch varName
        case 'NO3'; varLim = [0 18]; varTicks = 0:3:18;
        case 'SiO4'; varLim = [0 15]; varTicks = 0:3:15;
        case 'NH4'; varLim = [0 0.25]; varTicks = 0:0.05:0.25;
        case 'PO4'; varLim = [0 1.25]; varTicks = 0:0.25:1.25;
        case 'S1'; varLim = [0 0.5]; varTicks = 0:0.1:0.5;
        case 'S2'; varLim = [0 1.5]; varTicks = 0:0.3:1.5;
        case 'Z1'; varLim = [0 0.25]; 
        case 'Z2'; varLim = [0 0.25];
        case 'DN'; varLim = [0 0.5]; varTicks = 0:0.1:0.5; 
        case 'DSi'; varLim = [0 0.5]; varTicks = 0:0.1:0.5;
        case 'DOX'; varLim = [135 325];  varTicks = 135:38:325;
    end

    h = subplot(2,4, sn);
    h.Position = h.Position + [0 -0.027 0.006 0.09];
    pcolor(time, -depth, varTmp)
    shading interp
    cbar = colorbar('Location', 'eastoutside', 'Color','k', 'Ticks', varTicks);
    cbar.Position(3) = cbar.Position(3)/1.5;
    cbar.Position = cbar.Position +[0.048 0 0 0];
    colormap(cmap)
    caxis(varLim) %#ok<CAXIS>
    text(0.04, 0.98, varName, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
    set(gca, 'LineWidth', 1, 'TickDir', 'in', 'Layer', 'Top')
    box on; grid on
    xlim([0 100])
    xticks(0:15:100)
    switch varName
        case {'NO3', 'S2'}
            ylabel('Depth (m)', 'FontWeight', 'bold', 'FontSize', 11)
        otherwise
            yticklabels('')
    end
    switch varName
        case {'NO3', 'SiO4', 'PO4', 'S1'}
            xticklabels('')
        otherwise
            % xlabel('Days', 'FontWeight', 'bold', 'FontSize', 11)
    end
end
%% END