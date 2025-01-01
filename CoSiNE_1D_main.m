function [TA, DA] = CoSiNE_1D_main(Mobj, InitCnd)
% A 1D ecosystem model based on the CoSiNE model
%
%% Syntax
%
%
%% Description
%
%
%% Examples
% See CoSiNE_1D_demo.m
%
%% Input Arguments
% Mobj --- a data struct contaning the configuration info.
% InitCnd --- a data struct containing the initial conditions and
% atmospheric forcing data.
%
%% Output Arguments
% TA --- all state variables simulated from the model
%  DA --- temporal variables used for diagnosis.
% diffusivity --- diffusivity data (m^2/d)
%
%% Notes
% This function is adapted from the CoSiNE module of SCHSIM modeling
% system. Please refer to the Supporting Information of Wu et al. (2024,
% ERL), or the CoSiNE manual written by Dr.Zhengui Wang
% for more details.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marin Science in 2024.
% Last Updated on 2024-08-25.
% Email: wwu@vims.edu
%
% See also:

%% Load the BGC parameters from cosine.nml
global idelay ndelay ibgraze idapt alpha_corr zeptic iz2graze ico2s ispm spm0 ised  ...
          gmaxs gammas pis kno3s knh4s kpo4s kco2s ksio4 kns alphas betas aks ...
          betaz alphaz gammaz kez kgz rhoz ...
          ipo4 TR kox wss2 wsdn wsdsi si2n p2n o2no o2nh c2n gamman pco2a kmdn kmdsi ...
          iws NO3c ws1 ws2 ...
          iclam deltaZ kcex Nperclam Wclam Fclam nclam0 ...
          fS2 fDN fDSi rkS2 rkDN rkDSi mkS2 mkDN mkDSi %#ok<*NUSED,*GVMIS>

global o2flx co2flx nh4flx sio4flx po4flx flxS2 flxDN flxDSi PS2 PDN PDSi RS2 RDN RDSi ...
    NPS1 RPS1 MTS1 NPS2 RPS2 MTS2 EXZ1 MTZ1 EXZ2 MTZ2 ...
    GS1Z1 GS2Z2 GZ1Z2 GDNZ2 GTZ2 MIDN MIDSi Nit ...
    sLight SPM ak2 ak3 ADPT OXR pih1 pih2 pnh4s1 pnh4s2 bfNO3S1 bfNH4S1 fPO4S1 fCO2S1 ...
    bfS1 bfNO3S2 bfNH4S2 fNO3S2 fNH4S2 fPO4S2 fCO2S2 fSiO4S2 bfS2

global dt time dz depth nLevs

Config = read_cosine_nml(Mobj.nml_file);
parList = fieldnames(Config);
for iVar= 1:numel(parList)
    varName = parList{iVar};
    eval([varName,' = Config.(varName);'])
end
%% Basic Paramaters
% Latitude
siteLat = Mobj.latitude;

% Temporal resolution
nyears = Mobj.nyears;                     % number of years of run
dt = Mobj.dt;                                     % time step (days)
nTimes = round(365*nyears/dt);      % number of time steps
time = (1:nTimes)*dt;

% Vertical resolution
dz = Mobj.dz;                                    % number of levels
nLevs = Mobj.depmax/dz+1;                 % number of vertical layers
depth = (1:nLevs)*dz-dz;

switch ispm
    case 0
        SPM = ones(nLevs, 1)*spm0;
    case 1
        SPM = InitCnd.SPM(:);
end

varList = {'NO3', 'SiO4', 'NH4', 'PO4', 'S1', 'S2', 'Z1', 'Z2', 'DN', 'DSi', 'DOX'};
nVars = numel(varList);

T0 = nan(nLevs, nVars);
for iVar = 1:nVars
    varName = varList{iVar};
    varTmp = InitCnd.(varName);
    if isscalar(varTmp)
        InitCnd.(varName) = ones(nLevs, 1)*varTmp;
    elseif numel(varTmp)~=nLevs
        error(['the dimension of ', varName, 'is not arropriate!'])
    end
    T0(:, iVar) = varTmp(:);
end

%% Atmopheric forcing
srad = surface_light(siteLat, time);       % Incoming light
srad_0801 = surface_light(siteLat, 210+(0:dt:1-dt));

div_val = fix(length(srad)/length(srad_0801));
res_val = length(srad)/length(srad_0801)-div_val;
srad = [repmat(srad_0801(:), fix(div_val), 1); srad_0801(1:length(srad_0801)*res_val)'];

PARfrac = 0.46;
PAR = max(PARfrac*srad, 0); 
%% Sepcify output variables
varList2 = {'NPS1', 'RPS1', 'MTS1', 'NPS2', 'RPS2', 'MTS2', 'EXZ1', 'MTZ1', 'EXZ2', 'MTZ2', 'GS1Z1', 'GS2Z2', 'MIDN', 'MIDSi', 'Nit', ...
    'bfNO3S1', 'bfNH4S1', 'fPO4S1', 'fPO4S2', 'bfNO3S2', 'bfNH4S2', 'fSiO4S2', 'bfS1', 'bfS2', 'pih1', 'pih2', 'o2flx', 'flxS2', 'flxDN', 'flxDSi'};
nVars2 = numel(varList2);
for iVar = 1:nVars2
    varName = varList2{iVar};
    DA.(varName) = nan(nLevs, nTimes);
end
%% Model run
PS2 = zeros(3,1); PDN = zeros(3,1); PDSi = zeros(3,1);
RS2 = zeros(3,1); RDN = zeros(3,1); RDSi = zeros(3,1);

T = T0;
TA = nan*ones(nLevs, nVars,nTimes);
for iTime = 1:nTimes
    progressbar(iTime/nTimes)
    MLD  = round(interp1(InitCnd.mld_time, InitCnd.mld, mod(time(iTime),365)));          % mixed layer depth in meters
    BBL  = round(interp1(InitCnd.bbl_time, InitCnd.bbl, mod(time(iTime),365)));          % bottom boundary layer depth in meters
    temp_prof = interp1(InitCnd.temp_time, InitCnd.temp', mod(time(iTime), 365))';            % interpolate temp profile in time

    T = physics(T, dt, dz, MLD, BBL);  % Vertically diffusion
    T = biology(T, dz, dt, temp_prof, PAR(iTime)); % Biological model

    TA(:,:, iTime) = T;
    if nargout == 2
        for iVar = 1:nVars2
            varName = varList2{iVar};
            eval(['varTmp = ',varName,';']) %#ok<EVLEQ> 
            DA.(varName)(:, iTime) = varTmp;
        end
    end
end

end

function T = physics(T, dt, dz, MLD, BBL)
% Physical part of the model
n = size(T,1);                                               % number of levels

% Values of diffusivities in surface/bottom mixing layers
nu1 = 0.000003*60*60*24;       % upper layer diffusivities (Converted from m2/s to m2/d)
nu2 = 0.000001*60*60*24;       % middle layer diffusivities (Converted from m2/s to m2/d)
nu3 = 0.000003*60*60*24;       % lower layer diffusivities (Converted from m2/s to m2/d)

% Define vector of diffusivities for all (n+1) cell boundaries
Zcell = dz*[0:n]';                                           %#ok<*NBRAK1> % Define cell boundaries
Nu = nu1*ones(n+1,1); 
Nu(Zcell>MLD) = nu2;
Nu(Zcell>BBL) = nu3;  % bottom boundary layer

% Define Crank-Nicolson matrices 
M = diag(-[Nu(1:n)+Nu(2:end)])+diag(Nu(2:end-1),-1)+diag(Nu(2:end-1),1);
M(1,1) = -Nu(2);                                                    % no flux top BC
M(end,end) = -Nu(end-1);                                     % no flux bottom BC
M1 = eye(n)-dt/2/dz^2*M;
M2 = eye(n)+dt/2/dz^2*M; 

% Update the state variables
T = M1\(M2*T);
end

function T = biology(T, dz, dt, temp_prof, PAR)
% BGC model

global idelay ndelay ibgraze idapt alpha_corr zeptic iz2graze ico2s ispm spm0 ised  ...
          gmaxs gammas pis kno3s knh4s kpo4s kco2s ksio4 kns alphas betas aks ...
          betaz alphaz gammaz kez kgz rhoz ...
          ipo4 TR kox wss2 wsdn wsdsi si2n p2n o2no o2nh c2n gamman pco2a kmdn kmdsi ...
          iws NO3c ws1 ws2 ...
          iclam deltaZ kcex Nperclam Wclam Fclam nclam0

global o2flx co2flx nh4flx sio4flx po4flx flxS2 flxDN flxDSi ...
    NPS1 RPS1 MTS1 NPS2 RPS2 MTS2 EXZ1 MTZ1 EXZ2 MTZ2 ...
    GS1Z1 GS2Z2 GZ1Z2 GDNZ2 GTZ2 MIDN MIDSi Nit ...
    sLight SPM ak2 ak3 ADPT OXR pih1 pih2 pnh4s1 pnh4s2 bfNO3S1 bfNH4S1 fPO4S1 fCO2S1 ...
    bfS1 bfNO3S2 bfNH4S2 fNO3S2 fNH4S2 fPO4S2 fCO2S2 fSiO4S2 bfS2

global time depth nLevs

Tadjust = exp(0.069*(temp_prof-TR));
%-------------------------------------------------------------------------
% Biological state variables
iNO3 = 1; iSiO4 = 2; iNH4 = 3; iPO4 = 4;  % nutrients
iS1 = 5; iS2 = 6; iZ1 = 7; iZ2 = 8;           % phytplankton/zooplankton
iDN = 9; iDSi = 10;                                  % detritus
iDOX = 11;              % carbon-oxygen

rKe = (aks(1)+aks(2)*(T(:, iS1)+T(:, iS2))+aks(3).*SPM).*dz;
rKe(1) = (aks(1)+aks(2)*(T(1, iS1)+T(1, iS2))+aks(3).*SPM(1)).*dz/2;
sLight = PAR*cumprod(exp(-rKe));

% NH4 inhibition for S1 and S2
pnh4s1 = min(1.0, exp(-pis(1)*T(:, iNH4)));
pnh4s2 = min(1.0, exp(-pis(2)*T(:, iNH4)));

% PO4,CO2,and SiO4 limiation factors
fPO4S1 = T(:, iPO4)./(kpo4s(1)+T(:, iPO4));

fPO4S2 = T(:, iPO4)./(kpo4s(2)+T(:, iPO4));
fSiO4S2 = T(:, iSiO4)./(ksio4 +T(:, iSiO4));

% Nitrogen limitation factors
rtmp = 1+T(:, iNH4)/knh4s(1)+pnh4s1.*T(:, iNO3)/kno3s(1);
bfNO3S1 = pnh4s1.*T(:, iNO3)./(kno3s(1).*rtmp);
bfNH4S1 = T(:, iNH4)./(knh4s(1).*rtmp);

rtmp =1+T(:, iNH4)./knh4s(2)+pnh4s2.*T(:, iNO3)/kno3s(2);
bfNO3S2=pnh4s2.*T(:, iNO3)./(kno3s(2).*rtmp);
bfNH4S2=T(:, iNH4)./(knh4s(2).*rtmp);

% final limitation
bfS1 = min([bfNO3S1+bfNH4S1, fPO4S1], [], 2);  %*pih1
bfS2 = min([bfNO3S2+bfNH4S2, fSiO4S2, fPO4S2], [], 2); %*pih2

%adjustment for nitrogen limitation factors
fNO3S1= bfS1.*bfNO3S1./(bfNO3S1+bfNH4S1+1.0e-6);
fNH4S1= bfS1.*bfNH4S1./(bfNO3S1+bfNH4S1+1.0e-6);

fNO3S2 = bfS2.*bfNO3S2./(bfNO3S2+bfNH4S2+1.0e-6);
fNH4S2 = bfS2.*bfNH4S2./(bfNO3S2+bfNH4S2+1.0e-6);

%Zooplankton grazing
GS1Z1= betaz(1).*T(:, iZ1).*T(:, iS1)./(kgz(1)+T(:, iS1));

mS2i = T(:, iS2);
mZ1i = T(:, iZ1);
mDNi = T(:, iDN);
mZ2i = T(:, iZ2);

rhot = rhoz(1).*mS2i+rhoz(2).*mZ1i+rhoz(3).*mDNi;
rhop = rhoz(1).*mS2i.*mS2i+rhoz(2).*mZ1i.*mZ1i+rhoz(3).*mDNi.*mDNi;

GS2Z2 = betaz(2)*rhoz(1).*mS2i.*mS2i.*mZ2i./(kgz(2)*rhot+rhop);
GZ1Z2 = betaz(2)*rhoz(2).*mZ1i.*mZ1i.*mZ2i./(kgz(2)*rhot+rhop);
GDNZ2 = betaz(2)*rhoz(3).*mDNi.*mDNi.*mZ2i./(kgz(2)*rhot+rhop);

ind = (rhot<=0 & rhop<=0) | iz2graze==0;
GS2Z2(ind) = 0;
GDNZ2(ind) = 0;
GZ1Z2(ind) = 0;

GTZ2 = GDNZ2+GZ1Z2+GS2Z2;

% oxidation rate of organic matter
OXR = T(:, iDOX)./(kox+T(:, iDOX));  
OXR(end) = 1;  % bottom

%-------------------------------------------------------------------
% Light Availibity
%-------------------------------------------------------------------
ADPT = 1;
zr = -(depth(:)+dz/2);
if idapt==1
    ADPT = alpha_corr*(1.0-4.0*zr/zeptic);
end
pih1 = (1.0-exp(-sLight.*ADPT.*alphas(1)./gmaxs(1))).*exp(-betas(1).*sLight./gmaxs(1));
pih2 = (1.0-exp(-sLight.*ADPT.*alphas(2)./gmaxs(2))).*exp(-betas(2.)*sLight./gmaxs(2));

%% S1
NPS1= gmaxs(1).*fNO3S1.*pih1.*T(:, iS1); %Growth
RPS1= gmaxs(1).*max(kns(1).*T(:, iNH4)./(knh4s(1)+T(:, iNH4)), fNH4S1.*pih1).*T(:, iS1); % Growth, nighttime uptake
MTS1 = gammas(1)*T(:, iS1); % Mortality

Tc = NPS1+RPS1-GS1Z1-MTS1;
T(:, iS1) = T(:, iS1) + dt*Tadjust.*Tc;
%% S2
NPS2 = gmaxs(2).*fNO3S2.*pih2.*T(:, iS2); % Growth
RPS2 = gmaxs(2).*max([kns(2).*T(:, iNH4)./(knh4s(2)+T(:, iNH4)), fNH4S2.*pih2], [], 2).*T(:, iS2); %Growth, nighttime uptake
MTS2 = gammas(2).*T(:, iS2); % Mortality

Tc = NPS2+RPS2-GS2Z2-MTS2;
T(:, iS2) = T(:, iS2) + dt*Tadjust.*Tc;
%% Z1
EXZ1 = OXR.*kez(1).*T(:, iZ1); %excretion
MTZ1 = gammaz(1).*T(:, iZ1).*T(:, iZ1); %Mortality

Tc = alphaz(1)*GS1Z1-EXZ1-GZ1Z2-MTZ1;
T(:, iZ1) = T(:, iZ1)  + dt*Tadjust.*Tc;
%% Z2
EXZ2 = OXR.*kez(2).*T(:, iZ2); % excretion
MTZ2 = gammaz(2).*T(:, iZ2).*T(:, iZ2); %Mortality

if ibgraze==1 %mimic bottom grazing
    MTZ2(end) = bgraze*gammaz(2)*T(end, iZ2);
end

Tc = alphaz(2)*GTZ2-EXZ2-MTZ2;
T(:, iZ2) = T(:, iZ2) + dt*Tadjust.*Tc;
%% DN
MIDN = max(kmdn(1)*temp_prof+kmdn(2), 0.05).*OXR.*T(:, iDN); %remineralization, 1.5 to increase dissolution

Tc = (1-alphaz(1))*GS1Z1+(1-alphaz(2))*GTZ2-GDNZ2 +MTS1+MTS2+MTZ1+MTZ2-MIDN;
T(:, iDN) = T(:, iDN) + dt*Tadjust.*Tc;
%% DSi
MIDSi = max(kmdsi(1)*temp_prof+kmdsi(2), 0.01).*T(:, iDSi); %remineralization, 1.5 to increase dissolution

Tc = (GS2Z2+MTS2)*si2n-MIDSi;
T(:, iDSi) = T(:, iDSi) + dt*Tadjust.*Tc;
%% NO3
Nit = gamman.*OXR.*T(:, iNH4); %Nitrification

Tc = -NPS1-NPS2+Nit;
T(:, iNO3) = T(:, iNO3) + dt*Tadjust.*Tc;
%% NH4
Tc = -RPS1-RPS2+EXZ1+EXZ2-Nit+MIDN;
T(:, iNH4) = T(:, iNH4) + dt*Tadjust.*Tc;

%% SiO4
Tc = -(NPS2+RPS2)*si2n+MIDSi;
T(:, iSiO4) = T(:, iSiO4) + dt*Tadjust.*Tc;

%% PO4
Tc = (EXZ1+EXZ2+MIDN-NPS1-RPS1-NPS2-RPS2)*p2n;
if ipo4==1
    Tc = Tc + MIDSi*p2n/si2n;
end
T(:, iPO4) = T(:, iPO4) + dt*Tadjust.*Tc;
%% DOX
Tc = (NPS1+NPS2)*o2no+(RPS1+RPS2-EXZ1-EXZ2-MIDN)*o2nh-2.0*Nit;
T(:, iDOX) = T(:, iDOX) + dt*Tc;

%% Mimic the sitotaxis of diatoms
if iws==1
    wss2 = zeros(nLevs, 1);
    ind = T(:, iNO3) < NO3c;
    wss2(ind) =  max(0.1, ws1*exp(-ws2*T(ind, iNO3)));
end
%% Sinking of phytoplankton and detritus.
Tc = [0.001; T(1:end-1, iS2)];
T(:,iS2)   = T(:,iS2) - dt*wss2.*(T(:,iS2)-Tc)./dz;

Tc = [0.001; T(1:end-1,iDSi)];
T(:,iDSi)   = T(:,iDSi) - dt*wsdsi*(T(:,iDSi)-Tc)/dz;

Tc = [0.001; T(1:end-1,iDN)];
T(:,iDN)   = T(:,iDN) - dt*wsdn*(T(:,iDN)-Tc)/dz;
%% Sediment Flux
if ised==1
    drat = dz/max(dz, 0.4);
    rat = 1;
    if T(end,iS2)<=2.5
        rat = 0;
    end
    flxS2 = rat*drat*wss2*T(end, iS2);
    flxDN = drat*wsdn*T(end, iDN);
    flxDSi = drat*wsdsi*T(end, iDSi);

    [nh4flx,sio4flx,po4flx,co2flx,o2flx] = calc_cosine_sedflux(flxS2, flxDN, flxDSi);

    T(end, iNO3) = T(end, iNO3) + dt*nh4flx/dz;
    T(end, iSiO4) = T(end, iSiO4) + dt*sio4flx/dz;
    T(end, iPO4) = T(end, iPO4) + dt*po4flx/dz;
    T(end, iDOX) = T(end, iDOX) + dt*o2flx/dz;
end

%% Check
T(T<0) = 0;

end

function srad = surface_light(latitude, time)
%
% This function calculates the hourly radiation from the formula
% Io = SC * ECC .* C_theta * PARfrac * cloud
% Where SC = solar constant
%       ECC = eccentricity correction (between 0.9 and 1.1 approx.)
%       C_theta = cosine of the solar zenith angle, as function of
%                 latitude, time of the year and hour of the day
%       PARfrac = photosynthetically active radiation fraction 0.43 of total
%                 incoming light
%       cloud = fraction of light available after cloud attenuation
%

% Constants
SC = 1366.1;  % solar constant W m^{-2}
PHI = latitude*2*pi/360; % latitude angle
cloud = 0.3;

% Eccentricity correction
%  calculated using the Spencer's fourier series.
gam = 2*pi*time./365;  % time of the year in radians
ECC = 1.000110 + (0.034221 * cos(gam)) + (0.001280 * sin (gam)) + ...
    (0.000719*cos(2*gam) + (0.000077 * sin(2*gam))); % Eccentricity: correction factor of Earth's orbit (~0.9 - ~1.1)

% Calculate declination angle
delta = 0.006918 - 0.399912 * cos(gam) + 0.070257 * sin (gam) - ...
    0.006758 * cos(2* gam) + 0.000907 * sin (2*gam) - 0.002697 * cos(3*gam) +...
    0.001480*sin(3*gam);

% Calculate hour angle
%  The hour angle is equal to zero at local noon and increases in magnitude
%  by pi/12 or 15 degrees for every hour before and after noon.
hour = 24*mod(time,1);
h = pi/12 * (12 - hour);

% Calculate cosine of the solar zenith angle (C_theta)
C_theta = sin(PHI) * sin(delta)+ cos(PHI) * cos(delta) .* cos(h);

% Calculate hourly solar radiation
srad = SC * ECC.* C_theta * (1-cloud); % solar radiation after attenuation by clouds
night = srad <0;% Set night hours to zero
srad(night)=0;

end

function [nh4flx,sio4flx,po4flx,co2flx,o2flx] = calc_cosine_sedflux(flxS2,flxDN, flxDSi)
% Calculate the sediment oxygen demand
global fS2 fDN fDSi rkS2 rkDN rkDSi mkS2 mkDN mkDSi
global si2n p2n c2n o2nh ipo4
global PS2 PDN PDSi RS2 RDN RDSi JS2 JDN JDSi
global dt

% for each G class
for m = 1:3
    % diagensis flux
    JS2(m) = RS2(m) *PS2(m);     %#ok<*AGROW> %*drat
    JDN(m) = RDN(m) *PDN(m);   %*drat
    JDSi(m) = RDSi(m)*PDSi(m);  %*drat

    % sediment POM conc.
    PS2(m) = PS2(m) +dt*(fS2(m)*flxS2-JS2(m));
    PDN(m) = PDN(m) +dt*(fDN(m)*flxDN-JDN(m));
    PDSi(m) = PDSi(m)+dt*(fDSi(m)*flxDSi-JDSi(m));

    % sediment POM decay rate
    RS2(m) = RS2(m) +dt*(rkS2(m) -RS2(m)*fS2(m)*flxS2/max(1e-6, PS2(m)));
    RDN(m) = RDN(m) +dt*(rkDN(m) -RDN(m)*fDN(m)*flxDN/max(1e-6, PDN(m)));
    RDSi(m) = RDSi(m)+dt*(rkDSi(m)-RDSi(m)*fDSi(m)*flxDSi/max(1e-6, PDSi(m)));

    % check
    RS2(m) = max([min(RS2(m), mkS2(m)),0]);
    RDN(m) = max([min(RDN(m), mkDN(m)),0]);
    RDSi(m) = max([min(RDSi(m), mkDSi(m)),0]);
end

% sediment flux
nh4flx = sum(JS2)+sum(JDN);
sio4flx = si2n*sum(JS2)+sum(JDSi);
po4flx = p2n*nh4flx;
co2flx = c2n*nh4flx;
o2flx =- o2nh*nh4flx;

if ipo4==1
    po4flx = po4flx+p2n*sum(JDSi)/si2n;
end
end

function D = read_cosine_nml(filepath) %#ok<STOUT>
% Load the cosine.nml file
RD = importdata(filepath,'%/s',inf);
nLines = numel(RD);

clear D
for ii = 1:nLines
    tmp_line = strip(RD{ii});
    if strcmp(tmp_line(1), '!') || strcmp(tmp_line(1), '&')  || strcmp(tmp_line(1), '/')
    else
        ind_comt = strfind(tmp_line, '!')-1;
        if isempty(ind_comt)
            ind_comt = 10000;
        end
        ind1 = min(strfind(tmp_line, '='));
        ind2 = min(ind_comt, numel(tmp_line), 'omitnan');
        tmp_line = [tmp_line(1:ind1), '[', tmp_line(ind1+1:ind2), ']'];
        eval(['D.', tmp_line, ';']) %#ok<EVLDOT>
    end
end

end





















