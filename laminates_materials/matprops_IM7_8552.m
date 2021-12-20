%% Material properties for IM7/8552 CFRP
%% From Fortunato et al. 2021
% pp = ply property
%% SG props
pp.E_1 = 146.86e9;  % (pm 0.02)
pp.E_2 = 9.27e9;    % (pm 0.03)
pp.G_12 = 5.01e9;   % (pm 0.01)
pp.nu_12 = 0.34;     % (pm 0.01)
%% AVG props (between SG and DIC)
% pp.E_1 = 148.8e9;  % (pm 2.27)
% pp.E_2 = 9.19e9;    % (pm 0.10)
% pp.G_12 = 5.06e9;   % (pm 0.07)
% pp.nu_12 = 0.34;     % (pm 0.01)
% %  DIC 
% pp.E_1 = 150.75e9;  % (pm 0.5)
% pp.E_2 = 9.11e9;    % (pm 0.04)
% pp.G_12 = 5.126e9;   % (pm 0.15)
% pp.nu_12 = 0.35;     % (pm 0.01)
%%
pp.nu_21 = pp.E_2/pp.E_1*pp.nu_12; 
pp.rho = 1556.75;    % (pm 72.5)
pp.alpha1 = -0.5e-6; % **Ref 23
pp.alpha2 = 25.8e-6; % **Ref 23
% pp.alpha1 = -0.1e-6; % **Ref 24
% pp.alpha2 = 31e-6;   % **Ref 24
% pp.alpha1 = -0.3e-6; % ** Avg used in Table 7
% pp.alpha2 = 28.4e-6; % ** Table 7
pp.Cp = 857;         % **?
pp.Ce = pp.Cp;  % assumption made by Emery & Barton
pp.K1 = pp.alpha1/(pp.rho*pp.Cp);
pp.K2 = pp.alpha2/(pp.rho*pp.Cp);
% pp.K1 = -0.17e-12;        % pm 0.06 (1/Pa) (experimental)
% pp.K2 = 16.0e-12;         % pm 1.0 XXX these are corrected from the paper
pp.k=0.8412;  %  plot range top and bottom +/- 15%

%% Literature thermal properites for resin (from Literature**)
% rp = resin property
rp.rhor = 1153;
rp.Cpr = 1100;
rp.alphaR = 53.5e-6; 
rp.Er = 3.8e9;
rp.nur = 0.35;
rp.Kr = rp.alphaR/(rp.rhor*rp.Cpr);
% Temporary Test for old error
%rp.Kr = rp.alphaR/(pp.rho*pp.Cp)
%k = 0; % conduction coefficient (not sure of correct value)
%sprintf('Material Properties Loaded')