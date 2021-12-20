%% Test with experimental data

addpath('Code')

StrainOptions = 0;  % 0 - Strain Gagues; 1 - LSDIC

L=[];  % laminate structure (holds sample class)
RES=[]; % results structure (holds arrays)
tic
%% Laminate info file
LF = 'laminates_materials\IM7_90_0.m';

%% Experimental data file
DF1 = load('IJF_Experimental\CFRP_TSA_mean.mat');
DF2 = load('IJF_Experimental\StrainGauges_CFRP.mat'); 
DF3 = load('IJF_Experimental\LSDIC_CFRP_mean.mat'); 
disp("Loading Complete   "); toc

%% Experimental Data for [90/0]3s laminate
T0 =  DF1.CFRP.CP90.FR383.Tmean;
dT =  DF1.CFRP.CP90.FR383.dTmean;
depsx = DF2.SG.CP90.Exx;
depsy = DF2.SG.CP90.Eyy;
lsdic_y = DF3.LSDIC.CP90.Exx.vector.meanval; % DIC x and y are switched
lsdic_x = DF3.LSDIC.CP90.Eyy.vector.meanval;


freqs = [30.1 , 20.1, 10.1, 5.1, 3.1];
srrl_thick = 5.1e-6;

%% Cycle through TSA samples
tic
for i_s = 1:length(freqs)
    sprintf('Sample %i of %i' , [  i_s , length(freqs) ] )
    if StrainOptions
        L{i_s} = TSAsample(LF, T0(i_s), lsdic_x(i_s), lsdic_y(i_s), 0, freqs(i_s), srrl_thick );
    else
        L{i_s} = TSAsample(LF, T0(i_s), depsx(i_s), depsy(i_s), 0, freqs(i_s), srrl_thick );
    end
    L{i_s}.set_experimental(dT(i_s));
    L{i_s}.read_laminfo()
    L{i_s}.read_mpf()
    L{i_s}.calc_tau()
    % calculate predicted dToT from data
    L{i_s}.build_plies()
    L{i_s}.evaluate_dToTAQ()
    L{i_s}.evaluate_surfply()
    L{i_s}.evaluate_SRRL()
    L{i_s}.calc_globalCTE()
    L{i_s}.calc_Kxy()
    L{i_s}.evaluate_global()
    
    RES.freqs(i_s) = freqs(i_s);
    RES.dToT_exp(i_s) = L{i_s}.exp_dToT;
    RES.dToT_e20(i_s) = L{i_s}.dToT_eq20;
    RES.dToT_srrl(i_s) = L{i_s}.dToT_SRRL;
    RES.dToT_surfp(i_s) = L{i_s}.dToT_SURFP;
    RES.dToT_glob(i_s) = L{i_s}.dToT_GLOB;

    disp("Frequency Loop Iteration Complete"); 
    toc
end

% emmisivity
em = 0.95;
% Just Plot
figno = 13;
ff=figure(figno);
clf
plot(RES.freqs , RES.dToT_exp, 'r*')
hold on
plot(RES.freqs , RES.dToT_e20*em, 'ko')
plot(RES.freqs , RES.dToT_srrl*em, 'bs')
plot(RES.freqs , RES.dToT_surfp*em, 'bd')
plot(RES.freqs , RES.dToT_glob*em, 'bv')
% 
legend('Exp','eq20','SRRL','Surf Ply','Global','NumColumns',2)
xlabel('Freq (Hz)'); ylabel('\Delta T/T_0 (K/K)')