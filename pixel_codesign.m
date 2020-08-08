% polymorph pixel co-design
% CRT
% 05/19/2020 - 07/22/2020

%% setup
close all
clear all

%% starting assumptions and physical constants

% constants
k = 1.38e-23;
T = 37 + 273;   % body temp

% system constraints
V_DD = 1.5;
R_elec = 24e3;
pixel_dim = 152.5e-6;
area_limit = (pixel_dim)^2;

% amplitude range
% peak voltages for AP (300Hz - 10kHz) and a little LFP (50Hz - 10kHz)
V_AP = 2e-3;
V_LFP = 4e-3;

%% design targets

% except for here, as noted, all noise is squared voltage
noise_Vrms_target = 5e-6;
noise_Vsq_target = noise_Vrms_target^2;

A_target = 100;
A_OL_target = 1e3;
A_f = 1;

f_LP_AP = 10e3;
f_HP_AP = 300;

ENBW = f_LP_AP * pi/2;

% number of time constants for "complete" settling
settling_TCs = 5;

% target for OTA noise source other than M1
excess = 1.2;

% swing for A1
V_swing_target = V_LFP * A_target;

%% technology information

L_min = 130e-9;

% C_ox and C_P (C_gg) [F/m^2]
% C_ox for flicker noise
C_ox = 0.01;
% C_P for parasitics
M_input_area = 120e-6 * 600e-9;
C_gg = 252e-15;
C_P_N = C_gg / M_input_area;
C_P_P = 0.00633;

% flicker noise coefficients [J]
K_flicker_N = 1.24e-25;
K_flicker_P = 1.98e-25;

% mismatch coefficients
% [Vm] and [m]
A_VT = 10e-9;
A_beta = 0.013e-6;

% noise factor gamma (often quoted as 2/3)
gamma_W = 1;
gamma_S = 2/3;

% estimates for device area
L_diff = 200e-9;
L_guard = 1e-6;

%% pseudoR cutoff and noise

% possible pseudoR values from mc: 90GOhm to 5TOhm
pseudoR = logspace(11.25,11.77,100);
% assume C_F is minimum size
C_F = 10e-15;
f_pseudo = 1./(2.*pi.*pseudoR.*C_F);
% reduction factor to account for not all of kT/C
% appearing in the 300Hz+ band
reduc_factor = 1./(1 + f_HP_AP./f_pseudo);
noise_pseudoR = 2 .* k .* T ./ C_F .* reduc_factor;
% input-refer - note that this uses ideal gain
noise_pseudoR_in = noise_pseudoR ./ (A_target.^2);

% noise due to pseudoR as a function of Cin
figure
semilogx(pseudoR/1e9, sqrt(noise_pseudoR_in) .* 1e6);
xlabel('R_{pseud} (G\Omega)');
ylabel('V_{n-in-RMS} (\muV)');
title('pseudoR Noise vs. pseudoR');

%% closed-loop amp design: sweep Cin

% assumed pseudoR based on 10fF feedback cap
pseudoR = 500e9;

% sweep Cin
%C_in = 1e-12;
C_in = linspace(100e-15, 3e-12, 1000)';
C_F = C_in / A_target;

% integrated pseudoR noise for 300Hz sw cap HP cutoff
f_pseudo = 1./(2.*pi.*pseudoR.*C_F);
% reduction factor to account for not all of kT/C
% appearing in the 300Hz+ band
reduc_factor = 1./(1 + f_HP_AP./f_pseudo);
noise_pseudoR = 2 .* k .* T ./ C_F .* reduc_factor;
% input-refer - note that this uses ideal gain
noise_pseudoR_in = noise_pseudoR ./ (A_target.^2);

% ideal closed-loop gain including electrode attenuation
% at three freqs: 300, 1000, 10000 Hz
% assumes ideal opamp with no parasitics
A_elec_300 = C_in ./ sqrt(C_F.^2 + (C_F .* C_in .* R_elec .* (2*pi*300)).^2);
A_elec_1e3 = C_in ./ sqrt(C_F.^2 + (C_F .* C_in .* R_elec .* (2*pi*1e3)).^2);
A_elec_1e4 = C_in ./ sqrt(C_F.^2 + (C_F .* C_in .* R_elec .* (2*pi*1e4)).^2);

%% plots from Cin choice

% noise due to pseudoR as a function of Cin
figure
plot(C_in.*(1e12), sqrt(noise_pseudoR_in) .* 1e6);
xlabel('C_{in} (pF)');
ylabel('V_{n-in-RMS} (\muV)');
title('pseudoR Noise');

% input gain including electrode atten
% should be negligible
figure
plot(C_in.*(1e12), A_elec_300);
hold on
plot(C_in.*(1e12), A_elec_1e3);
plot(C_in.*(1e12), A_elec_1e4);
xlabel('C_{in} (pF)');
ylabel('A_{with-elec}');
axis([100e-3 3 98 102]);
legend('300 Hz', '1kHz', '10kHz');
title('LNA Gain Including Electrode Attenuation');

%% results from first design section

% baseline C_in for noise, assuming we meet target gain
C_in = 1e-12;
C_F = C_in / A_target;
% now adjust C_in to pre-warp the CL gain, gain error due to finite OL gain
C_in = C_F * A_target * (1 + 1/A_OL_target) / (1 - A_target/A_OL_target);
% integrated pseudoR noise for 300Hz sw cap HP cutoff
f_pseudo = 1./(2.*pi.*pseudoR.*C_F);
reduc_factor = 1./(1 + f_HP_AP./f_pseudo);
noise_pseudoR = 2 .* k .* T ./ C_F .* reduc_factor;

%% begin OTA design: sweep WL

% transistor area, from 10^-14 m^2 to 10^-9 m^2
% for reference, pixel is 2.33*10^-8 m^2
% 120um/600nm device is 7.2*10^-11 m^2both
% and a 130nm/130nm device is 1.69*10^-14 m^2
WL = logspace(-14,-8,1000);

C_P = WL .* C_P_N;

% resulting actual closed-loop gain
A_actual = C_in ./ (C_F + (C_F + C_in + C_P) ./ A_OL_target);
noise_pseudoR_in = noise_pseudoR ./ (A_actual.^2);

% noise budget
noise_OTA = noise_Vsq_target - noise_pseudoR_in;

% noise gain
A_noise = (C_P + C_in + C_F) ./ C_in;

% flicker noise
noise_flicker = 2 .* K_flicker_N ./ C_ox ./ WL .* (A_noise).^2 ...
    .* log(f_LP_AP / f_HP_AP);

noise_thermal = noise_OTA - noise_flicker;

% thermal noise
% NOTE: we only care about thermal noise in the signal band (to 10kHz)
% technically, ENBW is not correct here because we can do a software
% brickwall filter at 10kHz, and this will be valid since we have used
% appropriate anti-aliasing.
% thus, ENBW is used here only to provide a safety margin.
gm_req = 2 .* 4 .* k .* T .* gamma_W...
    ./ noise_thermal .* (A_noise).^2 .* ENBW .* excess;

% first approx of mismatch: assume input pair dominates
sigma_Vth_1 = A_VT ./ sqrt(WL);
V_os_3sigma = sigma_Vth_1 .* 3;

%% plots from flicker noise sizing choice

% noise breakdown
figure
semilogx(WL, (noise_Vrms_target*1e6).*ones(size(WL)));
hold on
semilogx(WL, sqrt(noise_pseudoR_in) .* 1e6);
semilogx(WL, sqrt(noise_OTA) .* 1e6);
semilogx(WL, sqrt(noise_flicker) .* 1e6);
semilogx(WL, sqrt(noise_thermal) .* 1e6);
semilogx(WL, sqrt(noise_flicker) ./ (A_noise).^2 .* 1e6);
semilogx(WL, sqrt(noise_thermal) ./ (A_noise).^2 .* 1e6);
xlabel('area_{diff} (m^2)');
ylabel('V_{n-RMS} (\muV)');
axis([1e-14 1e-8 0 10]);
legend('Budget', 'pseudoR', 'OTA', 'Flicker', 'Thermal',...
    'Flicker/A_{noise}', 'Thermal/A_{noise}');
title('Contributors to LNA Noise');

% mismatch
figure
semilogx(WL, V_os_3sigma .* 1e3);
xlabel('area_{diff} (m^2)');
ylabel('V_{OS} (mV)');
title('3\sigma V_{OS}');

% gm requirement
figure
semilogx(WL, gm_req.*(1e6));
xlabel('area_{diff} (m^2)');
ylabel('g_m (\muS)');
axis([1e-14 1e-8 0 500]);
title('g_m Requirement');

% closed-loop gain
figure
semilogx(WL, A_actual);
xlabel('area_{diff} (m^2)');
ylabel('A_{actual}');
title('Closed-Loop Gain due to C_P');

%% gm/id, gmro characterization

% info about number of characterization sweeps
num_lengths_N = 10;
num_lengths_P = 10;
points_per_length = 25;
L_N = [130e-9 200e-9 500e-9 1e-6 2e-6 4e-6 6e-6 8e-6 10e-6 12e-6]';
L_P = L_N;

% NMOS

% import the characterization data
gmid_n_130n = readmatrix('charac/charac_gmid_n_130n.csv');
gmid_n_200n = readmatrix('charac/charac_gmid_n_200n.csv');
gmid_n_500n = readmatrix('charac/charac_gmid_n_500n.csv');
gmid_n_1u = readmatrix('charac/charac_gmid_n_1u.csv');
gmid_n_2u = readmatrix('charac/charac_gmid_n_2u.csv');
gmid_n_4u = readmatrix('charac/charac_gmid_n_4u.csv');
gmid_n_6u = readmatrix('charac/charac_gmid_n_6u.csv');
gmid_n_8u = readmatrix('charac/charac_gmid_n_8u.csv');
gmid_n_10u = readmatrix('charac/charac_gmid_n_10u.csv');
gmid_n_12u = readmatrix('charac/charac_gmid_n_12u.csv');

gmro_n_130n = readmatrix('charac/charac_gmro_n_130n.csv');
gmro_n_200n = readmatrix('charac/charac_gmro_n_200n.csv');
gmro_n_500n = readmatrix('charac/charac_gmro_n_500n.csv');
gmro_n_1u = readmatrix('charac/charac_gmro_n_1u.csv');
gmro_n_2u = readmatrix('charac/charac_gmro_n_2u.csv');
gmro_n_4u = readmatrix('charac/charac_gmro_n_4u.csv');
gmro_n_6u = readmatrix('charac/charac_gmro_n_6u.csv');
gmro_n_8u = readmatrix('charac/charac_gmro_n_8u.csv');
gmro_n_10u = readmatrix('charac/charac_gmro_n_10u.csv');
gmro_n_12u = readmatrix('charac/charac_gmro_n_12u.csv');

cgg_n_130n = readmatrix('charac/charac_cgg_n_130n.csv');
cgg_n_200n = readmatrix('charac/charac_cgg_n_200n.csv');
cgg_n_500n = readmatrix('charac/charac_cgg_n_500n.csv');
cgg_n_1u = readmatrix('charac/charac_cgg_n_1u.csv');
cgg_n_2u = readmatrix('charac/charac_cgg_n_2u.csv');
cgg_n_4u = readmatrix('charac/charac_cgg_n_4u.csv');
cgg_n_6u = readmatrix('charac/charac_cgg_n_6u.csv');
cgg_n_8u = readmatrix('charac/charac_cgg_n_8u.csv');
cgg_n_10u = readmatrix('charac/charac_cgg_n_10u.csv');
cgg_n_12u = readmatrix('charac/charac_cgg_n_12u.csv');

% concatenate the results for easier looping
gmid_n = horzcat(gmid_n_130n, gmid_n_200n(:,2), gmid_n_500n(:,2),...
    gmid_n_1u(:,2), gmid_n_2u(:,2), gmid_n_4u(:,2), gmid_n_6u(:,2),...
    gmid_n_8u(:,2), gmid_n_10u(:,2), gmid_n_12u(:,2));
gmro_n = horzcat(gmro_n_130n, gmro_n_200n(:,2), gmro_n_500n(:,2),...
    gmro_n_1u(:,2), gmro_n_2u(:,2), gmro_n_4u(:,2), gmro_n_6u(:,2),...
    gmro_n_8u(:,2), gmro_n_10u(:,2), gmro_n_12u(:,2));
cgg_n = horzcat(cgg_n_130n, cgg_n_200n(:,2), cgg_n_500n(:,2),...
    cgg_n_1u(:,2), cgg_n_2u(:,2), cgg_n_4u(:,2), cgg_n_6u(:,2),...
    cgg_n_8u(:,2), cgg_n_10u(:,2), cgg_n_12u(:,2));

% plot the characterization results
figure
semilogx(gmid_n(:,1)*1e6, gmid_n(:,2));
hold on
semilogx(gmid_n(:,1)*1e6, gmid_n(:,3));
semilogx(gmid_n(:,1)*1e6, gmid_n(:,4));
semilogx(gmid_n(:,1)*1e6, gmid_n(:,5));
semilogx(gmid_n(:,1)*1e6, gmid_n(:,6));
semilogx(gmid_n(:,1)*1e6, gmid_n(:,7));
semilogx(gmid_n(:,1)*1e6, gmid_n(:,8));
semilogx(gmid_n(:,1)*1e6, gmid_n(:,9));
semilogx(gmid_n(:,1)*1e6, gmid_n(:,10));
semilogx(gmid_n(:,1)*1e6, gmid_n(:,11));
title('g_m/I_D vs. I_D for NMOS');
xlabel('I_D (\muA)');
ylabel('g_m/I_D (V^{-1})');
legend('L = 130nm', 'L = 200nm', 'L = 500nm', 'L = 1\mum', 'L = 2\mum',...
    'L = 4\mum', 'L = 6\mum', 'L = 8\mum', 'L = 10\mum', 'L = 12\mum');

figure
semilogx(gmro_n(:,1)*1e6, gmro_n(:,2));
hold on
semilogx(gmro_n(:,1)*1e6, gmro_n(:,3));
semilogx(gmro_n(:,1)*1e6, gmro_n(:,4));
semilogx(gmro_n(:,1)*1e6, gmro_n(:,5));
semilogx(gmro_n(:,1)*1e6, gmro_n(:,6));
semilogx(gmro_n(:,1)*1e6, gmro_n(:,7));
semilogx(gmro_n(:,1)*1e6, gmro_n(:,8));
semilogx(gmro_n(:,1)*1e6, gmro_n(:,9));
semilogx(gmro_n(:,1)*1e6, gmro_n(:,10));
semilogx(gmro_n(:,1)*1e6, gmro_n(:,11));
title('g_mr_o vs. I_D for NMOS');
xlabel('I_D (\muA)');
ylabel('g_mr_o');
legend('L = 130nm', 'L = 200nm', 'L = 500nm', 'L = 1\mum', 'L = 2\mum',...
    'L = 4\mum', 'L = 6\mum', 'L = 8\mum', 'L = 10\mum', 'L = 12\mum');

figure
semilogx(cgg_n(:,1)*1e6, cgg_n(:,2) / (130e-9 * 260e-9));
hold on
semilogx(cgg_n(:,1)*1e6, cgg_n(:,3) / (200e-9 * 400e-9));
semilogx(cgg_n(:,1)*1e6, cgg_n(:,4) / (500e-9 * 1000e-9));
semilogx(cgg_n(:,1)*1e6, cgg_n(:,5) / (1000e-9 * 2000e-9));
semilogx(cgg_n(:,1)*1e6, cgg_n(:,6) / (2000e-9 * 4000e-9));
semilogx(cgg_n(:,1)*1e6, cgg_n(:,7) / (4000e-9 * 8000e-9));
semilogx(cgg_n(:,1)*1e6, cgg_n(:,8) / (6000e-9 * 12000e-9));
semilogx(cgg_n(:,1)*1e6, cgg_n(:,9) / (8000e-9 * 16000e-9));
semilogx(cgg_n(:,1)*1e6, cgg_n(:,10) / (10000e-9 * 20000e-9));
semilogx(cgg_n(:,1)*1e6, cgg_n(:,11) / (12000e-9 * 24000e-9));
title('C_{gg} vs. I_D for NMOS');
xlabel('I_D (\muA)');
ylabel('C_{gg} (F/m^2)');
legend('L = 130nm', 'L = 200nm', 'L = 500nm', 'L = 1\mum', 'L = 2\mum',...
    'L = 4\mum', 'L = 6\mum', 'L = 8\mum', 'L = 10\mum', 'L = 12\mum');


% PMOS

% import the characterization data
gmid_p_130n = readmatrix('charac/charac_gmid_p_130n.csv');
gmid_p_200n = readmatrix('charac/charac_gmid_p_200n.csv');
gmid_p_500n = readmatrix('charac/charac_gmid_p_500n.csv');
gmid_p_1u = readmatrix('charac/charac_gmid_p_1u.csv');
gmid_p_2u = readmatrix('charac/charac_gmid_p_2u.csv');
gmid_p_4u = readmatrix('charac/charac_gmid_p_4u.csv');
gmid_p_6u = readmatrix('charac/charac_gmid_p_6u.csv');
gmid_p_8u = readmatrix('charac/charac_gmid_p_8u.csv');
gmid_p_10u = readmatrix('charac/charac_gmid_p_10u.csv');
gmid_p_12u = readmatrix('charac/charac_gmid_p_12u.csv');

gmro_p_130n = readmatrix('charac/charac_gmro_p_130n.csv');
gmro_p_200n = readmatrix('charac/charac_gmro_p_200n.csv');
gmro_p_500n = readmatrix('charac/charac_gmro_p_500n.csv');
gmro_p_1u = readmatrix('charac/charac_gmro_p_1u.csv');
gmro_p_2u = readmatrix('charac/charac_gmro_p_2u.csv');
gmro_p_4u = readmatrix('charac/charac_gmro_p_4u.csv');
gmro_p_6u = readmatrix('charac/charac_gmro_p_6u.csv');
gmro_p_8u = readmatrix('charac/charac_gmro_p_8u.csv');
gmro_p_10u = readmatrix('charac/charac_gmro_p_10u.csv');
gmro_p_12u = readmatrix('charac/charac_gmro_p_12u.csv');

% concatenate the results for easier looping
gmid_p = horzcat(gmid_p_130n, gmid_p_200n(:,2), gmid_p_500n(:,2),...
    gmid_p_1u(:,2), gmid_p_2u(:,2), gmid_p_4u(:,2), gmid_p_6u(:,2),...
    gmid_p_8u(:,2), gmid_p_10u(:,2), gmid_p_12u(:,2));
gmro_p = horzcat(gmro_p_130n, gmro_p_200n(:,2), gmro_p_500n(:,2),...
    gmro_p_1u(:,2), gmro_p_2u(:,2), gmro_p_4u(:,2), gmro_p_6u(:,2),...
    gmro_p_8u(:,2), gmro_p_10u(:,2), gmro_p_12u(:,2));

% plot the characterization results
figure
semilogx(gmid_p(:,1)*1e6, gmid_p(:,2));
hold on
semilogx(gmid_p(:,1)*1e6, gmid_p(:,3));
semilogx(gmid_p(:,1)*1e6, gmid_p(:,4));
semilogx(gmid_p(:,1)*1e6, gmid_p(:,5));
semilogx(gmid_p(:,1)*1e6, gmid_p(:,6));
semilogx(gmid_p(:,1)*1e6, gmid_p(:,7));
semilogx(gmid_p(:,1)*1e6, gmid_p(:,8));
semilogx(gmid_p(:,1)*1e6, gmid_p(:,9));
semilogx(gmid_p(:,1)*1e6, gmid_p(:,10));
semilogx(gmid_p(:,1)*1e6, gmid_p(:,11));
title('g_m/I_D vs. I_D for PMOS');
xlabel('I_D (\muA)');
ylabel('g_m/I_D (V^{-1})');
legend('L = 130nm', 'L = 200nm', 'L = 500nm', 'L = 1\mum', 'L = 2\mum',...
    'L = 4\mum', 'L = 6\mum', 'L = 8\mum', 'L = 10\mum', 'L = 12\mum');

figure
semilogx(gmro_p(:,1)*1e6, gmro_p(:,2));
hold on
semilogx(gmro_p(:,1)*1e6, gmro_p(:,3));
semilogx(gmro_p(:,1)*1e6, gmro_p(:,4));
semilogx(gmro_p(:,1)*1e6, gmro_p(:,5));
semilogx(gmro_p(:,1)*1e6, gmro_p(:,6));
semilogx(gmro_p(:,1)*1e6, gmro_p(:,7));
semilogx(gmro_p(:,1)*1e6, gmro_p(:,8));
semilogx(gmro_p(:,1)*1e6, gmro_p(:,9));
semilogx(gmro_p(:,1)*1e6, gmro_p(:,10));
semilogx(gmro_p(:,1)*1e6, gmro_p(:,11));
title('g_mr_o vs. I_D for PMOS');
xlabel('I_D (\muA)');
ylabel('g_mr_o');
legend('L = 130nm', 'L = 200nm', 'L = 500nm', 'L = 1\mum', 'L = 2\mum',...
    'L = 4\mum', 'L = 6\mum', 'L = 8\mum', 'L = 10\mum', 'L = 12\mum');

%% transistor-level OTA design: tunable parameters

% set minimum W
W_min = 200e-9;

% assumptions for WI or SI
gmid_WI = 25;
gmid_MI = 20;
gmid_SI = 10;

% M1 specs determined from optimization above
gm_min = 100e-6;
WL_min = 3.85e-11;

% ratio between main branch and output branch
current_scale = 16;
% ratio between load branch and one branch of CMFB amp
cmfb_scale = 1;

% partition of gains
A_OL_1_target = 320;
A_OL_2_target = A_OL_target / A_OL_1_target;
A_OL_CMFB_target = 20;

% degen resistor headroom consumption
V_R_degen = 0.2;
% estimated headroom for saturation for a single device (even SI)
V_hr_SI = 300e-3;
V_hr_WI = 250e-3;

%% swing requirement calculation

V_swing_telecas = V_DD - 1*V_hr_SI - 4*V_hr_WI - V_R_degen;
V_swing_foldcas = V_DD - 1*V_hr_SI - 3*V_hr_WI - V_R_degen;
V_swing_stage2 = V_DD - V_hr_SI - V_hr_WI;

%% amp/filter interface

% assume min size cap for C_RA
C_RA = 10e-15;

% biquad filter params
omega_0 = 2*pi * sqrt(f_HP_AP * f_LP_AP);
BW = 2*pi * (f_LP_AP - f_HP_AP);
Q = omega_0 / BW;

% input cap size
C_R1 = C_RA * A_f / Q;

% sweep possible switching frequencies
f_swcap = logspace(5,7,1000)';

f_BW_amp = f_swcap ./ 10;

% resistive approximation
R_1 = 1 ./ (C_R1 * f_swcap);

% one-stage version. settling will be ok, so look at gain degredation
% using 60dB gain for one-stage (80dB not possible)
rout_assumed = (1e3) / gm_min;
Rout_amp_required = (rout_assumed .* R_1) ./ (R_1 - rout_assumed);
% using reasonable assumptions for degen R and load device,
% calculate what gmro would be necessary for up and down cascode devices
R_degen_assumed = 50e3;
ro_load_assumed = 4e7;
Rout_branch = 2 * Rout_amp_required;
gmro_required_up = sqrt(Rout_branch ./ R_degen_assumed);
gmro_required_down = Rout_branch ./ ro_load_assumed;

figure
semilogx(f_swcap, gmro_required_up);
hold on
semilogx(f_swcap, gmro_required_down);
xlabel('f_{swcap} (Hz)');
ylabel('g_mr_o');
title('Single-Stage: Cascode g_mr_o Required vs. Sw Cap Frequency');
legend('Up Branch (fold)', 'Down Branch (load)');

% settling limit
T_settle_limit = 1 ./ (f_swcap * 2);

% assume switch resistance - settling through this
R_sw = 14e3;
T_settle_sw = R_sw * C_R1 * settling_TCs;

%% transistor-level OTA design: action

% M1: input pair (NMOS)

% vectors to keep track of results for comparison
W = zeros(num_lengths_N, 1);
WL = zeros(num_lengths_N, 1);
gmro = zeros(num_lengths_N, 1);
gmid = zeros(num_lengths_N, 1);

% try each length
for i = 1:num_lengths_N
    
    % column for gmid, gmro data
    col = 1 + i;
    
    % determine corresponding length
    L_charac = L_N(i);
    W_charac = L_charac * 2;
    
    for j = 1:points_per_length
        
        % go from SI to WI to get minimum required gmid
        point = points_per_length - j + 1;
        
        gmid(i) = gmid_n(point, col);
        
        if (gmid(i) > gmid_WI)
            % found required gmid
            
            % determine corresponding current density for this length
            ID_charac = gmid_n(point, 1);
            ID_density = ID_charac / W_charac;
            
            % determine actual ID required
            ID_req = gm_min / gmid(i);
            ID = 1e-6 * ceil(ID_req * 1e6);
            
            % determine actual W required
            W(i) = ID / ID_density;
            WL(i) = W(i) .* L_charac;
            
            % keep track of corresponding gmro
            gmro(i) = gmro_n(point, col);
            break
        end      
    end
end

% plot results for M1
figure
plot(L_N*1e6, W*1e6);
xlabel('L (\mum)');
ylabel('W (\mum)');
title('M1: W vs. L');
figure
plot(L_N*1e6, WL*1e12);
xlabel('L (\mum)');
ylabel('W ((\mum)^2)');
title('M1: WL vs. L');
figure
plot(L_N*1e6, gmro);
xlabel('L (\mum)');
ylabel('g_mr_o');
title('M1: g_mr_o vs. L');

% pick L that gives us required WL for parasitics/flicker
for i = 1:num_lengths_N
    WL_this = WL(i);
    
    if (WL_this > WL_min)
        % determine final dimensions
        W_M1 = W(i);
        L_M1 = L_N(i);
        % determine final small signal params
        gmid_M1 = gmid(i);
        gmro_M1 = gmro(i);
        break
    end
end

% compute some currents
ID_main = ID;
ID_load = ID_main / current_scale;
ID_both = ID_main + ID_load;

% compute M1 small signal params
gm_M1 = gmid_M1 * ID_main;
ro_M1 = gmro_M1 / gm_M1;

% size the degen resistors for reasonable voltage swing
R_degen = V_R_degen / ID_main;

% beginning to compute the gain

% from equation for gain of first stage...
gmro_factor = gm_M1 * ro_M1^2 * R_degen / (ro_M1 + R_degen)^2;
% ... we can determine a target for gmro of M3 and M5
gmro_35_target = A_OL_1_target / gmro_factor;

% size 3,5 based on this target

% M3: telecascode device (NMOS)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_MI, 1, ID_main,...
    num_lengths_N, points_per_length);
[L_M3, W_M3, gmid_M3, gmro_M3] = pick_L(L_N, W, gmid, gmro,...
    gmro_35_target, num_lengths_N, W_min);
% compute M3 small signal params
gm_M3 = gmid_M3 * ID_main;
ro_M3 = gmro_M3 / gm_M3;

% M5: current source device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_MI, 1, ID_main,...
    num_lengths_P, points_per_length);
[L_M5, W_M5, gmid_M5, gmro_M5] = pick_L(L_P, W, gmid, gmro,...
    gmro_35_target, num_lengths_P, W_min);
% compute M5 small signal params
gm_M5 = gmid_M5 * ID_main;
ro_M5 = gmro_M5 / gm_M5;

% calculate actual first stage gain

% transconductance
Gs3 = (gmro_M3 + 1) / (ro_M3 + gmro_M5 * R_degen);
Gm_reduc = (Gs3 * ro_M1) / (1 + Gs3 * ro_M1);
Gm_1 = gm_M1 * Gm_reduc;

% output resistance
Rout_1 = 1/(1/(gmro_M5 * R_degen) + 1/(gmro_M3 * ro_M1));

A_OL_1_actual = Gm_1 * Rout_1;

% compute gain target for second stage
A_OL_2_target_updated = max(A_OL_target / A_OL_1_actual, A_OL_2_target);

% compute gmro target based on necessity that gm7ro9 > A2 for later math
gm_M7_assumed = gmid_WI * ID_load;
gm_M9_assumed = gmid_SI * ID_load;
ro_9_target = A_OL_2_target_updated / gm_M7_assumed;
gmro_9_target = gm_M9_assumed * ro_9_target;

% M9: load device (NMOS)
% the only strong-inversion transistor (bc noise)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_SI, 0, ID_load,...
    num_lengths_N, points_per_length);
[L_M9, W_M9, gmid_M9, gmro_M9] = pick_L(L_N, W, gmid, gmro,...
    gmro_9_target, num_lengths_N, W_min);
% compute M11 small signal params
gm_M9 = gmid_M9 * ID_load;
ro_M9 = gmro_M9 / gm_M9;

% compute target gmro7 based on target gain
ro_M7_target = (A_OL_2_target_updated * ro_M9) /...
    (gm_M7_assumed * ro_M9 - A_OL_2_target_updated);
gmro_7_target = gm_M7_assumed * ro_M7_target;

% M7: second stage input device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_WI, 1, ID_load,...
    num_lengths_P, points_per_length);
[L_M7, W_M7, gmid_M7, gmro_M7] = pick_L(L_P, W, gmid, gmro,...
    gmro_7_target, num_lengths_P, W_min);
% compute M7 small signal params
gm_M7 = gmid_M7 * ID_load;
ro_M7 = gmro_M7 / gm_M7;

% calculate actual second stage gain
Rout_2 = (ro_M7 * ro_M9) / (ro_M7 + ro_M9);
A_OL_2_actual = gm_M7 * Rout_2;

A_OL_actual = A_OL_1_actual * A_OL_2_actual;

%% post-design OTA characterization

% noise - note this is just open-loop OTA noise for sanity check

% compute second stage noise
noise_final_A2 = 8*k*T*(gamma_W/gm_M7 + gamma_S*gm_M9/(gm_M7^2))*ENBW...
    + 2 * K_flicker_P / (W_M7 * L_M7 * C_ox) * log(f_LP_AP / f_HP_AP);

% compute first stage noise
noise_final_A1 = 8*k*T*(gamma_W/gm_M1 + 1/(R_degen * gm_M1^2))*ENBW...
    + 2 * K_flicker_N / (W_M1 * L_M1 * C_ox) * log(f_LP_AP / f_HP_AP);

% compute total input-referred OTA noise
noise_final_OTA = noise_final_A1 + (noise_final_A2 / (A_OL_1_actual^2));
noise_final_OTA_V = sqrt(noise_final_OTA);

% offset: it can be shown that all contributions are neglig except input
sigma_Vth_1_final = A_VT / sqrt(W_M1 * L_M1);
V_os_3sigma_final = sigma_Vth_1_final * 3;

% settling - need to satisfy swcap settling assumption
T_settle_OTA = Rout_2 * C_R1 * settling_TCs;

T_settle_OTA_vec = T_settle_OTA .* ones(size(f_swcap));
T_settle_sw_vec = T_settle_sw .* ones(size(f_swcap));

figure
semilogx(f_swcap, T_settle_limit*1e6);
hold on
semilogx(f_swcap, T_settle_sw_vec*1e6);
semilogx(f_swcap, T_settle_OTA_vec*1e6);
xlabel('f_{swcap} (Hz)');
ylabel('\tau (\mus)');
title('Amp/Filter Settling Time vs. Sw Cap Frequency');
legend('Limit', 'Switch Settling', 'OTA 2nd Stage Settling');

%% CMFB
ID_cmfb = ID_load / cmfb_scale;

% compute gmro target based on necessity that gm13ro15 > A for later math
gm_M13_assumed = gmid_WI * ID_cmfb;
gm_M15_assumed = gmid_SI * ID_cmfb;
ro_15_target = A_OL_CMFB_target / gm_M13_assumed;
gmro_15_target = gm_M15_assumed * ro_15_target;

% M15: load device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_SI, 0, ID_cmfb,...
    num_lengths_P, points_per_length);
[L_M15, W_M15, gmid_M15, gmro_M15] = pick_L(L_P, W, gmid, gmro,...
    gmro_15_target, num_lengths_P, W_min);
% compute M15 small signal params
gm_M15 = gmid_M15 * ID_cmfb;
ro_M15 = gmro_M15 / gm_M15;

% compute target gmro13 based on target gain
ro_M13_target = (A_OL_CMFB_target * ro_M15) /...
    (gm_M13_assumed * ro_M15 - A_OL_CMFB_target);
gmro_13_target = gm_M13_assumed * ro_M13_target;

% M13: input device (NMOS)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_WI, 1, ID_cmfb,...
    num_lengths_N, points_per_length);
[L_M13, W_M13, gmid_M13, gmro_M13] = pick_L(L_N, W, gmid, gmro,...
    gmro_13_target, num_lengths_N, W_min);
% compute M13 small signal params
gm_M13 = gmid_M13 * ID_cmfb;
ro_M13 = gmro_M13 / gm_M13;

% calculate actual CMFB amp gain
A_OL_CMFB_actual = (gmro_M13 * ro_M15) / (ro_M13 + ro_M15);

%% BPF A2 design: tunable parameters

% gain target
A_2_required = 300;
% A3_M1 specs determined from f_swcap optimization in other script
gm_min_A2 = 2e-6;

WL_min_A2 = 9e-12;

% partition gains
A_2_1_required = 30;

%% BPF A2 design: action

% A2_M1: input pair (NMOS)

% vectors to keep track of results for comparison
W = zeros(num_lengths_N, 1);
WL = zeros(num_lengths_N, 1);
gmro = zeros(num_lengths_N, 1);
gmid = zeros(num_lengths_N, 1);

% try each length
for i = 1:num_lengths_N
    
    % column for gmid, gmro data
    col = 1 + i;
    
    % determine corresponding length
    L_charac = L_N(i);
    W_charac = L_charac * 2;
    
    for j = 1:points_per_length
        
        % go from SI to WI to get minimum required gmid
        point = points_per_length - j + 1;
        
        gmid(i) = gmid_n(point, col);
        
        if (gmid(i) > gmid_WI)
            % found required gmid
            
            % determine corresponding current density for this length
            ID_charac = gmid_n(point, 1);
            ID_density = ID_charac / W_charac;
            
            % determine actual ID required
            ID_req = gm_min_A2 / gmid(i);
            ID = 0.01e-6 * ceil(ID_req * 100e6);
            
            % determine actual W required
            W(i) = ID / ID_density;
            WL(i) = W(i) .* L_charac;
            
            % keep track of corresponding gmro
            gmro(i) = gmro_n(point, col);
            break
        end      
    end
end

% pick L that gives us required WL for parasitics/flicker
for i = 1:num_lengths_N
    WL_this = WL(i);
    
    if (WL_this > WL_min_A2)
        % determine final dimensions
        W_A2_M1 = W(i);
        L_A2_M1 = L_N(i);
        % determine final small signal params
        gmid_A2_M1 = gmid(i);
        gmro_A2_M1 = gmro(i);
        break
    end
end

% corresponding current
ID_A2 = ID;

% compute M1 small signal params
gm_A2_M1 = gmid_A2_M1 * ID_A2;
ro_A2_M1 = gmro_A2_M1 / gm_A2_M1;

% need to calculate load device ro
ro_A2_M3_required = A_2_1_required * ro_A2_M1 /...
    (gmro_A2_M1 - A_2_1_required);

gm_A2_M3_estimate = gmid_SI * ID_A2;

gmro_A2_M3_target = ro_A2_M3_required * gm_A2_M3_estimate;

% A2_M3: load device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_SI, 0, ID_A2,...
    num_lengths_P, points_per_length);
[L_A2_M3, W_A2_M3, gmid_A2_M3, gmro_A2_M3] = pick_L(L_P, W, gmid, gmro,...
    gmro_A2_M3_target, num_lengths_P, W_min);
% compute M3 small signal params
gm_A2_M3 = gmid_A2_M3 * ID_A2;
ro_A2_M3 = gmro_A2_M3 / gm_A2_M3;

% calculate actual first-stage gain
A_2_1 = gmro_A2_M1 * ro_A2_M3 / (ro_A2_M1 + ro_A2_M3);
% calculate required second-stage gain
A_2_2_required = A_2_required / A_2_1;

% A2_M5: second stage input pair (PMOS)

% vectors to keep track of results for comparison
W = zeros(num_lengths_P, 1);
WL = zeros(num_lengths_P, 1);
gmro = zeros(num_lengths_P, 1);
gmid = zeros(num_lengths_P, 1);

% try each length
for i = 1:num_lengths_P
    
    % column for gmid, gmro data
    col = 1 + i;
    
    % determine corresponding length
    L_charac = L_P(i);
    W_charac = L_charac * 2;
    
    for j = 1:points_per_length
        
        % go from SI to WI to get minimum required gmid
        point = points_per_length - j + 1;
        
        gmid(i) = gmid_p(point, col);
        
        if (gmid(i) > gmid_WI)
            % found required gmid
            
            % determine corresponding current density for this length
            ID_charac = gmid_p(point, 1);
            ID_density = ID_charac / W_charac;
            
            % determine actual ID required
            ID_req = gm_min_A2 / gmid(i);
            ID = 0.01e-6 * ceil(ID_req * 100e6);
            
            % determine actual W required
            W(i) = ID / ID_density;
            WL(i) = W(i) .* L_charac;
            
            % keep track of corresponding gmro
            gmro(i) = gmro_p(point, col);
            break
        end      
    end
end

% pick L that gives us required WL for parasitics/flicker
for i = 1:num_lengths_P
    WL_this = WL(i);
    
    if (WL_this > WL_min_A2)
        % determine final dimensions
        W_A2_M5 = W(i);
        L_A2_M5 = L_P(i);
        % determine final small signal params
        gmid_A2_M5 = gmid(i);
        gmro_A2_M5 = gmro(i);
        break
    end
end

% corresponding current
ID_A2_2 = ID;

% compute M5 small signal params
gm_A2_M5 = gmid_A2_M5 * ID_A2_2;
ro_A2_M5 = gmro_A2_M5 / gm_A2_M5;

% need to calculate load device ro
ro_A2_M7_required = A_2_2_required * ro_A2_M5 /...
    (gmro_A2_M5 - A_2_2_required);

gm_A2_M7_estimate = gmid_SI * ID_A2_2;

gmro_A2_M7_target = ro_A2_M7_required * gm_A2_M7_estimate;

% A2_M7: load device (NMOS)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_WI, 1, ID_A2_2,...
    num_lengths_N, points_per_length);
[L_A2_M7, W_A2_M7, gmid_A2_M7, gmro_A2_M7] = pick_L(L_N, W, gmid, gmro,...
    gmro_A2_M7_target, num_lengths_N, W_min);
% compute M7 small signal params
gm_A2_M7 = gmid_A2_M7 * ID_A2_2;
ro_A2_M7 = gmro_A2_M7 / gm_A2_M7;

% calculate actual second-stage gain
A_2_2 = gmro_A2_M5 * ro_A2_M7 / (ro_A2_M5 + ro_A2_M7);
% calculate total gain
A_2_actual = A_2_1 * A_2_2;

% CMFB sizing?

%% BPF A3 design: tunable parameters

% gain target
A_3_required = 350;
% A3_M1 specs determined from f_swcap optimization in other script
gm_min_A3 = 2e-6;

WL_min_A3 = 9e-12;

gmro_cascode_target = 30;

%% BPF A3 design: action

% A3_M1: input pair (NMOS)

% vectors to keep track of results for comparison
W = zeros(num_lengths_N, 1);
WL = zeros(num_lengths_N, 1);
gmro = zeros(num_lengths_N, 1);
gmid = zeros(num_lengths_N, 1);

% try each length
for i = 1:num_lengths_N
    
    % column for gmid, gmro data
    col = 1 + i;
    
    % determine corresponding length
    L_charac = L_N(i);
    W_charac = L_charac * 2;
    
    for j = 1:points_per_length
        
        % go from SI to WI to get minimum required gmid
        point = points_per_length - j + 1;
        
        gmid(i) = gmid_n(point, col);
        
        if (gmid(i) > gmid_WI)
            % found required gmid
            
            % determine corresponding current density for this length
            ID_charac = gmid_n(point, 1);
            ID_density = ID_charac / W_charac;
            
            % determine actual ID required
            ID_req = gm_min_A3 / gmid(i);
            ID = 0.01e-6 * ceil(ID_req * 100e6);
            
            % determine actual W required
            W(i) = ID / ID_density;
            WL(i) = W(i) .* L_charac;
            
            % keep track of corresponding gmro
            gmro(i) = gmro_n(point, col);
            break
        end      
    end
end

% pick L that gives us required WL for parasitics/flicker
for i = 1:num_lengths_N
    WL_this = WL(i);
    
    if (WL_this > WL_min_A3)
        % determine final dimensions
        W_A3_M1 = W(i);
        L_A3_M1 = L_N(i);
        % determine final small signal params
        gmid_A3_M1 = gmid(i);
        gmro_A3_M1 = gmro(i);
        break
    end
end

% corresponding current
ID_A3 = ID;

% compute M1 small signal params
gm_A3_M1 = gmid_A3_M1 * ID_A3;
ro_A3_M1 = gmro_A3_M1 / gm_A3_M1;

% need to calculate required Rout and cascode gmro
Rout_required_A3 = A_3_required / gm_A3_M1;
Rdown_required_A3 = Rout_required_A3 * 2;

gmro_A3_M3_target = max(Rdown_required_A3 / ro_A3_M1, gmro_cascode_target);

% A3_M3: telecascode device (NMOS)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_WI, 1, ID_A3,...
    num_lengths_N, points_per_length);
[L_A3_M3, W_A3_M3, gmid_A3_M3, gmro_A3_M3] = pick_L(L_N, W, gmid, gmro,...
    gmro_A3_M3_target, num_lengths_N, W_min);
% compute M3 small signal params
gm_A3_M3 = gmid_A3_M3 * ID_A3;
ro_A3_M3 = gmro_A3_M3 / gm_A3_M3;

% calculate actual Rup required
Rdown_actual_A3 = gmro_A3_M3 * ro_A3_M1;
Rup_required_A3 = 1 / (1/Rout_required_A3 - 1/Rdown_actual_A3);

% assume Rup_required_A3 = ro_M7

% A3_M5: load cascode device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_MI, 1, ID_A3,...
    num_lengths_P, points_per_length);
[L_A3_M5, W_A3_M5, gmid_A3_M5, gmro_A3_M5] = pick_L(L_P, W, gmid, gmro,...
    gmro_cascode_target, num_lengths_P, W_min);
% compute M5 small signal params
gm_A3_M5 = gmid_A3_M5 * ID_A3;
ro_A3_M5 = gmro_A3_M5 / gm_A3_M5;

% size M7 based on estimated gm and ro required
gm_A3_M7_est = gmid_SI * ID_A3;
ro_A3_M7_required = Rup_required_A3 / gmro_A3_M5;
gmro_A3_M7_target = gm_A3_M7_est * ro_A3_M7_required;

% A3_M7: load device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_SI, 0, ID_A3,...
    num_lengths_P, points_per_length);
[L_A3_M7, W_A3_M7, gmid_A3_M7, gmro_A3_M7] = pick_L(L_P, W, gmid, gmro,...
    gmro_A3_M7_target, num_lengths_P, W_min);
% compute M7 small signal params
gm_A3_M7 = gmid_A3_M7 * ID_A3;
ro_A3_M7 = gmro_A3_M7 / gm_A3_M7;

Rup_actual_A3 = gmro_A3_M5 * ro_A3_M7;
Rout_actual_A3 = 1 / (1/Rup_actual_A3 + 1/Rdown_actual_A3);

A_3_actual = Rout_actual_A3 * gm_A3_M1;

% CMFB sizing?

%% SF1 design: tunable parameters

% gain target
A_SF1_required = 0.95;
% SF1_M1 specs determined from f_swcap optimization in other script
gm_min_SF1 = 1e-6;

WL_min_SF1 = 9e-12;

%% SF1 design: action

% SF1_M1: input device (PMOS)

% vectors to keep track of results for comparison
W = zeros(num_lengths_N, 1);
WL = zeros(num_lengths_N, 1);
gmro = zeros(num_lengths_N, 1);
gmid = zeros(num_lengths_N, 1);

% try each length
for i = 1:num_lengths_P
    
    % column for gmid, gmro data
    col = 1 + i;
    
    % determine corresponding length
    L_charac = L_P(i);
    W_charac = L_charac * 2;
    
    for j = 1:points_per_length
        
        % go from SI to WI to get minimum required gmid
        point = points_per_length - j + 1;
        
        gmid(i) = gmid_p(point, col);
        
        if (gmid(i) > gmid_WI)
            % found required gmid
            
            % determine corresponding current density for this length
            ID_charac = gmid_p(point, 1);
            ID_density = ID_charac / W_charac;
            
            % determine actual ID required
            ID_req = gm_min_SF1 / gmid(i);
            ID = 0.01e-6 * ceil(ID_req * 100e6);
            
            % determine actual W required
            W(i) = ID / ID_density;
            WL(i) = W(i) .* L_charac;
            
            % keep track of corresponding gmro
            gmro(i) = gmro_n(point, col);
            break
        end      
    end
end

% pick L that gives us required WL for parasitics/flicker
for i = 1:num_lengths_P
    WL_this = WL(i);
    
    if (WL_this > WL_min_SF1)
        % determine final dimensions
        W_SF1_M1 = W(i);
        L_SF1_M1 = L_P(i);
        % determine final small signal params
        gmid_SF1_M1 = gmid(i);
        gmro_SF1_M1 = gmro(i);
        break
    end
end

% corresponding current
ID_SF1 = ID;

% compute M1 small signal params
gm_SF1_M1 = gmid_SF1_M1 * ID_SF1;
ro_SF1_M1 = gmro_SF1_M1 / gm_SF1_M1;

% size M2 based on estimated gm and ro required
gm_SF1_M2_est = gmid_SI * ID_SF1;
ro_SF1_M2_required = A_SF1_required / gm_SF1_M1 / (1 - A_SF1_required);
gmro_SF1_M2_target = gm_SF1_M2_est * ro_SF1_M2_required;

% SF1_M2: load device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_SI, 0, ID_A3,...
    num_lengths_P, points_per_length);
[L_SF1_M2, W_SF1_M2, gmid_SF1_M2, gmro_SF1_M2] = pick_L(L_P, W, gmid, gmro,...
    gmro_SF1_M2_target, num_lengths_P, W_min);
% compute M2 small signal params
gm_SF1_M2 = gmid_SF1_M2 * ID_SF1;
ro_SF1_M2 = gmro_SF1_M2 / gm_SF1_M2;

A_SF1_actual = gm_SF1_M1 * ro_SF1_M2 / (1 + gm_SF1_M1 * ro_SF1_M2);

%% post-design SF characterization

% noise - note this is just open-loop OTA noise for sanity check

% compute first stage noise
noise_final_SF1 = 8*k*T*(gamma_W/gm_SF1_M1 +...
    gamma_S * gm_SF1_M2 / (gm_SF1_M1^2))*ENBW...
    + 2 * K_flicker_P / (W_SF1_M1 * L_SF1_M1 * C_ox) *...
    log(f_LP_AP / f_HP_AP);

% compute total input-referred OTA noise
noise_final_SF1_in = noise_final_SF1 / (A_target)^2;
noise_final_SF1_V = sqrt(noise_final_SF1_in);

% offset: it can be shown that all contributions are neglig except input
sigma_Vth_SF1_1 = A_VT / sqrt(W_SF1_M1 * L_SF1_M1);
V_os_3sigma_SF1 = sigma_Vth_SF1_1 * 3;

%% area estimate

% calculate A1 OTA area
area_A1 = 0;
area_A1 = area_A1 + get_area(L_M1, W_M1, L_guard, L_diff);
area_A1 = area_A1 + get_area(L_M3, W_M3, L_guard, L_diff);
area_A1 = area_A1 + get_area(L_M5, W_M5, L_guard, L_diff);
area_A1 = area_A1 + get_area(L_M7, W_M7, L_guard, L_diff);
area_A1 = area_A1 + get_area(L_M9, W_M9, L_guard, L_diff);
area_A1 = area_A1 + get_area(L_M13, W_M13, L_guard, L_diff);
area_A1 = area_A1 + get_area(L_M15, W_M15, L_guard, L_diff);

% calculate A3 OTA area
area_A3 = 0;
area_A3 = area_A3 + get_area(L_A3_M1, W_A3_M1, L_guard, L_diff);
area_A3 = area_A3 + get_area(L_A3_M3, W_A3_M3, L_guard, L_diff);
area_A3 = area_A3 + get_area(L_A3_M5, W_A3_M5, L_guard, L_diff);
area_A3 = area_A3 + get_area(L_A3_M7, W_A3_M7, L_guard, L_diff);
% assuming same CMFB amp as A1
area_A3 = area_A3 + get_area(L_M13, W_M13, L_guard, L_diff);
area_A3 = area_A3 + get_area(L_M15, W_M15, L_guard, L_diff);

% source follower area 
area_SF1 = 0;
area_SF1 = area_SF1 + get_area(L_SF1_M1, W_SF1_M1, L_guard, L_diff);
area_SF1 = area_SF1 + get_area(L_SF1_M2, W_SF1_M2, L_guard, L_diff);

% total area
area_total = area_A1 + 2*area_A3 + 2*area_SF1;

%% sizing helper functions

% W_vec: vector of widths (one for each length)
% gmid_vec: vector of actual gmid
% gmro_vec: vector of actual gmro

% L_vector: either N or P lengths
% gmid_matrix: either N or P characterization data
% gmro_matrix: either N or P characterization data
% gmid_target: min or max gmid
% inv_reg: 1 for WI i.e. min gmid, 0 for SI i.e. max gmid
% current: ID for this device
% num_lens: size of L_vector
% points_per_length: number of characterization points
function [W_vec, gmid_vec, gmro_vec] = size_M(L_vector, gmid_matrix,...
    gmro_matrix, gmid_target, inv_reg, current, num_lens, points_per_len)

    % vectors to keep track of results for comparison
    W_vec = zeros(num_lens, 1);
    gmro_vec = zeros(num_lens, 1);
    gmid_vec = zeros(num_lens, 1);

    % try each length
    for i = 1:num_lens

        % column for gmid, gmro data
        col = 1 + i;

        % determine corresponding length
        L_charac = L_vector(i);
        W_charac = L_charac * 2;
        
        % WI
        if (inv_reg == 1)
            for j = 1:points_per_len
                
                point = points_per_len - j + 1;
                gmid_vec(i) = gmid_matrix(point, col);
                
                if (gmid_vec(i) > gmid_target)
                    % found required gmid
                    % determine corresponding current density for this length
                    ID_charac = gmid_matrix(point, 1);
                    ID_density = ID_charac / W_charac;                 
                    % calculate W necessary from current for this device
                    W_vec(i) = current / ID_density;
                    % keep track of gmro
                    gmro_vec(i) = gmro_matrix(point, col);
                    break
                end
            end
        % SI
        else
            for j = 1:points_per_len
                
                point = j;
                gmid_vec(i) = gmid_matrix(point, col);
                
                if (gmid_vec(i) < gmid_target)
                    % found required gmid
                    % determine corresponding current density for this length
                    ID_charac = gmid_matrix(point, 1);
                    ID_density = ID_charac / W_charac;                 
                    % calculate W necessary from current for this device
                    W_vec(i) = current / ID_density;
                    % keep track of gmro
                    gmro_vec(i) = gmro_matrix(point, col);
                    break
                end
            end
        end
    end
end

% L_val: chosen L value
% W_val: chosen W value
% gmid_val: chosen gmid value
% gmro_val: chosen gmro value

% L_vector: either N or P lengths
% W_vector: corresponding W options
% gmid_vector: corresponding gmid
% gmro_vector: corresponding gmro
% gmro_target: min gmro
% num_lens: size of L_vector
function [L_val, W_val, gmid_val, gmro_val] = pick_L(L_vector, W_vector,...
    gmid_vector, gmro_vector, gmro_target, num_lens, min_W)

    for i = 1:num_lens
        gmro_this = gmro_vector(i);
        W_this = W_vector(i);
        
        % determine dimensions
        W_val = W_vector(i);
        L_val = L_vector(i);
        % determine small signal params
        gmid_val = gmid_vector(i);
        gmro_val = gmro_vector(i);

        if ((W_this > min_W) && (gmro_this > gmro_target))
            % if conditions met, these dimensions and params are final
            break
        end
    end
end

% area: silicon area estimate for pair of devices

% L: device L
% W: device W
% guard: estimate of additional guard-ring length
% diff: estimate of additional diffusion length
function area = get_area(L, W, guard, diff)

    area = (2*W + 2*guard) * (L + 2*guard + 2*diff);

end