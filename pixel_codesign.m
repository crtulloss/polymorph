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
V_LFP = 6e-3;

% assumed compensation cap for filter OTAs
C_load = 100e-15;
% mux line wiring capacitance
C_mux = 20e-15;

%% design targets

% frequencies
f_mux = 32 * 30e3;
f_LP_AP = 10e3;
f_HP_AP = 300;
% biquad filter params
omega_0 = 2*pi * sqrt(f_HP_AP * f_LP_AP);
BW = 2*pi * (f_LP_AP - f_HP_AP);
Q = omega_0 / BW;

% except for here, as noted, all noise is squared voltage
noise_Vrms_target = 5e-6;
noise_Vsq_target = noise_Vrms_target^2;
% target for OTA noise source other than M1
excess = 1.2;
% equivalent noise bandwidth
ENBW = f_LP_AP * pi/2;

% gains
A_target = 100;
A_OL_target = 1e3;
A_f = 1;
A_LP = 5/3;

% BPF nonidealities
gain_error = 0.01;
% stopband is defined as the point where 20dB rolloff stops and
% the low-freq zero kicks in
stopband_edge = 5;

% settling accuracy and associated time constant
settle_acc = 0.001;
settling_TCs = -log(settle_acc);
% fraction of period available for settling
settling_frac = 0.4;

% swing for A1
V_swing_target = V_LFP * A_target;

% BPF calcs: setup

% number of sweep points
num_points = 1000;
% filter design
f_swcap = logspace(5,7,num_points)';

% vector version of area limit for plots
area_limit_vec = area_limit * ones(num_points, 1);

%% technology information
% info for thicc oxide devices is given by metric3 instead of metric

L_min = 130e-9;
L_min3 = 350e-9;

% C_ox and C_P (C_gg) [F/m^2]
% C_ox for flicker noise
C_ox = 0.01;
C_ox3 = 0.005;
% C_P for parasitics
M_input_area = 120e-6 * 600e-9;
C_gg = 252e-15;
C_P_N = C_gg / M_input_area;
C_P_P = 0.00633;
C_P_N3 = 0.0017;
C_P_P3 = 0.0017;

% cap / area relationship for BPF calcs
cap_per_area_mim = 6.04e-15 / (2e-6)^2;    % [F/m^2]
cap_per_area_mos = 7.0e-15 / (1e-6)^2;
total_cap = area_limit * cap_per_area_mim;

% flicker noise coefficients [J]
K_flicker_N = 1.24e-25;
K_flicker_P = 1.98e-25;
K_flicker_N3 = 1.55e-25;
K_flicker_P3 = 4.27e-25;
K_ratio_3 = K_flicker_P3 / K_flicker_N3;

% mismatch coefficients
% [Vm] and [m]
A_VT = 10e-9;
A_beta = 0.013e-6;

% noise factor gamma (often quoted as 2/3)
gamma_W = 1;
gamma_S = 2/3;
gamma_W3 = 1.5;

% estimates for device area
L_diff = 200e-9;
L_guard = 1e-6;

% Resistance: ohms/sq
Rs = 1222;
WR_min = 400e-9;
% R mismatch coefficient [m]
A_R = 0.016e-6;



%% pseudoR cutoff and noise

% possible pseudoR values from mc: 90GOhm to 5TOhm
pseudoR = logspace(11.25,11.77,100);
% assume C_F is minimum size
C_F = 10e-15;
f_pseudo = 1./(2.*pi.*pseudoR.*C_F);

% reduction factor to account for not all of kT/C appearing in final BW
% we apply a HP filter in the SC stage, and only care about noise
% in the 300Hz+ band
ratio = f_HP_AP ./ f_pseudo;
reduc_factor = (2/pi./((ratio).^2 - 1)) .* ...
    (ratio*pi/4 - pi/2 + atan(ratio));
% resulting noise due to 2 pseudoR
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

% smallest pseudoR under mc
pseudoR = 178e9;

% sweep Cin
%C_in = 1e-12;
C_in = linspace(100e-15, 3e-12, 1000)';
C_F = C_in / A_target;

% integrated pseudoR noise for 300Hz sw cap HP cutoff
f_pseudo = 1./(2.*pi.*pseudoR.*C_F);
% reduction factor to account for not all of kT/C appearing in final BW
% we apply a HP filter in the SC stage, and only care about noise
% in the 300Hz+ band
ratio = f_HP_AP ./ f_pseudo;
reduc_factor = (2/pi./((ratio).^2 - 1)) .* ...
    (ratio*pi/4 - pi/2 + atan(ratio));
% resulting noise due to 2 pseudoR
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
ratio = f_HP_AP ./ f_pseudo;
reduc_factor = (2/pi./((ratio).^2 - 1)) .* ...
    (ratio*pi/4 - pi/2 + atan(ratio));
noise_pseudoR = 2 .* k .* T ./ C_F .* reduc_factor;



%% begin OTA design: sweep WL

% transistor area, from 10^-14 m^2 to 10^-9 m^2
% for reference, pixel is 2.33*10^-8 m^2
% 120um/600nm device is 7.2*10^-11 m^2both
% and a 130nm/130nm device is 1.69*10^-14 m^2
WL = logspace(-14,-8,1000);

C_P = WL .* C_P_N3;

% resulting actual closed-loop gain
A_actual = C_in ./ (C_F + (C_F + C_in + C_P) ./ A_OL_target);
noise_pseudoR_in = noise_pseudoR ./ (A_actual.^2);

% noise budget
noise_OTA = noise_Vsq_target - noise_pseudoR_in;

% noise gain
A_noise = (C_P + C_in + C_F) ./ C_in;

% flicker noise
noise_flicker = 2 .* K_flicker_N3 ./ C_ox3 ./ WL .* (A_noise).^2 ...
    .* log(f_LP_AP / f_HP_AP);

noise_thermal = noise_OTA - noise_flicker;

% thermal noise
% NOTE: we only care about thermal noise in the signal band (to 10kHz)
% technically, ENBW is not correct here because we can do a software
% brickwall filter at 10kHz, and this will be valid since we have used
% appropriate anti-aliasing.
% thus, ENBW is used here only to provide a safety margin.
gm_req = 2 .* 4 .* k .* T .* gamma_W3...
    ./ noise_thermal .* (A_noise).^2 .* ENBW .* excess;

% first approx of mismatch: assume input pair dominates
sigma_Vth_1 = A_VT ./ sqrt(WL);
V_os_3sigma = sigma_Vth_1 .* 3;

%% plots from WL sizing choice

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

%% begin OTA design: sweep WL (inverter-based)

% WL sweep:
% for the inverter-based option, this is the WL for the NMOS device
% in order to account for the difference in flicker noise coefficients,
% the PMOS device area will be scaled by K_flicker_P / K_flicker_N

C_P = WL .* (C_P_N3 + C_P_P3 * K_ratio_3);

% resulting actual closed-loop gain
A_actual = C_in ./ (C_F + (C_F + C_in + C_P) ./ A_OL_target);
noise_pseudoR_in = noise_pseudoR ./ (A_actual.^2);

% noise budget
noise_OTA = noise_Vsq_target - noise_pseudoR_in;

% noise gain
A_noise = (C_P + C_in + C_F) ./ C_in;

% for the inverter-based design, there are 4 devices contributing noise,
% but the contribution of each is reduced by a factor of 4 due to the
% double gm

% flicker noise
noise_flicker = K_flicker_N3 ./ C_ox3 ./ WL .* (A_noise).^2 ...
    .* log(f_LP_AP / f_HP_AP);

noise_thermal = noise_OTA - noise_flicker;

% thermal noise
% NOTE: we only care about thermal noise in the signal band (to 10kHz)
% technically, ENBW is not correct here because we can do a software
% brickwall filter at 10kHz, and this will be valid since we have used
% appropriate anti-aliasing.
% thus, ENBW is used here only to provide a safety margin.
gm_req = 4 .* k .* T .* gamma_W3...
    ./ noise_thermal .* (A_noise).^2 .* ENBW .* excess;

% first approx of mismatch: assume input pair dominates
sigma_Vth_1 = A_VT ./ sqrt(WL);
V_os_3sigma = sigma_Vth_1 .* 3;

%% plots from WL sizing choice (inverter-based)

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
title('Contributors to LNA Noise: Inverter-Based Option');

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
title('g_m Requirement: Inverter-Based Option');

% closed-loop gain
figure
semilogx(WL, A_actual);
xlabel('area_{diff} (m^2)');
ylabel('A_{actual}');
title('Closed-Loop Gain due to C_P: Inverter-Based Option');



%% gm/id, gmro characterization

% info about number of characterization sweeps
num_lengths_N = 10;
num_lengths_P = 10;
points_per_length = 25;
L_N = [130e-9 200e-9 500e-9 1e-6 2e-6 4e-6 6e-6 8e-6 10e-6 12e-6]';
L_P = L_N;

% for thick-oxide devices
num_lengths_N3 = 9;
num_lengths_P3 = 9;
L_N3 = [350e-9 500e-9 1e-6 2e-6 4e-6 6e-6 8e-6 10e-6 12e-6]';
L_P3 = L_N3;

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

gmid_N_350n = readmatrix('charac/charac_gmid_N_350n.csv');
gmid_N_500n = readmatrix('charac/charac_gmid_N_500n.csv');
gmid_N_1u = readmatrix('charac/charac_gmid_N_1u.csv');
gmid_N_2u = readmatrix('charac/charac_gmid_N_2u.csv');
gmid_N_4u = readmatrix('charac/charac_gmid_N_4u.csv');
gmid_N_6u = readmatrix('charac/charac_gmid_N_6u.csv');
gmid_N_8u = readmatrix('charac/charac_gmid_N_8u.csv');
gmid_N_10u = readmatrix('charac/charac_gmid_N_10u.csv');
gmid_N_12u = readmatrix('charac/charac_gmid_N_12u.csv');

gmro_N_350n = readmatrix('charac/charac_gmro_N_350n.csv');
gmro_N_500n = readmatrix('charac/charac_gmro_N_500n.csv');
gmro_N_1u = readmatrix('charac/charac_gmro_N_1u.csv');
gmro_N_2u = readmatrix('charac/charac_gmro_N_2u.csv');
gmro_N_4u = readmatrix('charac/charac_gmro_N_4u.csv');
gmro_N_6u = readmatrix('charac/charac_gmro_N_6u.csv');
gmro_N_8u = readmatrix('charac/charac_gmro_N_8u.csv');
gmro_N_10u = readmatrix('charac/charac_gmro_N_10u.csv');
gmro_N_12u = readmatrix('charac/charac_gmro_N_12u.csv');

cgg_N_350n = readmatrix('charac/charac_cgg_N_350n.csv');
cgg_N_500n = readmatrix('charac/charac_cgg_N_500n.csv');
cgg_N_1u = readmatrix('charac/charac_cgg_N_1u.csv');
cgg_N_2u = readmatrix('charac/charac_cgg_N_2u.csv');
cgg_N_4u = readmatrix('charac/charac_cgg_N_4u.csv');
cgg_N_6u = readmatrix('charac/charac_cgg_N_6u.csv');
cgg_N_8u = readmatrix('charac/charac_cgg_N_8u.csv');
cgg_N_10u = readmatrix('charac/charac_cgg_N_10u.csv');
cgg_N_12u = readmatrix('charac/charac_cgg_N_12u.csv');

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

gmid_N = horzcat(gmid_N_350n, gmid_N_500n(:,2),...
    gmid_N_1u(:,2), gmid_N_2u(:,2), gmid_N_4u(:,2), gmid_N_6u(:,2),...
    gmid_N_8u(:,2), gmid_N_10u(:,2), gmid_N_12u(:,2));
gmro_N = horzcat(gmro_N_350n, gmro_N_500n(:,2),...
    gmro_N_1u(:,2), gmro_N_2u(:,2), gmro_N_4u(:,2), gmro_N_6u(:,2),...
    gmro_N_8u(:,2), gmro_N_10u(:,2), gmro_N_12u(:,2));
cgg_N = horzcat(cgg_N_350n, cgg_N_500n(:,2),...
    cgg_N_1u(:,2), cgg_N_2u(:,2), cgg_N_4u(:,2), cgg_N_6u(:,2),...
    cgg_N_8u(:,2), cgg_N_10u(:,2), cgg_N_12u(:,2));

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

gmid_P_350n = readmatrix('charac/charac_gmid_P_350n.csv');
gmid_P_500n = readmatrix('charac/charac_gmid_P_500n.csv');
gmid_P_1u = readmatrix('charac/charac_gmid_P_1u.csv');
gmid_P_2u = readmatrix('charac/charac_gmid_P_2u.csv');
gmid_P_4u = readmatrix('charac/charac_gmid_P_4u.csv');
gmid_P_6u = readmatrix('charac/charac_gmid_P_6u.csv');
gmid_P_8u = readmatrix('charac/charac_gmid_P_8u.csv');
gmid_P_10u = readmatrix('charac/charac_gmid_P_10u.csv');
gmid_P_12u = readmatrix('charac/charac_gmid_P_12u.csv');

gmro_P_350n = readmatrix('charac/charac_gmro_P_350n.csv');
gmro_P_500n = readmatrix('charac/charac_gmro_P_500n.csv');
gmro_P_1u = readmatrix('charac/charac_gmro_P_1u.csv');
gmro_P_2u = readmatrix('charac/charac_gmro_P_2u.csv');
gmro_P_4u = readmatrix('charac/charac_gmro_P_4u.csv');
gmro_P_6u = readmatrix('charac/charac_gmro_P_6u.csv');
gmro_P_8u = readmatrix('charac/charac_gmro_P_8u.csv');
gmro_P_10u = readmatrix('charac/charac_gmro_P_10u.csv');
gmro_P_12u = readmatrix('charac/charac_gmro_P_12u.csv');

cgg_P_350n = readmatrix('charac/charac_cgg_P_350n.csv');
cgg_P_500n = readmatrix('charac/charac_cgg_P_500n.csv');
cgg_P_1u = readmatrix('charac/charac_cgg_P_1u.csv');
cgg_P_2u = readmatrix('charac/charac_cgg_P_2u.csv');
cgg_P_4u = readmatrix('charac/charac_cgg_P_4u.csv');
cgg_P_6u = readmatrix('charac/charac_cgg_P_6u.csv');
cgg_P_8u = readmatrix('charac/charac_cgg_P_8u.csv');
cgg_P_10u = readmatrix('charac/charac_cgg_P_10u.csv');
cgg_P_12u = readmatrix('charac/charac_cgg_P_12u.csv');

% concatenate the results for easier looping
gmid_p = horzcat(gmid_p_130n, gmid_p_200n(:,2), gmid_p_500n(:,2),...
    gmid_p_1u(:,2), gmid_p_2u(:,2), gmid_p_4u(:,2), gmid_p_6u(:,2),...
    gmid_p_8u(:,2), gmid_p_10u(:,2), gmid_p_12u(:,2));
gmro_p = horzcat(gmro_p_130n, gmro_p_200n(:,2), gmro_p_500n(:,2),...
    gmro_p_1u(:,2), gmro_p_2u(:,2), gmro_p_4u(:,2), gmro_p_6u(:,2),...
    gmro_p_8u(:,2), gmro_p_10u(:,2), gmro_p_12u(:,2));

gmid_P = horzcat(gmid_P_350n, gmid_P_500n(:,2),...
    gmid_P_1u(:,2), gmid_P_2u(:,2), gmid_P_4u(:,2), gmid_P_6u(:,2),...
    gmid_P_8u(:,2), gmid_P_10u(:,2), gmid_P_12u(:,2));
gmro_P = horzcat(gmro_P_350n, gmro_P_500n(:,2),...
    gmro_P_1u(:,2), gmro_P_2u(:,2), gmro_P_4u(:,2), gmro_P_6u(:,2),...
    gmro_P_8u(:,2), gmro_P_10u(:,2), gmro_P_12u(:,2));
cgg_P = horzcat(cgg_P_350n, cgg_P_500n(:,2),...
    cgg_P_1u(:,2), cgg_P_2u(:,2), cgg_P_4u(:,2), cgg_P_6u(:,2),...
    cgg_P_8u(:,2), cgg_P_10u(:,2), cgg_P_12u(:,2));

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
% PI: "pretty inverted", P is between M and S
gmid_WI = 25;
gmid_MI = 20;
gmid_PI = 15;
gmid_SI = 10;

% M1 specs determined from optimization above
gm_min = 150e-6;
WL_min = 6e-11;

% ratio between main branch and output branch
current_scale = 16;
% ratio between load branch and one branch of CMFB amp
cmfb_scale = 1;

% partition of gains
A_OL_CMFB_target = 20;

% assumption that (gmro)9ro11 << (gmro)7(gmro)5R,
% but this will not be possible while meeting other requirements
% to account for this, let gg_factor_1(gmro)9ro11 = (gmro)7(gmro)5R
gg_factor_1 = 0.5;
% assumption that (gmro)3ro1 >> (gmro)5R.
% let (gmro)3ro1 = gg_factor_2*(gmro)5R
gg_factor_2 = 5;

% only assumed gmro - it will be overwritten later so it's ok
gmro_9_target = 30;

% degen resistor headroom consumption
% chosen to achieve a tradeoff between Rout and swing
V_R_degen = 0.2;
% estimated headroom for saturation for a single device (even SI)
V_hr_SI = 300e-3;
V_hr_WI = 250e-3;

%% swing requirement calculation

V_swing_telecas = V_DD - 1*V_hr_SI - 4*V_hr_WI - V_R_degen;
V_swing_foldcas = V_DD - 1*V_hr_SI - 3*V_hr_WI - V_R_degen;
V_swing_stage2 = V_DD - V_hr_SI - V_hr_WI;

%% transistor-level OTA design: action

% M1: input pair (NMOS)

% vectors to keep track of results for comparison
W = zeros(num_lengths_N3, 1);
WL = zeros(num_lengths_N3, 1);
gmro = zeros(num_lengths_N3, 1);
gmid = zeros(num_lengths_N3, 1);

% try each length
for i = 1:num_lengths_N3
    
    % column for gmid, gmro data
    col = 1 + i;
    
    % determine corresponding length
    L_charac = L_N3(i);
    W_charac = L_charac * 2;
    
    for j = 1:points_per_length
        
        % go from SI to WI to get minimum required gmid
        point = points_per_length - j + 1;
        
        gmid(i) = gmid_N(point, col);
        
        if (gmid(i) > gmid_WI)
            % found required gmid
            
            % determine corresponding current density for this length
            ID_charac = gmid_N(point, 1);
            ID_density = ID_charac / W_charac;
            
            % determine actual ID required
            ID_req = gm_min / gmid(i);
            ID = 1e-6 * ceil(ID_req * 1e6);
            
            % determine actual W required
            W(i) = ID / ID_density;
            WL(i) = W(i) .* L_charac;
            
            % keep track of corresponding gmro
            gmro(i) = gmro_N(point, col);
            break
        end      
    end
end

% plot results for M1
figure
plot(L_N3*1e6, W*1e6);
xlabel('L (\mum)');
ylabel('W (\mum)');
title('M1: W vs. L');
figure
plot(L_N3*1e6, WL*1e12);
xlabel('L (\mum)');
ylabel('W ((\mum)^2)');
title('M1: WL vs. L');
figure
plot(L_N3*1e6, gmro);
xlabel('L (\mum)');
ylabel('g_mr_o');
title('M1: g_mr_o vs. L');

% pick L that gives us required WL for parasitics/flicker
for i = 1:num_lengths_N3
    WL_this = WL(i);
    
    if (WL_this > WL_min)
        % determine final dimensions
        W_M1 = W(i);
        L_M1 = L_N3(i);
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
R_degen = V_R_degen / ID_both;

% confirm resistor size is ok
W_R = WR_min;
L_R = R_degen / Rs * W_R;
sigma_R = A_R ./ sqrt(L_R .* W_R) * R_degen;

% M9 and M11 based on gain equation
gmro9ro11_target = A_OL_target * (1+1/gg_factor_1)^2 / gm_M1;

% calculte M11 requirements
gm_M11_assumed = gmid_SI * ID_load;
ro_M11_target = gmro9ro11_target / gmro_9_target;
gmro_11_target = gm_M11_assumed * ro_M11_target;

% M11: load device (NMOS)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_PI, 0, ID_load,...
    num_lengths_N, points_per_length);
[L_M11, W_M11, gmid_M11, gmro_M11] = pick_L(L_N, W, gmid, gmro,...
    gmro_11_target, num_lengths_N, W_min);
% compute M11 small signal params
gm_M11 = gmid_M11 * ID_load;
ro_M11 = gmro_M11 / gm_M11;

% re-calculate gmro9 target
gmro_9_target = gmro9ro11_target / ro_M11;

% M9: load cascode device (NMOS)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_MI, 1, ID_load,...
    num_lengths_N, points_per_length);
[L_M9, W_M9, gmid_M9, gmro_M9] = pick_L(L_N, W, gmid, gmro,...
    gmro_9_target, num_lengths_N, W_min);
% compute M9 small signal params
gm_M9 = gmid_M9 * ID_load;
ro_M9 = gmro_M9 / gm_M9;

% gmro_57 based on ratio gg_factor_1 (see above)
gmro_57_target = sqrt(gg_factor_1 * gmro_M9 * ro_M11 / R_degen);

% size 5,7 based on this target

% M5: current source device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_MI, 1, ID_both,...
    num_lengths_P, points_per_length);
[L_M5, W_M5, gmid_M5, gmro_M5] = pick_L(L_P, W, gmid, gmro,...
    gmro_57_target, num_lengths_P, W_min);
% compute M5 small signal params
gm_M5 = gmid_M5 * ID_main;
ro_M5 = gmro_M5 / gm_M5;

% M7: folded cascode device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_MI, 1, ID_load,...
    num_lengths_P, points_per_length);
[L_M7, W_M7, gmid_M7, gmro_M7] = pick_L(L_P, W, gmid, gmro,...
    gmro_57_target, num_lengths_P, W_min);
% compute M7 small signal params
gm_M7 = gmid_M7 * ID_load;
ro_M7 = gmro_M7 / gm_M7;

% size M3 such that (gmro)3ro1 >> (gmro)5R
gmro_3_target = gg_factor_2 * gmro_M5 * R_degen / ro_M1;

% M3: telecascode device (NMOS)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_MI, 1, ID_main,...
    num_lengths_N, points_per_length);
[L_M3, W_M3, gmid_M3, gmro_M3] = pick_L(L_N, W, gmid, gmro,...
    gmro_3_target, num_lengths_N, W_min);
% compute M3 small signal params
gm_M3 = gmid_M3 * ID_main;
ro_M3 = gmro_M3 / gm_M3;

% calculate the actual gain

% impedance looking into M3 drain
r_inbranch = gmro_M3 * ro_M1;
% impedance looking into M5 drain
r_top = gmro_M5 * R_degen;
% impedance looking into folding node
r_fold = 1 / (1/(r_top) + 1/(r_inbranch));
% impedance looking up from output
r_up = r_fold * gmro_M7;
% impedance looking down from output
r_down = gmro_M9 * ro_M11;

% output resistance
Rout_actual = 1/(1/r_down + 1/r_up);

% transconductance
Gs7 = (gm_M7 + 1/ro_M7) / (1 + gmro_M9 * ro_M11 / ro_M7);
Gd5 = 1 / (ro_M5 * (1 + R_degen / ro_M5 + gm_M5 * R_degen));
Gs3 = (gm_M3 + 1/ro_M3) / (1 + 1 / (ro_M3 * (Gs7 + Gd5)));

input_split = Gs3 * ro_M1 / (1 + Gs3 * ro_M1);
fold_split = Gs7 / (Gs7 + Gd5);

Gm_actual = gm_M1 * input_split * fold_split;

A_OL_actual = Rout_actual * Gm_actual;

%% post-design OTA characterization

% noise - note this is just open-loop OTA noise for sanity check

% compute first stage noise
noise_final_A1 = 8*k*T*(gamma_W/gm_M1 + 1/(R_degen * gm_M1^2))*ENBW...
    + 2 * K_flicker_N / (W_M1 * L_M1 * C_ox) * log(f_LP_AP / f_HP_AP);

% compute total input-referred OTA noise
noise_final_OTA = noise_final_A1;
noise_final_OTA_V = sqrt(noise_final_OTA);

% offset: it can be shown that all contributions are neglig except input
sigma_Vth_1_final = A_VT / sqrt(W_M1 * L_M1);
V_os_3sigma_final = sigma_Vth_1_final * 3;

% settling - need to satisfy swcap settling assumption

% assume min size cap for C_RA
C_RA = 10e-15;
% input cap size
C_R1 = C_RA * A_f / Q;

% settling limit
T_settle_limit = (1 ./ f_swcap) * settling_frac;

% assume switch resistance - settling through this
R_sw = 14e3;
T_settle_sw = R_sw * C_R1 * settling_TCs;

loop_gain_approx = A_OL_actual / A_target;
T_settle_OTA = (Rout_actual / loop_gain_approx) * C_R1 * settling_TCs;

T_settle_OTA_vec = T_settle_OTA .* ones(size(f_swcap));
T_settle_sw_vec = T_settle_sw .* ones(size(f_swcap));

figure
semilogx(f_swcap, T_settle_limit*1e6);
hold on
semilogx(f_swcap, T_settle_sw_vec*1e6);
semilogx(f_swcap, T_settle_OTA_vec*1e6);
xlabel('f_{swcap} (Hz)');
ylabel('n\tau (\mus)');
title('Amp/Filter Settling Time vs. Sw Cap Frequency');
legend('Limit', 'Switch Settling', 'OTA 2nd Stage Settling');

% fold node resistance - for comparison with output R to determine dom pole
R_fold = 1/(Gs7 + 1/r_inbranch + 1/r_top);

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



%% area contrib of first-stage CT filter capacitor and signal-path caps
area_sig = 2 .* C_in ./ cap_per_area_mim;

figure
semilogx(f_swcap, area_limit_vec*(1e12));
hold on
semilogx(f_swcap, area_sig * ones(num_points, 1)*(1e12));

f_LPCT = f_swcap ./ 10;
C_L = Gm_actual ./ (2*pi .* f_LPCT .* A_target);

area_CT = C_L ./ cap_per_area_mos + area_sig;

semilogx(f_swcap, area_CT*(1e12));

%% design based on CT passive second-order BPF

% C_R2 will be the smallest
C_R2 = 10e-15;
R_2 = 1 ./ (f_swcap .* C_R2);
% R_1 must be much smaller in order to gain unity gain in passband
C_R1 = C_R2 .* 10;
R_1 = 1 ./ (f_swcap .* C_R1);
% corresponding caps
C_1 = 1 ./ (2*pi*f_HP_AP.*R_1);
C_2 = 1 ./ (2*pi*f_LP_AP.*R_2);

% for differential operation, need to double shunt R and half shunt C
C_passive = 2*C_1 + C_R1/2 + 2*C_R2 + C_2/2;

% area
area_passive = C_passive ./ cap_per_area_mim;

% plots
semilogx(f_swcap, area_passive*(1e12));
semilogx(f_swcap, (area_passive+area_CT)*(1e12));

%% design based on CT first-order BPF (see Vol. 1 p. 185)

% swcap filter cap calculations
C_4 = 10e-15;
noise_inref = sqrt(k*T ./ C_4) ./ A_target;

C_2 = f_swcap .* C_4 ./ (2*pi*f_HP_AP);

C_1 = A_f .* C_2;

C_3 = (2*pi*f_LP_AP) .* C_1 ./ f_swcap;

C_firstorder = 2.*(C_1 + C_2 + C_3 + C_4);

% area
area_firstorder = C_firstorder ./ cap_per_area_mim;

% plots
semilogx(f_swcap, area_firstorder*(1e12));
semilogx(f_swcap, (area_firstorder+area_CT)*(1e12));

%% design based on CT biquad - TT

% swcap filter cap calculations
C_RA = 10e-15;

R_A = 1 ./ (C_RA .* f_swcap);
C_A = 1 ./ (omega_0 .* R_A);

R_B = R_A .* Q;
C_RB = 1 ./ (R_B .* f_swcap);

R_1 = R_B ./ A_f;
C_R1 = 1 ./ (R_1 .* f_swcap);

C_biquad = 2.*(2.*(C_A + C_RA) + C_RB + C_R1);

% area
area_biquad = C_biquad ./ cap_per_area_mim;

% plots
semilogx(f_swcap, area_biquad*(1e12));
semilogx(f_swcap, (area_biquad+area_CT)*(1e12));

%% design based on CT biquad - TT modified for high-Q

% swcap filter cap calculations
C_RAHQ = 10e-15;

R_AHQ = 1 ./ (C_RAHQ .* f_swcap);
C_AHQ = 1 ./ (omega_0 .* R_AHQ);

C_BHQ = C_AHQ ./ Q;

C_2HQ = A_f .* C_BHQ;

C_HQ = 2.*(2.*(C_AHQ + C_RAHQ) + C_BHQ + C_2HQ);

% area
area_HQ = C_HQ ./ cap_per_area_mim;

% plots
semilogx(f_swcap, area_HQ*(1e12));
semilogx(f_swcap, (area_HQ+area_CT)*(1e12));

%% design based on CT biquad - TT with specified LPF gain

% ok to overwrite cap values if we end up using this topology

% 2 degrees of freedom: C_RA and C_RAP

% swcap filter cap calculations
C_RA = 40e-15;
R_A = 1 ./ (C_RA .* f_swcap);

R_1 = R_A ./ A_LP;
C_R1 = 1 ./ (R_1 .* f_swcap);

R_B = R_1 .* A_f;
C_RB = 1 ./ (R_B .* f_swcap);

C_A = 1./ (R_B .* BW);

C_RAP = 10e-15;
R_AP = 1 ./ (C_RAP .* f_swcap);

C_AP = 1 ./ (R_A .* R_AP .* C_A .* omega_0.^2);

C_biquad = 2.*(C_A + C_RA + C_AP + C_RAP + C_RB + C_R1);

% area
area_biquad = C_biquad ./ cap_per_area_mim;

% plots
semilogx(f_swcap, area_biquad*(1e12));
semilogx(f_swcap, (area_biquad+area_CT)*(1e12));

%% plot settings
xlabel('f_{swcap} (Hz)');
ylabel('Capacitor Area (\mum)^2');
title('Capacitor Area vs. Sw Cap Frequency');
legend('Limit', 'LNA Signal Path', 'CT',...
    'Passive', 'CT+Passive',...
    'FirstOrder', 'CT+FirstOrder',...
    'TTBiquad', 'CT+TTBiquad',...
    'HQBiquad', 'CT+HQBiquad',...
    'TTBSwing', 'CT+TTBSwing');
axis([1e5, 1e7, 0, 25000]);

%% further analysis of TT topology

% noise contribution of switches at each capacitor
% noise in power!
% assume that this sampled noise will appear spread over the fs/2 BW
noise_C_RB = 2.*k.*T ./ (C_A + C_RB);
noise_C_RAP = 2.*k.*T .* C_RA ./ (C_AP .* C_RB);
noise_C_RA = 2.*k.*T .* C_RA ./ (C_RB .* (C_A + C_RB));

% total noise - assume neglig contrib at C_A,
% and assume C_R1 yields same contrib as C_RA
noise_total = noise_C_RB + noise_C_RAP + 2.*noise_C_RA;
noise_total_sqrt = sqrt(noise_total);

% plot as a function of fs
figure
semilogx(f_swcap, sqrt(noise_C_RB)*1e6/A_target);
hold on
semilogx(f_swcap, sqrt(noise_C_RAP)*1e6/A_target);
semilogx(f_swcap, sqrt(noise_C_RA)*1e6/A_target);
semilogx(f_swcap, noise_total_sqrt*1e6/A_target);
xlabel('f_{swcap} (Hz)');
ylabel('V_{n-in-rms} (\muV)');
title('Noise Approximation on Caps vs. Sw Cap Frequency');
legend('C_{RB}', 'C_{RA1}', 'C_{RA2},C_{R1}', 'Total');



%% Source Follower Required Specs

% first calculate settling time: 0.4 of clock period
settling_time = settling_frac ./ f_swcap;
% associated time constant required
time_const = settling_time / settling_TCs;

% calculate gm required for this time constant
gm_SF1_req = C_R1 ./ time_const;

% power calculation
I_SF1_req = 2 .* (gm_SF1_req / gmid_WI);

%% Mux Driver Required Specs (A2 or SF2)

% first calculate settling time: 0.4 of mux clock period
settling_time_mux = settling_frac ./ f_mux;
% associated time constant required
time_const_mux = settling_time_mux / settling_TCs;

% if we use A2 to drive the mux line directly
% calculate equivalent cap used in time constant formula
cap_eq_2_mux = ((C_load+C_mux) .* max(C_R1, C_RA) +...
    (C_load+C_mux) .* (C_RB + C_A) +...
    max(C_R1, C_RA) .* (C_RB + C_A)) ./ (C_RB + C_A);
% calculate required gm
Gm_2_req_mux = cap_eq_2_mux ./ time_const_mux;
I_2_req_mux = 2 .* (Gm_2_req_mux / gmid_WI);

% if we use a source follower
% calculate gm required for this time constant
gm_SF2_req = C_mux ./ time_const_mux;
% power calculation
I_SF2_req = 2 .* (gm_SF2_req / gmid_WI);

%% Opamp Required Specs

% gain specs
A_3_required = 1 ./ (R_AP .* C_AP .* 2 * pi * stopband_edge);
A_2_required = (1-gain_error).*(1 + A_f) ./ gain_error;

% swing specs
A2_swing = V_AP .* A_target .* A_f;
A3_swing = V_LFP .* A_target .* A_LP;

% settling specs
% required Gm for each amp
cap_eq_2 = (C_load .* max(C_R1, C_RA) + C_load .* (C_RB + C_A) +...
    max(C_R1, C_RA) .* (C_RB + C_A)) ./ (C_RB + C_A);
cap_eq_3 = (C_load .* C_RAP + C_load .* C_AP + C_RAP .* C_AP) ./ C_AP;
Gm_2_req = cap_eq_2 ./ time_const;
Gm_3_req = cap_eq_3 ./ time_const;

% power calculation from linear-charging assumption
I_2_req = 2 .* (Gm_2_req / gmid_WI);
I_3_req = 2 .* (Gm_3_req / gmid_WI);
P_swcap = V_DD .* (I_2_req + I_3_req + I_SF1_req);
P_swcap_A2mux = V_DD .* (I_2_req_mux + I_3_req + I_SF1_req);
P_swcap_SFmux = V_DD .* (I_2_req + I_3_req + I_SF1_req + I_SF2_req);

% now what if amp is slewing?
slew_cap_2 = 1 ./ (1./(C_RB + C_A) + 1./(C_R1 + C_RA));
slew_cap_3 = 1 ./ (1./(C_AP) + 1./(C_RAP));
% assume main signal at A2 output is AP
Iss_2_req = A2_swing ./ settling_time .* slew_cap_2;
% assume main signal at A3 output is LFP
Iss_3_req = A3_swing ./ settling_time .* slew_cap_3;

% another view of bandwidth requirements: GBW
gain_A2 = max(C_RA, C_R1) ./ (C_RB + C_A);
gain_A3 = C_RAP ./ C_AP;
% calculate required BW
BW_req = 1 ./ (time_const .* 2 .* pi);
% calculate required GBW
GBW_A2 = gain_A2 .* BW_req;
GBW_A3 = gain_A3 .* BW_req;

figure
semilogx(f_swcap, Gm_2_req*1e6);
hold on
semilogx(f_swcap, Gm_2_req_mux*1e6);
semilogx(f_swcap, Gm_3_req*1e6);
semilogx(f_swcap, gm_SF1_req*1e6);
semilogx(f_swcap, gm_SF2_req*1e6*ones(size(gm_SF1_req)));
xlabel('f_{swcap} (Hz)');
ylabel('G_m (\muS)');
title('Required Transconductance for Linear Cap Charging');
legend('A_2', 'A_2 driving mux', 'A_3', 'SF_1', 'SF_2');

figure
semilogx(f_swcap, Iss_2_req*1e6);
hold on
semilogx(f_swcap, Iss_3_req*1e6);
xlabel('f_{swcap} (Hz)');
ylabel('I_{SS} (\muA)');
title('Required Opamp Current for Slewing Cap Charging');
legend('A_2', 'A_3');

figure
semilogx(f_swcap, GBW_A2/1e3);
hold on
semilogx(f_swcap, GBW_A3/1e3);
xlabel('f_{swcap} (Hz)');
ylabel('GBW (kHz)');
title('Required GBW for Closed-Loop Integrator Gain');
legend('A_2', 'A_3');

figure
yyaxis left
semilogx(f_swcap, P_swcap*1e6);
hold on
semilogx(f_swcap, P_swcap_A2mux*1e6);
semilogx(f_swcap, P_swcap_SFmux*1e6);
ylabel('P (\muW)');
yyaxis right
semilogx(f_swcap, (area_biquad+area_CT)*(1e12));
ylabel('A (\mum)^2');
xlabel('f_{swcap} (Hz)');
title('Tradeoff Between Power and Area');
legend('P: no mux', 'P: A2 drives mux', 'P: SF drives mux', 'A');



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

% TODO MT2 from maximum expected current. plug in gm/id in MI, use max I
% to find W/L, use min L

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



%% appendix: confirming biquad frequency response (ideal)

% uses TT design with specified LPF gain
f = logspace(1,6,1000)';
jomega = 1i * f .* 2 .* pi;
num = jomega .* BW;
den = jomega.^2 + jomega .* BW + omega_0.^2;
TF = num./den;
TF_dB = 20.*log10(abs(TF));

num_LP = omega_0.^2 .* A_LP;
TF_LP = num_LP ./ den;
TF_LP_dB = 20.*log10(abs(TF_LP));

figure
semilogx(f, TF_dB);
axis([1e1, 1e6, -60, 20]);
hold on
semilogx(f, TF_LP_dB);
xlabel('f (Hz)');
ylabel('|H(f)| (dB)');
title('Desired biquad frequency response');
legend('BPF Output', 'LPF Output');

%% appendix: validation of z-domain transfer function

% uses TT design with specified LPF gain

% need to use only one switched-cap frequency for this
f_swcap = 480e3;

% swcap filter cap calculations
C_RA = 40e-15;
R_A = 1 ./ (C_RA .* f_swcap);

R_1 = R_A ./ A_LP;
C_R1 = 1 ./ (R_1 .* f_swcap);

R_B = R_1 .* A_f;
C_RB = 1 ./ (R_B .* f_swcap);

C_A = 1./ (R_B .* BW);

C_RAP = 10e-15;
R_AP = 1 ./ (C_RAP .* f_swcap);

C_AP = 1 ./ (R_A .* R_AP .* C_A .* omega_0.^2);

% validation
f = logspace(-3,6,1000)';
jomega = 1i .* f .* 2 .* pi;
num = -jomega .* C_R1 .* f_swcap ./ C_A;
den = (jomega).^2 + jomega .* C_RB .* f_swcap ./ C_A +...
    C_RA .* C_RAP .* (f_swcap).^2 ./ C_A ./ C_AP;
TF = num./den;
TF_dB = 20.*log10(abs(TF));

num_LP = -C_R1 .* C_RAP .* (f_swcap).^2 ./ C_A ./ C_AP;
TF_LP = num_LP ./ den;
TF_LP_dB = 20.*log10(abs(TF_LP));

z = exp(jomega ./ f_swcap);
diff = 1 - z.^-1;

num = -C_R1 .* C_AP .* diff;
den = C_A.*C_AP .* (diff).^2 + C_RB .* C_AP .* diff + C_RA.*C_RAP .* z.^-1;
num_LP = -C_R1 .* C_RAP;

ZTF = num./den;
ZTF_dB = 20.*log10(abs(ZTF));

ZTF_LP = num_LP ./ den;
ZTF_LP_dB = 20.*log10(abs(ZTF_LP));

figure
semilogx(f, TF_dB);
axis([1e-3, f_swcap/2, -60, 20]);
hold on
semilogx(f, ZTF_dB);
semilogx(f, TF_LP_dB);
semilogx(f, ZTF_LP_dB);
xlabel('f (Hz)');
ylabel('|H(f)| (dB)');
legend('BPF CT approx', 'BPF Actual DT', 'LPF CT approx', 'LPF Actual DT');
title('Comparison between CT approx and Actual DT (z-domain)');



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