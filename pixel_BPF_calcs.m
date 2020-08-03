% polymorph active swcap BPF calculations
% CRT
% 03/25/2020 - 07/22/2020

%% housekeeping
close all
clear all

%% design parameters

% basic
V_DD = 1.5;
gmid_WI = 25;

% info from pixel system-level design
C_in = 1e-12;

% number of sweep points
num_points = 1000;

% filter design
f_swcap = logspace(5,7,num_points)';
%f_swcap = 480e3;
f_HP = 300;
f_LP = 10e3;

% mux freq
f_mux = 32 * 30e3;

% first stage gain
A_1 = 100;
A_OL1 = 10000;
% first stage transconductance
g_m1 = 230e-6;
R_O1 = A_OL1 / g_m1;

% filter stage gain
A_f = 1;
% parasitic LPF gain (output of A3)
A_LP = 5/3;

% expected signal amplitudes
sig_AP = 2e-3;
sig_LFP = 3e-3;

% amplitude range
V_pp = 4e-3;
V_p = V_pp / 2;
V_p_single = V_p / 2;

% nonidealities
gain_error = 0.01;
% stopband is defined as the point where 20dB rolloff stops and
% the low-freq zero kicks in
stopband_edge = 5;

% settling accuracy and associated time constant
settle_acc = 0.001;
settling_TCs = -log(settle_acc);
% assumed compensation cap for OTAs
C_load = 100e-15;
% mux line wiring capacitance
C_mux = 20e-15;
% fraction of period available for settling
settling_frac = 0.4;

%% physical constants
k = 1.38e-23;
T = 300;

%% technology information
cap_per_area = 6.04e-15 / (2*2);    % [F/(um)^2]
cap_per_area_mos = 7.0e-15;
total_area = 152.5^2;
total_cap = total_area * cap_per_area;

area_limit = total_area * ones(num_points, 1);

%% area contrib of LNA signal-path caps
area_sig = 2 .* C_in ./ cap_per_area;

figure
semilogx(f_swcap, area_limit);
hold on
semilogx(f_swcap, area_sig * ones(num_points, 1));

%% area contrib of first-stage CT filter capacitor and signal-path caps
f_LPCT = f_swcap ./ 10;
C_L = g_m1 ./ (2*pi .* f_LPCT .* A_1);

area_CT = C_L ./ cap_per_area_mos + area_sig;

semilogx(f_swcap, area_CT);

%% design based on CT first-order BPF (see Vol. 1 p. 185)

% swcap filter cap calculations
C_4 = 10e-15;
noise_inref = sqrt(k*T ./ C_4) ./ A_1;

C_2 = f_swcap .* C_4 ./ (2*pi*f_HP);

C_1 = A_f .* C_2;

C_3 = (2*pi*f_LP) .* C_1 ./ f_swcap;

C_firstorder = 2.*(C_1 + C_2 + C_3 + C_4);

% area
area_firstorder = C_firstorder ./ cap_per_area;

% plots
semilogx(f_swcap, area_firstorder);
semilogx(f_swcap, area_firstorder+area_CT);

%% design based on CT biquad - TT

% biquad filter params
omega_0 = 2*pi*sqrt(f_HP .* f_LP);
BW = 2*pi*(f_LP - f_HP);
Q = omega_0 ./ BW;

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
area_biquad = C_biquad ./ cap_per_area;

% plots
semilogx(f_swcap, area_biquad);
semilogx(f_swcap, area_biquad+area_CT);

%% design based on CT biquad - TT modified for high-Q

% swcap filter cap calculations
C_RAHQ = 10e-15;

R_AHQ = 1 ./ (C_RAHQ .* f_swcap);
C_AHQ = 1 ./ (omega_0 .* R_AHQ);

C_BHQ = C_AHQ ./ Q;

C_2HQ = A_f .* C_BHQ;

C_HQ = 2.*(2.*(C_AHQ + C_RAHQ) + C_BHQ + C_2HQ);

% area
area_HQ = C_HQ ./ cap_per_area;

% plots
semilogx(f_swcap, area_HQ);
semilogx(f_swcap, area_HQ+area_CT);

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
area_biquad = C_biquad ./ cap_per_area;

% plots
semilogx(f_swcap, area_biquad);
semilogx(f_swcap, area_biquad+area_CT);

%% plot settings
xlabel('f_{swcap} (Hz)');
ylabel('Capacitor Area (\mum)^2');
title('Capacitor Area vs. Sw Cap Frequency');
legend('Limit', 'LNA Signal Path', 'CT',...
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
semilogx(f_swcap, sqrt(noise_C_RB)*1e6/A_1);
hold on
semilogx(f_swcap, sqrt(noise_C_RAP)*1e6/A_1);
semilogx(f_swcap, sqrt(noise_C_RA)*1e6/A_1);
semilogx(f_swcap, noise_total_sqrt*1e6/A_1);
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
A2_swing = sig_AP .* A_1 .* A_f;
A3_swing = sig_LFP .* A_1 .* A_LP;

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
semilogx(f_swcap, area_biquad+area_CT);
ylabel('A (\mum)^2');
xlabel('f_{swcap} (Hz)');
title('Tradeoff Between Power and Area');
legend('P: no mux', 'P: A2 drives mux', 'P: SF drives mux', 'A');

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