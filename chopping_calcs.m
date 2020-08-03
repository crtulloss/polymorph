% polymorph AFE chopping calculations
% CRT
% 01/14/2020

%% wang paper - current-mode capacitively coupled chopping amp

% specifications
gain_stage1 = 100;
gain_stage2 = 1;
omega_hp1 = 2 * pi * 100;
omega_hp2 = 2 * pi * 0.1;
omega_lp_amp1 = 2 * pi * 100e6;
omega_lp_amp2 = 2 * pi * 100e6;

% component values
C_a = 9e-15;
C_in = C_a * gain_stage1;
C_b = 100e-15;
C_x = C_b * gain_stage2;
R_a = 1 / (omega_hp1 * C_a);
R_b = 1 / (omega_hp2 * C_b);
A_O_1 = 1e6;
A_O_2 = 1e6;
omega_0 = 2 * pi * 0;

% node frequencies
G_a = 1 / R_a;
omega_in1 = G_a/C_in - 1i * omega_0;
G_b = 1 / R_b;
omega_in2 = G_b/C_x + 1i * omega_0;

% frequency input
f = logspace(-3,15, 1000)';
omega = f * 2 * pi;
jomega = omega * 1i;

% amplifier transfer functions
A_1_s = A_O_1 ./ (1 + jomega./omega_lp_amp1);
A_1_s_dB = 20 * log10(abs(A_1_s));
A_2_s = A_O_2 ./ (1 + jomega./omega_lp_amp2);
A_2_s_dB = 20 * log10(abs(A_2_s));

close all
figure
semilogx(f, A_1_s_dB);
hold on
semilogx(f, A_2_s_dB);
xlabel('f (Hz)');
ylabel('|A(j\omega)| (dB)');
legend('A_1', 'A_2');
title('Amplifier Frequency Response');
axis([1e-3 1e10 50 130]);

% loop gain for first stage
LG_1_s = A_1_s .* C_a ./ C_in .* (jomega + omega_hp1) ./ (jomega + omega_in1);
LG_1_s_dB = 20 * log10(abs(LG_1_s));
% forward gain for first stage
G_1_s = A_1_s ./ C_in ./ (jomega + omega_in1);
G_1_s_dB = 20 * log10(abs(G_1_s));
% feedback factor for first stage
FB_1_s = C_a .* (jomega + omega_hp1);
FB_1_s_dB = 20 * log10(abs(FB_1_s));
% inverse feedback factor - ideal closed-loop gain if large loop gain
H_CL_1_s = 1 ./ FB_1_s;
H_CL_1_s_dB = 20 * log10(abs(H_CL_1_s));
% actual closed-loop gain
A_CL_1_s = G_1_s ./ (1 - G_1_s .* FB_1_s);
A_CL_1_s_dB = 20 * log10(abs(A_CL_1_s));

% note
% LG is unitless
% forward gain, inverse fb, and actual CL gain are current to voltage
% feedback factor is voltage to current

figure
semilogx(f, LG_1_s_dB);
hold on
semilogx(f, G_1_s_dB);
semilogx(f, FB_1_s_dB);
semilogx(f, H_CL_1_s_dB);
semilogx(f, A_CL_1_s_dB);
xlabel('f (Hz)');
ylabel('|T(j\omega)| (dB)');
legend('Loop Gain', 'Forward Gain', 'FB Factor',...
    'Closed-Loop Gain, LG Large', 'Actual Closed-Loop Gain');
title('Loop Gain and Closed-Loop Gain, Stage 1');

% transfer function from various voltage sources to V_int
% computed using ideal closed-loop gain and actual closed-loop gain
% for each, VI_signal_stage_s is the voltage-to-current conversion
% for that signal. approx_TF_signal_stage_s is the approx overall TF
% for that signal and stage, assuming large loop gain, and
% actual_TF_signal_stage_s is the actual overall TF for that signal, stage

% VIN(s-s_0) - should be correct (i.e. 40dB) around 3.84MHz
VI_VIN_1_s = C_in .* (jomega - 1i*omega_0);
VI_VIN_1_s_dB = 20 * log10(abs(VI_VIN_1_s));
approx_TF_VIN_1_s = VI_VIN_1_s .* H_CL_1_s;
approx_TF_VIN_1_s_dB = 20 * log10(abs(approx_TF_VIN_1_s));
actual_TF_VIN_1_s = VI_VIN_1_s .* A_CL_1_s;
actual_TF_VIN_1_s_dB = 20 * log10(abs(actual_TF_VIN_1_s));

figure
semilogx(f, H_CL_1_s_dB);
hold on
semilogx(f, A_CL_1_s_dB);
semilogx(f, VI_VIN_1_s_dB);
semilogx(f, approx_TF_VIN_1_s_dB);
semilogx(f, actual_TF_VIN_1_s_dB);
xlabel('f (Hz)');
ylabel('|T(j\omega)| (dB)');
title('V_{IN}(s-s_0) -> V_{INT}(s) Transfer Function Computation');
legend('Approx Closed-Loop ItoV Gain', 'Actual Closed-Loop ItoV Gain',...
    'V_{IN} VtoI Gain', 'Approx VtoV Gain', 'Actual VtoV Gain');