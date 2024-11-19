% This file contains an example of use of the codes for applying the
% algorithm for tidal volume (TV) estimation from ECG-derived respiration
% (EDR) signals described in J. Lazaro et al., "Tracking Tidal Volume from
% Holter and Wearable Armband Electrocardiogram Monitoring" IEEE J Biomed
% Health Inform, 2024, DOI: 10.1109/JBHI.2024.3383232
%
% First, artificial EDR data is synthesized. Subsequently, the amplitude of
% those EDR signals (AEDR signals) is computed by the functions
% "compute_AEDR_*.m". Then, the models are trained with the function
% "train_linear_model.m".
%
% In addition, new EDR data is synthesized in order to apply the
% previoussly trained models with the function "apply_TV_model".
%
% Created by Jesus Lazaro <jlazarop@unizar.es> in 2024


clear; close all;


%% Algorithm parameters:
N_EDR = 3; %Number of EDR signals to synthesize
fs = 4; %Sampling rate of the EDR signals
fresp = 0.3; %Respiratory frequency of the synthesized EDR signals


%% Algorithm parameters:
fc = 0.05;


%% Synthesize EDR signals:
EDR_t = (0:1/fs:60).';
phase = 2*pi*rand(1,N_EDR);
EDR_sig = 0.5*(cos(fresp*2*pi*EDR_t + phase));

a=0.1; b=30; c=15; %Gaussian parameters for TV simulation
TV = 0.1+a*exp(-(EDR_t-b).^2/(c^2)); % TV simulation

EDR_sig = EDR_sig.*repmat(TV, 1, 3); %Add TV influence


%% Compute AEDR signals:
for k_EDR=1:N_EDR
    eval(['[AEDR' num2str(k_EDR) 'peak, AEDR' num2str(k_EDR) 'peak_t] = compute_AEDR_peak(EDR_sig(:,' num2str(k_EDR) '), fs, fc);']); %Peak-to-peak amplitude-estimation method
    eval(['[AEDR' num2str(k_EDR) 'ana, AEDR' num2str(k_EDR) 'ana_t] = compute_AEDR_ana(EDR_sig(:,' num2str(k_EDR) '), fs, fc);']); %Analytic amplitude-estimation method
    eval(['[AEDR' num2str(k_EDR) 'rms, AEDR' num2str(k_EDR) 'rms_t] = compute_AEDR_rms(EDR_sig(:,' num2str(k_EDR) '), fs, fc);']); %RMS amplitude-estimation method
end


%% Train TV model based on peak-to-peak AEDR signals:
peak_t_ini = -inf;
peak_t_end = inf;
for k_EDR=1:N_EDR
    eval(['peak_t_ini = max(peak_t_ini, AEDR' num2str(k_EDR) 'peak_t(1));'])
    eval(['peak_t_end = min(peak_t_end, AEDR' num2str(k_EDR) 'peak_t(end));'])
end
AEDRpeak_t = peak_t_ini:1/fs:peak_t_end;
AEDRsignals = nan(length(AEDRpeak_t), N_EDR);
for k_EDR=1:N_EDR
    eval(['aux_ind = AEDR' num2str(k_EDR) 'peak_t>=peak_t_ini & AEDR' num2str(k_EDR) 'peak_t<=peak_t_end;']);
    eval(['AEDRsignals(:, k_EDR) = AEDR' num2str(k_EDR) 'peak(aux_ind).'';']);
end

[bpeak, Xpeak] = train_linear_model(AEDRsignals, TV(EDR_t>=peak_t_ini & EDR_t<=peak_t_end));

TV_est_peak = bpeak.'*Xpeak.';


%% Train TV model based on analytic AEDR signals:
AEDRsignals = [nan(length(TV), N_EDR)];
for k_EDR=1:N_EDR
    eval(['AEDRsignals(:, k_EDR) = AEDR' num2str(k_EDR) 'ana.'';']);
end

[bana, Xana] = train_linear_model(AEDRsignals, TV);

TV_est_ana = bana.'*Xana.';


%% Train TV model based on RMS AEDR signals:
AEDRsignals = [nan(length(TV), N_EDR)];
for k_EDR=1:N_EDR
    eval(['AEDRsignals(:, k_EDR) = AEDR' num2str(k_EDR) 'rms.'';']);
end


[brms, Xrms] = train_linear_model(AEDRsignals, TV);

TV_est_rms = brms.'*Xrms.';



%% Figure:
figure; hold on;
ax(1) = subplot(2,1,1);
plot(EDR_t, EDR_sig);
ylabel('EDR signals (EDR units)');
ax(2) = subplot(2,1,2); hold on;
plot(EDR_t, TV, 'k');
plot(AEDRpeak_t, TV_est_peak, 'm')
plot(EDR_t, TV_est_ana, 'b');
plot(EDR_t, TV_est_rms, 'r');
xlabel('Time (s)');
ylabel('(TV units)');
legend({'TV reference', 'TV estimated (peak-to-peak)', 'TV estimated (analytic)', 'TV estimated (RMS)'});
linkaxes(ax, 'x');


%% How to estimate TV from new data with the trained models:
% Synthetise new EDR data:
newEDR_t = (0:1/fs:60).';
phase = 2*pi*rand(1,N_EDR);
newEDR_sig = 0.5*(cos(fresp*2*pi*newEDR_t + phase));

a=0.2; b=15; c=20; %Gaussian parameters for TV simulation
newTV = 0.1+a*exp(-(newEDR_t-b).^2/(c^2)); % TV simulation

newEDR_sig = newEDR_sig.*repmat(newTV, 1, 3); %Add TV influence

% Apply desired model, e.g., analytic model (coefs in 'bana'):
plotflag = true;
[newTV_est, newTV_est_t] = apply_TV_model(newEDR_sig, fs, bana, fc, 'ana', plotflag);