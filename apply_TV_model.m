function [TV_est, TV_est_t] = apply_TV_model(EDRsignals, fs, b, fc, AEDRmethod, plotflag)
% APPLY_TV_model  Apply a trained model for TV estimation from EDR signals.
%                 Described in J. Lazaro et al., "Tracking Tidal Volume
%                 from Holter and Wearable Armband Electrocardiogram
%                 Monitoring", IEEE J Biomed Health Inform, 2024,
%                 DOI: 10.1109/JBHI.2024.3383232
%
% Created by Jesus Lazaro <jlazarop@unizar.es> in 2024
%--------
%   Sintax: [TV_est, TV_est_t] = apply_TV_model(EDRsignals, fs, b, fc, AEDRmethod, plotflag)
%   In:   EDRsignals = matrix with ECG-derived respiration signals in columns
%         fs = sampling rate of EDRsignals (Hz)
%         b = coeficients of the trained model
%         fc = cutoff frequency for lowpass filtering (Hz) [Default: 0.05]
%         AEDRmethod = 'peak' for peak-to-peak AEDR method (default)
%                      'ana' for analytic AEDR method
%                      'rms' for RMS AEDR method
%         plotflag = if 1, plots a figure with PPG and SSF [Default: 0]
%
%   Out:  TV_est = Estimated tidal volem
%         TV_est_t = time vector for TV_est

    if nargin<3
        error('EDR signals, sampling rate, and coefficients of the trained model need to be provided');
    end
    
    if nargin<3
        fc=0.05;
    end

    if nargin<4
        AEDRmethod = 'peak';
    end
    
    if nargin<5
        plotflag = false;
    end

    if ~strcmp(AEDRmethod, 'peak') && ~strcmp(AEDRmethod, 'ana') && ~strcmp(AEDRmethod, 'rms')
        error('AEDRmethod must be either ''peak'', ''ana'', or ''rms''');
    end
    
    N_EDR = size(EDRsignals, 2);

    %% Compute AEDR signals:
    for k_EDR=1:N_EDR
        eval(['[AEDR' num2str(k_EDR) AEDRmethod ', AEDR' num2str(k_EDR) AEDRmethod '_t] = compute_AEDR_' AEDRmethod '(EDRsignals(:,' num2str(k_EDR) '), fs, fc);']); %amplitude-estimation method
    end
    

    %% Generate X matrix:
    if strcmp(AEDRmethod, 'peak') %peak-to-peak model
        peak_t_ini = -inf;
        peak_t_end = inf;
        for k_EDR=1:N_EDR
            eval(['peak_t_ini = max(peak_t_ini, AEDR' num2str(k_EDR) 'peak_t(1));'])
            eval(['peak_t_end = min(peak_t_end, AEDR' num2str(k_EDR) 'peak_t(end));'])
        end
        TV_est_t = peak_t_ini:1/fs:peak_t_end;
        X = [ones(length(TV_est_t), 1), nan(length(TV_est_t), N_EDR)];
        for k_EDR=1:N_EDR
            eval(['aux_ind = AEDR' num2str(k_EDR) 'peak_t>=peak_t_ini & AEDR' num2str(k_EDR) 'peak_t<=peak_t_end;']);
            eval(['X(:, k_EDR+1) = AEDR' num2str(k_EDR) 'peak(aux_ind).'';']);
        end
    else % non peak-to-peak models
        eval(['TV_est_t = AEDR1' AEDRmethod '_t']);
        X = [ones(length(TV_est_t), 1), nan(length(TV_est_t), N_EDR)];
        for k_EDR=1:N_EDR
            eval(['X(:, k_EDR+1) = AEDR' num2str(k_EDR) AEDRmethod '.'';']);
        end
    end
    

    %% Apply model:
    TV_est = b.'*X.';
    

    %% Figure:
    if plotflag       
        EDR_t = 0:1/fs:(length(EDRsignals)-1)/fs;

        figure;
        ax(1) = subplot(2,1,1);
        plot(EDR_t, EDRsignals);
        ylabel('EDR signals (EDR untis)');
        ax(2) = subplot(2,1,2); hold on;
        plot(TV_est_t, TV_est, 'r');
        xlabel('Time (s)');
        ylabel('TV estimated (TV units)');
        linkaxes(ax, 'x');
    end



end