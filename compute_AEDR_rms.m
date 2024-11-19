function [AEDR_sig, AEDR_t] = compute_AEDR_rms(EDR_sig, fs, fc, plotflag)
%COMPUTE_AEDR_ANA   Compute RMS amplitude-estimation method
%                   described in J. Lazaro et al., "Tracking Tidal Volume
%                   from Holter and Wearable Armband Electrocardiogram
%                   Monitoring", IEEE J Biomed Health Inform, 2024,
%                   DOI: 10.1109/JBHI.2024.3383232
%
% Created by Jesus Lazaro <jlazarop@unizar.es> in 2024
%--------
%   Sintax: [AEDR] = compute_AEDR_peak(EDR, fs, fc)
%   In:   EDR_sig = ECG-derived respiration signal
%         fs = sampling rate (Hz)
%         fc = cutoff frequency for lowpass filtering (Hz) [Default: 0.05]
%         plotflag = if 1, plots a figure with PPG and SSF [Default: 0]
%
%   Out:  AEDR_signal = AEDR signal
%         AEDR_t = time vector for AEDR signal

    if nargin<2
        error('EDR signal and sampling rate need to be provided');
    end
    
    if nargin<3
        fc=0.05;
    end
    
    if nargin<4
        plotflag = false;
    end

    AEDR_t = 1:length(EDR_sig);
    
    [Me, me] = envelope(EDR_sig, round(10*fs), 'rms');
    
    AEDR_sig = Me - me;
    AEDR_t = (AEDR_t-1)/fs;

    if fc~=inf
        if plotflag
            AEDR_nofilt = AEDR_sig;
        end
        
        [bb, aa] = butter(3, fc*2/fs, 'low');
        AEDR_sig = filtfilt(bb, aa, AEDR_sig);
    end
    
    if plotflag
        EDR_t = 0:1/fs:(length(EDR_sig)-1)/fs;
        figure;
        ax(1) = subplot(2,1,1); hold on;
        plot(EDR_t, EDR_sig, 'k');
        plot(EDR_t, me, 'b');
        plot(EDR_t, Me, 'r');
        legend({'EDR signal', 'Upper envelope', 'Lower envelope'});
        ax(2) = subplot(2,1,2); hold on;
        plot(AEDR_t, AEDR_nofilt, 'b');
        plot(AEDR_t, AEDR_sig, 'k');
        legend({'AEDR unfiltered', 'AEDR filtered'});
        xlabel('Time (s)');
        linkaxes(ax, 'x');
    end
    
end