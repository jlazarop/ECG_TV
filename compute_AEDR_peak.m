function [AEDR_sig, AEDR_t] = compute_AEDR_peak(EDR_sig, fs, fc, plotflag)
%COMPUTE_AEDR_PEAK  Compute peak-to-peak amplitude-estimation method
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
    
    if nargin<3
        plotflag = false;
    end
    
    if nargin<4
        plotflag = false;
    end
    
    if size(EDR_sig,1)>size(EDR_sig,2)
        EDR_sig = EDR_sig.';
    end
    
    [Mpeaks, Mlocs] = findpeaks(EDR_sig, 'MinPeakDistance', round(1*fs)-1);
    [mpeaks, mlocs] = findpeaks(-EDR_sig, 'MinPeakDistance', round(1*fs)-1);
    mpeaks = -mpeaks;
    
    %remove 2 consecutive mins
    Mlocs = [1, Mlocs, length(EDR_sig)]; 
    mask = true(size(mlocs));
    for kk=1:(length(Mlocs)-1)
        aux = mlocs>Mlocs(kk) & mlocs<Mlocs(kk+1);
        if(nnz(aux)>1)
            aux2 = double(aux); aux2(aux2==0)=inf;
            [~, auxmax] = min(mpeaks+aux2);
            aux(auxmax) = false;
            mask = mask & ~aux;
        end
    end
    Mlocs = Mlocs(2:(end-1));
    mlocs = mlocs(mask);
   
    %remove 2 consecutive maxs:
    mlocs = [1, mlocs, length(EDR_sig)]; 
    mask = true(size(Mlocs));
    for kk=1:(length(mlocs)-1)
        aux = Mlocs>mlocs(kk) & Mlocs<mlocs(kk+1);
        if(nnz(aux)>1)
            aux2 = double(aux); aux2(aux2==0)=-inf;
            [~, auxmax] = max(Mpeaks+aux2);
            aux(auxmax) = false;
            mask = mask & ~aux;
        end
    end
    mlocs = mlocs(2:(end-1));
    Mlocs=Mlocs(mask);

    mpeaks = mpeaks(mlocs>Mlocs(1));
    mlocs = mlocs(mlocs>Mlocs(1));

    aux_AEDR_sig = Mpeaks(1:length(mpeaks)) - mpeaks;
    aux_t = (Mlocs(1:length(mpeaks))+mlocs)/2;
    aux_t = (aux_t-1)/fs;
    
    AEDR_t = 0:1/fs:aux_t(end);
    AEDR_t = AEDR_t(AEDR_t>=aux_t(1));
    
    AEDR_sig = spline(aux_t, aux_AEDR_sig, AEDR_t);
%     AEDR_t = (AEDR_t-1)/fs;
    
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
        plot((mlocs-1)/fs, EDR_sig(mlocs), 'b*');
        plot((Mlocs-1)/fs, EDR_sig(Mlocs), 'r*');
        legend({'EDR signal', 'Min', 'Max'});
        ax(2) = subplot(2,1,2); hold on;
        plot(AEDR_t, AEDR_nofilt, 'b');
        plot(aux_t, aux_AEDR_sig, 'b*');
        plot(AEDR_t, AEDR_sig, 'k');
        legend({'AEDR unfiltered', 'AEDR filtered'});
        xlabel('Time (s)');
        linkaxes(ax, 'x');
    end


end