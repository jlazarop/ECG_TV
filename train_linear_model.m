function [b, X] = train_linear_model(AEDR_signals, TVref)
%TRAIN_LINEAR_MODEL Train linear model for TV estimation as described in
%                   J. Lazaro et al., "Tracking Tidal Volume from Holter
%                   and Wearable Armband Electrocardiogram Monitoring"
%                   IEEE J Biomed Health Inform, 2024,
%                   DOI: 10.1109/JBHI.2024.3383232
%
% Created by Jesus Lazaro <jlazarop@unizar.es> in 2024
%--------
%   Sintax: [b] = train_linear_model(EDR_signals)
%   In:   AEDR_signals = AEDR signals in columns
%         TVref = Reference TV
%
%   Out:  b = Coefficients of the obtained linear model
%         X = Matrix containing a first column of ones followed by
%             the input AEDR_signals. The estimation of TV from the
%             training data can be obtained by b.'*X.';

if nargin<2
    error('Not enough input arguments');
end

N_EDR = size(AEDR_signals, 2); % Number of EDR signals;

lb = [-inf zeros(1,N_EDR)];
ub = ones(1,N_EDR+1);
Aeq = ones(1,N_EDR+1);
beq = 1;

X = [ones(size(AEDR_signals, 1), 1), AEDR_signals];

b = lsqlin(X,TVref,[],[],[],[],lb,[]);