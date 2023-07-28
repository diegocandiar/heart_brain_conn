function struct_conn = compute_connectivity(data_eeg, pks_indx)

%% Load HRV
time = data_eeg.time{1};
t_heartbeats = time(pks_indx);
IBI = diff(t_heartbeats);
t_IBI = time(pks_indx(1:length(IBI)));
wind = 15;
[CSI_out, CVI_out, t_out] = compute_CSI_CVI(IBI, t_IBI, wind);

%% frequency analysis
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.pad          = 'nextpow2';
cfg.foi          = 0:0.5:45;                         
cfg.t_ftimwin    = ones(length(cfg.foi),1).*2;  
cfg.toi          = '50%';  % 2s 50p,                  
freq = ft_freqanalysis(cfg, data_eeg);

time2 = freq.time;
time3 = time2(2:length(time2)-2);
t1 = max([time3(1) t_out(1)]);
t2 = min([time3(end) t_out(end)]);
time3 = time3(time3 >=t1 & time3 <=t2);

cfg = [];
%cfg.frequency = [1 45];
cfg.latency = [time3(1) time3(end)];
freq = ft_selectdata(cfg, freq);

CSI = interp1(t_out,CSI_out,time3,'spline');
CVI = interp1(t_out,CVI_out,time3,'spline');

%% Separate TF
freq_epoch = freq; freq_epoch = rmfield(freq_epoch,{'powspctrm','freq'});
freq_delta = freq_epoch;
freq_theta = freq_epoch;
freq_alpha = freq_epoch;
freq_beta = freq_epoch;
freq_gamma = freq_epoch;

f1a = 1; f1b = 4; f1 = freq.freq >= f1a & freq.freq <= f1b;
f2a = 4; f2b = 8; f2 = freq.freq >= f2a & freq.freq <= f2b;
f3a = 8; f3b = 12; f3 = freq.freq >= f3a & freq.freq <= f3b;
f4a = 12; f4b = 30; f4 = freq.freq >= f4a & freq.freq <= f4b;
f5a = 30; f5b = 45; f5 = freq.freq >= f5a & freq.freq <= f5b;

freq_delta.trial{1} = squeeze(trapz(freq.powspctrm(:,f1,:),2));  
freq_theta.trial{1} = squeeze(trapz(freq.powspctrm(:,f2,:),2));    
freq_alpha.trial{1} = squeeze(trapz(freq.powspctrm(:,f3,:),2));  
freq_beta.trial{1} = squeeze(trapz(freq.powspctrm(:,f4,:),2));   
freq_gamma.trial{1} = squeeze(trapz(freq.powspctrm(:,f5,:),2));
freq_delta.dimord = 'chan_time';
freq_theta.dimord = 'chan_time';
freq_alpha.dimord = 'chan_time';
freq_beta.dimord = 'chan_time';
freq_gamma.dimord = 'chan_time';

%% Prepare inputs
win_RR = 15; % seconds
time = time3;
timearx = time(floor(win_RR/2) +1 : end-ceil(win_RR/2));

%% run ARX 

Fs = 1;

% delta
EEG_comp = [freq_delta.trial{1} ; CSI ; CVI]; 
[conn_delta, mrkv_delta] = arx_model(EEG_comp, Fs, timearx, wind); 

%theta
EEG_comp = [freq_theta.trial{1} ; CSI ; CVI]; 
[conn_theta, mrkv_theta] = arx_model(EEG_comp, Fs, timearx, wind); 

%alpha
EEG_comp = [freq_alpha.trial{1} ; CSI ; CVI];  
[conn_alpha, mrkv_alpha] = arx_model(EEG_comp, Fs, timearx, wind); 

%beta
EEG_comp = [freq_beta.trial{1} ; CSI ; CVI];  
[conn_beta, mrkv_beta] = arx_model(EEG_comp, Fs, timearx, wind); 

%gamma
EEG_comp = [freq_gamma.trial{1} ; CSI ; CVI]; 
[conn_gamma, mrkv_gamma] = arx_model(EEG_comp, Fs, timearx, wind); 

% conn(ch1, ch2, t): ch1 <- ch2

%% check outputs
struct_conn = struct;

struct_conn.time = time;
struct_conn.timearx = timearx;

struct_conn.CSI = CSI;
struct_conn.CVI = CVI;

struct_conn.freq_delta = freq_delta;
struct_conn.freq_theta = freq_theta;
struct_conn.freq_alpha = freq_alpha;
struct_conn.freq_beta = freq_beta;
struct_conn.freq_gamma = freq_gamma;

struct_conn.conn_delta = conn_delta;
struct_conn.conn_theta = conn_theta;
struct_conn.conn_alpha = conn_alpha;
struct_conn.conn_beta = conn_beta;
struct_conn.conn_gamma = conn_gamma;

struct_conn.mrkv_delta = mrkv_delta;
struct_conn.mrkv_theta = mrkv_theta;
struct_conn.mrkv_alpha = mrkv_alpha;
struct_conn.mrkv_beta = mrkv_beta;
struct_conn.mrkv_gamma = mrkv_gamma;

end