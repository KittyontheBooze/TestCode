%clc; clear;
%close all;

s_azim = 180;
s_elev = 0;

n_azim = 270;
n_elev = 0;

ITD_threshold = 3;

corr_onoff = 0;
corr_threshold = 0.85;
ITD_onoff = 0;

coh_onoff = 1;
coh_threshold = 0.8;

SNR = 0;

%% Source load
root_source = 'sound_source';
source_name = sprintf('speech_01');
source = load_wav(root_source, strcat(source_name, '.wav'));

%% Noise Loading
noise_name = sprintf('dog_01');
noise = load_wav(root_source, strcat(noise_name, '.wav'));

%% Source
x_s = source;
N = length(x_s);

sound(x_s,48000);

%% Noise Length Adjusting
x_n_temp = noise(fs*1.5:fs*1.5+N-1);

%% Volume Normalizing
insp_length = 1024;
num_insp = fix((N-insp_length)/ceil(insp_length*0.5))+1;

p_s_mx = 0;    % Extacting instantaneous power of source signal
for i = 1:num_insp
    start = (i-1)*ceil(insp_length*0.5) + 1;
    p_s = rms(x_s(start:start+insp_length-1));
    if p_s > p_s_mx
        p_s_mx = p_s;
    end
end

p_n_mx = 0;    % Extacting instantaneous power of noise signal
for i = 1:num_insp
    start = (i-1)*ceil(insp_length*0.5) + 1;
    p_n = rms(x_n_temp(start:start+insp_length-1));
    if p_n > p_n_mx
        p_n_mx = p_n;
    end
end

x_n = (p_s_mx/p_n_mx)*exp(-SNR/20).*x_n_temp;

num_test = fix((N-nfbin)/ceil(0.5*nfbin))+1;
%t = (0:num_test-1)*(ceil(0.5*nfbin)/48000).';

%% Source Generation

s_azim_index = s_azim/5 + 1;
s_elev_index = s_elev/5 + 9;

x_sL = conv(hrtf_db(1).azim(s_azim_index).elev(s_elev_index).hrtf_L, x_s);
x_sL = x_sL(1:N);
x_sR = conv(hrtf_db(1).azim(s_azim_index).elev(s_elev_index).hrtf_R, x_s);
x_sR = x_sR(1:N);

%% Noise Generation

n_azim_index = n_azim/5 + 1;
n_elev_index = n_elev/5 + 9;

x_nL = conv(hrtf_db(1).azim(n_azim_index).elev(n_elev_index).hrtf_L, x_n);
x_nL = x_nL(1:N);
x_nR = conv(hrtf_db(1).azim(n_azim_index).elev(n_elev_index).hrtf_R, x_n);
x_nR = x_nR(1:N);
%}
%% Making Label
label = zeros(num_test,1);
SNR_temporal = zeros(size(label));
for i = 1 : num_test
    %fprintf('%d\n', i)
    start = (i-1)*ceil(0.5*nfbin) + 1;%seg_ind = 1:num_test
    x_sLi = x_sL(start:start+nfbin-1);
    x_sRi = x_sR(start:start+nfbin-1);
    x_nLi = x_nL(start:start+nfbin-1);
    x_nRi = x_nR(start:start+nfbin-1);
    
    X_sLi = fft(x_sLi);
    X_sRi = fft(x_sRi);
    X_nLi = fft(x_nLi);
    X_nRi = fft(x_nRi);
    
    S_sLi = X_sLi.*conj(X_sLi);
    S_sRi = X_sRi.*conj(X_sRi);
    S_nLi = X_nLi.*conj(X_nLi);
    S_nRi = X_nRi.*conj(X_nRi);
    
    P_sLi = sum(S_sLi);
    P_sRi = sum(S_sRi);
    P_nLi = sum(S_nLi);
    P_nRi = sum(S_nRi);
    
    P_s = (P_sLi+P_sRi)/2;
    P_n = (P_nLi+P_nRi)/2;
    
    if P_s >= P_n
        label(i) = 1;
        SNR_temporal(i) = 20*log(P_s/P_n);
    else
        label(i) = 0;
        SNR_temporal(i) = 20*log(P_n/P_s);
    end
end
comparison = [label SNR_temporal];