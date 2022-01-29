%clc; clear;
%close all;
%
%% Setting

fs = 48000;
dt = 1/fs;

nfbin = 4096;
f = ((1:nfbin/2-1)*(fs/nfbin)).';

cut_lo = 100;
cut_up = 5000;
[ind_lo, ind_up] = find_edge(f, cut_lo, cut_up);

%% Filter Definition

fid = fopen('lpf_16k.mat','r');
if fid == -1
    lpf_16k = designfilt('lowpassfir', 'PassbandFrequency', 1.6e3, 'StopbandFrequency', 1.85e3, 'PassbandRipple', 0.2, 'StopbandAttenuation', 60, 'SampleRate', 48000);
    save('lpf_16k.mat', 'lpf_16k')
else
    load('lpf_16k.mat');
    fclose(fid);
end

%% Source
N = 94208;
rng(1)
x_s = 0.01*randn(N,1);
rng(2)
x_nL = 0.01*randn(N,1);
rng(3)
x_nR = 0.01*randn(N,1);


%% Test over All Azimuth and Elevation
num_test = fix((N-nfbin)/ceil(0.5*nfbin))+1;
t = (0:num_test-1)*(ceil(0.5*nfbin)/48000).';

%% Averaging range of Coherence Estimation

[cxy, f_cxy] = mscohere(x_nL(1:nfbin), x_nR(1:nfbin), 1024, 0.5*1024, [], 48000);
[ind_low, ind_high] = find_edge(f, 500, 1600);

%% Test Start
xcorr_wolpf = zeros(N_azims, N_elevs, num_test);
lag_wolpf = zeros(size(xcorr_wolpf));
xcorr_lpf = zeros(size(xcorr_wolpf));
lag_lpf = zeros(size(xcorr_wolpf));
avgcoh = zeros(size(xcorr_wolpf));
for a = 1:N_azims
    for b = 1:N_elevs   
        %% Source Generation
        source_location = [(a-1)*5, (b-9)*5];
        fprintf('(%d, %d)\n', source_location(1), source_location(2));
        
        x_sL = conv(hrtf_db(1).azim(a).elev(b).hrtf_L, x_s);
        x_sL = x_sL(1:N);
        x_sR = conv(hrtf_db(1).azim(a).elev(b).hrtf_R, x_s);
        x_sR = x_sR(1:N);
        
        %% Test Start
        
        for i = 1 : num_test
            fprintf('%d\n', i)
            start = (i-1)*ceil(0.5*nfbin) + 1;%seg_ind = 1:num_test
            x_L1 = x_sL(start:start+nfbin-1);
            x_R1 = x_sR(start:start+nfbin-1);
            
            x_L2 = x_nL(start:start+nfbin-1);
            x_R2 = x_nR(start:start+nfbin-1);
            %% Noise Generation
            p_s = (rms(x_sL)+rms(x_sR))/2;
            p_nL = rms(x_nL);
            p_nR = rms(x_nR);
            
            amp = p_s*exp(-(45-i)/20);
            x_L2 = amp/p_nL*x_L2;
            x_R2 = amp/p_nR*x_R2;
            %% Mixing
            
            x_L = x_L1 + x_L2;
            x_R = x_R1 + x_R2;
            
            y_L = filter(lpf_16k, x_L);
            y_R = filter(lpf_16k, x_R);
            
            [corrs_temp,lags_temp] = xcorr(x_L, x_R, 35+ITD_threshold, 'normalized');
            [mx_temp, ind_temp] = max(corrs_temp);
            lag_temp = lags_temp(ind_temp);
            xcorr_wolpf(a, b, i) = mx_temp;
            lag_wolpf(a, b, i) = lag_temp;
            
            [corrs_temp,lags_temp] = xcorr(y_L, y_R, 35+ITD_threshold, 'normalized');
            [mx_temp, ind_temp] = max(corrs_temp);
            lag_temp = lags_temp(ind_temp);
            xcorr_lpf(a, b, i) = mx_temp;
            lag_lpf(a, b, i) = lag_temp;
            
            [cxy, f_cxy] = mscohere(x_L, x_R, 1024, 0.5*1024, [], 48000);
            avgcoh(a, b, i) = mean(cxy(ind_low:ind_high));
        end
    end
end
%}
%{
save('xcorr_wolpf.mat','xcorr_wolpf')
save('lag_wolpf.mat','lag_wolpf')
save('xcorr_lpf.mat','xcorr_lpf')
save('lag_lpf.mat','lag_lpf')
save('avgcoh.mat','avgcoh')
%}
%{
xcorr_wolpf_plot = permute(xcorr_wolpf,[1, 3, 2]);
xcorr_lpf_plot = permute(xcorr_lpf,[1, 3, 2]);
avgcoh_plot = permute(avgcoh,[1, 3, 2]);
%}
for l = 1:27
figure
imagesc(avgcoh_plot(:,:,l), [0 1]);
end

