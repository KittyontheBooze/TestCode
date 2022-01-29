%clc; clear;
%close all;
%{
%% Setting

fs = 48000;
dt = 1/fs;

nfbin = 4096;
f = ((1:nfbin/2-1)*(fs/nfbin)).';

n_azim = 270;
n_elev = 0;

%% Filter Definition

fid = fopen('lpf_16k.mat','r');
if fid == -1
    lpf_16k = designfilt('lowpassfir', 'PassbandFrequency', 1.6e3, 'StopbandFrequency', 1.85e3, 'PassbandRipple', 0.2, 'StopbandAttenuation', 60, 'SampleRate', 48000);
    save('lpf_16k.mat', 'lpf_16k')
else
    load('lpf_16k.mat');
    fclose(fid);
end
for exp_case = 1:3
    % Sound load
    
    root_source = 'sound_source';
    
    switch exp_case
        case 1
            x_s = load_wav(root_source, 'man_01.wav');
            x_n = load_wav(root_source, 'man_02.wav');
        case 2
            x_s = load_wav(root_source, 'woman_01.wav');
            x_n = load_wav(root_source, 'woman_02.wav');
        case 3
            x_s = load_wav(root_source, 'man_01.wav');
            x_n = load_wav(root_source, 'woman_02.wav');
        case 4
            x_s = load_wav(root_source, 'woman_01.wav');
            x_n = load_wav(root_source, 'man_02.wav');
    end
    
    N = fs*5;
    x_s = x_s(1:N);
    x_n = x_n(1:N);
    
    %% Noise Generation
    noise_location = [n_azim, n_elev];
    n_azim_index = n_azim/5 + 1;
    n_elev_index = n_elev/5 + 9;
    
    x_nL = conv(hrtf_db(1).azim(n_azim_index).elev(n_elev_index).hrtf_L, x_n);
    x_nL = x_nL(1:N);
    x_nR = conv(hrtf_db(1).azim(n_azim_index).elev(n_elev_index).hrtf_R, x_n);
    x_nR = x_nR(1:N);
    
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
    
    save("xcorr_wolpf_case"+num2str(exp_case)+".mat",'xcorr_wolpf')
    save("lag_wolpf_case"+num2str(exp_case)+".mat",'lag_wolpf')
    save("xcorr_lpf_case"+num2str(exp_case)+".mat",'xcorr_lpf')
    save("lag_lpf_case"+num2str(exp_case)+".mat",'lag_lpf')
    save("avgcoh_case"+num2str(exp_case)+".mat",'avgcoh')
    
end

%{
xcorr_wolpf_plot = permute(xcorr_wolpf,[1, 3, 2]);
xcorr_lpf_plot = permute(xcorr_lpf,[1, 3, 2]);
avgcoh_plot = permute(avgcoh,[1, 3, 2]);

for l = 1:27
figure
imagesc(avgcoh_plot(:,:,l), [0 1]);
end
%}