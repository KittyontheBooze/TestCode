%    lpf_16k = designfilt('lowpassfir', 'PassbandFrequency', 1.6e3, 'StopbandFrequency', 2.1e3, 'PassbandRipple', 0.01, 'StopbandAttenuation', 80, 'SampleRate', 48000);
%    save('lpf_16k.mat', 'lpf_16k')

fid = fopen('lpf_16k.mat','r');
if fid == -1
    lpf_16k = designfilt('lowpassfir', 'PassbandFrequency', 1.6e3, 'StopbandFrequency', 1.85e3, 'PassbandRipple', 0.2, 'StopbandAttenuation', 60, 'SampleRate', 48000);
    save('lpf_16k.mat', 'lpf_16k')
else
    load('lpf_16k.mat');
    fclose(fid);
end

x_s = load_wav(root_source, 'man_01.wav');
x_s = x_s(1:N);
x_n = load_wav(root_source, 'man_02.wav');
x_n = x_n(1:N);
x_sL = filter(lpf_16k, x_s);
x_nL = filter(lpf_16k, x_n);

figure
spectrogram(x_sL,kaiser(nfbin,7), ceil(0.95*nfbin), nfbin, fs, 'yaxis')

figure
spectrogram(x_nL,kaiser(nfbin,7), ceil(0.95*nfbin), nfbin, fs, 'yaxis')

%fvtool(lpf_16k, 'Fs', 48000)

%freqz(lpf_16k)
%{
figure
plot(lpf_16k.Coefficients)

isminphase(lpf_16k)

figure
a = 1;
b = lpf_16k.Coefficients;
%freqz(b,a,'whole')

%min_phase_filter(b,a)

%
%zplane(b,a)

fs = 48000;
dt = 1/fs;

ITD_threshold = 3;

nfbin = 4096;
f = ((1:nfbin/2-1)*(fs/nfbin)).';

fid = fopen(strcat('hrtf_expdb_nfbin=', num2str(nfbin), '.mat'),'r');
if fid == -1
    root_hrtf = 'HRTF_HATS';
    hrtf_db = load_expDB(root_hrtf,nfbin);
    save(strcat('hrtf_expdb_nfbin=', num2str(nfbin), '.mat'), 'hrtf_db')
else
    load(strcat('hrtf_expdb_nfbin=', num2str(nfbin), '.mat'));
    fclose(fid);
end

n_azim = 270;
n_elev = 0;

root_source = 'sound_source';
x_s = load_wav(root_source, 'man_01.wav');
x_n = load_wav(root_source, 'man_02.wav');

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
xcorr_lpf = zeros(72, num_test);
lag_lpf = zeros(size(xcorr_wolpf));
for x = 1:72
    for y = 9
        %% Source Generation
        source_location = [(x-1)*5, (y-9)*5];
        fprintf('(%d, %d)\n', source_location(1), source_location(2));
        
        x_sL = conv(hrtf_db(1).azim(x).elev(y).hrtf_L, x_s);
        x_sL = x_sL(1:N);
        x_sR = conv(hrtf_db(1).azim(x).elev(y).hrtf_R, x_s);
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
            
            y_L = filter(b, a, x_L);
            y_R = filter(b, a, x_R);
            
            [corrs_temp,lags_temp] = xcorr(y_L(250:end), y_R(250:end), 35+ITD_threshold, 'normalized');
            [mx_temp, ind_temp] = max(corrs_temp);
            lag_temp = lags_temp(ind_temp);
            xcorr_lpf(x,i) = mx_temp;
            lag_lpf(x,i) = lag_temp;
        end
    end
end
%{
load("angdist_tensor_case1_(270,0).mat",'angdist_tensor')

angdist_tensor_deg = real(acos(angdist_tensor)/pi*180);
ER20_tensor = angdist_tensor_deg;
for i = 1:size(ER20_tensor,1)
    for j = 1:size(ER20_tensor, 2)
        for k = 1:size(ER20_tensor,3)
            if ER20_tensor(i,j,k) <= 20
                ER20_tensor(i,j,k) = 1;
            else
                ER20_tensor(i,j,k) = 0;
            end
        end
    end
end
ER20_tensor = permute(ER20_tensor,[3, 2, 1]);

accuracy_ER20 = sum(ER20_tensor(1,9,:))/numel(ER20_tensor(1,9,:));
mask_xcorr_lpf = xcorr_lpf;
for i = numel(xcorr_lpf)
    if mask_xcorr_lpf <=0.75
        mask_xcorr_lpf(i) = 1;
    else
        mask_xcorr_lpf(i) = 0;
    end
end
acc_ER20_xcorr_lpf = sum(ER20_tensor(:,9,1).*mask_xcorr_lpf)/sum(mask_xcorr_lpf);
%}
figure
imagesc(xcorr_lpf)
%}