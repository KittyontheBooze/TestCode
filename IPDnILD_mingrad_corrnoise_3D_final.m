clc; clear;
close all;

s_azim = 0;
s_elev = 0;

n_azim = 0;
n_elev = 0;

ITD_threshold = 3;
corr_threshold = 0.90;

corr_onoff = 1;
%% Setting

fs = 48000;
dt = 1/fs;

nfbin = 4096;
f = ((1:nfbin/2-1)*(fs/nfbin)).';

cut_lo = 100;
cut_up = 5000;
[ind_lo, ind_up] = find_edge(f, cut_lo, cut_up);

azims = (0:5:355).';
N_azims = length(azims);
elevs = (-40:5:90).';
N_elevs = length(elevs);

fid = fopen(strcat('hrtf_expdb_nfbin=', num2str(nfbin), '.mat'),'r');
if fid == -1
    root_hrtf = 'HRTF_HATS';
    hrtf_db = load_expDB(root_hrtf,nfbin);
    save(strcat('hrtf_expdb_nfbin=', num2str(nfbin), '.mat'), 'hrtf_db')
else
    load(strcat('hrtf_expdb_nfbin=', num2str(nfbin), '.mat'));
    fclose(fid);
end

b = zeros(nfbin/2-1,1);
for i = 1:(nfbin/2-1)
    if f(i) < 750
        b(i) = 1;
    elseif 750 <= f(i) && f(i) < 1500
        b(i) = -1/750*(f(i)-1500);
    else
        b(i) = 0;
    end
end

%% Filter Definition

fid = fopen('lpf_16k.mat','r');
if fid == -1
    lpf_16k = designfilt('lowpassfir', 'PassbandFrequency', 1.6e3, 'StopbandFrequency', 1.85e3, 'PassbandRipple', 0.2, 'StopbandAttenuation', 60, 'SampleRate', 48000);
    save('lpf_16k.mat', 'lpf_16k')
else
    load('lpf_16k.mat');
    fclose(fid);
end

%% Criteria Loading

mode = 'minimum';
gamma = 1;

if strcmp('minimum',mode)
    fid = fopen(strcat('mag_crit_nfbin=', num2str(nfbin), '_min_gamma=',num2str(gamma),'.mat'),'r');
    if fid == -1
        root_hrtf = 'HRTF_HATS';
        [mag_crit, pha_crit] = load_crits(root_hrtf,mode,gamma,nfbin);
        save(strcat('mag_crit_nfbin=', num2str(nfbin), '_min_gamma=',num2str(gamma),'.mat'), 'mag_crit');
        save(strcat('pha_crit_nfbin=', num2str(nfbin), '_min_gamma=',num2str(gamma),'.mat'), 'pha_crit');
    else
        load(strcat('mag_crit_nfbin=', num2str(nfbin), '_min_gamma=',num2str(gamma),'.mat'));
        fclose(fid);
    end
    fid = fopen(strcat('pha_crit_nfbin=', num2str(nfbin), '_min_gamma=',num2str(gamma),'.mat'),'r');
    if fid == -1
        root_hrtf = 'HRTF_HATS';
        [mag_crit, pha_crit] = load_crits(root_hrtf,mode,gamma,nfbin);
        save(strcat('mag_crit_nfbin=', num2str(nfbin), '_min_gamma=',num2str(gamma),'.mat'), 'mag_crit');
        save(strcat('pha_crit_nfbin=', num2str(nfbin), '_min_gamma=',num2str(gamma),'.mat'), 'pha_crit');
    else
        load(strcat('pha_crit_nfbin=', num2str(nfbin), '_min_gamma=',num2str(gamma),'.mat'));
        fclose(fid);
    end
elseif strcmp('maximum',mode)
    fid = fopen(strcat('mag_crit_nfbin=', num2str(nfbin), '_max_gamma=',num2str(gamma),'.mat'),'r');
    if fid == -1
        root_hrtf = 'HRTF_HATS';
        [mag_crit, pha_crit] = load_crits(root_hrtf,mode,gamma,nfbin);
        save(strcat('mag_crit_nfbin=', num2str(nfbin), '_max_gamma=',num2str(gamma),'.mat'), 'mag_crit');
        save(strcat('pha_crit_nfbin=', num2str(nfbin), '_max_gamma=',num2str(gamma),'.mat'), 'pha_crit');
    else
        load(strcat('mag_crit_nfbin=', num2str(nfbin), '_max_gamma=',num2str(gamma),'.mat'));
        fclose(fid);
    end
    fid = fopen(strcat('pha_crit_nfbin=', num2str(nfbin), '_max_gamma=',num2str(gamma),'.mat'),'r');
    if fid == -1
        root_hrtf = 'HRTF_HATS';
        [mag_crit, pha_crit] = load_crits(root_hrtf,mode,gamma,nfbin);
        save(strcat('mag_crit_nfbin=', num2str(nfbin), '_max_gamma=',num2str(gamma),'.mat'), 'mag_crit');
        save(strcat('pha_crit_nfbin=', num2str(nfbin), '_max_gamma=',num2str(gamma),'.mat'), 'pha_crit');
    else
        load(strcat('pha_crit_nfbin=', num2str(nfbin), '_max_gamma=',num2str(gamma),'.mat'));
        fclose(fid);
    end
else
    error('Improper mode selection');
end

%% Source load
root_source = 'sound_source';
source_name = sprintf('man_01');
source = load_wav(root_source, strcat(source_name, '.wav'));

%% Noise Loading
noise_name = sprintf('man_02');
noise = load_wav(root_source, strcat(noise_name, '.wav'));

%% Source
N = fs*5;
x = source(1:N);

sound(x,48000);

%% Noise Length Adjusting
x_n = noise(1:N);

%% Volume Normalizing
insp_length = 1024;
num_insp = fix((N-insp_length)/ceil(insp_length*0.5))+1;

p_s_mx = 0;    % Extacting instantaneous power of source signal
for i = 1:num_insp
    start = (i-1)*ceil(insp_length*0.5) + 1;
    p_s = rms(x(start:start+insp_length-1));
    if p_s > p_s_mx
        p_s_mx = p_s;
    end
end

p_n_mx = 0;    % Extacting instantaneous power of noise signal
for i = 1:num_insp
    start = (i-1)*ceil(insp_length*0.5) + 1;
    p_n = rms(x_n(start:start+insp_length-1));
    if p_n > p_n_mx
        p_n_mx = p_n;
    end
end

x_n = (p_s_mx/p_n_mx).*x_n;

%% Source Generation

s_azim_index = s_azim/5 + 1;
s_elev_index = s_elev/5 + 9;

x_sL = conv(hrtf_db(1).azim(s_azim_index).elev(s_elev_index).hrtf_L, x);
x_sL = x_sL(1:N);
x_sR = conv(hrtf_db(1).azim(s_azim_index).elev(s_elev_index).hrtf_R, x);
x_sR = x_sR(1:N);

%% Noise Generation

n_azim_index = n_azim/5 + 1;
n_elev_index = n_elev/5 + 9;

x_nL = conv(hrtf_db(1).azim(n_azim_index).elev(n_elev_index).hrtf_L, x_n);
x_nL = x_nL(1:N);
x_nR = conv(hrtf_db(1).azim(n_azim_index).elev(n_elev_index).hrtf_R, x_n);
x_nR = x_nR(1:N);

%% Mixing

x_L = x_sL + x_nL;
x_R = x_sR + x_nR;

%% Test Start

num_test = fix((N-nfbin)/ceil(0.5*nfbin))+1;
t = (0:num_test-1)*(ceil(0.5*nfbin)/48000).';

cor_lr_mx = zeros(num_test,1);
lag_lr_mx = zeros(num_test,1);
infer_rawdata = zeros(N_azims, N_elevs, num_test);
infer_result_raw = zeros(num_test,2);
tic;
for i = 1 : num_test
    fprintf('%d\n', i)
    start = (i-1)*ceil(0.5*nfbin) + 1;%seg_ind = 1:num_test
    xs_L = x_L(start:start+nfbin-1);
    xs_R = x_R(start:start+nfbin-1);
    
    y_L = filter(lpf_16k, xs_L);
    y_R = filter(lpf_16k, xs_R);
    
    [corrs,lags] = xcorr(y_L, y_R, 35+ITD_threshold, 'normalized');
    [mx, ind_corr] = max(corrs);
    cor_lr_mx(i) = mx;
    lag_lr_mx(i) = lags(ind_corr);
    %
    if corr_onoff && mx < corr_threshold
        fprintf('mx<%.2f\n', corr_threshold)
        infer_rawdata(:,:,i) = NaN;
        infer_result_raw(i,:) = NaN;
        continue;
    end
    %
    X_L = fft(xs_L);
    X_R = fft(xs_R);
    X_LR = X_L(1:nfbin/2-1).*conj(X_R(1:nfbin/2-1));
    
    x_IPD = unwrap(angle(X_LR));
    x_ILD = 20*log10(abs(X_L(1:nfbin/2-1))./abs(X_R(1:nfbin/2-1)));
    
    max_ELV = zeros(N_azims, N_elevs);
    for j = 1 : N_azims
        for k = 1 : N_elevs
            if abs(hrtf_db(1).azim(j).elev(k).ITD_raw - lags(ind_corr)) <= ITD_threshold
                h_IPD = hrtf_db(1).azim(j).elev(k).IPD(1:nfbin/2-1);
                h_ILD = hrtf_db(1).azim(j).elev(k).ILD(1:nfbin/2-1);
                
                ELV = 0;
                for l = ind_lo:ind_up
                    if abs(h_ILD(l) - x_ILD(l)) < mag_crit(j, k, l)
                        D_m = 1;
                    else
                        D_m = 0;
                    end
                    
                    if abs(h_IPD(l) - x_IPD(l)) < pha_crit(j, k, l)
                        D_p = 1;
                    else
                        D_p = 0;
                    end
                    
                    D = D_p*b(l) + D_m*(1-b(l));
                    ELV = ELV + (abs(X_L(l))+abs(X_R(l)))/2*D;
                end
                infer_rawdata(j,k,i) = ELV;
            else
                infer_rawdata(j,k,i) = NaN;
            end
        end
    end
    isnan_ELV = isnan(infer_rawdata(:,:,i));
    if ~all(all(isnan_ELV==1))
        max_ELV = max(max(infer_rawdata(:,:,i)));
        [infer_result_raw(i,1), infer_result_raw(i,2)] = find(infer_rawdata(:,:,i)==max_ELV);
    else
        fprintf('singular\n')
        infer_result_raw(i,:) = NaN;
    end
end
time_elapse = toc;
infer_result(:,1) = (infer_result_raw(:,1) - 1)*5;
infer_result(:,2) = (infer_result_raw(:,2) - 9)*5;

source_location = [s_azim s_elev];
source_location_rad = source_location*pi/180;
source_location_cart = [cos(source_location_rad(1))*cos(source_location_rad(2)), sin(source_location_rad(1))*cos(source_location_rad(2)), sin(source_location_rad(2))];
infer_result_rad = infer_result*pi()/180;
infer_location_cart = [cos(infer_result_rad(:,1)).*cos(infer_result_rad(:,2)), sin(infer_result_rad(:,1)).*cos(infer_result_rad(:,2)), sin(infer_result_rad(:,2))];
inner = sum(source_location_cart.*infer_location_cart,2);

accuracy_tf = sum(all(infer_result==source_location,2))/sum(all(~isnan(infer_result),2));
accuracy_cos10 = sum(inner>cos(10*pi/180))/sum(~isnan(inner));
accuracy_cos15 = sum(inner>cos(15*pi/180))/sum(~isnan(inner));

isnan_infer = zeros(1,2);
infer_result_plot = zeros(num_test,2);
for i = 1:num_test
    isnan_infer = isnan(infer_result(i,:));
   if all(all(isnan_infer==1))
      if i == 1
          infer_result_plot(i,:) = 0;
      else
          infer_result_plot(i,:) = infer_result_plot(i-1,:);
      end
   else
       infer_result_plot(i,:) = infer_result(i,:);
   end
end

all_data_table = [accuracy_tf accuracy_cos10 accuracy_cos15 time_elapse];

figure
scatter(infer_result_plot(:,1), infer_result_plot(:,2), 80,'LineWidth',4); %,'ro-'
set(gca, 'FontSize',17)
xlim([0 355])
ylim([-40 90])
xlabel('Azimuth Angle [deg]','fontsize',21); ylabel('Elevation Angle [deg]','fontsize',21);
title({"Sound Source Localization in 3D Space", "[Environmental Noise, Source Location : ("+num2str(s_azim)+" ,"+num2str(s_elev)+")]"},'fontsize',32);

set(gcf, 'Position', [0 0 1024 768])
saveas(gcf, "plots/Sound Source Localization in 3D Space [Environmental Noise, Source Location _ ("+num2str(s_azim)+" ,"+num2str(s_elev)+")].png")
saveas(gcf, "plots/fig files/Sound Source Localization in 3D Space [Environmental Noise, Source Location _ ("+num2str(s_azim)+" ,"+num2str(s_elev)+")].fig")

% GIF 내보내기
h1 = figure;
% This ensures that getframe() returns a consistent size.
axis tight manual
filename = strcat("plots/Sound Source Localization in 3D Space [Environmental Noise, Source Location _ ("+num2str(s_azim)+" ,"+num2str(s_elev)+")].gif");

for n = 1:num_test
    scatter(infer_result_plot(n,1), infer_result_plot(n,2), 80,'LineWidth',4); %,'ro-'
    set(gca, 'FontSize',17)
    xlim([0 355])
    ylim([-40 90])
    xlabel('Azimuth Angle [deg]','fontsize',21); ylabel('Elevation Angle [deg]','fontsize',21);
    title({"Sound Source Localization in 3D Space", "[Environmental Noise, Source Location : ("+num2str(s_azim)+" ,"+num2str(s_elev)+")]"},'fontsize',32);
    
    set(gcf, 'Position', [0 0 1024 768])
    
    % Capture the plot as an image
    frame = getframe(h1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append');
    end
end
%GIF 내보내기 끝