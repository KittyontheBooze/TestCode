%clc; clear;
close all;

%{

n_azim = 270;
n_elev = 0;

ITD_threshold = 3;
corr_threshold = 0.75;

corr_onoff = 0;
SNR = 0;
%
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

beta = zeros(nfbin/2-1,1);
for i = 1:(nfbin/2-1)
    if f(i) < 750
        beta(i) = 1;
    elseif 750 <= f(i) && f(i) < 1500
        beta(i) = -1/750*(f(i)-1500);
    else
        beta(i) = 0;
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

a_ipd = zeros(nfbin/2-1,1);
a_ild = zeros(nfbin/2-1,1);
fid = fopen(strcat('a_ipd_nfbin=', num2str(nfbin),'.mat'),'r');
if fid == -1
    root_hrtf = 'HRTF_HATS';
    [a_ipd, a_ild] = load_crits(root_hrtf,nfbin);
    save(strcat('a_ipd_nfbin=', num2str(nfbin),'.mat'), 'a_ipd');
    save(strcat('a_ild_nfbin=', num2str(nfbin),'.mat'), 'a_ild');
else
    load(strcat('a_ipd_nfbin=', num2str(nfbin),'.mat'));
    fclose(fid);
end
fid = fopen(strcat('a_ild_nfbin=', num2str(nfbin),'.mat'),'r');
if fid == -1
    root_hrtf = 'HRTF_HATS';
    [a_ipd, a_ild] = load_crits(root_hrtf,nfbin);
    save(strcat('a_ipd_nfbin=', num2str(nfbin),'.mat'), 'a_ipd');
    save(strcat('a_ild_nfbin=', num2str(nfbin),'.mat'), 'a_ild');
else
    load(strcat('a_ild_nfbin=', num2str(nfbin),'.mat'));
    fclose(fid);
end

for exp_case = 3
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
    
    TF_tensor = zeros(N_azims, N_elevs, num_test);
    ER15_tensor = zeros(size(TF_tensor));
    MaxELV_tensor = zeros(size(TF_tensor));
    angdist_tensor = zeros(size(TF_tensor));
    label_tensor = zeros(size(TF_tensor));
    tempSNR_tensor = zeros(size(TF_tensor));
    P_s_tensor = zeros(size(TF_tensor));
    P_n_tensor = zeros(size(TF_tensor));
    infer_azim = zeros(size(TF_tensor));
    infer_elev = zeros(size(TF_tensor));
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
                %% Amplitude Comparison
                X_L1 = fft(x_L1);
                X_R1 = fft(x_R1);
                X_L2 = fft(x_L2);
                X_R2 = fft(x_R2);
                
                S_L1 = X_L1.*conj(X_L1);
                S_R1 = X_R1.*conj(X_R1);
                S_L2 = X_L2.*conj(X_L2);
                S_R2 = X_R2.*conj(X_R2);
                
                P_L1 = sum(S_L1);
                P_R1 = sum(S_R1);
                P_L2 = sum(S_L2);
                P_R2 = sum(S_R2);
                
                P_s = (P_L1+P_R1)/2;
                P_s_tensor(a, b, i) = P_s;
                P_n = (P_L2+P_R2)/2;
                P_n_tensor(a, b, i) = P_n;
                
                if P_s >= P_n
                    label = 1;
                    label_tensor(a, b, i) = label;
                    tempSNR_tensor(a, b, i) = 20*log(P_s/P_n);
                else
                    label = 0;
                    label_tensor(a, b, i) = label;
                    tempSNR_tensor(a, b, i) = 20*log(P_n/P_s);
                end
                %% Mixing
                
                x_L = x_L1 + x_L2;
                x_R = x_R1 + x_R2;
                
                X_L = fft(x_L);
                X_R = fft(x_R);
                X_LR = X_L(1:nfbin/2-1).*conj(X_R(1:nfbin/2-1));
                
                x_IPD = unwrap(angle(X_LR));
                x_ILD = 20*log10(abs(X_L(1:nfbin/2-1))./abs(X_R(1:nfbin/2-1)));
                
                ELV_map = zeros(N_azims, N_elevs);
                for j = 1 : N_azims
                    for k = 1 : N_elevs
                        h_IPD = hrtf_db(1).azim(j).elev(k).IPD(1:nfbin/2-1);
                        h_ILD = hrtf_db(1).azim(j).elev(k).ILD(1:nfbin/2-1);
                        
                        ELV = 0;
                        for l = ind_lo:ind_up
                            if abs(h_ILD(l) - x_ILD(l)) < a_ild(l)
                                D_m = 1;
                            else
                                D_m = 0;
                            end
                            
                            if abs(h_IPD(l) - x_IPD(l)) < a_ipd(l)
                                D_p = 1;
                            else
                                D_p = 0;
                            end
                            
                            D = D_p*beta(l) + D_m*(1-beta(l));
                            ELV = ELV + (abs(X_L(l))+abs(X_R(l)))/2*D;
                        end
                        ELV_map(j,k) = ELV;
                    end
                end
                
                MaxELV = max(max(ELV_map));
                [temp1, temp2] = find(ELV_map==MaxELV);
                MaxELV_tensor(a, b, i) = MaxELV;
                %% Infer Result Choice
                if length(temp1) == 1  || length(temp2)  == 1
                    infer_ind1 = temp1;
                    infer_ind2 = temp2;
                else % When the solution of predition has multiple solution
                    if i ~= 1
                        prev_ind1 = infer_azim(a, b, i-1);
                        prev_ind2 = infer_elev(a, b, i-1);
                    else
                        if label == 1
                            prev_ind1 = a;
                            prev_ind2 = b;
                        else
                            prev_ind1 = n_azim_index;
                            prev_ind2 = n_elev_index;
                        end
                    end
                    prev_loc = [(prev_ind1-1)*5, (prev_ind2-9)*5];
                    prev_loc_rad = prev_loc*pi/180;
                    prev_loc_cart = [cos(prev_loc_rad(1))*cos(prev_loc_rad(2)), sin(prev_loc_rad(1))*cos(prev_loc_rad(2)), sin(prev_loc_rad(2))];
                    
                    infer_loc = [(temp1-1)*5, (temp2-9)*5];
                    infer_loc_rad = infer_loc.*pi/180;
                    infer_loc_cart = [cos(infer_loc_rad(:,1)).*cos(infer_loc_rad(:,2)), sin(infer_loc_rad(:,1)).*cos(infer_loc_rad(:,2)), sin(infer_loc_rad(:,2))];
                    
                    inner = sum(prev_loc_cart.*infer_loc_cart,2);
                    
                    [~, temp_ind] = max(inner);
                    
                    infer_ind1 = temp1(temp_ind);
                    infer_ind2 = temp2(temp_ind);
                end
                
                infer_azim(a, b, i) = infer_ind1;
                infer_elev(a, b, i) = infer_ind2;
                
                if infer_ind1 == a && infer_ind2 == b
                    TF_tensor(a, b, i) = 1;
                else
                    TF_tensor(a, b, i) = 0;
                end
                
                infer_loc = [(infer_ind1-1)*5, (infer_ind2-9)*5];
                infer_loc_rad = infer_loc*pi()/180;
                infer_loc_cart = [cos(infer_loc_rad(1)).*cos(infer_loc_rad(2)), sin(infer_loc_rad(1)).*cos(infer_loc_rad(2)), sin(infer_loc_rad(2))];
                
                if label == 1
                    source_location_rad = source_location*pi/180;
                else
                    source_location_rad = noise_location*pi/180;
                end
                source_location_cart = [cos(source_location_rad(1))*cos(source_location_rad(2)), sin(source_location_rad(1))*cos(source_location_rad(2)), sin(source_location_rad(2))];
                inner = sum(source_location_cart.*infer_loc_cart,2);
                
                angdist_tensor(a, b, i) = inner;
                
                if inner >= cos(15*pi/180)
                    ER15_tensor(a, b, i) = 1;
                else
                    ER15_tensor(a, b, i) = 0;
                end
                
            end
        end
    end
    
    save("TF_tensor_case"+num2str(exp_case)+".mat",'TF_tensor')
    save("ER15_tensor_case"+num2str(exp_case)+".mat",'ER15_tensor')
    save("MaxELV_tensor_case"+num2str(exp_case)+".mat",'MaxELV_tensor')
    save("angdist_tensor_case"+num2str(exp_case)+".mat",'angdist_tensor')
    save("label_tensor_case"+num2str(exp_case)+".mat",'label_tensor')
    save("P_s_tensor_case"+num2str(exp_case)+".mat",'P_s_tensor')
    save("P_n_tensor_case"+num2str(exp_case)+".mat",'P_n_tensor')
    save("tempSNR_tensor_case"+num2str(exp_case)+".mat",'tempSNR_tensor')
    save("infer_azim_case"+num2str(exp_case)+".mat",'infer_azim')
    save("infer_elev_case"+num2str(exp_case)+".mat",'infer_elev')
end
%}
%
TF_tensor_plot = permute(TF_tensor,[1, 3, 2]);
ER15_tensor_plot = permute(ER15_tensor,[1, 3, 2]);
angdist_tensor_plot = permute(angdist_tensor,[1, 3, 2]);
label_tensor_plot = permute(label_tensor,[1, 3, 2]);
P_s_tensor_plot = permute(P_s_tensor,[1, 3, 2]);
tempSNR_tensor_plot = permute(tempSNR_tensor,[1, 3, 2]);
xcorr_lpf_plot = permute(xcorr_lpf,[1, 3, 2]);

figure
imagesc(TF_tensor_plot(:,:,9));
colormap gray

figure
imagesc(ER15_tensor_plot(:,:,9));
colormap gray

figure
imagesc(angdist_tensor_plot(:,:,9), [cos(30*pi/180), 1]);
colormap gray

figure
imagesc(label_tensor_plot(:,:,9));
colormap gray

figure
imagesc(P_s_tensor_plot(:,:,9));
colormap gray

figure
imagesc(tempSNR_tensor_plot(:,:,9));
colormap gray

figure
imagesc(tempSNR_tensor_plot(:,:,9)>40);
colormap gray

figure
imagesc(xcorr_lpf_plot(:,:,9), [0.75 1]);
colormap gray
%{
len_SNRdata = sum(sum(sum(~isinf(tempSNR_tensor))));
SNRdata = zeros(len_SNRdata,1);
ind = 1;
for d = 1:numel(tempSNR_tensor)
    if ~isinf(tempSNR_tensor(d))
        SNRdata(ind) = tempSNR_tensor(d);
        ind = ind + 1;
    end
end

[nlogL, AVAR] = gamlike(gamfit(SNRdata), SNRdata)

acc_ER15_elev_0 = sum(sum(ER15_tensor_plot(:,:,9)))/numel(ER15_tensor_plot(:,:,9));

temp = tempSNR_tensor_plot(:,:,9);
%}

%{
figure
imagesc(ax,ax,infer_rawdump(:,:,k).', [0 c_thres]); axis xy
set(gca,'FontSize',19)
colormap jet
xlabel('Source Angle [deg.]','fontsize',24); ylabel('HRTF Angle [deg.]','fontsize',24);
if m == 1
    title({"ELV Map for "+string+"Estimation (Proposed)", "(Noiseless, t = "+num2str(t(k),'%0.2f')+"s)"},'fontsize',27);
elseif m == 2
    title({"ELV Map for "+string+"Estimation (Proposed)", "(Environmental Noise, Dog Barking at (270,0)), t = "+num2str(t(k),'%0.2f')+"s"},'fontsize',27);
else
    title({"ELV Map for "+string+"Estimation (Proposed)", "(Reverberant, RT60 = 0."+num2str(rt60(rt))+"s), t = "+num2str(t(k),'%0.2f')+"s"},'fontsize',27);
end
grid on

set(gcf, 'Position', [0 0 1024 1024])
%}