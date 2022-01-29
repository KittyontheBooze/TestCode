fs = 48000;
dt = 1/fs;

nfbin = 4096;
f = ((1:nfbin/2-1)*(fs/nfbin)).';

cut_lo = 100;
cut_up = 5000;
[ind_lo, ind_up] = find_edge(f, cut_lo, cut_up);

azim = 0:5:355;
elev = -40:5:90;
N_azims = length(azim);
N_elevs = length(elev);

s_azim = 0;
s_elev = 0;
s_azim_index = s_azim/5 + 1;
s_elev_index = s_elev/5 + 9;

n_azim = 270;
n_elev = 0;
n_azim_index = n_azim/5 + 1;
n_elev_index = n_elev/5 + 9;

label_tot = [];
%% Source load
root_source = 'sound_source';
source_name = "speech_"+num2str(num_file,'%02d');
source = load_wav(root_source, strcat(source_name, '.wav'));

%% Noise Loading
noise_name = sprintf('dog_01');
noise = load_wav(root_source, strcat(noise_name, '.wav'));

%% Source
N = length(source);

x_s = source;
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

x_n = (p_s_mx/p_n_mx).*x_n_temp;

%% Absolute spectrum inspection
num_test = fix((N-nfbin)/ceil(0.5*nfbin))+1;
t = (0:num_test-1)*(ceil(0.5*nfbin)/48000).';

label = zeros(num_test,1);
for i = 1 : num_test
    %fprintf('%d\n', i)
    start = (i-1)*ceil(0.5*nfbin) + 1;%seg_ind = 1:num_test
    x_si = x_s(start:start+nfbin-1);
    x_ni = x_n(start:start+nfbin-1);
    
    X_si = fft(x_si);
    X_ni = fft(x_ni);
    
    S_si = X_si.*conj(X_si);
    S_ni = X_ni.*conj(X_ni);
    
    P_si = sum(S_si);
    P_ni = sum(S_ni);
    
    if P_si >= P_ni
        label(i) = 1;
    else
        label(i) = 0;
    end
end
label_tot = [label_tot;label];

temp = load('experimental results/original_corrnoise_azim_infer_rawdump_allexp_nfbin_4096.mat');
infer_rawdump = cell2mat(table2array(cell2table(struct2cell(temp))));

infer_max = zeros(size(infer_rawdump));
ind_table_tot = zeros(N_azims,size(infer_rawdump,3));
tf_table_tot = zeros(size(ind_table_tot));
for i = 1:size(infer_rawdump,3)
    for j = 1:N_azims
        [maximum, ind_temp] = max(infer_rawdump(j,:,i));
        ind_table_tot(j,i) = ind_temp;
        infer_max(j,ind_temp,i) = maximum;
        
        if label_tot(i) == 1 %해당 Frame의 답이 Speech의 위치일 경우
            if ind_temp == j
                tf_table_tot(j,i)=1;
            else
                tf_table_tot(j,i)=0;
            end
        else %해당 Frame의 답이 Dog Barking의 위치일 경우
            if ind_temp == n_azim_index
                tf_table_tot(j,i)=1;
            else
                tf_table_tot(j,i)=0;
            end
        end
    end
end
%
position_table = (ind_table_tot - 1) * 5;
source_azimuth = zeros(size(tf_table_tot));
source_elevation = s_elev;
for i = 1:N_azims
    for j = 1:size(infer_rawdump,3)
        if label_tot(j) == 1 %해당 Frame의 답이 Speech의 위치일 경우
            source_azimuth(i,j) = (i - 1) * 5;
        else %해당 Frame의 답이 Dog Barking의 위치일 경우
            source_azimuth(i,j) = n_azim;
        end
    end
end

position_table_rad = position_table*pi/180;
source_azimuth_rad = source_azimuth*pi/180;
source_elevation_rad = source_elevation*pi/180;

cos_table_tot  = zeros(size(tf_table_tot));
for i = 1:N_azims
    for j = 1:size(infer_rawdump,3)
        source_location_cart = [cos(source_azimuth_rad(i,j))*cos(source_elevation_rad), sin(source_azimuth_rad(i,j))*cos(source_elevation_rad), sin(source_elevation_rad)];
        infer_location_cart = [cos(position_table_rad(i,j))*cos(source_elevation_rad), sin(position_table_rad(i,j))*cos(source_elevation_rad), sin(source_elevation_rad)];
        cos_table_tot(i,j) = sum(source_location_cart.*infer_location_cart);
    end
end

accuracy_cos10 = sum(sum(cos_table_tot>cos(10*pi/180)))/numel(cos_table_tot);
accuracy_cos15 = sum(sum(cos_table_tot>cos(15*pi/180)))/numel(cos_table_tot);

accuracy_total = sum(sum(tf_table_tot))/numel(tf_table_tot);
accuracy_azims = sum(tf_table_tot,2)/size(tf_table_tot, 2);

save('original_corrnoise_azim_tf_table_tot.mat', 'tf_table_tot');
save('original_corrnoise_azim_cos_table_tot.mat', 'cos_table_tot');

