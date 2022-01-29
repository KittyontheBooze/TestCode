function hrtf_db = load_expDB(root_hrtf, nfbin)

if nargin == 1, nfbin = 512; end
if nfbin < 512, error('filter length must be equal to or larger than 512'); end

fs = 48000;
dt = 1/fs;
upsample_coef = 1;

fid = fopen('lpf_16k.mat','r');
if fid == -1
    lpf_16k = designfilt('lowpassfir', 'PassbandFrequency', 1.6e3, 'StopbandFrequency', 1.85e3, 'PassbandRipple', 0.2, 'StopbandAttenuation', 60, 'SampleRate', 48000);
    save('lpf_16k.mat', 'lpf_16k')
else
    load('lpf_16k.mat');
    fclose(fid);
end

azims = linspace(0,355,72);
elevs = linspace(-40,90,27);

hrtf_db.azim.elev.hrtf_L = zeros(512,1);
hrtf_db.azim.elev.hrtf_R = zeros(512,1);
hrtf_db.azim.elev.ITD = 0;
hrtf_db.azim.elev.ITD_raw = 0;
hrtf_db.azim.elev.IPD = zeros(nfbin,1);
hrtf_db.azim.elev.ILD = zeros(nfbin,1);

for i = 1:length(azims)
    for j = 1:length(elevs)
        if elevs(j) >= 0
            path_name = sprintf('%s/a%03d/a%03de+%02d.txt', root_hrtf,azims(i),azims(i),elevs(j));
        else
            path_name = sprintf('%s/a%03d/a%03de-%02d.txt', root_hrtf,azims(i),azims(i),abs(elevs(j)));
        end
        
        fid = fopen(path_name,'r');
        if fid == -1, error('cannot open file : %s',path_name); end
        
        temp = (fscanf(fid,'%f %f',[2 inf])).';
        hrtf_db(1).azim(i).elev(j).hrtf_L = temp(:,1);
        hrtf_db(1).azim(i).elev(j).hrtf_R = temp(:,2);
        
        temp_L = filter(lpf_16k, temp(:,1));
        temp_R = filter(lpf_16k, temp(:,2));
        
        x_L = resample(temp_L, upsample_coef, 1);
        x_R = resample(temp_R, upsample_coef, 1);
        
        [corrs,lags] = xcorr(x_L,x_R);
        [~, ITD_index] = max(corrs);
        hrtf_db(1).azim(i).elev(j).ITD_raw = lags(ITD_index);
        hrtf_db(1).azim(i).elev(j).ITD = lags(ITD_index)*(dt/upsample_coef);
        
        h_L = [temp(:,1); zeros(nfbin - 512,1)];
        h_R = [temp(:,2); zeros(nfbin - 512,1)];
        
        H_L = fft(h_L);
        H_R = fft(h_R);
        
        cross_spec = H_L.*conj(H_R);
        hrtf_db(1).azim(i).elev(j).IPD = unwrap(angle(cross_spec));
        hrtf_db(1).azim(i).elev(j).ILD = 20*log10(abs(H_L)./abs(H_R));
        
        fclose(fid);
        
        j = j + 1;
    end
    i = i + 1;
end

end