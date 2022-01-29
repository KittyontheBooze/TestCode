%% RIR Loading
[h, fs] = audioread('rt60_81_L.wav');
%{
fid = fopen('hrtf_db.mat','r');
if fid == -1
    root_hrtf = 'HRTF_HATS';
    hrtf_db = load_DB(root_hrtf);
    save('hrtf_db.mat', 'hrtf_db')
else
    load('hrtf_db.mat');
    fclose(fid);
end

h = hrtf_db(1).azim(1).elev(9).hrtf_L;
fs = 48000;
%}
%{
airpar.fs = 48e3;
airpar.rir_type = 1;
airpar.room = 2;
airpar.channel = 1;
airpar.head = 0;
airpar.rir_no = 1;
airpar.azimuth = 0;

[h,air_info] = load_air(airpar);
fs = airpar.fs;
%}
N = length(h);
dt = 1/fs;
t = (1:N)/fs;

%% Energy Decay Curve
h_sqflip = flip(h.^2);

edc_temp = zeros(N,1);
cumm = 0;
for i = 1:N
    cumm = cumm + h_sqflip(i);
    edc_temp(i) = cumm;
end
edc = 10*log10(flip(edc_temp)./edc_temp(end));

figure
plot(t,edc,'LineWidth', 1.5)
set(gca,'FontSize',19)
colormap jet
xlabel('Time [sec]','fontsize',27); ylabel('Energy Decay [db]','fontsize',27);
title({"Energy Decay Curve from Aachen RIR data", "(Meeting Room, RT60 = 0.81s)"},'fontsize',32);
grid on%RT60은 이 플롯의 Free decay 구간의 Slope의 역수 * 60

set(gcf, 'Position', [0 0 1024 512])
saveas(gcf, strcat('plots/Energy Decay Curve from Aachen RIR data (Meeting Room, RT60 = 0.81s).png'))
saveas(gcf, strcat('plots/fig files/Energy Decay Curve from Aachen RIR data (Meeting Room, RT60 = 0.81s).fig'))
%
%% Averaged Slope

avg_len = 2^13;
num_bin = fix((N-avg_len)/ceil(0.05*avg_len))+1;
avgslope = zeros(num_bin,1);
for j = 1 : num_bin
    start = (j-1)*ceil(0.05*avg_len) + 1;
    frame = edc(start:start+avg_len-1);
    avgslope(j) = (frame(end)-frame(1))./(dt*avg_len);
end

RT60_instant = 60./avgslope;
t_RT60 = (0:num_bin-1)*(ceil(0.05*avg_len)*dt).';

figure
yyaxis left
plot(t,edc); hold on
yyaxis right
plot(t_RT60,RT60_instant)
%

