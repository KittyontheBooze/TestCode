clc; clear;
close all;

nfbin = 4096;
fs = 48000;

fid = fopen(strcat('hrtf_expdb_nfbin=',num2str(nfbin),'.mat'),'r');
if fid == -1
    root_hrtf = 'HRTF_HATS';
    hrtf_db = load_expDB(root_hrtf,nfbin);
    save(strcat('hrtf_expdb_nfbin=',num2str(nfbin),'.mat'), 'hrtf_db')
else
    load(strcat('hrtf_expdb_nfbin=',num2str(nfbin),'.mat'));
    fclose(fid);
end

f = (1:nfbin)*(fs/nfbin);

itd = zeros(27,72);

[azim, elev] = meshgrid(0:5:355,-40:5:90);

for i = 1:72
    for j = 1:27
        itd(j,i) = hrtf_db(1).azim(i).elev(j).ITD_raw;
    end
end

itd_in_s = itd/fs*1000;

figure
plot(azim, itd_in_s(9,:),'LineWidth', 1.5); hold on;
plot(azim, itd_in_s(12,:),'LineWidth', 1.5)
plot(azim, itd_in_s(15,:),'LineWidth', 1.5)
plot(azim, itd_in_s(18,:),'LineWidth', 1.5)
plot(azim, itd_in_s(21,:),'LineWidth', 1.5)
plot(azim, itd_in_s(24,:),'LineWidth', 1.5)
plot(azim, itd_in_s(27,:),'LineWidth', 1.5)
legend('ITD at elev=0°','ITD at elev=15°','ITD at elev=30°','ITD at elev=45°','ITD at elev=60°','ITD at elev=75°','ITD at elev=90°','fontsize',13);
set(gca,'FontSize',19,'XTick', [0 45 90 135 180 225 270 315 ])
xlabel('Azimuth Angle [deg]','fontsize',21);
ylabel('ITD [ms]','fontsize',21);
xlim([0 355])
grid on

set(gcf, 'Position', [0 0 1024 640])
saveas(gcf, strcat('plots/ITD Distribution.png'))

figure
surf(azim,elev,itd);
xlabel('Azimuth Angle [deg.]','fontsize',21);
ylabel('Elevation Angle [deg.]','fontsize',21);
zlabel('ITD [s]','fontsize',21);
title('ITD Distribution','fontsize',28);
xlim([0 355])
ylim([-40 90])
grid on

set(gcf, 'Position', [0 0 1024 1024])
