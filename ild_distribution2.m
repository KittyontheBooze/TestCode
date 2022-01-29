clc; clear;
close all;

nfbin = 24000;

fid = fopen(strcat('hrtf_expdb_nfbin=',num2str(nfbin),'.mat'),'r');
if fid == -1
    root_hrtf = 'HRTF_HATS';
    hrtf_db = load_expDB(root_hrtf,nfbin);
    save(strcat('hrtf_expdb_nfbin=',num2str(nfbin),'.mat'), 'hrtf_db')
else
    load(strcat('hrtf_expdb_nfbin=',num2str(nfbin),'.mat'));
    fclose(fid);
end

f = (1:nfbin)*(48000/nfbin);

ind = find_indices(f, 250, 500, 1000, 2000, 3000,5000, 7000);
ild_ind1 = zeros(27, 72);
ild_ind2 = zeros(27, 72);
ild_ind3 = zeros(27, 72);
ild_ind4 = zeros(27, 72);
ild_ind5 = zeros(27, 72);
ild_ind6 = zeros(27, 72);
ild_ind7 = zeros(27, 72);

[azim, elev] = meshgrid(0:5:355,-40:5:90);

for i = 1:72
    for j = 1:27
        ild_ind1(j, i) = hrtf_db(1).azim(i).elev(j).ILD(ind(1));
        ild_ind2(j, i) = hrtf_db(1).azim(i).elev(j).ILD(ind(2));
        ild_ind3(j, i) = hrtf_db(1).azim(i).elev(j).ILD(ind(3));
        ild_ind4(j, i) = hrtf_db(1).azim(i).elev(j).ILD(ind(4));
        ild_ind5(j, i) = hrtf_db(1).azim(i).elev(j).ILD(ind(5));
        ild_ind6(j, i) = hrtf_db(1).azim(i).elev(j).ILD(ind(6));
        ild_ind7(j, i) = hrtf_db(1).azim(i).elev(j).ILD(ind(7));
    end
end

figure
surf(azim,elev,ild_ind1);
xlabel('Azimuth Angle [deg.]','fontsize',21);
ylabel('Elevation Angle [deg.]','fontsize',21);
zlabel('ILD [db]','fontsize',21);
title('ILD at 250 Hz [db]','fontsize',24);
xlim([0 355])
ylim([-40 90])
grid on

set(gcf, 'Position', [0 0 768 768])
saveas(gcf, strcat('plots/ILD at 250 Hz [db].png'))

figure
surf(azim,elev,ild_ind2);
set(gca,'FontSize',19)
xlabel('Azimuth Angle [deg.]','fontsize',21);
ylabel('Elevation Angle [deg.]','fontsize',21);
zlabel('ILD [db]','fontsize',21);
title('ILD at 500 Hz [db]','fontsize',24);
xlim([0 355])
ylim([-40 90])
grid on

set(gcf, 'Position', [0 0 768 768])
saveas(gcf, strcat('plots/ILD at 500 Hz [db].png'))

figure
surf(azim,elev,ild_ind3);
set(gca,'FontSize',19)
xlabel('Azimuth Angle [deg.]','fontsize',21);
ylabel('Elevation Angle [deg.]','fontsize',21);
zlabel('ILD [db]','fontsize',21);
title('ILD at 1000 Hz [db]','fontsize',24);
xlim([0 355])
ylim([-40 90])
grid on

set(gcf, 'Position', [0 0 768 768])
saveas(gcf, strcat('plots/ILD at 1000 Hz [db].png'))

figure
surf(azim,elev,ild_ind4);
set(gca,'FontSize',19)
xlabel('Azimuth Angle [deg.]','fontsize',21);
ylabel('Elevation Angle [deg.]','fontsize',21);
zlabel('ILD [db]','fontsize',21);
title('ILD at 2000 Hz [db]','fontsize', 24);
xlim([0 355])
ylim([-40 90])
grid on

set(gcf, 'Position', [0 0 768 768])
saveas(gcf, strcat('plots/ILD at 2000 Hz [db].png'))

figure
surf(azim,elev,ild_ind5);
set(gca,'FontSize',19)
xlabel('Azimuth Angle [deg.]','fontsize',21);
ylabel('Elevation Angle [deg.]','fontsize',21);
zlabel('ILD [db]','fontsize',21);
title('ILD at 3000 Hz [db]','fontsize',24);
xlim([0 355])
ylim([-40 90])
grid on

set(gcf, 'Position', [0 0 768 768])
saveas(gcf, strcat('plots/ILD at 3000 Hz [db].png'))

figure
surf(azim,elev,ild_ind6);
set(gca,'FontSize',19)
xlabel('Azimuth Angle [deg.]','fontsize',21);
ylabel('Elevation Angle [deg.]','fontsize',21);
zlabel('ILD [db]','fontsize',21);
title('ILD at 5000 Hz [db]','fontsize',24);
xlim([0 355])
ylim([-40 90])
grid on

set(gcf, 'Position', [0 0 768 768])
saveas(gcf, strcat('plots/ILD at 5000 Hz [db].png'))

figure
surf(azim,elev,ild_ind7);
set(gca,'FontSize',19)
xlabel('Azimuth Angle [deg.]','fontsize',21);
ylabel('Elevation Angle [deg.]','fontsize',21);
zlabel('ILD [db]','fontsize',21);
title('ILD at 7000 Hz [db]','fontsize',24);
xlim([0 355])

set(gcf, 'Position', [0 0 1024 640])
%saveas(gcf, strcat('plots/ILD Distribution.png'))

figure
plot(azim, ild_ind3(9,:),azim, ild_ind4(9,:),azim, ild_ind5(9,:), azim, ild_ind6(9,:), azim, ild_ind7(9,:), 'LineWidth', 2);
set(gca,'FontSize',19,'XTick', [0 45 90 135 180 225 270 315 ])
xlabel('Azimuth Angle [deg]','fontsize',21);
ylabel('ILD [db]','fontsize',21);
legend('ILD at 1000 Hz', 'ILD at 2000 Hz', 'ILD at 3000 Hz', 'ILD at 5000 Hz', 'ILD at 7000 Hz', 'Location', 'southeast', 'FontSize', 14)
xlim([0 355])
grid on

set(gcf, 'Position', [0 0 1024 640])
saveas(gcf, strcat('plots/ILD Distribution.png'))

xlim([0 355])

set(gcf, 'Position', [0 0 1024 640])
%saveas(gcf, strcat('plots/ITD Distribution.png'))