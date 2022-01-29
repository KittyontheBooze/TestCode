clc; clear;
close all;

nfbin = 1024;

fid = fopen('hrtf_expdb.mat','r');
if fid == -1
    root_hrtf = 'HRTF_HATS';
    hrtf_db = load_expDB(root_hrtf,nfbin);
    save('hrtf_expdb.mat', 'hrtf_db')
else
    load('hrtf_expdb.mat');
    fclose(fid);
end

f = (1:nfbin)*(48000/nfbin);

ind = find_indices(f, 250, 500, 1000, 2000, 3000);
ild_ind1 = zeros(72*27,1);
ild_ind2 = zeros(72*27,1);
ild_ind3 = zeros(72*27,1);
ild_ind4 = zeros(72*27,1);
ild_ind5 = zeros(72*27,1);
azim = zeros(72*27,1);
elev = zeros(72*27,1);

for i = 1:72
    for j = 1:27
        ild_ind1((i-1)*27+j) = hrtf_db(1).azim(i).elev(j).ILD(ind(1));
        ild_ind2((i-1)*27+j) = hrtf_db(1).azim(i).elev(j).ILD(ind(2));
        ild_ind3((i-1)*27+j) = hrtf_db(1).azim(i).elev(j).ILD(ind(3));
        ild_ind4((i-1)*27+j) = hrtf_db(1).azim(i).elev(j).ILD(ind(4));
        ild_ind5((i-1)*27+j) = hrtf_db(1).azim(i).elev(j).ILD(ind(5));
        azim((i-1)*27+j) = (i-1) * 5;
        elev((i-1)*27+j) = (j-9) * 5;
    end
end

figure
scatter3(azim,elev,ild_ind1);
xlabel('Azimuth Angle [deg.]','fontsize',12);
ylabel('Elevation Angle [deg.]','fontsize',12);
zlabel('ILD at 250 Hz [db]','fontsize',12);
title('ILD at 250 Hz [db]','fontsize',14);
grid on

set(gcf, 'Position', [0 0 1024 1024])

figure
scatter3(azim,elev,ild_ind2);
xlabel('Azimuth Angle [deg.]','fontsize',12);
ylabel('Elevation Angle [deg.]','fontsize',12);
zlabel('ILD at 500 Hz [db]','fontsize',12);
title('ILD at 500 Hz [db]n','fontsize',14);
grid on

set(gcf, 'Position', [0 0 1024 1024])

figure
scatter3(azim,elev,ild_ind3);
xlabel('Azimuth Angle [deg.]','fontsize',12);
ylabel('Elevation Angle [deg.]','fontsize',12);
zlabel('ILD at 1000 Hz [db]','fontsize',12);
title('ILD at 1000 Hz [db]','fontsize',14);
grid on

set(gcf, 'Position', [0 0 1024 1024])

figure
scatter3(azim,elev,ild_ind4);
xlabel('Azimuth Angle [deg.]','fontsize',12);
ylabel('Elevation Angle [deg.]','fontsize',12);
zlabel('ILD at 2000 Hz [db]','fontsize',12);
title('ILD at 2000 Hz [db]','fontsize',14);
grid on

set(gcf, 'Position', [0 0 1024 1024])

figure
scatter3(azim,elev,ild_ind5);
xlabel('Azimuth Angle [deg.]','fontsize',12);
ylabel('Elevation Angle [deg.]','fontsize',12);
zlabel('ILD at 3000 Hz [db]','fontsize',12);
title('ILD at 3000 Hz [db]','fontsize',14);
grid on

set(gcf, 'Position', [0 0 1024 1024])

ind_1 = 0.05*abs(max(max(ild_ind1))-min(min(ild_ind1)));
ind_2 = 0.05*abs(max(max(ild_ind2))-min(min(ild_ind2)));
ind_3 = 0.05*abs(max(max(ild_ind3))-min(min(ild_ind3)));
ind_4 = 0.05*abs(max(max(ild_ind4))-min(min(ild_ind4)));
ind_5 = 0.05*abs(max(max(ild_ind5))-min(min(ild_ind5)));