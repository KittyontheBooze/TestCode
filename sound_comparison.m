fs = 48000;

root_source = 'sound_source';
x_s = load_wav(root_source, 'man_01.wav');
x_n = load_wav(root_source, 'man_02.wav');

N = fs*5;
x_s = x_s(1:N);
x_n = x_n(1:N);

t = (1:N)/fs;

figure
subplot(211)
plot(t, x_s, 'LineWidth', 1.5)
ylabel('Amplitude [au]','fontsize',17);

subplot(212)
plot(t, x_n, 'LineWidth', 1.5)
xlabel('Time [sec]','fontsize',17); ylabel('Amplitude [au]','fontsize',17);
saveas(gcf, "plots/fig files/Comparison of 2 Sound Sources.fig")


set(gcf, 'Position', [0 0 1024 640])
saveas(gcf, "plots/Comparison of 2 Sound Sources.png")


%{
figure

subplot(311)
sub(1) = subplot(3,1,1);
imagesc(t,ax,tf_table_ori,[0 1]); axis xy
colormap(sub(1),gray)
set(gca,'YDir','normal','FontSize',13)
if ae == 1
    set(gca,'YDir','normal','YTick', [0 90 180 270 355], 'FontSize',13)
    ylabel('Azimuth Angle [deg]','fontsize',17);
else
    set(gca,'YDir','normal','YTick', [-40 0 45 90], 'FontSize',13)
    ylabel('Elevation Angle [deg]','fontsize',17);
end
title("Unmasked T/F Table for SSL Algorithm",'fontsize',24);

subplot(312)
sub(2) = subplot(3,1,2);
imagesc(t,ax,isnan(tf_table_ori_mask_cor_wolpf),[0 1]); axis xy
colormap(sub(2),mycolormap)
set(gca,'YDir','normal','FontSize',13)
if ae == 1
    set(gca,'YDir','normal','YTick', [0 90 180 270 355], 'FontSize',13)
    ylabel('Azimuth Angle [deg]','fontsize',17);
else
    set(gca,'YDir','normal','YTick', [-40 0 45 90], 'FontSize',13)
    ylabel('Elevation Angle [deg]','fontsize',17);
end
title("Mask with Table with ICC without LPF",'fontsize',24);

subplot(313)
sub(3) = subplot(3,1,3);
imagesc(t,ax,isnan(tf_table_ori_mask_cor_lpf),[0 1]); axis xy
colormap(sub(3),mycolormap)
set(gca,'YDir','normal','YTick', [-40 0 45 90], 'FontSize',13)
if ae == 1
    set(gca,'YDir','normal','YTick', [0 90 180 270 355], 'FontSize',13)
    xlabel('Time [sec]','fontsize',17); ylabel('Azimuth Angle [deg]','fontsize',17);
else
    set(gca,'YDir','normal','YTick', [-40 0 45 90], 'FontSize',13)
    xlabel('Time [sec]','fontsize',17); ylabel('Elevation Angle [deg]','fontsize',17);
end
title("Mask with Table with ICC with LPF 1.6kHz",'fontsize',24);
grid on

set(gcf, 'Position', [0 0 1024 1024])
saveas(gcf, "plots/Masking Region Comparison ICC with and without(Reverberant, RT60 = 0."+num2str(rt60_bin(rt))+", "+string+").png")
saveas(gcf, "plots/fig files/Masking Region Comparison ICC with and without(Reverberant, RT60 = 0."+num2str(rt60_bin(rt))+", "+string+").fig")
%}