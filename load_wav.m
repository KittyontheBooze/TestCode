function [s] = load_wav(root, file_name)
%
% function [s] = load_wav(root, dir_ch, file_name)
% Return sound source to test in 48kHz sampling rate
%	root is root directory.
%	dir_ch is directory character, '/' (unix) or ':' (mac).
%	file_name is the name of the sound file to read
%
path_temp = sprintf('%s/%s',root, file_name);
[s_temp,Fs] = audioread(path_temp);

if Fs ~= 48000
    s = resample(s_temp, 48000, Fs);
else
    s = s_temp;
end
end


