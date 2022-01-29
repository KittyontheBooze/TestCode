function [a_ipd, a_ild] = load_crits(root_hrtf,nfbin)

fid = fopen(strcat('hrtf_expdb_nfbin=', num2str(nfbin), '.mat'),'r');
if fid == -1
    hrtf_db = load_expDB(root_hrtf,nfbin);
    save(strcat('hrtf_expdb_nfbin=', num2str(nfbin), '.mat'), 'hrtf_db')
else
    load(strcat('hrtf_expdb_nfbin=', num2str(nfbin), '.mat'));
    fclose(fid);
end

temp_ipd = zeros(72,27);
temp_ild = zeros(72,27);
a_ipd = zeros(nfbin/2-1,1);
a_ild = zeros(nfbin/2-1,1);
for i = 1:(nfbin/2-1)
    for k = 1:72
        for l = 1:27
            temp_ipd(i,k) = hrtf_db(1).azim(k).elev(l).IPD(i);
            temp_ild(i,k) = hrtf_db(1).azim(k).elev(l).ILD(i);
        end
    end
    a_ipd(i) = 0.05*(max(max(temp_ipd))-min(min(temp_ipd)));
    a_ild(i) = 0.05*(max(max(temp_ild))-min(min(temp_ild)));
end

end