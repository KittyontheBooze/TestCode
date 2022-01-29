function hrtf_db = load_expandedDB(root_hrtf)

fs = 48000;
dt = 1/fs;

azims = linspace(0,355,72);
elevs = linspace(-40,90,27);

hrtf_db.azim.elev.hrtf_L = zeros(512,1);
hrtf_db.azim.elev.hrtf_R = zeros(512,1);
hrtf_db.azim.elev.ITD = 0;

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
        
        fclose(fid);
        
        j = j + 1;
    end
    i = i + 1;
end

end