function [index_low, index_high] = find_edge(f, f_low, f_high)

flag = 0;
for i = 1:length(f)
    if flag == 0
        if f(i) > f_low
            index_low = i;
            flag = 1;
            i = i + 1;
            continue;
        end
    else
        if f(i) > f_high
            index_high = i - 1;
            break;
        end
    end
    i = i + 1;
end

end