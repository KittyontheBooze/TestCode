function f_index = find_indices(f_bin, varargin)

f_index = [];

for i =  1 :nargin-1
    val = f_bin-varargin{i};
    [~, index] = min(abs(val));
    f_index = [f_index index];
end

end