function array_3d = load_c_3Ddata(filename,size)

temp = sprintf('../output/%s',filename);
data = load(temp);
array_3d = reshape(data', size); % MATLAB reads column-major, so transpose and reshape
array_3d = permute(array_3d, [3, 2, 1]); % Permute to get original dimensions (n x m x l)
end