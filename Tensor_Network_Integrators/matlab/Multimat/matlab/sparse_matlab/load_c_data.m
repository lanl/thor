function fout = load_c_data(filename)

temp = sprintf('../output/%s',filename);
fout = load(temp);
end