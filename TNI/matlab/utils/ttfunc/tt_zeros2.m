function  T = tt_zeros2(C)
%this function will create a zero tensor train with the dimensions
%specificed in the cell array C

T = cellfun(@(c) zeros(1,c),C,'UniformOutput',false);
T = cell2core(tt_tensor, T);
end