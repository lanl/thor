function  T = tt_ones2(C)
%this function will create an all-one tensor train with the dimensions
%specificed in the cell array C

T = cellfun(@(c) ones(1,c),C,'UniformOutput',false);
T = cell2core(tt_tensor, T);
end