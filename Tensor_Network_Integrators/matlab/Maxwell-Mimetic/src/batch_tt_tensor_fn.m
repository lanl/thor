function [varargout] = batch_tt_tensor_fn(C)
%C is a cell array
for i = 1:nargout
  varargout{i} = tt_tensor(C{i});
end