function [difmat,intp,intp_noBC] = generate_1D_op_mat_fn(nx,dx,spsfmt)

if nargin<3;spsfmt=0;end %set default sparse format

%positive mu and eta
difmatp = 1/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));
intpmatp = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1));
intpmatp_noBC = intpmatp;
intpmatp_noBC(1,:) = 0;

%negative mu and eta
difmatm = 1/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
intpmatm = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),1));
intpmatm_noBC = intpmatm;
intpmatm_noBC(end,:) = 0;

difmat ={difmatm,difmatp};
intp = {intpmatm,intpmatp};
intp_noBC = {intpmatm_noBC, intpmatp_noBC};

if spsfmt
  for i = 1:numel(difmat)
    difmat{i} = sparse(difmat{i});
    intp{i} = sparse(intp{i});
    intp_noBC{i} = sparse(intp_noBC{i});
  end
end

end