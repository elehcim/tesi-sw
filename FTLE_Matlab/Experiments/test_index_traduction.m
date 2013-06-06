clear all
n_vec=[50 3 50 3]
di=floor(rand*(prod(n_vec)));
[i,j,k,l]=di2ijkl(di,n_vec);
fprintf('   di \t = %i -> i,j,k,l = %i %i %i %i\n',di,i,j,k,l)
di_check=ijkl2di(i,j,k,l,n_vec);
fprintf('Verifica:\ndi_check = %i\n',di_check)