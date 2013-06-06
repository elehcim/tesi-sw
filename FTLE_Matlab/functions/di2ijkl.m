function [i,j,k,l]=di2ijkl(di,n_vec)
denom=[n_vec 1];
leng=length(n_vec);
for i=1:leng
	% metodo divisione intera
	denominat=prod(denom(end-(i-1):end));
	%floor((di-1)/denominat)
	ris(i)=mod((floor((di-1)/denominat)),n_vec(end-(i-1)))+1;
end
i=ris(4);
j=ris(3);
k=ris(2);
l=ris(1);
%TODO modify in a way that returns a vector
