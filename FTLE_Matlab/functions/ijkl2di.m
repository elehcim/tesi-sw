function  di= ijkl2di(i,j,k,l,n_vec)
di=n_vec(2)*n_vec(3)*n_vec(4)*(i-1)+...
	n_vec(3)*n_vec(4)*(j-1)+n_vec(4)*(k-1)+(l-1) +1;

%TODO generalizzare
% Fare una cosa del genere (sicuramente vanno aggiustati gli estremi degli
% indici)
%
% function  di= ijkl2di(ind_vec,n_vec)
% di=1;
% for i=1:length(ind_vec)
% 	di=di+prod(n_vec(end-i:end))*(ind_vec(i)-1);
% end
