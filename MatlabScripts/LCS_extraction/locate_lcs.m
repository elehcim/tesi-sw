function lcs=locate_lcs(v, limit)
% Suppress ftle less than a percentage of the maximum one.
if limit <=0
	error('limit must be positive');
elseif limit >= 1
	error('limit must be <= 1');
else
	v_max=max(v(:));
	v(v<(limit*v_max))=0;
	lcs=v;
end