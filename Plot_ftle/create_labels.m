function labels = create_labels( my_vis_var )
labels=cell(size(my_vis_var));
for i = 1:length(my_vis_var)
	switch my_vis_var{i}
		case 'nx'
			lab='$x$ (non-dim)';
		case 'ny'
			lab='$y$ (non-dim)';
		case 'nvx'
			lab='$\dot{x}$ (non-dim)';
		case 'nvy'
			lab='$\dot{y}$ (non-dim)';
		case 'ne'
			lab='$e$ (non-dim)';
	end
	labels{i}=lab;
end