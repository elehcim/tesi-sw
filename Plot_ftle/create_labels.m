function labels = create_labels( my_vis_var )
labels=cell(size(my_vis_var));
for i = 1:length(my_vis_var)
labels{i}=my_vis_var{i}(2:end); % remove 'n' from the name
end