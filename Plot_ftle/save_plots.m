function save_plots(handle,fig_name,output_dir,n1,n2,t0)
% TODO salvare roba 3d
% Salvare solo FTLE o solo LCS a seconda
if n1==n2
	f_save_name=sprintf('FTLE_%i_%.2f.fig',n1,t0);
	g_save_name=sprintf('LCS_%i_%.2f.fig',n1,t0);
else
	f_save_name=sprintf('FTLE_%ix%i_%.2f.fig',n1,n2,t0);
	g_save_name=sprintf('LCS_%ix%i_%.2f.fig',n1,n2,t0);
end
hgsave(ftle_fig_handle,[output_dir,f_save_name])
%hgsave(lcs_fig_handle,[output_dir,g_save_name])
end