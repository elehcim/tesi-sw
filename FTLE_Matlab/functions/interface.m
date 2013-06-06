function [x_T, y_T, vx_T, e_T] = interface(filename,nx, ny, nvx, ne)
% this function returns as vectors the results of the integration given as 
% an input the file "filename" as an input and the size of the grid

%TODO Ripulire il codice
% nx=n_vec(1);
% ny=n_vec(2);
% nvx=n_vec(3);
% ne=n_vec(4);

%TODO Inizializzare i vettori x_T y_T ...

%% Apro e leggo il file
fid=fopen(filename,'r');
[vett,cont]=fscanf(fid,'%g',[5 inf]);

%% Riordino i dati in colonna
%data=reshape(vett,[],5);
% size(data)
data=vett';
%% Elimino la colonnna relativa a vy
% data(:,4)=[];
% size(data)
% vett=reshape(data,[],1);
% length(vett)
% size(vett)

%% Memorizzo i dati
index=1;
for i=1:nx
	for j=1:ny
		for k=1:nvx
			for l=1:ne
				%ijkl2di(i,j,k,l,[nx,ny,nvx,ne]);
				x_T(i,j,k,l)=data(index,1);
				y_T(i,j,k,l)=data(index,2);
				vx_T(i,j,k,l)=data(index,3);
				e_T(i,j,k,l)=data(index,5);
				index=index+1;
			end
		end
	end
end
fclose(fid);