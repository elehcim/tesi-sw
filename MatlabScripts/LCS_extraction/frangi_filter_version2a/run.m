%compile needed mex file
%mex eig3volume.c

%load('ExampleVolumeStent');
load ml
V=ftle_matrix;
% Frangi Filter the stent volume
options.BlackWhite=true;
options.FrangiScaleRange=[1 5];
Vfiltered=FrangiFilter3D(V,options);
sliceomatic(Vfiltered)
% figure,
% subplot(2,2,1), imshow(squeeze(max(V,[],2)),[])
% subplot(2,2,2), imshow(squeeze(max(Vfiltered,[],2)),[])
% subplot(2,2,3), imshow(V(:,:,100),[])
% subplot(2,2,4), imshow(Vfiltered(:,:,100),[])