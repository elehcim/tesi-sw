function Files_path = folder
% Returns the path to the "Nostro software" directory
% This function returns the path to the Nostro software directory taking
% into account the OS of the user.
if ispc % per Windows user...
	userdir = getenv('USERPROFILE');
	Files_path = [userdir '\Dropbox\Tesi\Software\SoftwareTesi\'];
	%Files_path = [userdir '\Dropbox\Tesi\Software\Nostro software\'];
else % per Unix user
	userdir = getenv('HOME');
	Files_path = [userdir '/Dropbox/Tesi/Software/SoftwareTesi/'];
	%Files_path = [userdir '/Dropbox/Tesi/Software/Nostro software/'];
end