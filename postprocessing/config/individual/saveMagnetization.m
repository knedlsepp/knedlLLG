function saveMagnetization(mesh,path,name,mhj,j,variant) %#ok<*INUSL>
%SAVEMAGNETIZATION    Saves the magnetization.
%   See also:
%   PLOTANDSAVEPOSTPROCESSOR
%   PLOTPOSTPROCESSOR
%   SAVEPOSTPROCESSOR
%   ENERGYPLOTPOSTPROCESSOR
%
%   Author: Josef Kemetmueller - 16.12.2013
VTKQUIVERURL = 'http://www.mathworks.com/matlabcentral/fileexchange/33656-vtkpipe';
filename = sprintf('%s/m_%s_%d',path,name,j);
switch variant
    case 'txt'
        dlmwrite([filename,'.txt'], mhj);
    case 'mat'
        save(filename,'mhj');
    otherwise
        error('Unknown option %s', variant);
end
try
	coo = 1e9*mesh.coordinates; % Rescale from nanometers
	fileout = [filename, '.vtu'];
	vtkquiver(coo(:,1), coo(:,2), coo(:,3), ...
	          mhj(:,1), mhj(:,2), mhj(:,3), 'magnetization', fileout);
	%%
	bdN = mesh.bd.nodesOf3DMesh;
	fileout = [filename, '_bd.vtu'];
	vtkquiver(coo(bdN,1), coo(bdN,2), coo(bdN,3), mhj(bdN,1), mhj(bdN,2), mhj(bdN,3), 'magnetization', fileout);
catch
	disp(['To save the magnetization as vtu-file, download vtkquiver from here: ',VTKQUIVERURL]);
end
