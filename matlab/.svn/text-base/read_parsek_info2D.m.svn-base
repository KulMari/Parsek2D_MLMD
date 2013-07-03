%read_parsek_info2D(results_dir) - function to read information about *.hdf files saved in a Parsek2D simulation
% results_dir is a string variable with the path of the directory where *.hdf files are.
%
% Example:
% read_parsek_info2D('/home/results')
%
% This functions returns the Dataset saved in the *.hdf files, and some information about the simulation parameters.
%
% author: Enrico Camporeale
% e-mail: e.camporeale@qmul.ac.uk
% date:   01/01/07 

function read_parsek_info(results_dir)

if isempty(results_dir)
    results_dir=char(cellstr(pwd))
end
processor_name=[results_dir, '/resultsproc0.hdf'];
setting_name=[results_dir, '/resultssettings.hdf'];

info=hdf5info(processor_name);
list=[];
nGroups=size(info.GroupHierarchy.Groups,2);
    
for i=1:nGroups
    
    nnGroups=size(info.GroupHierarchy.Groups(i).Groups,2);

    for ii=1:nnGroups
      if strfind(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/species_')
      nnnGroups=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups,2);
        for iii=1:nnnGroups
             list=[list;{info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name}];
            
        end
      
      else
          list=[list;{info.GroupHierarchy.Groups(i).Groups(ii).Name}];
      end
    end
    

end

disp(['Files in directory ' results_dir ' contain the folowing Datasets:']);
disp(' ')
disp(list)
disp('The simulation is run with the following parameters:')
disp(['Lx = ' num2str( hdf5read(setting_name,'/collective/Lx')) '   Ly = ' num2str( hdf5read(setting_name,'/collective/Ly'))])
disp(['Nx = ' num2str( hdf5read(setting_name,'/collective/Nxc')) '   Ny = ' num2str( hdf5read(setting_name,'/collective/Nyc'))])
disp(['dx = ' num2str( hdf5read(setting_name,'/collective/Dx')) '   dy = ' num2str( hdf5read(setting_name,'/collective/Dy'))])
disp(['dt = ' num2str( hdf5read(setting_name,'/collective/Dt')) '   Th = ' num2str( hdf5read(setting_name,'/collective/Th'))   ])
disp(['# cycles = ' num2str( hdf5read(setting_name,'/collective/Ncycles')) '   # species = ' num2str( hdf5read(setting_name,'/collective/Ns'))   ])
disp(' ');disp('%%%%%%%%%%%%%%%%%%%%')

for i=0:hdf5read(setting_name,'/collective/Ns')-1
disp(['Species  ' num2str(i) ' -->  # particles = ' num2str(hdf5read(setting_name,['/collective/species_',num2str(i),'/Np'])) ...
   '   qom = ' num2str(hdf5read(setting_name,['/collective/species_',num2str(i),'/qom']))])
end
disp('%%%%%%%%%%%%%%%%%%%%');disp(' ')

disp(['The processor topology is ' num2str( hdf5read(setting_name,'/topology/XLEN')) 'x' num2str( hdf5read(setting_name,'/topology/YLEN'))])

%Np(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/Np']);
