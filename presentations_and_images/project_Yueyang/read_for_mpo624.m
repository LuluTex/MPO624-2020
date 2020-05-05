
% bsub -J mpo624 -P ome -o pt.o%J -e pt.e%J -W 0:15 -q general -n 1 -R "rusage[mem=20000]" matlab -r read_for_mpo624

addpath(genpath('/nethome/yxl1496/HYCOM'));
addpath(genpath('/nethome/yxl1496/mytoolbox'));

clear
hycom_domain = 'GSH';
read_HYCOM_grid;

day2sec=86400;            % seconds per day
deg2rad=pi/180;           % rad per degree
coeff=9.8/(4*pi/day2sec); % g/f0 where f0 is the planetary vort at pole

%% The times and layers you want to read
dt_ori = 0.5;
[day_s, day_interv] = deal(370, 0.5);
layers = [1 12 20 25]; % 4.77m, ~250m, ~1100m deep. [Surf, below ML, deep layer]
nlayers = length(layers);

day_e = day_s + day_interv;
tuv = day_s:dt_ori:day_e;
ntuv = length(tuv);

%% load the layer thicknesses at each time 

% thicknesses of each layer that lays above and the layer
layers_thcks = zeros(JDM,IDM,ntuv,nlayers);
[mixlay_depth,ssh_al,uGeo,vGeo] = deal(zeros(JDM,IDM,ntuv));
[uGS,vGS] = deal(zeros(JDM,IDM,ntuv,nlayers));

k_dm = 4;

% vars to be loaded
THCK_INDEX = 1; 
MLD_INDEX = 1;
UVEL_INDEX = 1;
VVEL_INDEX = 1;
SSH_INDEX = 1;

iic = 2:IDM; jjc = 2:JDM;

% loop for time first
nind = 1;
for ndy = tuv
    
    % day & year: 366 days for yr8 and 133 for yr9
    if ndy <= 366.5
        nday = floor(ndy);
        nyear = 8;
    else
        nday = floor(ndy-366);
        nyear = 9;
    end
    % hour: '00' midnight, '12' noon
    hour_num = num2str(mod(ndy,1)*24,'%2.2i');
    
    % read thck from top down
    for ilay = 1:nlayers
        % load thckness of current layer
        klay = layers(ilay);
        read_arch_GSHR
        layers_thcks(:,:,nind,ilay) = thck;
        mixlay_depth(:,:,nind) = mld;
        ssh_al(:,:,nind) = ssh;
        uGS(:,:,nind,ilay) = uvel;
        vGS(:,:,nind,ilay) = vvel;
        uGeo(jjc,iic,nind) = -(ssh(jjc,iic).*scpx(jjc,iic)-ssh(jjc-1,iic).*scpx(jjc-1,iic))...
            ./(scvx(jjc,iic).*scvy(jjc,iic))*coeff./sin(deg2rad*ulat(jjc,iic));
        vGeo(jjc,iic,nind) = (ssh(jjc,iic).*scpy(jjc,iic)-ssh(jjc,iic-1).*scpy(jjc,iic-1))...
            ./(scux(jjc,iic).*scuy(jjc,iic))*coeff./sin(deg2rad*vlat(jjc,iic));
        fclose('all');
    end
    nind = nind + 1;
end

% layers_depths = cumsum(layers_thcks,k_dm,'omitnan') - layers_thcks / 2;

%% save the data
otptdir = '/nethome/yxl1496/HYCOM/';
otptfile = [otptdir 'modeldata_mpo624.nc'];

nccreate(otptfile,'layers_thcks','Dimensions',{'dim1',JDM,'dim2',IDM,'dim3',ntuv,'dim4',nlayers});
nccreate(otptfile,'mixlay_depth','Dimensions',{'dim1',JDM,'dim2',IDM,'dim3',ntuv});
nccreate(otptfile,'ssh','Dimensions',{'dim1',JDM,'dim2',IDM,'dim3',ntuv});
nccreate(otptfile,'uGeo','Dimensions',{'dim1',JDM,'dim2',IDM,'dim3',ntuv});
nccreate(otptfile,'vGeo','Dimensions',{'dim1',JDM,'dim2',IDM,'dim3',ntuv});
nccreate(otptfile,'uGS','Dimensions',{'dim1',JDM,'dim2',IDM,'dim3',ntuv,'dim4',nlayers});
nccreate(otptfile,'vGS','Dimensions',{'dim1',JDM,'dim2',IDM,'dim3',ntuv,'dim4',nlayers});

ncwrite(otptfile,'layers_thcks',layers_thcks);
ncwrite(otptfile,'mixlay_depth',mixlay_depth);
ncwrite(otptfile,'ssh',ssh_al);
ncwrite(otptfile,'uGS',uGS);
ncwrite(otptfile,'vGS',vGS);
ncwrite(otptfile,'uGeo',uGeo);
ncwrite(otptfile,'vGeo',vGeo);
