clear

% grid
hycom_domain = 'GSH';
read_HYCOM_grid;

% vel m/s
thcks = ncread('modeldata_mpo624.nc','layers_thcks');
mld = ncread('modeldata_mpo624.nc','mixlay_depth');
ssh = ncread('modeldata_mpo624.nc','ssh');
u = ncread('modeldata_mpo624.nc','uGS');
v = ncread('modeldata_mpo624.nc','vGS');
ugeo = ncread('modeldata_mpo624.nc','uGeo');
vgeo = ncread('modeldata_mpo624.nc','vGeo');

nt = size(u,3);
nz = size(u,4);

KE = 0.5 * (u.^2 + v.^2);
KEgeo = 0.5 * (ugeo.^2 + vgeo.^2);
KEnogeo = 0.5 * ((u-ugeo).^2 + (v-vgeo).^2);

%-------------------------------------------
%  calc mean and eddy fields
%-------------------------------------------
[um,vm] = deal(zeros(size(u)));
dx_mod = 0.02; % [deg]
smdeg_half = round(1 ./ dx_mod);
for iz = 1:4
    % coarsen data in space
    um(:,:,:,iz) = smcs_UV_GSH(squeeze(u(:,:,:,iz)),0,0,1,smdeg_half,scux.*scuy);
    vm(:,:,:,iz) = smcs_UV_GSH(squeeze(v(:,:,:,iz)),0,0,1,smdeg_half,scvx.*scvy);
end
[ue,ve] = deal(u-um,v-vm);

KEm = 0.5 * (um.^2 + vm.^2);

%%  Calc the wavenumber power spectral density (psd)  
% Perform 1-D Fourier decomposition along each latitude line

it = 1;
iz = 1;

% field to to
f_2d = squeeze(KE(:,:,it,iz));
% figure;plot_field_model(sqrt(f_2d),plon1d,plat1d,'ice');colorbar

% Define the unit time/length. Here [m per unit len] is determined by model.
unitdist = 1 * 111e3;

% Sample rate. # of samples per unit time/length. This MUST be consistent
% with the unit time/length above.
Fs = 48;

% # of FFT points used to calculate the PSD estimate
NFFT = 512;

% zonal wavenumber psd
Pk = zeros(JDM,NFFT/2+1);
% metidional wavenumber psd
Pl = zeros(NFFT/2+1,IDM);

% loop over all latitudes
for js = 1:1:JDM 
    % for a specific lat
    f_1d = f_2d(js,:);
    f_1d(f_1d==0 | isnan(f_1d)) = [];
    if ~isempty(f_1d)
        % Fk: vector of normalized frequencies at which the PSD is estimated
        [Pk(js,:),Fk] = pwelch(f_1d,[],[],NFFT,Fs,'onesided');  
    else
        Pk(js,:) = NaN;
    end
end

% Fk [cycles / unit_len] to [cycles / km]
Fk_new = Fk / unitdist *1e3;

% Pk [m2/s2 / (cyc/unit_len)] to [m2/s2 / (cyc/m)]
Pk_new = Pk * unitdist;

%% plot
[jjc,iic] = deal(1:JDM,1:IDM);

figure;plot_field_model(sqrt(f_2d),plon1d,plat1d,'ice');caxis([0,2]);colorbar
title('Surface speed [m/s]')
hold on
m_rectangle(plon1d(1),plat1d(jjc(1)),plon1d(end)-plon1d(1),plat1d(jjc(end))-plat1d(jjc(1)),0,...
    'EdgeColor','r','LineWidth',3);

%%
figure
loglog(Fk_new,nanmean(Pk_new(jjc,:),1),'LineWidth',2.2)
xlabel('Wavenumber $[km^{-1}]$','Interpreter','latex','FontSize',14)
ylabel('KE spectral density $[m^3 s^{-2}]$','Interpreter','latex','FontSize',14)
% reference line
hold on 
x_3 = 8e-3:1e-3:1e-1;
y_3 = x_3.^ (-3) * 10^(-2.5);
loglog(x_3,y_3,'LineWidth',1)
text(x_3(20),y_3(10),'k^{-3}','FontSize',12)
hold on 
x_53 = 1e-3:1e-4:7e-3;
y_53 = x_53.^ (-5/3) * 10^(-.2);
loglog(x_53,y_53,'LineWidth',1)
text(x_53(20),y_53(10),'k^{-5/3}','FontSize',12)
hold on 
x_2 = 1e-2:1e-3:1e-1;
y_2 = x_2 .^ (-2) * 10^(-1);
% loglog(x_2,y_2,'LineWidth',1)

ax = gca;
ax.TickLength = [.015, .015];
ax.LineWidth = 1.5;
ax.XLim = [5e-4 3e-1];
ax.XGrid = 'off';
ax.XMinorGrid = 'off';
ax.XMinorTick = 'on';
ax.YLim = [1e-2 1e5]; %5e-4 4e4
ax.YGrid = 'off';
ax.YMinorGrid = 'off';
ax.YMinorTick = 'off';

axis square

% set(gcf,'PaperPositionMode','auto');
% print(gcf,['psd_all_Z',num2str(iz)],'-dpng','-r400')

%%



