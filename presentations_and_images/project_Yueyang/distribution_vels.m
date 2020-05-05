
clear

% grid
hycom_domain = 'GSH';
read_HYCOM_grid;
% clearvars -except JDM IDM plon1d plat1d plon plat

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

% figure;plot_field_model(um(:,:,1,1),plon1d,plat1d,'balance');caxis([-2 2]);colorbar

%%
%-------------------------------------------
%  velocity grouped into 2deg bins & normalize them
%-------------------------------------------
[dx,dy] = deal(100); % 100 grid points = ~2deg for GSH model
[nx,ny] = deal(ceil(IDM/dx),ceil(JDM/dy));

% vels normalized over data points within each bin
[u_norm,v_norm,um_norm,vm_norm,ue_norm,ve_norm] = deal(zeros(size(u)));

% 
for iz = 1:nz
    for it = 1:nt
        for i = 1:nx
            for j = 1:ny
                
                % indices of grid points within the bin
                [jj, ii] = deal((j-1)*dy+1 : j*dy, (i-1)*dx+1 : i*dx);
                if jj(end)>JDM
                    jj = jj(1):JDM;
                end
                if ii(end)>IDM
                    ii = ii(1):IDM;
                end
                
                u_norm(jj,ii,it,iz) = (u(jj,ii,it,iz) - nanmean(u(jj,ii,it,iz),'all')) / nanstd(u(jj,ii,it,iz),0,'all');
                v_norm(jj,ii,it,iz) = (v(jj,ii,it,iz) - nanmean(v(jj,ii,it,iz),'all')) / nanstd(v(jj,ii,it,iz),0,'all');
                
                um_norm(jj,ii,it,iz) = (um(jj,ii,it,iz) - nanmean(um(jj,ii,it,iz),'all')) / nanstd(um(jj,ii,it,iz),0,'all');
                vm_norm(jj,ii,it,iz) = (vm(jj,ii,it,iz) - nanmean(vm(jj,ii,it,iz),'all')) / nanstd(vm(jj,ii,it,iz),0,'all');
                
                ue_norm(jj,ii,it,iz) = (ue(jj,ii,it,iz) - nanmean(ue(jj,ii,it,iz),'all')) / nanstd(ue(jj,ii,it,iz),0,'all');
                ve_norm(jj,ii,it,iz) = (ve(jj,ii,it,iz) - nanmean(ve(jj,ii,it,iz),'all')) / nanstd(ve(jj,ii,it,iz),0,'all');
                
            end
        end
    end
end

%-------------------------------------------
%   energetic regions  
%-------------------------------------------
% either the magnitude of normalized zonal or meridional vel exceeded 3.5
% ifenergetic = sqrt(u_norm.^2 + v_norm.^2) > 3.5 ;
% ifenergetic = abs(u_norm) > 4 | abs(v_norm) > 4;

% figure
% m_proj('Equidistant Cylindrical','lon',[min(plon1d),max(plon1d)],'lat',[min(plat1d),max(plat1d)]);
% m_line(plon(ifenergetic(:,:,it,iz)),plat(ifenergetic(:,:,it,iz)),'marker','.','color','r','linest','none','markersize',3.5);
% m_grid('linestyle','none','tickdir','out','xtick',280:5:310,'xticklabels',[],...
%     'ytick',30:5:45,'yticklabels',[],'linewidth',1.2);
% m_coast('patch',[.8 .8 .8]);


%% set edges & binning the vel 

it = 1;
iz = 1;
[jjc,iic] = deal(1:JDM,1:IDM);
% figure;plot_field_model(u(jjc,iic,it,3),plon1d(iic),plat1d(jjc),'balance');cmocean('balance','pivot',0);caxis([-2 2]);colorbar

% set manually
f_u = u_norm(jjc,iic,it,iz);
f_v = v(jjc,iic,it,iz);

[vel_min,vel_max] = deal(-3,3);
[ymin, ymax] = deal(5e-5,1e1);
dvel_bin = 0.05;

%-------------------------------------------
%      edges
%-------------------------------------------
edges = vel_min:dvel_bin:vel_max;

% field to do [col vector]
f_u(isnan(f_u)) = [];
f_u = f_u';
f_v(isnan(f_v)) = [];
f_v = f_v';

%-------------------------------------------
%   calc pdfs 
%-------------------------------------------
[histValues_u,pdf_norm_u,pdf_exp_u,edges_center] = calc_pdfs(f_u,edges);
[histValues_v,pdf_norm_v,pdf_exp_v,~] = calc_pdfs(f_v,edges);

% statistics for the shape of distribution
[kurt_u,kurt_v] = deal(kurtosis(f_u),kurtosis(f_v));


%% plot

figure 
[ha, ~] = tight_subplot(1,2,[.05 .02],[.10 .05],[.05 .02]);
axes(ha(1))
% data
semilogy(edges_center,histValues_u,'o','Color','k','LineWidth',1.2,'MarkerSize',8)
xlabel('u [m/s]','Interpreter','latex','FontSize',16)
% normal dist curve
hold on 
plot(edges_center,pdf_norm_u,'b-','LineWidth',1.0)
% exponential dist curve
hold on
plot(edges_center,pdf_exp_u,'r-','LineWidth',1.2)

axes(ha(2))
semilogy(edges_center,histValues_v,'o','Color','k','LineWidth',1.2,'MarkerSize',8)
xlabel('v [m/s]','Interpreter','latex','FontSize',16)
% normal dist curve
hold on 
plot(edges_center,pdf_norm_v,'b-','LineWidth',1.0)
% exponential dist curve
hold on
plot(edges_center,pdf_exp_v,'r-','LineWidth',1.2)

set(ha(2),'YTickLabel','')

ax = ha(1);
ax.TickLength = [.015, .015];
ax.LineWidth = 1.5;
ax.XLim = [vel_min vel_max];
ax.YLim = [ymin, ymax];
ax.YMinorTick = 'off';
ax = ha(2);
ax.TickLength = [.015, .015];
ax.LineWidth = 1.5;
ax.XLim = [vel_min vel_max];
ax.YLim = [ymin, ymax];
ax.YMinorTick = 'off';

legend('data', 'norm', 'exp','fontsize',12)
% legend(ha(2),'data', 'norm','fontsize',12)
if iz == 1
    title(ha(2),'surface','fontsize',14,'Position',[-3.2,12,0]) % [-7.5,1.1,0]
elseif iz == 4
    title(ha(2),'deep (>1000m)','fontsize',14,'Position',[-1.8,12,0])
end

% set(gcf,'PaperPositionMode','auto');
% print(gcf,['pdf_uv_Z',num2str(iz)],'-dpng','-r400')

%%

function [histValues,pdf_norm,pdf_exp,edges_center] = calc_pdfs(f,edges)
% A local function to reduce the lines of code...

% edge centers
edges_center = 0.5 * (edges(1:end-1) + edges(2:end)); 
Nc = length(edges_center);

%   histogram of the data
[histValues,~] = histcounts(f,edges,'Normalization','pdf');

% normal fitting
pd_norm = fitdist(f,'Normal');
pdf_norm = pdf(pd_norm,edges_center);

% exponential fitting
pd_exp_r = fitdist(f(f>0),'Exponential');
pdf_exp_r = pdf(pd_exp_r,edges_center(Nc/2+1:Nc)) / 2;
pd_exp_l = fitdist(abs(f(f<0)),'Exponential');
pdf_exp_l = pdf(pd_exp_l,abs(edges_center(1:Nc/2))) / 2;

% two sided
pdf_exp = [pdf_exp_l pdf_exp_r];

end

