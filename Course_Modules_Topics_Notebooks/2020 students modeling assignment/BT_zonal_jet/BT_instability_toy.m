% 
% Study barotropic instability of zonal jet on the f-plane.
% The flow is a zonal eastward jet initially and evolves according to a set
% of shallow water equations (SWE). No external forcing and dissipation are
% imposed.
% 
% The shallow water model has a zonally periodic boundary condition.
% 
%%
clear
load('BT_inst.mat')

[nj,ni,~] = size(u);

%% 
% zonal mean field
um = repmat(mean(u,2),[1 ni 1]);
vm = repmat(mean(v,2),[1 ni 1]);
zetam = repmat(mean(zeta,2),[1 ni 1]);

% eddy field
ue = u - um;
ve = v - vm;
zetae = zeta - zetam;

% Reynolds stress, i.e., eddy momentum flux: <u'v'>
ueve = repmat(mean(ue .* ve,2),[1 ni 1]);
conv_ueve = zeros(size(u));
conv_ueve(2:end-1,:,:) = - (ueve(3:end,:,:) - ueve(1:end-2,:,:)) ./ dy / 2;
% : <u'u'>
ueue = repmat(mean(ue .* ue,2),[1 ni 1]);

% gradient of mean vorticity
dzmdy = zeros(size(u));
dzmdy(2:end-1,:,:) =  (zetam(3:end,:,:) - zetam(1:end-2,:,:)) ./ dy / 2;

% gradient of mean zonal vel
dumdy = zeros(size(u));
dumdy(2:end-1,:,:) =  (um(3:end,:,:) - um(1:end-2,:,:)) ./ dy / 2;

% Scales for the mean quantities
U0 = max(um(:,:,1),[],'all');
Z0 = max(zetam(:,:,1),[],'all');
ZY0 = max(dzmdy(:,:,1),[],'all');

% time derivative 
um_dt = zeros(size(u));
um_dt(:,:,2:end-1) = (um(:,:,3:end) - um(:,:,1:end-2)) ./ dt_save / 2;

ueve_dt = zeros(size(u));
ueve_dt(:,:,2:end-1) = (ueve(:,:,3:end) - ueve(:,:,1:end-2)) ./ dt_save / 2;

ueue_dt = zeros(size(u));
ueue_dt(:,:,2:end-1) = (ueue(:,:,3:end) - ueue(:,:,1:end-2)) ./ dt_save / 2;

%% Zonal mean vel  vs.  mean PV (RV)
figure
% 'it' is hour
for it = 50:50
    
    clf
    subplot(121)
    plot(um(:,1,it)./U0,0:nj-1,'k','linewidth',1.5)
    xlabel('$\bar{u}/U_{0max}$','Interpreter','latex','FontSize',14);
    ylabel('y (1000s of km)','Interpreter','latex','FontSize',14);
        
    text(0, nj, ['Time = ' num2str(t_save(it)./3600) ' hours'],...
        'verticalalignment','bottom','fontsize',12);
    set(gca,'xlim',[-0.05 1],'ylim',[0 nj-1])
    
    subplot(122)
    plot(zetam(:,1,it)./Z0,0:nj-1,'k','linewidth',1.5)
    xlabel('$\bar{\zeta}/\bar{\zeta}_{0max}$','Interpreter','latex','FontSize',14);
    set(gca,'xlim',[-1 1],'ylim',[0 nj-1])
    
    drawnow
    pause(0.2);
        
end

%% <u'v'> vs.  d<u>/dy
figure

for it = 51:51
    
%     subplot(11)
%     plot(ueve_dt(:,1,it),0:nj-1,'k','linewidth',1.5)
%     xlabel('$\overline{u^{\prime}v^{\prime}}$','Interpreter','latex','FontSize',14);
    
    subplot(131)
    plot(ueve(:,1,it),0:nj-1,'k','linewidth',1.5)
    xlabel('$\overline{u^{\prime}v^{\prime}}$','Interpreter','latex','FontSize',14);
    set(gca,'xlim',[-100 100],'ylim',[0 nj-1])
    text(0, nj, ['Time = ' num2str(t_save(it)./3600) ' hours'],...
        'verticalalignment','bottom','fontsize',12);
    
    subplot(132)
    plot(-dumdy(:,1,it),0:nj-1,'k','linewidth',1.5)
    xlabel('$-\partial_{y}\bar{u}$','Interpreter','latex','FontSize',14);
    set(gca,'xlim',[-5e-4 5e-4],'ylim',[0 nj-1])
    
    subplot(133)
    plot(-ueve(:,1,it)./dumdy(:,1,it),0:nj-1,'k','linewidth',1.5)
    xlabel('$- \overline{u^{\prime}v^{\prime}} / \partial_{y}\bar{u}$','Interpreter','latex','FontSize',14);
    set(gca,'xlim',[-5e5 5e5],'ylim',[0 nj-1])
    
end

%% dum/dt vs. div of flux

figure
for it = 51:5:300
    
    clf
    
    % d<u>/dt
    subplot(121)
    plot(um_dt(:,1,it),0:nj-1,'k','linewidth',1.4)
    % conv of eddy flux  - d<u'v'>/dy
    hold on
    plot(conv_ueve(:,1,it),0:nj-1,'k--','linewidth',1.4)
    %
    text(0, nj, ['Time = ' num2str(t_save(it)./3600) ' hours'],...
        'verticalalignment','bottom','fontsize',12);
    xlabel('$m/s^{2} $','Interpreter','latex','FontSize',14);
    ylabel('y (1000s of km)','Interpreter','latex','FontSize',14);
    legend('$\partial_{t}\overline{u}$','$ - \nabla \cdot \overline{u^{\prime}v^{\prime}} $',...
        'Interpreter','latex','FontSize',10);
    set(gca,'xlim',[-5e-4 5e-4],'ylim',[0 nj-1])

    drawnow
    pause(0.5);
        
end


%% d<u'u'>/dt vs. d<u>dy * <u'v'>

figure
for it = 51:51
    
    clf
    % d<u'u'>/dt
    subplot(121)
    plot(ueue_dt(:,1,it),0:nj-1,'k','linewidth',1.5)
    xlabel('$\partial_{t}\overline{{u^{\prime}}^{2}}$','Interpreter','latex','FontSize',14);
    ylabel('y (1000s of km)','Interpreter','latex','FontSize',14);
        
    text(0, nj, ['Time = ' num2str(t_save(it)./3600) ' hours'],...
        'verticalalignment','bottom','fontsize',12);
    
    % KE consersion: d<u>dy * <u'v'>
    subplot(122)
    plot(-dumdy(:,1,it).*ueve(:,1,it),0:nj-1,'k','linewidth',1.5)
    xlabel('$ - \partial_{y}\bar{u} \cdot\overline{u^{\prime}v^{\prime}} $','Interpreter','latex','FontSize',14);
        
%     set(gca,'xlim',[-6e-4 6e-4],'ylim',[0 nj-1])
    
    drawnow
    pause(0.5);
        
end

%%
outpath = '/Users/yueyang/Documents/MATLAB/Mo_CFD/barotrop_inst_simulation/BT_instability_toymodel/images/';
set(gcf,'units','inches');
pos = get(gcf,'position');
pos([3 4]) = [10.5 5];
set(gcf,'position',pos)

% Axis units are thousands of kilometers (x and y are in metres)
x_sc = xp(1,:).*1e-6;
y_sc = yp(:,1)'.*1e-6;
height_scale = 0.001;

% Interval between arrows in the velocity vector plot
interval = 6;

% Set this to "true" to save each frame as a png file
plot_frames = true;


% Loop through the frames of the animation
for it = 1%:4:ntsave-140
    clf
    
    % height and velocity field at this instant
    h2d = eta(:,:,it);
    u2d = u(:,:,it);
    v2d = v(:,:,it);
    zeta2d = zeta(:,:,it);
    
    %-----------------------------
    % Plot the height & velocity field
    %-----------------------------
    % Plot height field
    h1 = axes('Position',[0.06 0.56 0.74 0.4]);
    handle = pcolor(x_sc,y_sc,h2d*height_scale);%
    set(handle,'EdgeColor', 'none');
    cmocean('balance')
    caxis([9.5 10.5])
    % Plot the velocity vectors
    hold on;
    quiver(x_sc(3:interval:end),y_sc(3:interval:end), u2d(3:interval:end,3:interval:end),...
        v2d(3:interval:end,3:interval:end),'Color','w');
    ylabel('Y (1000s of km)');
    title(['$ \bf{\eta \, (km) \; and \; \vec{u} \,(m/s) }$'],'Interpreter','latex');
    cb = colorbar;
    set(cb,'position',[0.92 0.56 0.03 0.4]);
    set(gca,'tickdir','out','xticklabel',[]);
    % plot zonal mean u
    h2 = axes('Position',[0.82 0.56 0.08 0.4]);
    plot(um(:,1,it)./U0,0:nj-1,'k','linewidth',1.2)
    title(['$\bar{u}/U_{0max}$'],'Interpreter','latex');
    set(gca,'tickdir','out','xlim',[0 1],'ylim',[0 nj-1],'ytick',0:10:nj-1,'yticklabel',[]);
    
    T = text(0, nj, ['Time = ' num2str(t_save(it)./3600) ' hours'],...
        'verticalalignment','bottom','fontsize',14);
    set(T,'position',[-9.5 -9.24 0])
    
    %-----------------------------
    % Compute the vorticity
    %-----------------------------
    % Plot the vorticity
    h3 = axes('Position',[0.06 0.07 0.74 0.4]);
    handle = pcolor(x_sc, y_sc, zeta2d);
    set(handle,'EdgeColor', 'none');
    cmocean('curl')
    caxis([-3 3].*1e-4);
    ylabel('Y (1000s of km)');
    title(['$ \bf{\zeta \, (s^{-1})}$'],'Interpreter','latex');
    cb = colorbar;
    set(cb,'position',[0.92 0.07 0.03 0.4]);
    set(gca,'tickdir','out')
    %     plot mean vort
    h4 = axes('Position',[0.82 0.07 0.08 0.4]);
    plot(zetam(:,1,it)./Z0,0:nj-1,'k','linewidth',1.2)
    title('$\bar{\zeta}/\bar{\zeta}_{0max}$','Interpreter','latex');
    set(gca,'tickdir','out','xlim',[-1 1],'ylim',[0 nj-1],'ytick',0:10:nj-1,'yticklabel',[]);
    
    % Now render this frame
    drawnow
        
    % To make an animation we can save the frames as a
    % sequence of images
    
%     set(gcf,'PaperPositionMode','auto');
%     print(gcf,[outpath 'pic',num2str(it,'%03d')],'-dpng','-r200')
%     clf
end
