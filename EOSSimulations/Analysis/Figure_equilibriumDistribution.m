clc;clear;

warning off
format long
sub_label_x = -0.15;
sub_label_y = 0.97 ;

y_label_x = -0.16;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5 ;
x_label_y = -0.15;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 2.5;
yl_in = 2.25;
dx_l_in = 0.6 ;dx_r_in = 0.1;
dy_b_in = 0.45;dy_t_in = 0.2;
xwidth = 2*(xl_in+dx_l_in)+dx_r_in;
ywidth = 1*(yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 1 ;
figbox_linewidth = 1;
ms = 4;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];
pos12 = [(dx_l+xl)+dx_l,dy_b,xl,yl];

% ghermez, banafsh, abi, zard, sabz, abi roshan
COLOR = {[0.9,0,0],'b'};
MARKER = {'o','s'};
numberSimulations = 7;
H = 40;
wetFractions = [0, 0.5];
positions = [pos11;pos12];
surfaces = {'uncharged','charged'};
labels = {'(a)','(b)'};
for j = 1:2
    surface = surfaces{j};
    subplot('position',positions(j,:))

    for i = 1:length(wetFractions)
        LOC = ['../Equilibrium_' surface 'Surface-longrun-wetFraction=' num2str(wetFractions(i)) '/'];
        hold on

        k = 1;
        DATA_sim = importdata([LOC 'densityplus' num2str(k) '.dat']);
        m = length(DATA_sim);
        DATA_SIM = zeros(200,numberSimulations);
        z = (DATA_sim(:,1)+1)/2;
        % LEGEND = cell(1,numberSimulations);
        for k = 1:numberSimulations
            DATA_sim = importdata([LOC 'densityplus' num2str(k) '.dat']);
            DATA_SIM(:,k) = DATA_sim(:,2);
        end
        DATA_SIM_MEAN = mean(DATA_SIM,2);
        DATA_SIM_STD = std(DATA_SIM,1,2);
        [ z_c ] = coarsen_x( z, H );
        DATA_SIM_MEAN_coarse = spline(z,DATA_SIM_MEAN,z_c);
        DATA_SIM_STD_coarse = spline(z,DATA_SIM_STD,z_c);
        errorbar(z_c,DATA_SIM_MEAN_coarse,DATA_SIM_STD_coarse,'Linestyle','none','Color','k','Marker',MARKER{i},'linewidth',linewidth/2,'markerfacecolor',COLOR{i},'MarkerEdgeColor','k','MarkerSize',ms);
        drawnow
    end
    xlim([0,1])
%     ylim([0,1])
    XLIM = xlim;YLIM = ylim;

    xlab = '$z/H$';
    x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
    y_dim = x_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
    clear text
    text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
    ylab = '$\displaystyle\frac{n_+}{N}$';
    x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
    y_dim = y_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
    clear text
    text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

    if j == 1
        LEGEND = {'$0$','$0.5$'};
        for i = 1:length(wetFractions)
            hand(i) = plot(10,10,'Linestyle','none','Color','k','Marker',MARKER{i},'linewidth',linewidth/2,'markerfacecolor',COLOR{i},'MarkerEdgeColor','k','MarkerSize',ms);
        end
        leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.45*xl,dy_b+0.2*yl,0.05,0.05]);
        leg.ItemTokenSize = [16,1];
        leg.NumColumns = 2;
        legend boxoff
        title(leg,'$w$','interpreter','latex')
    end

    INFO_X = 0.5;
    INFO_Y = 1.05;
    info = [surface];
    x_dim = INFO_X*(XLIM(2)-XLIM(1))+XLIM(1);
    y_dim = INFO_Y*(YLIM(2)-YLIM(1))+YLIM(1);
    clear text
    text(x_dim,y_dim,info,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

    box on
    label = labels{j};
    x_dim = sub_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
    y_dim = sub_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
    clear text
    text(x_dim,y_dim,label,'HorizontalAlignment','center','Interpreter','latex','fontsize',man_fontsize)

    xtick = [0 0.2 0.4 0.6 0.8 1];
    xticklabels = {'$0$','$0.2$','$0.4$','$0.6$','$0.8$','$1$'};
    set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','xtick',xtick,'xticklabels',xticklabels);
    hold off
end

NAME = ['Figure_equilibriumDistribution.pdf'];
print(NAME,'-dpdf','-painters')