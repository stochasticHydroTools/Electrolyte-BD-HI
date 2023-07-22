clc;clear;

warning off
format long
sub_label_x = -0.15;
sub_label_y = 0.97 ;

y_label_x = -0.1;
y_label_y = 0.5 ;

x_label_x = 0.5 ;
x_label_y = -0.1;

xl_in = 3.5;
yl_in = 3.5;
dx_l_in = 0.6;dx_r_in = 0.1;
dy_b_in = 0.5;dy_t_in = 0.2;
xwidth = 1*(xl_in+dx_l_in)+dx_r_in;
ywidth = 1*(yl_in+dy_b_in)+dy_t_in;

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 1 ;
figbox_linewidth = 1;
ms = 5;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];

% ghermez, banafsh, abi, zard, sabz, abi roshan
COLOR = {[0.9,0,0],[0,0.7,0],'b'};
MARKER = {'o','s','^'};
numberSimulations = 32;
H = 8.29876;
wetFractions = [0,0.5,1];
dt = 0.1;
results_DIR = '../results/Lxy=33.195a-H=8.29876a-(a=2.5e-10)/w-Effect/';

subplot('position',pos11)
LEGEND = {'$0$','$0.5$','$1$'};
for i = 1:length(wetFractions)
    LOC = [results_DIR 'Leimkuhler-Equilibrium_chargedSurface-longrun-wetFraction=' num2str(wetFractions(i)) '-dt=' num2str(dt) '/'];
    hold on

    DATA_sim = importdata([LOC 'densityplus1.dat']);
    m = length(DATA_sim);
    DATA_SIM = zeros(200,numberSimulations);
    z = (DATA_sim(:,1)+1)/2;
    for k = 1:numberSimulations
        DATA_sim = importdata([LOC 'densityplus' num2str(k) '.dat']);
        DATA_SIM(:,k) = DATA_sim(:,2);
    end
    DATA_SIM_MEAN = mean(DATA_SIM,2);
    DATA_SIM_STD = std(DATA_SIM,1,2);
    [ z_c ] = coarsen_x( z, 9 );
    DATA_SIM_MEAN_coarse = spline(z,DATA_SIM_MEAN,z_c);
    DATA_SIM_STD_coarse = spline(z,DATA_SIM_STD,z_c);
    errorbar(z_c,DATA_SIM_MEAN_coarse,DATA_SIM_STD_coarse,'Linestyle','none','Color','k','Marker',MARKER{i},'linewidth',linewidth/2,'markerfacecolor',COLOR{i},'MarkerEdgeColor','k','MarkerSize',ms);
    hand(i) = plot(10,10,'Linestyle','none','Color','k','Marker',MARKER{i},'linewidth',linewidth/2,'markerfacecolor',COLOR{i},'MarkerEdgeColor','k','MarkerSize',ms);
    drawnow
end
xlim([0,1])
ylim([0,2.5])
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


LEGEND = cell(1,3);
for i = 1:3
    LEGEND{i} = ['$' num2str(wetFractions(i)) '$'];
end
for i = 1:length(wetFractions)
    
end
leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.45*xl,dy_b+0.85*yl,0.05,0.05]);
leg.ItemTokenSize = [16,1];
leg.NumColumns = 1;
legend boxoff
title(leg,'$w$','interpreter','latex')

box on
xtick = [0 0.2 0.4 0.6 0.8 1];
xticklabels = {'$0$','$0.2$','$0.4$','$0.6$','$0.8$','$1$'};
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','xtick',xtick,'xticklabels',xticklabels);
hold off

NAME = ['Figure_equilibriumDistribution.pdf'];
print(NAME,'-dpdf','-vector')