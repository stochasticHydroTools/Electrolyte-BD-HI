clc;clear;

warning off
format long
sub_label_x = -0.15;
sub_label_y = 0.97 ;

y_label_x = -0.13;
y_label_y = 0.5  ;

x_label_x = 0.5  ;
x_label_y = -0.13;

xl_in = 3;
yl_in = 3;
dx_l_in = 0.6;dx_r_in = 0.1;
dy_b_in = 0.5;dy_t_in = 0.25;
xwidth = 1*(xl_in+dx_l_in)+dx_r_in;
ywidth = 1*(yl_in+dy_b_in)+dy_t_in;

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 11;
fontsize     = 10;
small_fontsize = 9;
linewidth    = 1 ;
figbox_linewidth = 1;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];
COLOR = {'w',[0,0.7,0],'b',[0.9,0,0],'k',[0.4940 0.1840 0.5560]};
MARKER = {'^','d','s','o'};
ms = [5,4,4,4];

results_DIR = '../results/Lxy=33.195a-H=8.29876a-(a=2.5e-10)/dt-Effect/';

subplot('position',pos11);
hold on

numberSimulations = 32;
H = 8.29876;
wetFraction = 1;
dts = [0.05,0.15,0.15,0.2];
methods = {'Leimkuhler','EulerMaruyama','Leimkuhler','Leimkuhler'};
LEGEND = {'ref','EM $\Delta t/\tau_D=0.15$','LK $\Delta t/\tau_D=0.15$','LK $\Delta t/\tau_D=0.2$'};
for j = 1:length(dts)
    dt = dts(j);
    method = methods{j};

    LOC = [results_DIR method, '-Equilibrium_chargedSurface-longrun-wetFraction=' num2str(wetFraction) '-dt=' num2str(dt) '/'];

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
    errorbar(z_c,DATA_SIM_MEAN_coarse,DATA_SIM_STD_coarse,'Linestyle','none','Color','k','Marker',MARKER{j},'linewidth',linewidth/2,'markerfacecolor',COLOR{j},'MarkerEdgeColor','k','MarkerSize',ms(j));
    hand(j) = plot(10,10,'Linestyle','none','Color','k','Marker',MARKER{j},'linewidth',linewidth/2,'markerfacecolor',COLOR{j},'MarkerEdgeColor','k','MarkerSize',ms(j));
    drawnow
end

xlim([0,1])
ylim([0,2.5])
XLIM = xlim;YLIM = ylim;
plot([1/H,1/H],[-10,10],'--k','linewidth',0.5)
plot([1-1/H,1-1/H],[-10,10],'--k','linewidth',0.5)

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

%     leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[i*dx_l+(i-1)*xl+0.45*xl,dy_b+0.75*yl,0.05,0.05]);
leg = legend(hand,LEGEND,'interpreter','latex','location','north');
leg.ItemTokenSize = [16,1];
leg.NumColumns = 1;
legend boxoff
% title(leg,'$\Delta t/\tau_D$','interpreter','latex')

box on
xtick = [0 0.2 0.4 0.6 0.8 1];
xticklabels = {'$0$','$0.2$','$0.4$','$0.6$','$0.8$','$1$'};
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k','xtick',xtick,'xticklabels',xticklabels);
hold off

NAME = ['Figure_equilibriumDistribution_UpdateMethod.pdf'];
print(NAME,'-dpdf','-vector')