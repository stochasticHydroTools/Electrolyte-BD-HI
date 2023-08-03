clc;clear;

warning off
format long
sub_label_x = -0.15;
sub_label_y = 0.97 ;

y_label_x = -0.1;
y_label_y = 0.5 ;
y_label_x_inset = -0.27;
y_label_y_inset = 0.5  ;

x_label_x = 0.5 ;
x_label_y = -0.13;
x_label_x_inset = 0.5 ;
x_label_y_inset = -0.3;

xl_in = 4;
yl_in = 3;
dx_l_in = 0.55;dx_r_in = 0.1;
dy_b_in = 0.55;dy_t_in = 0.35;
xwidth = 1*(xl_in+dx_l_in)+dx_r_in;
ywidth = 1*(yl_in+dy_b_in+dy_t_in);

figure('color','white','Units', 'inches','Position', [1, 1, xwidth, ywidth], 'PaperUnits', 'inches', 'PaperSize', [xwidth, ywidth])
hold on
man_fontsize = 12;
fontsize     = 11;
small_fontsize = 10;
linewidth    = 1.5 ;
figbox_linewidth = 1;
ms = 7;
markersize = ones(1,10)*ms;

dx_l = dx_l_in/xwidth;dx_r = dx_r_in/xwidth;
dy_b = dy_b_in/ywidth;dy_t = dy_t_in/ywidth;

xl = xl_in/xwidth;
yl = yl_in/ywidth;

pos11 = [dx_l,dy_b,xl,yl];

% ghermez, banafsh, abi, zard, sabz, abi roshan
COLOR = {'b',[0.9,0,0],[0 0.7 0],[0.85,0.33,0.1],[0 0.3 0.7],[0.49,0.18,0.56],[0.93,0.69,0.13],[0,0.45,0.74],[0.47,0.67,0.19],[0.3,0.75,0.93],'r','b','k'};
MARKER = {'s','o','^','x','V','d','>','<','p','h','+'};
markersize(1) = ms+1;

%%
subplot('position',pos11)
hold on

wetFractions = [0,0.5,1];
LEGEND = {'$0$','$\frac{1}{2}$','$1$'};
slopes = zeros(1,length(wetFractions));
for i = 1:length(wetFractions)
    DATA = importdata(['../Current-vs-E-SI-chargedSurface_wetFraction=' num2str(wetFractions(i)) '.dat']);
    p = polyfit(DATA(:,1),DATA(:,2),1);
    slopes(i) = p(1);
    range = [0,1.1*DATA(end,1)];
    plot(range,p(1)*range+p(2),'--','COLOR',COLOR{i},'Linewidth',linewidth)
    hand(i) = plot(DATA(:,1),DATA(:,2),'Linestyle','none','Color','k','Marker',MARKER{i},'linewidth',linewidth/2,'markerfacecolor',COLOR{i},'MarkerEdgeColor',COLOR{i},'MarkerSize',markersize(i));
    % errorbar(DATA(:,1),DATA(:,2),DATA(:,3),'Linestyle','none','Color','k','Marker',MARKER{i},'linewidth',linewidth/2,'markerfacecolor',COLOR{i},'MarkerEdgeColor','k','MarkerSize',markersize(i));
end

xlim(range)
ylim([0,1e9])
XLIM = xlim;YLIM = ylim;

xlab = '$E\;(\mathrm{V/m})$';
x_dim = x_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = x_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,xlab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
ylab = '$I\;(\mathrm{A/m^2})$';
x_dim = y_label_x*(XLIM(2)-XLIM(1))+XLIM(1);
y_dim = y_label_y*(YLIM(2)-YLIM(1))+YLIM(1);
clear text
text(x_dim,y_dim,ylab,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize,'Rotation',90)

leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.1*xl,dy_b+0.8*yl,0.05,0.05]);
leg.ItemTokenSize = [16,1];
leg.NumColumns = 1;
title(leg,'$w$','interpreter','latex')
legend boxoff

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k');
hold off

NAME = 'Figure_IvsE.pdf';
print(NAME,'-dpdf','-vector')

