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
dy_b_in = 0.5 ;dy_t_in = 0.2;
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
COLOR = {'b',[0.9,0,0],[0.85,0.33,0.1],[0 0.7 0],[0 0.3 0.7],[0.49,0.18,0.56],[0.93,0.69,0.13],[0,0.45,0.74],[0.47,0.67,0.19],[0.3,0.75,0.93],'r','b','k'};
MARKER = {'s','o','^','x','V','d','>','<','p','h','+'};
markersize(1) = ms+1;

%%
subplot('position',pos11)
hold on

charge = '100';
DATA_uncharged = importdata('../Current-vs-E-SI-unchargedSurface.dat');
DATA_charged = importdata(['../Current-vs-E-SI-chargedSurface-' charge 'mV.dat']);
range = [0,2*DATA_charged(end,1)];

p_uncharged = polyfit(DATA_uncharged(:,1),DATA_uncharged(:,2),1);
p_charged = polyfit(DATA_charged(:,1),DATA_charged(:,2),1);

k = 1;plot(range,p_uncharged(1)*range+p_uncharged(2),'--','COLOR',COLOR{k},'Linewidth',linewidth)
k = 2;plot(range,p_charged(1)*range+p_charged(2),'--','COLOR',COLOR{k},'Linewidth',linewidth)

% k = 1;hand(k) = errorbar(DATA_uncharged(:,1),DATA_uncharged(:,2),DATA_uncharged(:,3),'Linestyle','none','Color',COLOR{k},'Marker',MARKER{k},'linewidth',linewidth/2,'markerfacecolor',COLOR{k},'MarkerEdgeColor',COLOR{k},'MarkerSize',markersize(k));
% k = 2;hand(k) = errorbar(DATA_charged(:,1),DATA_charged(:,2),DATA_charged(:,3),'Linestyle','none','Color',COLOR{k},'Marker',MARKER{k},'linewidth',linewidth/2,'markerfacecolor',COLOR{k},'MarkerEdgeColor',COLOR{k},'MarkerSize',markersize(k));

k = 1;hand(k) = plot(DATA_uncharged(:,1),DATA_uncharged(:,2),'Linestyle','none','Color',COLOR{k},'Marker',MARKER{k},'linewidth',linewidth/2,'markerfacecolor',COLOR{k},'MarkerEdgeColor',COLOR{k},'MarkerSize',markersize(k));
k = 2;hand(k) = plot(DATA_charged(:,1),DATA_charged(:,2),'Linestyle','none','Color',COLOR{k},'Marker',MARKER{k},'linewidth',linewidth/2,'markerfacecolor',COLOR{k},'MarkerEdgeColor',COLOR{k},'MarkerSize',markersize(k));

xlim([0,3e7])
ylim([0,2e7])
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

LEGEND = cell(2,2);
LEGEND{1,1} = ['uncharged,'];
LEGEND{2,1} = ['charged,'];
LEGEND{1,2} = ['$\sigma=' num2str(p_uncharged(1),2) '\;\mathrm{S/m}$'];
LEGEND{2,2} = ['$\sigma=' num2str(p_charged(1),2) '\;\mathrm{S/m}$'];
% leg = legend(hand,LEGEND,'interpreter','latex','units','normalized','position',[dx_l+0.3*xl,dy_b+0.85*yl,0.05,0.05]);
% leg.ItemTokenSize = [16,1];
% leg.NumColumns = 1;
% legend boxoff

% making a legend 
dy = 0.07*(YLIM(2)-YLIM(1));dx = 0.05*(XLIM(2)-XLIM(1));
x = 0.07*(XLIM(2)-XLIM(1))+XLIM(1);
y = 0.9*(YLIM(2)-YLIM(1))+YLIM(1);
k = 1;plot(x,y,'Linestyle','none','Color',COLOR{k},'Marker',MARKER{k},'linewidth',linewidth/2,'markerfacecolor',COLOR{k},'MarkerEdgeColor',COLOR{k},'MarkerSize',markersize(k))
k = 2;plot(x,y-dy,'Linestyle','none','Color',COLOR{k},'Marker',MARKER{k},'linewidth',linewidth/2,'markerfacecolor',COLOR{k},'MarkerEdgeColor',COLOR{k},'MarkerSize',markersize(k))

clear text
text(x+dx,y,LEGEND{1,1},'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
clear text
text(x+5.5*dx,y,LEGEND{1,2},'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

clear text
text(x+dx,y-dy,LEGEND{2,1},'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)
clear text
text(x+5.5*dx,y-dy,LEGEND{2,2},'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','latex','fontsize',man_fontsize)

box on
set(gca,'linewidth',figbox_linewidth,'fontsize',fontsize,'TickLabelInterpreter','latex','xcolor','k','ycolor','k');
hold off

NAME = ['Figure_IvsE.pdf'];
print(NAME,'-dpdf','-painters')

