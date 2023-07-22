function [ x_c ] = coarsen_x_stretched( x, H )
i_rich_l = find(x>=1/H,1);
i_rich_h = find(x>=4/H,1);

n_l = 2;
n_m = 7;
n_h = 10;
x_l = [linspace(0.25/H,x(i_rich_l),n_l)];
x_m = [linspace(x(i_rich_l),x(i_rich_h),n_m)];
x_h = [linspace(x(i_rich_h),0.5,n_h)];
x_c = [x_l,x_m(2:end),x_h(2:end)];
x_c = [x_c,1-x_c(end-1:-1:1)];