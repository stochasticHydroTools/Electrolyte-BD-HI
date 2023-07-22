function [ x_c ] = coarsen_x( x, c )
IM = floor(length(x)/c);
x_c = linspace(0,0.5,IM);
x_c = [x_c,1-x_c(end-1:-1:1)];