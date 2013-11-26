clear all;
close all;

graphics_toolkit("gnuplot");

function y = fcn(x)
y= 1.2771*((cos((263*pi)/2000)^2*cos(x)^2)/(1-cos((263*pi)/2000)^2*sin(x)^2)^(3/2)-(cos((263*pi)/2000)^2*sin(x)^2)/(1-cos((263*pi)/2000)^2*sin(x)^2)^(3/2)+(3*cos((263*pi)/2000)^4*sin(x)^2*cos(x)^2)/(1-cos((263*pi)/2000)^2*sin(x)^2)^(5/2));
end

xx=0:0.01:2*pi;
y=arrayfun(@fcn,xx);

plot(xx,y);
print('d2.png','-dpng');
