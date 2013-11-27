clear all;
close all;

graphics_toolkit("gnuplot");

a = dlmread("../Prgm/salida.csv","\t");

figure
hold on

xlabel('Tiempo [Segundos]')
ylabel('Theta [Radianes]')
axis([0 30 -2.5 2.5])
grid minor

plot(a(2001:end,1),a(2001:end,2));
print('../informe/salida.png','-dpng');
