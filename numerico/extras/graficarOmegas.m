1;

close all;
clear all;
graphics_toolkit("gnuplot");

om = dlmread("../w-optimo.csv","\t",2,0);
figure;
plot (om(:,1), om(:,2));

xlabel('Omega')
ylabel('Iteraciones hasta convergencia')
grid minor

print('omega.png','-dpng');

