clear all;
close all;

graphics_toolkit("gnuplot");

a = dlmread("../Prgm/salida.csv","\t");

figure
hold on

xlabel('Tiempo [Segundos]')
ylabel('Theta [Radianes]')
%axis([0 30 -2.5 2.5])
grid minor

CNT=rows(a)/3

plot(a(1:CNT,1),a(1:CNT,2),'r');
axis([0 100 -200 200])
%plot(a(CNT+1:2*CNT,1),a(CNT+1:2*CNT,2),'g');
plot(a(1,1),a(1,1),'g');
%plot(a(2*CNT+1:end,1),a(2*CNT+1:end,2),'b');
plot(a(1,1),a(1,1),'b');
legend({"Euler","RK2","RK4"},'Location','SouthEast')	
print('../informe/salida_2.png','-dpng');

