1;

close all;
clear all;
graphics_toolkit("gnuplot");

figure;
datos = dlmread("dia4.csv","\t",0,0);
polar ((360-datos(:,1)+90)*pi/180, datos(:,2));
print('dia4.png','-dpng');