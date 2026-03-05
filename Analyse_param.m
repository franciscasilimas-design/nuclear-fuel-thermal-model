clc
clear
close all



% Puissance volumique (MW/m^3)
H = [100 200 300 400 500];

% Température centrale modèle linéaire (K)
Tc_lin = [864.6 1147.9 1431.4 1715.4 1999.4];

% Température centrale modèle non-linéaire (K)
Tc_nl = [821.5 1088.5 1384.5 1709.9 2069.1];

% Tracé
figure

plot(H,Tc_lin,'-o','LineWidth',2)
hold on

plot(H,Tc_nl,'-s','LineWidth',2)

grid on

xlabel('Puissance volumique H (MW/m^3)','FontSize',12)
ylabel('Température centrale T_c (K)','FontSize',12)

title('Température centrale en fonction de la puissance volumique','FontSize',14)

legend('Modèle linéaire','Modèle non-linéaire','Location','northwest')
set(gca,'FontSize',12)
xlim([0 550])
ylim([800 2100])
