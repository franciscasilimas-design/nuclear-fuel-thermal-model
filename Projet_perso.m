clc
clear
close all

R1 = 4e-3;      % Rayon combustible (m)
gap = 80e-6;    % Epaisseur gap (m)
R2 = R1 + gap;
R3 = R2 + 0.6e-3; % Epaisseur gaine (m)

k1 = 2.5;       % UO2 (W/m.K)
k2 = 0.15;      % Helium (W/m.K)
k3 = 16;        % Zircaloy (W/m.K)

qppp = 5e8;     % Production volumique (W/m^3)
h = 15000;      % Coefficient convection (W/m^2.K)
Tinf = 580;     % Température eau (K)

%% Constantes analytiques

A2 = - qppp*R1^2/(2*k2);
A3 = - qppp*R1^2/(2*k3);

B3 = Tinf + (qppp*R1^2)/(2*h*R3) - A3*log(R3);

B2 = B3 + (A3 - A2)*log(R2);

C2 = - qppp*R1^2/(2*k2)*log(R1) + B2 + qppp*R1^2/(4*k1);

%% Discrétisation

r1 = linspace(0,R1,200);
r2 = linspace(R1,R2,100);
r3 = linspace(R2,R3,150);

%% Température

T1 = -qppp/(4*k1)*r1.^2 + C2;
T2 = A2*log(r2) + B2;
T3 = A3*log(r3) + B3;

%% Flux radial

q1 = -k1*(-qppp/(2*k1)*r1);
q2 = -k2*(A2./r2);
q3 = -k3*(A3./r3);

figure
plot(r1,T1,'LineWidth',2)
hold on
plot(r2,T2,'LineWidth',2)
plot(r3,T3,'LineWidth',2)
xlabel('Rayon (m)')
ylabel('Température (K)')
legend('Combustible','Gap He','Gaine')
title('Profil radial de température')
grid on

figure
plot(r1,q1,'LineWidth',2)
hold on
plot(r2,q2,'LineWidth',2)
plot(r3,q3,'LineWidth',2)
xlabel('Rayon (m)')
ylabel('Flux radial (W/m^2)')
legend('Combustible','Gap He','Gaine')
title('Profil radial du flux thermique')
grid on

Q_generated = qppp*pi*R1^2;
Q_surface = q3(end)*2*pi*R3;

disp([Q_generated Q_surface])


h_values = linspace(5000,40000,50);
Tc_h = zeros(size(h_values));

for i=1:length(h_values)
    h = h_values(i);
    B3 = Tinf + (qppp*R1^2)/(2*h*R3) - A3*log(R3);
    B2 = B3 + (A3 - A2)*log(R2);
    C2 = - qppp*R1^2/(2*k2)*log(R1) + B2 + qppp*R1^2/(4*k1);
    Tc_h(i) = C2;
end

figure
plot(h_values,Tc_h,'LineWidth',2)
xlabel('Coefficient h (W/m^2.K)')
ylabel('Température centrale (K)')
legend('Température centrale')
title('Influence du coefficient de convection sur la Température centrale')
grid on


gap_values = linspace(20e-6,145e-6,50);
Tc_gap = zeros(size(gap_values));

for i=1:length(gap_values)
    R2 = R1 + gap_values(i);
    A2 = - qppp*R1^2/(2*k2);
    A3 = - qppp*R1^2/(2*k3);
    B3 = Tinf + (qppp*R1^2)/(2*h*R3) - A3*log(R3);
    B2 = B3 + (A3 - A2)*log(R2);
    C2 = - qppp*R1^2/(2*k2)*log(R1) + B2 + qppp*R1^2/(4*k1);
    Tc_gap(i) = C2;
end

figure
plot(gap_values*1e6,Tc_gap,'LineWidth',2)
xlabel('Epaisseur gap (µm)')
ylabel('Température centrale (K)')
legend('Température centrale')
title('Influence du gap sur la Température centrale')
grid on



%% Paramètres non-linéaires UO2

a = 1/6;
b = 1/6000;

%% Température d'interface issue du modèle linéaire
Ts = T1(end);   % température à r = R1

%% Profil non-linéaire (Kirchhoff)

Tc_linear = T1(1);


T_nl = ( (a + b*Ts).*exp( (b*qppp/4)*(R1^2 - r1.^2) ) - a )/b;

Tc_nonlinear = T_nl(1);

k_center = 1/(a + b*Tc_nonlinear);



figure
plot(r1,T1,'LineWidth',2)
hold on
plot(r1,T_nl,'--','LineWidth',2)
xlabel('Rayon (m)')
ylabel('Température (K)')
legend('Linéaire','Non-linéaire k(T)')
title('Comparaison profil température combustible')
grid on

%Tc_linear = T1(1);
%Tc_nonlinear = T_nl(1);

%% RESOLUTION NUMERIQUE NON-LINEAIRE

N = 200;
r = linspace(0,R1,N)';
dr = r(2)-r(1);

Ts = T1(end); % température interface
T = T1';      % initialisation avec solution linéaire

tol = 1e-6;
err = 1;
iter = 0;

while err > tol && iter < 1000

    T_old = T;

    % Conductivité locale
    k = 1./(a + b*T);

    A = zeros(N,N);
    bvec = zeros(N,1);

    % Condition centre (symétrie)
    A(1,1) = 1;
    A(1,2) = -1;
    bvec(1) = 0;

    % Noeuds internes
    for i = 2:N-1

        rp = r(i) + dr/2;
        rm = r(i) - dr/2;

        kp = (k(i)+k(i+1))/2;
        km = (k(i)+k(i-1))/2;

        A(i,i-1) = rm*km/dr^2;
        A(i,i)   = - (rp*kp + rm*km)/dr^2;
        A(i,i+1) = rp*kp/dr^2;

        bvec(i) = - r(i)*qppp;
    end

    % Condition interface
    A(N,N) = 1;
    bvec(N) = Ts;

    % Résolution système linéaire
    T = A\bvec;

    err = max(abs(T - T_old));
    iter = iter + 1;
end




figure
plot(r1,T_nl,'LineWidth',2)
hold on
plot(r,T,'--','LineWidth',2)
xlabel('Rayon (m)')
ylabel('Température (K)')
legend('Kirchhoff analytique','Numérique DF')
title('Validation solution non-linéaire')
grid on





disp(['Convergence en ', num2str(iter), ' iterations'])
disp(['Tc linéaire = ', num2str(Tc_linear)])
disp(['Tc non-linéaire = ', num2str(Tc_nonlinear)])
disp(['k_center non-linéaire = ', num2str(k_center)])





% Puissance volumique (MW/m^3)
%H = [100 200 300 400 500];

% Température centrale modèle linéaire (K)
%Tc_lin = [864.6 1147.9 1431.4 1715.4 1999.4];

% Température centrale modèle non-linéaire (K)
%Tc_nl = [821.5 1088.5 1384.5 1709.9 2069.1];

% Tracé
%figure

%plot(H,Tc_lin,'-o','LineWidth',2)
%hold on

%plot(H,Tc_nl,'-s','LineWidth',2)
%grid on

%xlabel('Puissance volumique H (MW/m^3)','FontSize',12)
%ylabel('Température centrale T_c (K)','FontSize',12)

%title('Température centrale en fonction de la puissance volumique','FontSize',14)

%legend('Modèle linéaire','Modèle non-linéaire','Location','northwest')


