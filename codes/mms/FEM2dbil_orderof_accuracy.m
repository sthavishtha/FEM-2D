%%
%% Script: ML_FEM_OrderOfAccuracy.m
%%
%% Compute order of accuracy using L2 error measure
%% Comment/uncomment array data below to reproduce the respective plot

clear all, close all, clc

%
% Data obtained from running 
% FEM2d_diff.m/FEM2Dbil_elasticity.m
% with the following mesh sequence
% no_elements = {4, 8, 16, 32, 64, 128}
%

% FEM2d_diff.m
data = [
2.5000e-01,  8.716565e-02;
1.2500e-01,  2.184438e-02; 
6.2500e-02,  5.464083e-03;
3.1250e-02,  1.366203e-03;
1.5625e-02,  3.415622e-04;
7.8125e-03,  8.539125e-05;
];

% FEM2Dbil_elasticity.m
% data = [
% 2.5000e-01, 1.195967e-02;
% 1.2500e-01,	2.970692e-03;
% 6.2500e-02,	7.410421e-04;
% 3.1250e-02,	1.851371e-04;
% 1.5625e-02, 4.627524e-05;
% ];

h = data(:,1);
E = data(:,2);

logh = log10(h);
logE = log10(E);

% plot log(h) vs log(E)
plot(logh,logE,'-ok'); 
xlabel('$log_{10}(h)$','Interpreter','latex','FontSize',20); 
ylabel('$log_{10}(E)$','Interpreter','latex','FontSize',20); 
grid on; 
drawnow;
set(findall(gcf,'type','axes'),'FontSize',18,'LineWidth',2)
set(findall(gcf,'type','line'),'LineWidth',2)
set(findall(gcf,'type','text'),'Interpreter','latex','FontSize',20)
print(gcf,'-dpng','-r300', '2Ddiffeqn_conv_bil');
% print(gcf,'-dpng','-r300', '2Delastic_conv_bil');

% perform linear regression and obtain the correlation coefficient
p = polyfit(logh,logE,1);
R = corrcoef(logh,logE);

fprintf(1,'Order of accuracy:       %1.4f \n', p(1) );
fprintf(1,'Correlation coefficient: %1.4f \n', R(1,2) );




