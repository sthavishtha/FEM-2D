%%
%% Script: ML_FEM_OrderOfAccuracy.m
%%
%% Compute order of accuracy using L2 error measure
%%

clear all, close all, clc

% Data obtained from running 
% FEM2Dbiquad_elasticity.m
% for biquadratic elements 
% 3x3 GL quadrature

data = [
2.5000e-01,	1.263989e-03;
1.2500e-01,	1.593697e-04;
6.2500e-02,	1.843803e-05;
3.1250e-02,	2.182951e-06;
1.5625e-02,	2.674489e-07;
];

h = data(:,1);
E = data(:,2);

logh = log10(h);
logE = log10(E);

% plot log(h) vs log(E)
plot(logh,logE,'-ok');
hold on
xlabel('$log_{10}(h)$','Interpreter','latex','FontSize',20); 
ylabel('$log_{10}(E)$','Interpreter','latex','FontSize',20); 
grid on; 
drawnow;
set(findall(gcf,'type','axes'),'FontSize',18,'LineWidth',2)
set(findall(gcf,'type','line'),'LineWidth',2)
set(findall(gcf,'type','text'),'Interpreter','latex','FontSize',20)
set(findall(gcf,'tag','legend'),'Interpreter','latex','FontSize',18,'LineWidth',1)
print(gcf,'-dpng','-r300', '2Delastic_conv_biquad');

% perform linear regression and obtain the correlation coefficient
p = polyfit(logh,logE,1);
R = corrcoef(logh,logE);

fprintf(1,'Order of accuracy 1:       %1.4f \n', p(1) );
fprintf(1,'Correlation coefficient 1: %1.4f \n', R(1,2) );




