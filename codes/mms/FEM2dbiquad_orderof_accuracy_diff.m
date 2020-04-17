%%
%% Script: ML_FEM_OrderOfAccuracy.m
%%
%% Compute order of accuracy using L2 error measure
%%

clear all, close all, clc

% Data obtained from running 
% FEM2d_diff.m 
% for biquadratic elements 
% 2x2 GL quadrature, 3x3 GL quadrature

data1 = [
1.2500e-01,	2.717701e-03;
6.2500e-02	3.422732e-04;
3.1250e-02,	4.288628e-05;
1.5625e-02,	5.364138e-06;
7.8125e-03,	6.706233e-07;
];

data2 = [
1.2500e-01,	2.703869e-03;
6.2500e-02,	3.383269e-04;
3.1250e-02,	4.231138e-05;
1.5625e-02,	5.289730e-06;
7.8125e-03,	6.612605e-07;
];

h1 = data1(:,1);
E1 = data1(:,2);
h2 = data2(:,1);
E2 = data2(:,2);

logh1 = log10(h1);
logE1 = log10(E1);
logh2 = log10(h2);
logE2 = log10(E2);

% plot log(h) vs log(E)
plot(logh1,logE1,'-<k'); 
hold on
plot(logh2,logE2,'->r'); 
xlabel('$log_{10}(h)$','Interpreter','latex','FontSize',20); 
ylabel('$log_{10}(E)$','Interpreter','latex','FontSize',20); 
legend({'$2 \times 2 \ quadrature$','$3 \times 3 \ quadrature$'},...
    'Interpreter','latex','Location','southeast');
grid on; 
drawnow;
set(findall(gcf,'type','axes'),'FontSize',18,'LineWidth',2)
set(findall(gcf,'type','line'),'LineWidth',2)
set(findall(gcf,'type','text'),'Interpreter','latex','FontSize',20)
set(findall(gcf,'tag','legend'),'Interpreter','latex','FontSize',18,'LineWidth',1)
print(gcf,'-dpng','-r300', '2Ddiffeqn_conv_biquad');

% perform linear regression and obtain the correlation coefficient
p1 = polyfit(logh1,logE1,1);
R1 = corrcoef(logh1,logE1);
p2 = polyfit(logh2,logE2,1);
R2 = corrcoef(logh2,logE2);

fprintf(1,'Order of accuracy 1:       %1.4f \n', p1(1) );
fprintf(1,'Correlation coefficient 1: %1.4f \n', R1(1,2) );
fprintf(1,'Order of accuracy 2:       %1.4f \n', p2(1) );
fprintf(1,'Correlation coefficient 2: %1.4f \n', R2(1,2) );




