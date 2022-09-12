tic
clc
clear all;

Ep=7e5;%Young's modulus of platelet
Em=1960;%Young's modulus of matrix
num=0.4;%Poisson's ratio of matrix
sigmcrit=243.05;%Normal strength of matrix
sigpcrit=1400;%Normal strength of platelet
taumcrit=140;%Shear strength of matrix



Gm=Em/(2*(1+num));%Shear modulus of matrix
sigtotau=sigpcrit/taumcrit;
gbe=Gm/Ep;
gam_m_crit=taumcrit/Gm;
rcount=0;
scount=0;
A=zeros();
vf=80;
m=0:0.01:1;
B=zeros();
Em_norm=Em/Ep;
sigm_norm=sigmcrit/sigpcrit;
wm_norm=(sigmcrit^2/(2*Em))/(sigpcrit*gam_m_crit);

[Eibar,sigibar,Uibar,rhopt0,bet0,rho0,C,Ec,sigc,Uc,fc]=deal(zeros(size([])));
for i=1:numel(m)
   n = 0:0.01:1;
      for j=1:numel(n)
          k=0:0.01:1;
          for q=1:numel(k)
              [sum]=deal(zeros(size(k)));
              sum(i,j,q)=k(q)+n(j)+m(i); %checks whether the sum m+n+k is 1
              if (sum(i,j,q)==1)
                  rcount=rcount+1;
                  %f=zeros(rcount,3);
                  A(rcount,1)=m(i);
                  A(rcount,2)=n(j);
                  A(rcount,3)=k(q);
                  %1/1+fc used for maximizing fc
                  fc=@(x)1/(1+(((x(1)*x(2))/(((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*((x(4)*(3*x(4) - 4))/(3*(x(4) - 1)^2) - (Ep*x(3)^2*x(4)^2*taumcrit^2*(2*num + 2)*(3*x(2) - 3))/(Em*x(1)*x(2)*sigpcrit^2*((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*(3*x(4) - 3)*(x(3) - 1)^2*(x(4) - 1)^2))))^m(i)*((x(1)*x(2)*(x(3) - 1)*(x(4) - 1))/(x(3)*x(4)))^n(j)*((Em*x(1)*x(2)*sigpcrit*((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*((x(4)*(3*x(4) - 4))/(3*(x(4) - 1)^2) - (Ep*x(3)^2*x(4)^2*taumcrit^2*(2*num + 2)*(3*x(2) - 3))/(Em*x(1)*x(2)*sigpcrit^2*((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*(3*x(4) - 3)*(x(3) - 1)^2*(x(4) - 1)^2))*(x(3) - 1)^2*(x(4) - 1)^2)/(2*Ep*x(3)^2*x(4)^2*taumcrit*(2*num + 2)))^k(q))); 
                  x0 = [(vf/100)-0.01,(vf/100)-0.01,2,2];
                  lb = x0;%[0.59,0.69,2,2];
                  ub = [vf/100,vf/100,20,20];
                  A = [];
                  b = [];
                  Aeq = [];
                  beq = [];
                  x = fmincon(fc,x0,A,b,Aeq,beq,lb,ub);
                  fm=(Em_norm^m(i))*(sigm_norm^n(j))*(wm_norm^k(q));
                  frat=@(x)(((x(1)*x(2))/(((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - ...(Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*((x(4)*(3*x(4) - 4))/(3*(x(4) - 1)^2) - (Ep*x(3)^2*x(4)^2*taumcrit^2*(2*num + 2)*(3*x(2) - 3))/(Em*x(1)*x(2)*sigpcrit^2*((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*(3*x(4) - 3)*(x(3) - 1)^2*(x(4) - 1)^2))))^m(i)*...
                      ((x(1)*x(2)*(x(3) - 1)*(x(4) - 1))/(x(3)*x(4)))^n(j)*...
                      ((Em*x(1)*x(2)*sigpcrit*((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*((x(4)*(3*x(4) - 4))/(3*(x(4) - 1)^2) - (Ep*x(3)^2*x(4)^2*taumcrit^2*(2*num + 2)*(3*x(2) - 3))/(Em*x(1)*x(2)*sigpcrit^2*((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*(3*x(4) - 3)*(x(3) - 1)^2*(x(4) - 1)^2))*(x(3) - 1)^2*(x(4) - 1)^2)/(2*Ep*x(3)^2*x(4)^2*taumcrit*(2*num + 2)))^k(q)))))/fm; 
                  Ec_norm=@(x)((x(1)*x(2))/(((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*((x(4)*(3*x(4) - 4))/(3*(x(4) - 1)^2) - (Ep*x(3)^2*x(4)^2*taumcrit^2*(2*num + 2)*(3*x(2) - 3))/(Em*x(1)*x(2)*sigpcrit^2*((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*(x(3) - 1)^2))*(3*x(4) - 3)*(x(3) - 1)^2*(x(4) - 1)^2))));
                  Sigc_norm=@(x)((x(1)*x(2)*(x(3) - 1)*(x(4) - 1))/(x(3)*x(4)));
                  wc_norm=@(x)((Em*x(1)*x(2)*sigpcrit*((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2)...
                      - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*...
                      (x(3) - 1)^2))*((x(4)*(3*x(4) - 4))/(3*(x(4) - 1)^2) - (Ep*x(3)^2*x(4)^2*taumcrit^2*(2*num + 2)*...
                      (3*x(2) - 3))/(Em*x(1)*x(2)*sigpcrit^2*((x(3)*(3*x(3) - 4))/(3*(x(3) - 1)^2) ...
                      - (Ep*x(3)^2*taumcrit^2*(2*num + 2)*(3*x(1) - 3))/(Em*x(1)*sigpcrit^2*(3*x(3) - 3)*...
                      (x(3) - 1)^2))*(3*x(4) - 3)*(x(3) - 1)^2*(x(4) - 1)^2))*(x(3) - 1)^2*...
                      (x(4) - 1)^2)/(2*Ep*x(3)^2*x(4)^2*taumcrit*(2*num + 2)));
 %frat=@(Ec_norm,Sigc_norm,wc_norm)(Ec_norm^m(i) * Sigc_norm^n(j) * wc_norm^k(q));
                      scount=scount+1;
                      B(scount,1)=m(i);
                      B(scount,2)=n(j);
                      B(scount,3)=k(q);
                      B(scount,4)=x(1);
                      B(scount,5)=x(2);
                      B(scount,6)=x(3);
                      B(scount,7)=x(4);
                      B(scount,8)=frat(x);
                      B(scount,9)=Ec_norm(x);
                      B(scount,10)=Sigc_norm(x);
                      B(scount,11)=wc_norm(x);

              end
          end
      end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','fc/fi');
%subplot(1,3,1)
% Plot the data
% First set the colormap (can't be done afterwards)
colormap(jet)
[hg,htick,hcb]=tersurf(B(:,1),B(:,2),B(:,3),B(:,8));
% Add the labels
hlabels=terlabel('m','n','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following modifications are not serious, just to illustrate how to
% use the handles:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--  Change the color of the grid lines
set(hg(:,3),'color','r')
set(hg(:,2),'color','b')
set(hg(:,1),'color','m')

%--  Modify the labels
set(hlabels,'fontsize',20)
set(hlabels(3),'color','r')
set(hlabels(2),'color','b')
set(hlabels(1),'color','m')
%--  Modify the tick labels
set(htick(:,3),'color','r','linewidth',3,'fontsize',14)
set(htick(:,2),'color','b','linewidth',3,'fontsize',14)
set(htick(:,1),'color','m','linewidth',3,'fontsize',14)
set(gca,'fontsize',16)
%--  Change the colorbar
set(hcb,'xcolor','w','ycolor','w')
%--  Modify the figure color
set(gcf,'color',[0 0 0.3])
set(gcf,'color','w')

%-- Change some defaults
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
h = colorbar;
ylabel(h, '$$\tilde{f_c}$$/$$\tilde{f_i}$$', 'Interpreter', 'LaTeX','fontsize',22,'color','b')
set(gca,'fontsize',16)
caxis([0, 25]);
filename1 = ['fcbyfi',num2str(vf),'.fig'];
saveas(gcf,filename1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('Name','\phi_1');
% % Plot the data
% % First set the colormap (can't be done afterwards)
% colormap(jet)
% [hg,htick,hcb]=tersurf(B(:,1),B(:,2),B(:,3),B(:,4));
% % Add the labels
% hlabels=terlabel('m','n','k');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % The following modifications are not serious, just to illustrate how to
% % use the handles:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %--  Change the color of the grid lines
% set(hg(:,3),'color','r')
% set(hg(:,2),'color','b')
% set(hg(:,1),'color','m')
% 
% %--  Modify the labels
% set(hlabels,'fontsize',20)
% set(hlabels(3),'color','r')
% set(hlabels(2),'color','b')
% set(hlabels(1),'color','m')
% %--  Modify the tick labels
% set(htick(:,3),'color','r','linewidth',3,'fontsize',14)
% set(htick(:,2),'color','b','linewidth',3,'fontsize',14)
% set(htick(:,1),'color','m','linewidth',3,'fontsize',14)
% set(gca,'fontsize',16)
% %--  Change the colorbar
% set(hcb,'xcolor','w','ycolor','w')
% %--  Modify the figure color
% set(gcf,'color',[0 0 0.3])
% set(gcf,'color','w')
% 
% %-- Change some defaults
% set(gcf,'paperpositionmode','auto','inverthardcopy','off')
% h = colorbar;
% ylabel(h, '\phi_1','fontsize',22,'color','b')
% set(gca,'fontsize',16)
% caxis([0, 0.9]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('Name','\phi_2');
% % Plot the data
% % First set the colormap (can't be done afterwards)
% colormap(jet)
% [hg,htick,hcb]=tersurf(B(:,1),B(:,2),B(:,3),B(:,5));
% % Add the labels
% hlabels=terlabel('m','n','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following modifications are not serious, just to illustrate how to
% use the handles:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %--  Change the color of the grid lines
% set(hg(:,3),'color','r')
% set(hg(:,2),'color','b')
% set(hg(:,1),'color','m')
% 
% %--  Modify the labels
% set(hlabels,'fontsize',20)
% set(hlabels(3),'color','r')
% set(hlabels(2),'color','b')
% set(hlabels(1),'color','m')
% %--  Modify the tick labels
% set(htick(:,3),'color','r','linewidth',3,'fontsize',14)
% set(htick(:,2),'color','b','linewidth',3,'fontsize',14)
% set(htick(:,1),'color','m','linewidth',3,'fontsize',14)
% set(gca,'fontsize',16)
% %--  Change the colorbar
% set(hcb,'xcolor','w','ycolor','w')
% %--  Modify the figure color
% set(gcf,'color',[0 0 0.3])
% set(gcf,'color','w')
% 
% %-- Change some defaults
% set(gcf,'paperpositionmode','auto','inverthardcopy','off')
% h = colorbar;
% ylabel(h, '\phi_2','fontsize',22,'color','b')
% set(gca,'fontsize',16)
% caxis([0, 0.9]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','n_1');
%subplot(1,3,2)
% Plot the data
% First set the colormap (can't be done afterwards)
colormap(jet)
[hg,htick,hcb]=tersurf(B(:,1),B(:,2),B(:,3),B(:,6));
% Add the labels
hlabels=terlabel('m','n','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following modifications are not serious, just to illustrate how to
% use the handles:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--  Change the color of the grid lines
set(hg(:,3),'color','r')
set(hg(:,2),'color','b')
set(hg(:,1),'color','m')

%--  Modify the labels
set(hlabels,'fontsize',20)
set(hlabels(3),'color','r')
set(hlabels(2),'color','b')
set(hlabels(1),'color','m')
%--  Modify the tick labels
set(htick(:,3),'color','r','linewidth',3,'fontsize',14)
set(htick(:,2),'color','b','linewidth',3,'fontsize',14)
set(htick(:,1),'color','m','linewidth',3,'fontsize',14)
set(gca,'fontsize',16)
%--  Change the colorbar
set(hcb,'xcolor','w','ycolor','w')
%--  Modify the figure color
set(gcf,'color',[0 0 0.3])
set(gcf,'color','w')

%-- Change some defaults
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
h = colorbar;
ylabel(h, 'n_1','fontsize',22,'color','b')
set(gca,'fontsize',16)
caxis([2, 20]);
filename2 = ['n1_',num2str(vf),'.fig'];
saveas(gcf,filename2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','n_2');
%subplot(1,3,3)
% Plot the data
% First set the colormap (can't be done afterwards)
colormap(jet)
[hg,htick,hcb]=tersurf(B(:,1),B(:,2),B(:,3),B(:,7));
% Add the labels
hlabels=terlabel('m','n','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following modifications are not serious, just to illustrate how to
% use the handles:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--  Change the color of the grid lines
set(hg(:,3),'color','r')
set(hg(:,2),'color','b')
set(hg(:,1),'color','m')

%--  Modify the labels
set(hlabels,'fontsize',20)
set(hlabels(3),'color','r')
set(hlabels(2),'color','b')
set(hlabels(1),'color','m')
%--  Modify the tick labels
set(htick(:,3),'color','r','linewidth',3,'fontsize',14)
set(htick(:,2),'color','b','linewidth',3,'fontsize',14)
set(htick(:,1),'color','m','linewidth',3,'fontsize',14)
set(gca,'fontsize',16)
%--  Change the colorbar
set(hcb,'xcolor','w','ycolor','w')
%--  Modify the figure color
set(gcf,'color',[0 0 0.3])
set(gcf,'color','w')

%-- Change some defaults
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
h = colorbar;
ylabel(h, 'n_2','fontsize',22,'color','b')
set(gca,'fontsize',16)
caxis([2, 20]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename3 = ['n2_',num2str(vf),'.fig'];
saveas(gcf,filename3);
filename = 'diffMP_ternary_data_04Apr2022.xlsx';
sheetname='Zhang_MP';
writematrix(B,filename,'Sheet',sheetname,'Range','A2')
Btitle={'m','n','k','phi1','phi2','n1','n2','fcrat','Ec_norm','Sigc_norm','wc_norm'};
writecell(Btitle,filename,'Sheet',sheetname,'Range','A1')
toc











