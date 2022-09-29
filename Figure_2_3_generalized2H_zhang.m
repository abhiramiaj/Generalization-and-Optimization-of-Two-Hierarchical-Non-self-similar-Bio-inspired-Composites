%%Generalized two hierarchical model Zhang
%%%%------------2HSS-----------%%%%
%ps--volume fraction at first level of hierarchy
%ps2--volume fraction at second level of hierarchy
%rs--aspect ratio rho at first level of hierarchy
%rs2--aspect ratio rho at second level of hierarchy
%n1--n at first level of hierarchy
%n2--n at second level of hierarchy
%as--alpha at first level of hierarchy
%ass--alpha at second level of hierarchy
clc
clear all

figure
subplot(1,2,1)
%Strength ratio 2H Regular staggered composites with Stairwise staggered
%composites as platelets
%for all values of n (4 CASES)
% syms y r
%rho=aspect ratio; 
%subplot(1,2,1);
n1=2:1:10; %Number of platelets in a period at the first level of hierarchy
ps1=0.5;%Volume fraction at the first level of hierarchy
ps2=0.5;%Volume fraction at the second level of hierarchy
gbe=1./1000;%Shear modulus of matrix to Young's modulus of platetet ratio
rs1=10;%Platelet aspect ratio at the first level of hierarchy
rs2=10;%Platelet aspect ratio at the second level of hierarchy
z=0.5;%overlap ratio
Ep=7e9;%Young's modulus of platelet
sigtotau=10;%Ratio between normal strength of platelet and shear strength of matrix
for n2=[2,3,4,5,6,7,8,9,10] %Number of platelets in a period at the second level of hierarchy
as=(ps1.*rs1.*rs1.*gbe)./(3*(1-ps1));
D=(((n1.*((3*n1)-4))./(3*((n1-1).^2)))+(((n1.*n1)./(3*(n1-1).*as))));
a_ss=(ps2.*rs2.*rs2.*gbe*D)./(3*(1-ps2)*ps1);
Erat2hss=(D.*(((n2.*((3*n2)-4))./(3*((n2-1).^2)))+(((n2.*n2)./(3*(n2-1).*a_ss))))).^-1
Ess=Erat2hss*Ep*ps1*ps2;
rs1crit2hss=(n1-1)*sigtotau;
plot (n1,Erat2hss,'linewidth', 1.5);
hold on;
end
%  title('SW within Regular')
xlabel ('n_1', 'fontsize', 18,'fontweight','bold')
ylabel ('E_{SS}/Ep \phi_{1} \phi_{2}', 'fontsize', 18,'fontweight','bold')
legend ('n_2 = 2', 'n_2 = 3', 'n_2 = 4','n_2 = 5','n_2 = 6','n_2 = 7','n_2 = 8','n_2 = 9','n_2 = 10','fontsize', 12);
pbaspect([5 4 1])
set(gca,'fontsize',14)
set(gcf,'color','w')

subplot(1,2,2)

for n2=[2 3 4 5 6 7 8 9 10]
    for i=1:numel(n1)
        as(i)=(ps1.*rs1.*rs1.*gbe)./(3*(1-ps1));
        D(i)=(((n1(i).*((3*n1(i))-4))./(3*((n1(i)-1).^2)))+(((n1(i).*n1(i))./(3*(n1(i)-1).*as(i)))));
        a_ss(i)=(ps2.*rs2.*rs2.*gbe*D(i))./(3*(1-ps2)*ps1);
        Erat2hss(i)=(D(i).*(((n2.*((3*n2)-4))./(3*((n2-1).^2)))+(((n2.*n2)./(3*(n2-1).*a_ss(i)))))).^-1
        Ess(i)=Erat2hss(i)*Ep*ps1*ps2;
        rs1crit2hss(i)=(n1(i)-1)*sigtotau;
        if rs1<=rs1crit2hss(i)
            rs2crit2hss(i)=((n2-1)*(n1(i)-1)*ps1*rs1*sigtotau)./(n1(i)*rs1crit2hss(i));
        else
            rs2crit2hss(i)=((n2-1)*ps1*(n1(i)-1)*sigtotau)./n1(i);
        end 
 
   
   %case 1 n 2
    if rs2<=rs2crit2hss(i)
        sigrat2hss(i) = rs1./(n2.*ps1.*sigtotau);
    %case 3
    elseif rs2>rs2crit2hss(i) && rs1<=rs1crit2hss(i)
        sigrat2hss(i)=((n2-1).*rs1)./(n1(i).*n2.*sigtotau);
    %case 4
    else
        sigrat2hss(i)=((n1(i)-1).*(n2-1))./(n1(i).*n2);
    end
    wrat2hss(i)=sigrat2hss(i).^2/Erat2hss(i);
    end
    plot (n1,sigrat2hss,'linewidth',1.5)
    hold on;
end

%  title('SW within Regular')
xlabel ('n_1', 'fontsize', 18,'fontweight','bold')
ylabel ('\sigma_{SS}/\sigma_p \phi_{1} \phi_{2}', 'fontsize', 18,'fontweight','bold')
legend ('n_2 = 2', 'n_2 = 3', 'n_2 = 4','n_2 = 5','n_2 = 6','n_2 = 7','n_2 = 8','n_2 = 9','n_2 = 10','fontsize', 12);
pbaspect([5 4 1])
set(gca,'fontsize',14)
set(gcf,'color','w')
figure
subplot(1,2,1)
for n2=[2,3,4,5,6,7,8,9,10]
    for i=1:numel(n1)
        as(i)=(ps1.*rs1.*rs1.*gbe)./(3*(1-ps1));
        D(i)=(((n1(i).*((3*n1(i))-4))./(3*((n1(i)-1).^2)))+(((n1(i).*n1(i))./(3*(n1(i)-1).*as(i)))));
        a_ss(i)=(ps2.*rs2.*rs2.*gbe*D(i))./(3*(1-ps2)*ps1);
        Erat2hss(i)=(D(i).*(((n2.*((3*n2)-4))./(3*((n2-1).^2)))+(((n2.*n2)./(3*(n2-1).*a_ss(i)))))).^-1;
        Ess(i)=Erat2hss(i)*Ep*ps1*ps2;
        rs1crit2hss(i)=(n1(i)-1)*sigtotau;
        if rs1<=rs1crit2hss(i)
            rs2crit2hss(i)=((n2-1)*(n1(i)-1)*ps1*rs1*sigtotau)./(n1(i)*rs1crit2hss(i));
        else
            rs2crit2hss(i)=((n2-1)*ps1*(n1(i)-1)*sigtotau)./n1(i);
        end 
        %case 1 n 2
        if rs2<=rs2crit2hss(i)
            sigrat2hss(i) = rs1./(n2.*ps1.*sigtotau);
        %case 3
        elseif rs2>rs2crit2hss(i) && rs1<=rs1crit2hss(i)
            sigrat2hss(i)=((n2-1).*rs1)./(n1(i).*n2.*sigtotau);
    %case 4
        else
            sigrat2hss(i)=((n1(i)-1).*(n2-1))./(n1(i).*n2);
        end
        erat2hss(i)=sigrat2hss(i)/Erat2hss(i);
    end
    plot (n1,erat2hss,'linewidth',1.5);
    hold on;
end


%  title('SW within Regular')
xlabel ('n_1', 'fontsize', 18,'fontweight','bold')
ylabel ('\epsilon_{SS}/\epsilon_p', 'fontsize', 18,'fontweight','bold')
legend ('n_2 = 2', 'n_2 = 3', 'n_2 = 4','n_2 = 5','n_2 = 6','n_2 = 7','n_2 = 8','n_2 = 9','n_2 = 10','fontsize', 12);
pbaspect([5 4 1])
set(gca,'fontsize',14)
set(gcf,'color','w')

subplot(1,2,2)
for n2=[2,3,4,5,6,7,8,9,10]
    for i=1:numel(n1)
        as(i)=(ps1.*rs1.*rs1.*gbe)./(3*(1-ps1));
        D(i)=(((n1(i).*((3*n1(i))-4))./(3*((n1(i)-1).^2)))+(((n1(i).*n1(i))./(3*(n1(i)-1).*as(i)))));
        a_ss(i)=(ps2.*rs2.*rs2.*gbe*D(i))./(3*(1-ps2)*ps1);
        Erat2hss(i)=(D(i).*(((n2.*((3*n2)-4))./(3*((n2-1).^2)))+(((n2.*n2)./(3*(n2-1).*a_ss(i)))))).^-1;
        Ess(i)=Erat2hss(i)*Ep*ps1*ps2;
        rs1crit2hss(i)=(n1(i)-1)*sigtotau;
        if rs1<=rs1crit2hss(i)
            rs2crit2hss(i)=((n2-1)*(n1(i)-1)*ps1*rs1*sigtotau)./(n1(i)*rs1crit2hss(i));
        else
            rs2crit2hss(i)=((n2-1)*ps1*(n1(i)-1)*sigtotau)./n1(i);
        end 
    %case 1 n 2
    if rs2<=rs2crit2hss(i)
        sigrat2hss(i) = rs1./(n2.*ps1.*sigtotau);
   
    %case 3
    elseif rs2>rs2crit2hss(i) && rs1<=rs1crit2hss(i)
        sigrat2hss(i)=((n2-1).*rs1)./(n1(i).*n2.*sigtotau);
    %case 4
    else
        sigrat2hss(i)=((n1(i)-1).*(n2-1))./(n1(i).*n2);
    end
    wrat2hss(i)=sigrat2hss(i).^2/Erat2hss(i);
    end
    plot (n1,wrat2hss,'linewidth',1.5)
    hold on;
end
xlabel ('n_1', 'fontsize', 18,'fontweight','bold')
ylabel ('w_{SS}/w_p \phi_{1} \phi_{2}', 'fontsize', 18,'fontweight','bold')
legend ('n_2 = 2', 'n_2 = 3', 'n_2 = 4','n_2 = 5','n_2 = 6','n_2 = 7','n_2 = 8','n_2 = 9','n_2 = 10','fontsize', 12);
pbaspect([5 4 1])
set(gca,'fontsize',14)
set(gcf,'color','w')

