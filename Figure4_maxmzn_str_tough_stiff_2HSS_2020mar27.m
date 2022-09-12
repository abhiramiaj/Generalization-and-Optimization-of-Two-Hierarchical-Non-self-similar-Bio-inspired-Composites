
figure
clc
Ep=7e12;
Em=2.086e10;
num=0.49;
Gm=Em/(2*(1+num));
gbe=Gm/Ep;
sigtotau=10;
n1=2:1:10;  %swm 1H
    %n2=10;
ps1=0.5; %swm 1H
ps2=0.5; %swm 2H
%z2=0.5; %rsm 2H
% for n2=2:1:10
%     
%for n2 =[2 3 4 5 6 7 8 9 10]
    figure
    for n2 =2
for m=1:numel(n1)
%%%%------------2HSS-----------%%%%
rs1=1:1:120;
[as,D,rs2crit2hss,A,Z,Y]=deal(zeros(size(rs1)));
   for i=1:numel(rs1)
   rs2 = 1:1:120;
   [ass,Erat2hss,E2hss,sigrat2hss,wrat2hss]=deal(zeros(size(rs2)));
   for j=1:numel(rs2)
as(i)=(ps1.*rs1(i).*rs1(i).*gbe)./(3*(1-ps1));
D(i)=(((n1(m).*((3*n1(m))-4))./(3*((n1(m)-1).^2)))+(((n1(m).*n1(m))./(3*(n1(m)-1).*as(i)))));
ass(j)=(ps2.*rs2(j).*rs2(j).*gbe.*D(i))./(ps1.*(1-ps2).*3);
Erat2hss(j)=(D(i).*(((n2.*((3*n2)-4))./(3*((n2-1).^2)))+(((n2.*n2)./(3*(n2-1).*ass(j)))))).^-1;
E2hss(j)=Erat2hss(j)*Ep*ps1*ps2;
rs1crit2hss=(n1(m)-1)*sigtotau;
if rs1(i)<=rs1crit2hss
       rs2crit2hss(i)=((n2-1)*ps1*rs1(i))./(n1(m));
        else
        rs2crit2hss(i)=((n2-1)*ps1*(n1(m)-1)*sigtotau)./n1(m);
end
    %case 1 n 2
    if rs2(j)<=rs2crit2hss(i)
    sigrat2hss(j) = rs2(j)./(n2*ps1*sigtotau);
%     wrat2hrs= rr2(j)*rr2(j)./(2*2*2*ps2*ps2*sigtotau*sigtotau*E2hrs)
    %case 3
    elseif rs2(j)>rs2crit2hss(j) && rs1(i)<=rs1crit2hss
   sigrat2hss(j)=((n2-1)*rs1(i))./(n1(m).*n2*sigtotau);
%    wrat2hrs=(rs1(i)*rs1(i))./(2*2*n1(m)*n1(m)*sigtotau*sigtotau*E2hrs)
    %case 4
    else
        sigrat2hss(j)=((n1(m)-1)*(n2-1))./(n1(m).*n2);
%         wrat2hrs=((n1(m)-1).^2)/(2*2*2*n1(m)*n1(m)*E2hrs)
    end
        wrat2hss(j)=sigrat2hss(j).^2/Erat2hss(j);
A(i,j)=wrat2hss(j);
Z(i,j)=Erat2hss(j);
Y(i,j)=sigrat2hss(j);
   end
   end
 maxtough(m)=max(A(:))
   [rs1_maxtough, rr2_maxtough] = find(ismember(A, max(A(:))))
   maxstrength(m)=max(Y(:))
 [rs1_maxstrength, rr2_maxstrength] = find(ismember(Y, max(Y(:))));
   Q=rs1_maxstrength;
   R=rr2_maxstrength;
   strengthatmaxtough(m)=Y(rs1_maxtough,rr2_maxtough)
  stiffnessatmaxtough(m)=Z(rs1_maxtough,rr2_maxtough)
end
%yyaxis left
%plot (n1,maxtough,'linewidth',1.5)
% title('2hss maximum toughness and corresponding strength vs n')
plot(n1,maxtough,'d', 'MarkerSize',8,'Color','Black','MarkerFaceColor','Black')
addaxis(n1,strengthatmaxtough,'h','MarkerSize',8,'Color','Blue','MarkerFaceColor','Blue')
%set(gca,'fontsize',16)
addaxis(n1,stiffnessatmaxtough,'^','MarkerSize',8,'Color','Magenta','MarkerFaceColor','Magenta')
legend({'w^{SS, max}_{critical}/\phi w^{p}_{critical}','\sigma^{SS}_{critical}/\phi\sigma^{p}_{critical}',['E^{SS}_{critical}/\phi\sigma^{p}_{critical}; n_2 = ' num2str(n2)]},'FontSize',14,'fontweight','bold')
addaxislabel(1,'w^{SS, max}_{critical}/\phi w^{p}_{critical}');
addaxislabel(2,'\sigma^{SS}_{critical}/\phi\sigma^{p}_{critical}');
addaxislabel(3,'E^{SS}_{critical}/\phi\sigma^{p}_{critical}');
xlabel('n_1','FontSize',18,'fontweight','bold')
pbaspect([5 4 1])
%legend(['n_1 =' num2str(n1)])
set(gcf,'color','w')
end
%For text formatting, edit addaxislabel file

%yyaxis right
%plot(n1,strengthatmaxtough,'linewidth',1.5)



%set(gca,'fontsize',16)
% end