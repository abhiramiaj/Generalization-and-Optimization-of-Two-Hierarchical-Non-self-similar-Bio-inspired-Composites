clear all
clc
syms Ep sigpcrit taumcrit Ep Em num ps1 ps2 n1 n2 m n k
%%%%%%%%%%MP1%%%%%%%%%%

Gm=Em/(2*(1+num));
sigtotau=sigpcrit/taumcrit;
gbe=Gm/Ep;
gam_m_crit=taumcrit/Gm;
 
%%%%------------2HSS-----------%%%%
%ps1--volume fraction at first level of hierarchy
%ps2--volume fraction at second level of hierarchy
%rs1--aspect ratio rho at first level of hierarchy
%rs2--aspect ratio rho at second level of hierarchy
%n1--n at first level of hierarchy
%n2--n at second level of hierarchy
%as--alpha at first level of hierarchy
%ass--alpha at second level of hierarchy
rs1crit2hss=(n1-1)*sigtotau;
rs1=rs1crit2hss;
rs2crit2hss=((n2-1)*ps1*(n1-1)*sigtotau)./n1;
rs2=rs2crit2hss;
as=(ps1.*rs1.*rs1.*gbe)./(3*(1-ps1));
D=(((n1.*((3*n1)-4))./(3*((n1-1).^2)))+(((n1.*n1)./(3*(n1-1).*as))));
ass=(ps2.*rs2.*rs2.*gbe*D)./(3*(1-ps2)*ps1);
Erat2hss=(D.*(((n2.*((3*n2)-4))./(3*((n2-1).^2)))+(((n2.*n2)./(3*(n2-1).*ass))))).^-1;
Ess=Erat2hss*Ep*ps1*ps2;
E_norm=Ess/Ep;

sigrat2hss=((n1-1).*(n2-1))./(n1*n2); %case 4
sig_norm=sigrat2hss*ps1*ps2;

wrat2hss=sigrat2hss.^2/Erat2hss;
wpcrit=sigpcrit^2/(2*Ep);
w_norm=wrat2hss*wpcrit*ps1*ps2/(sigpcrit*gam_m_crit)

%fc=(Erat2hss^m)*(sigrat2hss^n)*(wrat2hss^k);%%%%%%%%%This is the fitness function
fc=(E_norm^m)*(sig_norm^n)*(w_norm^k);%%%%%%%%%This is the modified fitness function 
%refering barthelat2014 stiffness and strength normalized with the platelet
%stiffness and strength, respectively; toughness normalized with the
%toughness of a hypothetical material with normal strength of platelet and
%shear strain of matrix
fc=char(fc);
fc=strrep(fc, 'ps1','x(1)'); %first modification, replacing
%ps1 with x(1)
fc=strrep(fc, 'ps2','x(2)');
%second modification, replacing
%ps2 with x(2)
fc=strrep(fc, 'n1','x(3)');
%third modification, replacing
%n1 with x(3)
fc=strrep(fc, 'n2','x(4)');
%fourth modification, replacing
%n2 with x(4)
fc=strrep(fc, '^m','^m(i)');
fc=strrep(fc, '^n','^n(j)');
fc=strrep(fc, '^k','^k(q)');
%modifications to change m,n,k to m(i),n(j) and k(q)
% Ep=1200;
% Em=0.67; 
% num=0.4;
% sigmcrit=1;
% sigpcrit=40;
% taumcrit=0.576;