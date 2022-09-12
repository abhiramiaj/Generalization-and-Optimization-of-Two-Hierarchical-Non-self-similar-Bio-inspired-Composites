
function [a,b] = Hessiancritptsolve()
clear all
options = optimset('Display','off');
fun=@critpt;
x0=[0.01,0.01];
[a,b]=fsolve(fun,x0,options);
end

%x1 tauibar, x2 phi
function F = critpt(x)
F=[(33*2^(13/25)*x(2)*(-x(1)*(x(2) - 1))^(53/100)*(-x(2)*(6229007507042843*x(1) - 4503599627370496))^(7/50)*((10141204801825835211973625643008*3^(1/2)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)^2*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + ((6229007507042843*x(1) + 4503599627370496)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1)))^2 + 3*coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))^2 - 4))/(8*x(1)*(6229007507042843*x(1) - 4503599627370496)) + (1125899906842624*3^(1/2)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1))))*(3*x(2) - 3))/(3*(6229007507042843*x(1) - 4503599627370496)*(x(2) - 1)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))))/(51200*((2251799813685248*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + 1)^2*(x(2)/((2251799813685248*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + 1))^(67/100)) - (53*2^(13/25)*(x(2)/4 - 1/4)*(x(2)/((2251799813685248*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + 1))^(33/100)*(-x(2)*(6229007507042843*x(1) - 4503599627370496))^(7/50))/(12800*(-x(1)*(x(2) - 1))^(47/100)) - (43603052549299901*2^(13/25)*x(2)*(x(2)/((2251799813685248*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + 1))^(33/100)*(-x(1)*(x(2) - 1))^(53/100))/(25600*(-x(2)*(6229007507042843*x(1) - 4503599627370496))^(43/50)),(33*2^(13/25)*(1/((2251799813685248*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + 1) - (x(2)*((coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1)))^2 + 3*coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))^2 - 4)/(8*x(2)*(x(2) - 1)) + (1125899906842624*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1))))*(3*x(2) - 3))/(3*x(2)*(6229007507042843*x(1) - 4503599627370496)*(x(2) - 1)^2*(-(x(1)*x(2))/(x(2) - 1))^(1/2))))/((2251799813685248*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + 1)^2)*(-x(1)*(x(2) - 1))^(53/100)*(-x(2)*(6229007507042843*x(1) - 4503599627370496))^(7/50))/(51200*(x(2)/((2251799813685248*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + 1))^(67/100)) - (53*2^(13/25)*x(1)*(x(2)/((2251799813685248*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + 1))^(33/100)*(-x(2)*(6229007507042843*x(1) - 4503599627370496))^(7/50))/(51200*(-x(1)*(x(2) - 1))^(47/100)) - (7*2^(13/25)*(6229007507042843*x(1) - 4503599627370496)*(x(2)/((2251799813685248*3^(1/2)*x(1)*(coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(27021597764222976*x(1))) + coth((3^(1/2)*(6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2))/(9007199254740992*x(1)))))/((6229007507042843*x(1) - 4503599627370496)*(-(x(1)*x(2))/(x(2) - 1))^(1/2)) + 1))^(33/100)*(-x(1)*(x(2) - 1))^(53/100))/(25600*(-x(2)*(6229007507042843*x(1) - 4503599627370496))^(43/50))];
end