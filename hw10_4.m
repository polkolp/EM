clear all; clc; close all;
%%
%The program is created by Professor Agrawal
% Effective indices of TE modes of a symmetric dielectric waveguide
waveln = 1.55; d_core = 0.8; % both in same units
n1 = 3.5; n2 = 1.5; % core and cladding ref. indices
k0a = (2*pi/waveln)*d_core/2; % k0*a (dimensionless)
R = k0a*sqrt(n1^2 - n2^2); % R parameter sets the number of modes
u = dslab(R); % Newton¡¯s method for finding eigenvalues
neff = sqrt(n1^2 - (u./k0a).^2); % effective mode indices

%% plot modes
a = d_core/2;
kc = (2*pi/waveln)*(sqrt(n1^2-neff.^2));
ac = (2*pi/waveln)*(sqrt(neff.^2-n2^2));
x = linspace(-d_core, d_core,500);
mode = zeros(length(x));

for i = 1: length(neff)
    figure();
    if mod(i,2) == 1
        mode(x >= a) = cos(kc(i)*a)*exp(ac(i)*a)*exp(-ac(i)*x(x >= a));
        mode(x <= -a) = cos(kc(i)*a)*exp(ac(i)*a)*exp(ac(i)*x(x <= -a));
        mode(x >= -a & x <= a) = cos(kc(i)*x(x >= -a & x <= a));
        plot(x,mode);
        xlabel('thickness [nm]');
        ylabel('E/E_1');
        str = sprintf('neff = %f',neff(i));
        title(str);
        y1=get(gca,'ylim');
        hold on
        plot([-a -a],y1)
        plot([a a],y1)
        hold off
    else
        mode(x >= a) = sin(kc(i)*a)*exp(ac(i)*a)*exp(-ac(i)*x(x >= a));
        mode(x <= -a) = -sin(kc(i)*a)*exp(ac(i)*a)*exp(ac(i)*x(x <= -a));
        mode(x >= -a & x <= a) = sin(kc(i)*x(x >= -a & x <= a));
        plot(x,mode);
        xlabel('thickness [nm]');
        ylabel('E/E_1');
        str = sprintf('neff = %f',neff(i));
        title(str);
        y1=get(gca,'ylim');
        hold on
        plot([-a -a],y1)
        plot([a a],y1)
        hold off
    end
end
        
            
            