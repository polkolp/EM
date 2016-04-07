clear all; clc; close all;
%%
waveln = 1.55; d_core = 0.8; % both in same units
nc = 1; nf = 3.5; ns = 1.5; % core and cladding ref. indices
neff = dwga3(d_core,waveln,ns,nf,nc);

%% plot modes
a = d_core/2;
kf = (2*pi/waveln)*(sqrt(nf^2-neff.^2));
ac = (2*pi/waveln)*(sqrt(neff.^2-nc^2));
as = (2*pi/waveln)*(sqrt(neff.^2-ns^2));
x = linspace(-d_core, d_core,500);
mode = zeros(length(x));

for i = 1: length(neff)
    figure();
    if mod(i,2) == 1
        mode(x >= a) = cos(kf(i)*a)*exp(ac(i)*a)*exp(-ac(i)*x(x >= a));
        mode(x <= -a) = cos(kf(i)*a)*exp(as(i)*a)*exp(as(i)*x(x <= -a));
        mode(x >= -a & x <= a) = cos(kf(i)*x(x >= -a & x <= a));
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
        mode(x >= a) = sin(kf(i)*a)*exp(ac(i)*a)*exp(-ac(i)*x(x >= a));
        mode(x <= -a) = -sin(kf(i)*a)*exp(as(i)*a)*exp(as(i)*x(x <= -a));
        mode(x >= -a & x <= a) = sin(kf(i)*x(x >= -a & x <= a));
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
        
            
            