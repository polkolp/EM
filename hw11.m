clear all; clc; close all;
%%
la0 = 650;
ef = 1.5^2;
ec = -19.6+0.44*1i; es = ec;
a = 100;
be0 = 1.1*sqrt(ef);

[be,Err] = pwga(la0,ef,ec,es,a,be0)

%%
ef = -19.6+0.44*1i; 
ec = 1.5^2; es = ec;

[be,Err] = pwga(la0,ef,ec,es,a,be0)