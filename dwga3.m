function neff = dwga3(d,lam0,ns,nf,nc,mode,r,tol)
% Effecive indices of TE or TM modes in asymmetric 3-layer dielectric waveguide
% d, lam0  = thickness and lambda (same units)
% ns,nf,nc = n values of substrate, core, cover      
% mode     = 'TE' or 'TM' (default TE)
% r        = relaxation parameter (default r=0.5)
% tol      = error tolerance (default tol=1e-6)

if nargin <5, help dwga3; return; end
if nargin<=7, tol = 1e-6; end
if nargin<=6, r = 0.5; end
if nargin<=5, mode = 'TE'; end

k0a = (2*pi/lam0)*d/2;      % k0*d/2
R = k0a*sqrt(nf^2-ns^2); d = (ns^2-nc^2)/(nf^2-ns^2);  % V & asymmetry
if strcmpi(mode,'TE'), ps = 1; pc = 1;
    else  ps =(nf/ns)^2; pc =(nf/nc)^2; 
end
   
M = floor((2*R-atan(pc*sqrt(d)))/pi); m = (0:M)'; % M+1 modes                            % mode indices
u = R*ones(M+1,1); v = zeros(M+1,1); w = R*sqrt(d)*ones(M+1,1); % initialize
Nit = 1;
while 1
   unew = r*(m*pi/2 + atan(ps*v./u)/2 + atan(pc*w./u)/2) + (1-r)*u;
   if norm(unew-u) <= tol, break; end
   Nit = Nit + 1; u = unew; 
   v = sqrt(R^2 - u.^2); w = sqrt(R^2*d + v.^2);
   if Nit > 1000, 
       fprintf('failed to converge with r = 0.5: try a smaller r'); break; 
   end
end
neff = sqrt(nf^2 - (u/k0a).^2);           % Eff. index beta/k0



