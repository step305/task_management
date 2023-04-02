%%
function p= make_covariance_ellipses(x,P)
% compute ellipses for plotting state covariances
N= 50;
inc= 2*pi/N;
phi= 0:inc:2*pi;

p = make_ellipse(x(1:2), P(1:2,1:2), 2, phi);

end