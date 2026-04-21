function R=R_coh(theta,sigma,lambda);
% Corrected 12.december 2013
R=exp(-2*(2*pi*sigma*sin(theta)./lambda).^2);
%R=exp(-2*(2*pi*sigma*sin(theta)./lambda.^2));




