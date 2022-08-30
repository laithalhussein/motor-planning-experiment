function confidence_ellipse_95(covmatrix,mu,num,plot_color,line_wid)
% CONFIDENCE_ELLIPSE: This function plots the error ellipse associated with
% standard errors of a 2-D distribution.  It requires the covariance
% matrix (covmatrix), the mean of the distribution (mu), and the number of
% subjects (num).  Optional inputs are Color (in RGB format) and LineWidth.

n=100;
theta=(0:n)/n*2*pi;

if nargin < 4
    plot_color = [0 0 0];
    line_wid = 0.5;
elseif nargin < 5
    line_wid = 0.5;
end
    
    


x=[cos(theta);sin(theta)];

% Compute quantile for the desired percentile
k = sqrt(qchisq(1-exp(-1/2),2)); % 2 is the number of dimensions (degrees of freedom)

[eigvec, eigval] = eig(covmatrix);
output = 1.96*x'*sqrt(eigval)*eigvec';
xvals=k*output(:,1);
yvals=k*output(:,2);


plot(mu(1)+xvals,mu(2)+yvals,'Color',plot_color,'LineWidth',line_wid)

function x=qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.
if nargin<2
  n=1;
end

s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;

for ii=1:14
  dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
  x(s) = x(s)+dx;
  if all(abs(dx) < 1e-6)
    break;
  end
end

function F=pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.
if nargin<2
  n=1;
end
F=zeros(size(x));

if rem(n,2) == 0
  s = x>0;
  k = 0;
  for jj = 0:n/2-1;
    k = k + (x(s)/2).^jj/factorial(jj);
  end
  F(s) = 1-exp(-x(s)/2).*k;
else
  for ii=1:numel(x)
    if x(ii) > 0
      F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
    else
      F(ii) = 0;
    end
  end
end

function f=dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
if nargin<2
  n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));