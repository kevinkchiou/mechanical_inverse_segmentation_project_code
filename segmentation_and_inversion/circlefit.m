function [r x y]=circlefit(v,override)
%computes radius and center of best fitting circle to a set of data v, with
%v(:,1)=x(:), v(:,2)=y(:) the points to fit.  will probably converge faster
%if points are in sequential ordering.  entered as circlefit(v).

if nargin<2
    override=0;
end

sz=size(v);
if(override==0)
    if(max(sz)<7)
        r=[];x=[];y=[];
        return;
    end
end

%pick 3 points to create seed circle
%r1=[v(1,1) v(1,2)];
%r2=[v(ceil(len/2),1) v(ceil(len/2),2)];
%r3=[v(len,1) v(len,2)];
%[r x y]=seedcircle(r1,r2,r3);

%pick one of these.
%[x y r]=circfit(v(:,1),v(:,2));
if length(v)==3
    [x y r]=circfit(v(:,1),v(:,2));
elseif length(v)>3
    par=PrattSVD(v);x=par(1);y=par(2);r=par(3);
end

%plotcircle(r,x,y,v);%function for testing results


function [r x y]=seedcircle(r1,r2,r3)
%creates seed circle uniquely given the three points r1,r2,r3. returns
%radius, and center of circle.

%rotation by pi/2 matrix
rot=[[0 -1];[1 0];];

%midpoints of chords
g1=(r2-r1)*0.5+r1;
g2=(r3-r2)*0.5+r2;

m1=(r2-r1)*rot;
m2=(r3-r2)*rot;

m1xm2=crosstwo(m1,m2);
if(m1xm2==0)
    error('we may have a problem...\n');
end

t2=crosstwo(m1,g1-g2)/crosstwo(m1,m2);

pos=m2*t2+g2;


rad1=sqrt((r1-pos)*(r1-pos)');
rad2=sqrt((r2-pos)*(r2-pos)');
rad3=sqrt((r3-pos)*(r3-pos)');
eps=0.0000000001*(rad1+rad2+rad3)/3;
if(abs(rad1-rad2)<eps && abs(rad2-rad3)<eps && abs(rad1-rad3)<eps)
    r=(rad1+rad2+rad3)/3;
end
x=pos(1);y=pos(2);

function plotcircle(r0,x0,y0,v)
%function for testing the accuracy of circlefit().

%{
t=linspace(0,2*pi,200);
x=r0*cos(t)+x0;
y=r0*sin(t)+y0;
plot(x,y);hold on;
%}
plot(v(:,1),v(:,2),'rs');

function a=crosstwo(x1,x2)
%two dimensional cross product

a=(x1(1)*x2(2)-x1(2)*x2(1));

%Kasa method. biased to small circles for small arc data
function   [xc,yc,R,a] = circfit(x,y)
%
%   [xc yx R] = circfit(x,y)
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991, 
   x=x(:); y=y(:);
   a=[x y ones(size(x))]\[-(x.^2+y.^2)];
   xc = -.5*a(1);
   yc = -.5*a(2);
   R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));

%this function is stolen from Nikolai Chernov's website. Pratt svd method.
function Par = PrattSVD(XY)

%--------------------------------------------------------------------------
%  
%     Algebraic circle fit by Pratt
%      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
%      Computer Graphics, Vol. 21, pages 145-152 (1987)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this is a version optimized for stability, not for speed
%
%--------------------------------------------------------------------------

centroid = mean(XY);   % the centroid of the data set

[U,S,V]=svd([(XY(:,1)-centroid(1)).^2+(XY(:,2)-centroid(2)).^2,...
        XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)],0);

if (S(4,4)/S(1,1) < 1e-12)   %  singular case
    A = V(:,4);
    disp('Pratt singular case');
else                         %  regular case
    W=V*S;
    Binv = [0 0 0 -0.5; 0 1 0 0; 0 0 1 0; -0.5 0 0 0];
    [E,D] = eig(W'*Binv*W);
    [Dsort,ID] = sort(diag(D));
    A = E(:,ID(2));
    for i=1:4
        S(i,i)=1/S(i,i); 
    end
    A = V*S*A;
end

Par = -(A(2:3))'/A(1)/2 + centroid;
Par = [Par , sqrt(A(2)^2+A(3)^2-4*A(1)*A(4))/abs(A(1))/2];


