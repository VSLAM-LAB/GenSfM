function plotc(x0,r,n)

tt = linspace(0,2*pi,n);

plot(x0(1)+r*cos(tt),x0(2)+r*sin(tt),'r-')
