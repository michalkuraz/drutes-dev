a=8.66424
n=3.44613
m=1-1/n
ts=1
tr=0

set xrange [-10:0]

C(x) = a*m*n*(-tr + ts)*(-(a*x))**(-1 + n)*(1 + (-(a*x))**n)**(-1 - m)

th(x) = 1/(1+(-a*x)**n)**m
K(x) = (1 - (-(a*x))**(m*n)/(1 + (-(a*x))**n)**m)**2*(1 + (-(a*x))**n)**(m/2.0)

set grid 
set y2tics

plot th(x) axes x1y1, K(x) axes x1y1, C(x) axes x1y2
