#!/usr/local/bin/python

import numpy as np
from numpy import array as arr
from numpy import zeros,sin,cos,exp
from sys import argv

From:float  = 0
To:float    = 1
alpha0:float = 1
alpha1:float = 0
beta0:float	= 0
beta1:float	= 1
A:float 	= 1
B:float 	= 2

Pi:list[float] = []
Qi:list[float] = []
Fi:list[float] = []

N: int

def _H (_N: int) -> float:
	return (To - From) / (_N-1)

def p (n: float) -> float:
	 return 2
def q ( n: float) -> float:
	return 0
def f ( n: float) -> float:
	return 4 * exp(n) * ( sin(n) + cos(n) )


def ASolution(x: float) -> float:
    C1 = exp(3)*(8*sin(1) + 4*cos(1))/10 - exp(2)
    C = 2*exp(1)/5 - C1
    return exp(x)*(6*sin(x) - 2*cos(x)) / 5+ C1*exp(-2*x) + C

def Dy( xi: int,  Y: list[float],  H: float, N: int) -> float:
    if (xi == 0):
        return (-Y[2] + 4*Y[1] - 3*Y[0]) / (2*H)
    elif (xi == N):
        return (3*Y[N] - 4*Y[N-1] + Y[N-2]) / (2*H)
    elif (xi < N and xi > 0):
        return ( Y[xi+1] - Y[xi-1] ) / ( 2* H)
    else:
        return 0

def DDy( xi: int, Y: list[float],  H: float, N: int) -> float:
    if (xi == 0):
        return ( -Dy(2,Y,H,N ) + 4*Dy(1,Y,H,N ) - 3 *Dy(0,Y,H,N ) ) / (2*H)
    elif (xi == N):
        return ( Dy(N,Y,H,N ) - 4*Dy(N-1,Y,H,N ) + 3 *Dy(N-2,Y,H,N ) ) / (2*H)
    elif (xi < N and xi > 0):
        return ( Dy(xi+1, Y, H,N) - Dy(xi-1, Y, H,N)) / ( 2* H)
    else:
        return 0;	

def lval( xi: int,  Y: list[float],  H: float, N: int) -> float:
    return DDy(xi,Y,H,N) + Dy(xi,Y,H,N)*Pi[xi] + Y[xi]*Qi[xi];

def iter(N: int) -> arr:
	H = _H(N)
	X = [ x * _H(N) + From for x in range(N) ]
	a: arr = zeros(N)
	b: arr = zeros(N+1)
	c: arr = zeros(N)
	d: arr = zeros(N+1)
	Ai: arr = zeros(N+1)
	Bi: arr = zeros(N+1)
	Yi: arr = zeros(N+1)
	for i in range(1, N):
		a[i-1] = 1 - (H * Pi[i])/2
		b[i] = H*H*Qi[i] - 2
		c[i] = 1+ (H*Pi[i])/2
		d[i] = H*H*Fi[i]
	b[0] = H*alpha0 - alpha1
	c[0] = alpha1
	d[0] = A*H
	a[N-1] = -beta1
	b[N] = H*beta0 + beta1
	d[N] = B*H
	Ai[0] = -c[0]/b[0]
	Bi[0] = d[0] / b[0]
	for i in range(1,N):
		Ai[i] = -c[i] / (b[i] + a[i-1]*Ai[i-1])
		Bi[i] = (d[i] - a[i-1]*Bi[i-1] ) / ( b[i] + a[i-1]*Ai[i-1] )
	Ai[N] = 0
	Bi[N] = (d[N] - a[N-1]*Bi[N-1] ) / ( b[N] + a[N-1]*Ai[N-1] )
	Yi[N] = Bi[N]
	i = N-1
	while i >= 0:
		Yi[i] = Ai[i]*Yi[i+1] + Bi[i]
		i = i - 1
	return Yi

def regen_grid(_N: int) -> list[float]:
	global Pi, Qi, Fi, grid, N
	N = _N
	grid = [ x * _H(_N) + From for x in range(_N) ]
	Pi = []
	Qi = []
	Fi = []
	for x in grid:
		Pi.append(p(x))
		Qi.append(q(x))
		Fi.append(f(x))
	return grid

def save_csv(_N: int, filename: str):
	N = _N
	regen_grid(N)
	Yi = []
	Yi = iter(N)

	with open(filename, 'w') as ff:
		for i, x in enumerate(grid):
			ff.write(
				f"{x:+02.6f}\t"+
				f"{Yi[i]:+02.6f}\t"+
				f"{ASolution(grid[i]):+02.6f}\t"+
				f"{Fi[i]:+02.6f}\t"+
				f"{lval(i,Yi,_H(N),N):+02.6f}\n"
			)

if __name__ == '__main__' :
	if (len(argv) > 1):
		save_csv( int(argv[1]), 'data.csv' )
	else:
		save_csv( 2, 'data.csv' )
