// C variant

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <stdint.h>

const double
	From = 0,
	To = 1,
	alpha0 = 1,
	alpha1 = 0,
	beta0 = 0,
	beta1 = 1,
	A = 1,
	B = 2;
uint32_t N = 2;

double p (double n) { return 2;}
double q (double n) { return 0;}
double f (double n) { return 4 * exp(n) * ( sin(n) + cos(n) );}

int32_t iter(double);
int32_t grid(double, double**);
double	rval(uint32_t, double);

double* x_iter; // dynamic
double* p_iter;
double* q_iter;
double* f_iter;
double* y_iter;

void main(void)
{
	double H = (To - From) / N;
	printf("H=%lf\tFrom=%lf\tTo=%lf\nAlloc()\n",H,From,To);
	p_iter = malloc(sizeof(double) * N);
	q_iter = malloc(N*sizeof(double));
	f_iter = malloc(N*sizeof(double));
	y_iter = malloc((N+1)*sizeof(double));

	grid(H, &x_iter);
	for (uint32_t x = 0; x < N; x++)
	{
		p_iter[x] = p(x_iter[x]);
		q_iter[x] = q(x_iter[x]);
		f_iter[x] = f(x_iter[x]);
	}
	printf("Iter()\n");
	iter(H);
	print("Saving\n");
	FILE* fd = fopen("data.csv", "w+");
	for (uint32_t i = 0; i < N; i++)
		fprintf(fd, "%2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\n",
			x_iter[i],
			y_iter[i],
			f_iter[i],
			rval(i, H) 
		);
	printf("Cleanup\n");
	fclose(fd);
	free(p_iter);
	free(q_iter);
	free(f_iter);
	free(y_iter);
	free(x_iter);
}

int32_t grid(double H, double** buf)
{
	*buf = malloc(sizeof(double) * N);
	for (uint32_t i = 0; i < N; i ++)
		(*buf)[i] = From + H * i;
	return N;
}

int32_t iter(double H)
{
	double
		*a = malloc(N),
		*b = malloc(N+1),
		*c = malloc(N),
		*d = malloc(N+1),
		*A_iter = malloc(N+1),
		*B_iter = malloc(N+1);
	
	b[0] = H*alpha0 - alpha1;
	c[0] = alpha1;
	d[0] = A*H;
	for (uint32_t i = 1; i < N - 1; i ++)
	{
		a[i-1] = 1 - (p_iter[i] * H) / 2;
		b[i] = H*H*q_iter[i] - 2;
		c[i] = 1 + H * p_iter[i] / 2;
		d[i] = H*H*f_iter[i];
	}
	a[N-1] = -beta1;
	b[N] = H*beta0 + beta1;
	d[N] = B*H;
	A_iter[0] = -c[0] / b[0];
	B_iter[0] = d[0] / b[0];
	for (uint32_t i = 1; i < N - 1; i ++)
	{
		A_iter[i] = -c[i] / (b[i] + a[i-1]*A_iter[i-1]);
		B_iter[i] = (d[i] - a[i-1]*B_iter[i-1]) / (b[i] + a[i-1]*A_iter[i-1]);
	}
	A_iter[N] = 0;
	B_iter[N] = (d[N] - a[N-1]*B_iter[N-1]) / (b[N] + a[N-1]*A_iter[N-1]);
	y_iter[N] = B_iter[N];
	for (uint32_t i = N - 1; i > 0; i --)
	{
		y_iter[i] = A_iter[i]*y_iter[i+1] + B_iter[i-1];
	}
	free(a);
	free(b);
	free(c);
	free(d);
	free(A_iter);
	Free(B_iter);
	return N;
}

double rval(uint32_t i, double H)
{
	if (i == 0)
		return  \
		( - ( y_iter[2-1] + y_iter[2+1]) \
		+ 4 * ( y_iter[1-1] + y_iter[1+1])  \
		- 3 * (-y_iter[2] + 4 * y_iter[1]  \
		- 3 * y_iter[0]) ) / (4*H*H);
	else if (i == 1)
		return ( ( -y_iter[2] + 4 * y_iter[1] - 3 * y_iter[0]) + (y_iter[i+1-1] + y_iter[i+1+1]) ) / (4*H*H) ;
	else if (i == N-1)
		return ( ( y_iter[i-1-1] + y_iter[i-1+1]) +( 3 * y_iter[N] - 4 * y_iter[N-1] + y_iter[N-2]) ) / (4*H*H) ;
	else if (i == N)
		return ( \
		3*	( 3 * y_iter[N] - 4 * y_iter[N-1] + y_iter[N-2]) \
		-4*	(y_iter[N-1-1] + y_iter[N-1+1]) \
		+	(y_iter[N-2-1] + y_iter[N-2+1]) ) / (4*H*H);
	else
		return ( (y_iter[i-1-1] + y_iter[i-1+1]) + (y_iter[i+1-1] + y_iter[i+1+1]) ) / (4*H*H) ;
}