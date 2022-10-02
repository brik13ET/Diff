import csv
import matplotlib.pyplot as plt
import sys

def draw_pic(filename: str, out: str):
	with open(filename, newline='\n') as f:
		ff = csv.reader(f,delimiter='\t')
		X = []
		Y = []
		F = []
		R = []
		S = []
		E = []
		for row in ff:
			X.append(float(row[0]))
			Y.append(float(row[1]))
			F.append(float(row[2]))
			R.append(float(row[3]))
			S.append(float(row[4]))
			E.append(float(row[5]))
		plt.switch_backend('Agg')
		plt.rcParams["figure.figsize"] = (15,24)
		fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1)
		ax1.plot(X, F, label='Дано')
		ax1.plot(X, R, label='Решение')
		ax2.plot(X, Y, label='Прогонка')
		ax2.plot(X, S, label='Решение')
		ax3.plot(X, [ abs(F[i] - R[i]) for i, _ in enumerate(X) ], label='Delta Рунге')
		
		ax4.plot(X, [ abs(S[i] - Y[i]) for i, _ in enumerate(X) ], label='Delta Абсолютное')

		ax1.legend()
		ax2.legend()
		ax3.legend()
		fig.savefig(out)



if __name__ == '__main__':
    draw_pic(sys.argv[1], 'fig.svg')