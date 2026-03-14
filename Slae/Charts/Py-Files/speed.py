import matplotlib.pyplot as plt
import numpy as np
alpha = []
t_dense = []
t_sparse = []
with open("build/data.txt", "r") as file:
    for line in file:
        s = line.split(" ")
        alpha.append(float(s[0]))
        t_dense.append(float(s[1]) * 1e6)
        t_sparse.append(float(s[2]) * 1e6)
x = np.array(alpha)
list = [t_dense, t_sparse]
labels = ['Dense Matrix', 'Sparse Matrix']

plt.figure(figsize=(12, 8))
for i, (data, labelss) in enumerate(zip(list, labels)):
    plt.plot(x, data, label=labelss)
    plt.grid(True)
    plt.title("w = Ax для разреженной и плотной квадратных матриц")
    plt.xlabel('alpha = n/N (доля ненулевых элементов), N = 1000 (размер матрицы)')
    plt.ylabel('average time, mus')
    plt.legend()
plt.show()