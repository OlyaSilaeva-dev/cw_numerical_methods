import numpy as np
import matplotlib.pyplot as plt

# Чтение данных из файла
data = np.loadtxt('output.txt')
x = np.linspace(-1, 1, len(data))

rho = data[:, 0]
u = data[:, 1]
p = data[:, 2]

# Построение графиков
plt.figure(figsize=(10, 6))
# plt.subplot(3, 1, 1)
plt.plot(x, rho, label="Плотность")
plt.title('Плотность')
plt.xlabel('x')
plt.ylabel('rho')

# plt.subplot(3, 1, 2)
plt.plot(x, u, label="Скорость", color='g')
# plt.title('Скорость')
# plt.xlabel('x')
plt.ylabel('u')

# plt.subplot(3, 1, 3)
plt.plot(x, p, label="Давление", color='r')
# plt.title('Давление')
# plt.xlabel('x')
plt.ylabel('p')

plt.tight_layout()
plt.show()
