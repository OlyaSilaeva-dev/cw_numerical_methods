import numpy as np
import matplotlib.pyplot as plt

# Чтение коэффициентов D из файла
with open('output.txt', 'r') as file:
    coefficients = np.array([float(value) for value in file.read().split()])

# Убедимся, что коэффициентов ровно 5
if len(coefficients) != 5:
    raise ValueError("Ожидается 5 коэффициентов в файле: D_L, D_L*, U, D_R*, D_R")

# Распределяем коэффициенты
D_L, D_L_star, U, D_R_star, D_R = coefficients

def positive_only(y):
    return np.where(y > 0, y, np.nan)
# Создание диапазона x
x = np.linspace(-1, 1, 500)

# Построение графика
plt.figure(figsize=(10, 6))

# Построение прямых для каждого коэффициента
plt.plot(x, positive_only(x / D_L), label='D_L')
plt.plot(x, positive_only(x / D_L_star), label='D_L*')
plt.plot(x, positive_only(x / U), label='U', linestyle='dashed')
plt.plot(x, positive_only(x / D_R_star), label='D_R*')
plt.plot(x, positive_only(x / D_R), label='D_R')

# Настройки графика
plt.xlabel('x')
plt.ylabel('t')
plt.title('Прямые t = D * x для каждого коэффициента D')
plt.legend()
plt.grid()

# Показать график
plt.savefig('plot.png')
