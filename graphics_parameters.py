import matplotlib.pyplot as plt
import pandas as pd

# Чтение данных из файла
filename = 'output_values.txt'  # Укажите путь к вашему файлу

# Читаем файл с помощью pandas
data = pd.read_csv(filename, sep=" ", header=None, names=["rho", "u", "p", "gamma", "rho_0", "c_0"])

# Построение графиков на одном графике
plt.figure(figsize=(10, 6))

# Рисуем каждую переменную разным цветом
plt.plot(data.index, data['rho'], label='rho', color='blue')
plt.plot(data.index, data['u'], label='u', color='green')
plt.plot(data.index, data['p'], label='p', color='red')
plt.plot(data.index, data['gamma'], label='gamma', color='orange')
# plt.plot(data.index, data['rho_0'], label='rho_0', color='purple')
# plt.plot(data.index, data['c_0'], label='c_0', color='brown')

# Добавляем заголовок и метки осей
plt.title('Parameters vs Index')
plt.xlabel('X')
plt.ylabel('Value')

# Добавляем легенду
plt.legend()

# Показать график
plt.tight_layout()
plt.savefig('graphics_parameters.png')
