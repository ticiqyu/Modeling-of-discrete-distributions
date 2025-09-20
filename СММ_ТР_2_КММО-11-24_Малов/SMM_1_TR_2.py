import numpy as np
import matplotlib.pyplot as plt
import math
import random
import scipy.stats
#8 вариант
lamb = 0.81


def f_exp(x):
    return lamb*math.exp(-(lamb*x))
def my_exp(x):
    return 1-math.exp(-lamb*x)
def reverse_exp(y):
    return round(-math.log(1-y)/(lamb),8)

def generate_ksi_exp(size = 200):
    ksi = []
    for _ in range(size):
        alpha = random.random()
        ksi.append(round(reverse_exp(alpha),5))
    return ksi

random.seed(10)
ksi_exp = generate_ksi_exp()
ksi_exp = np.array(ksi_exp)
output_ksi = "\n".join(map('{:.5f}'.format, ksi_exp))
print(f"Сгенерированная выборка САМНР: \n {output_ksi}")
ksi_exp.sort()
print(ksi_exp)
output_ksi = "\n".join(map('{:.5f}'.format, ksi_exp))
print(f"\n Отсортированная выборка САМНР: \n {output_ksi}\n")
ksi_exp_E_x = sum(ksi_exp)/200
print(f"Матожидание выборки САМНР:  {ksi_exp_E_x}")
squares = [(x-ksi_exp_E_x)**2 for x in ksi_exp]
ksi_exp_D_x = sum(squares)/199
print(f"Дисперсия выборки САМНР: {ksi_exp_D_x }")


np.random.seed(12)
data_exp = np.random.exponential(1/lamb,200)
output_data = "\n".join(map('{:.6f}'.format, data_exp))
print(f"Сгенерированная выборка numpy: \n {output_data}")
data_exp.sort()
output_data = "\n".join(map('{:.6f}'.format, data_exp))
print(f"\n Отсортированная выборка numpy: \n {output_data}\n")
data_exp_E_x = sum(data_exp)/200
print(f"Матожидание выборки САМНР: {data_exp_E_x}")
squares = [(x-data_exp_E_x)**2 for x in data_exp]
data_exp_D_x = sum(squares)/199
print(f"Дисперсия выборки САМНР: {data_exp_D_x}")

import numpy as np
import math

def my_exp(x, lamb=0.81):
    return 1 - math.exp(-lamb * x)

def compute_intervals(data, m):
    a_0 = 0
    a_m = max(data)
    interval = (a_m - a_0) / m
    a_i = [round(a_0 + i * interval, 5) for i in range(1, m + 1)]
    return a_i

def compute_frequencies(data, a_i):
    n_i = np.zeros(len(a_i))
    for x in data:
        i = 0
        while i < len(a_i) - 1 and x > a_i[i]:
            i += 1
        n_i[i] += 1
    w_i = [round(n_i[i] / len(data), 5) for i in range(len(n_i))]
    return n_i, w_i

def compute_probabilities(a_i, func):
    prob = []
    prev = 0
    for x in a_i:
        current = func(x)
        prob.append(round(current - prev, 5))
        prev = current
    return prob

# Пример данных (замените на ваши ksi_exp и data_exp)

m = 8

print("KSI_EXP")
a_i_ksi = compute_intervals(ksi_exp, m)
n_i_ksi, w_i_ksi = compute_frequencies(ksi_exp, a_i_ksi)
print("Интервалы:", a_i_ksi)
print("Частоты n_i:", n_i_ksi)
print("Относительные частоты w_i:", w_i_ksi, "\n")

print("DATA_EXP")
a_i_data = compute_intervals(data_exp, m)
n_i_data, w_i_data = compute_frequencies(data_exp, a_i_data)
print("Интервалы:", a_i_data)
print("Частоты n_i:", n_i_data)
print("Относительные частоты w_i:", w_i_data, "\n")

prob_1 = compute_probabilities(a_i_ksi, my_exp)
print("prob_1 =", prob_1)
print("Сумма prob_1 =", sum(prob_1))

prob_2 = compute_probabilities(a_i_data, my_exp)
exp_1 = [round(my_exp(x), 5) for x in a_i_data]
print("exp_1 =", exp_1)
print("prob_2 =", prob_2)
print("Сумма prob_2 =", sum(prob_2))



ksi_exp_dif = [round(w_i_ksi[i]-prob_1[i],5) for i in range(len(w_i_ksi))]
print(ksi_exp_dif)
data_exp_dif = [round(w_i_data[i]-prob_2[i],5) for i in range(len(w_i_data))]
print(data_exp_dif)

hi_ksi = [round(200*ksi_exp_dif[i]**2/prob_1[i],5) for i in range(len(ksi_exp_dif))]
print(hi_ksi)
hi_data = [round(200*data_exp_dif[i]**2/prob_2[i],5) for i in range(len(data_exp_dif))]
print(hi_data)

print(f"Хи-квадрат САМНР {sum(hi_ksi)}")
print(f"Хи-квадрат numpy {sum(hi_data)}")
scipy.stats.chi2.ppf(0.95,m-1)

plt.bar([x for x in range(len(w_i_ksi))],w_i_ksi,color = 'red',label = 'Относительные частоты от обратной функции')
plt.plot([x for x in range(len(prob_1))],prob_1,label = 'Теоретическое распределение плотностей')

plt.title("Гистограмма отн. частот обратной функции")
plt.xlabel("i")
plt.ylabel("w_i_ksi")
plt.grid()
plt.legend()
plt.show()

plt.bar([x for x in range(len(w_i_ksi))],w_i_ksi,color = 'Yellow',label = 'Относительные частоты от библ. функции')
plt.plot([x for x in range(len(prob_2))],prob_2,label = 'Теоретическое распределение плотностей')

plt.title("Гистограмма отн. частот np.random.exponential")
plt.xlabel("i")
plt.ylabel("w_i_data")
plt.grid()
plt.legend()
plt.show()

for a in a_i_ksi:
    print(round(my_exp(a),5))

for a in a_i_data:
    print(round(my_exp(a),5))

