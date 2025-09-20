import math
import random
import numpy as np
import scipy.stats as sps
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import chi2
n =15
p =0.564
q = 1 - p
lamb = 2.1

def frequancy(data):
    counts = Counter(data)
    n_arr = [counts.get(i, 0) for i in range(min(data), max(data) + 1)]
    return n_arr
def relative_frequancy(n_arr,size = 200):
    return [count/size for count in n_arr]
def element_is_zero1(frequancy1):
    for k in range(len(frequancy1)):
        if frequancy1[k] == 0:
            return True
    return False
def element_is_zero2(arr1, arr2):  # Переименовал для ясности
  """Проверяет, есть ли нулевые элементы хотя бы в одном из массивов."""
  if not arr1 or not arr2:
    return True # Считаем, что если один из массивов пустой, условие выполняется
  
  min_len = min(len(arr1), len(arr2))
  for i in range(min_len):
    if arr1[i] == 0 or arr2[i] == 0:
      return True
  return False
def pad_array(arr, min_val, max_val, n):
   
    for _ in range(min_val):
        arr.insert(0, 0)
    for _ in range(n - max_val -1): # -1 потому что множество не может содержать дубликаты
        arr.append(0)  # append быстрее, чем insert(len(arr), ...)
    return arr

def generate_ksi_binom(n, p, size=200):
    ksi = []
    P = np.zeros(n + 1)
    # Вычисление вероятностей для биномиального распределения
    P[0] = (1 - p) ** n  # Вероятность 0 успехов
    for i in range(1, n + 1):
        P[i] = P[i - 1] * (n - (i - 1)) * p / (i * (1 - p))
    for _ in range(size):
        k = 0
        R = P[0]
        alpha = random.random()
        Q = alpha - R
        while Q>0:
            R = R*(n-k)*p/((k+1)*(1-p))
            Q=Q-R
            k=k+1
        ksi.append(k)
    return ksi, P
s = 0
S_arr_binom = []
n_arr1_binom = []
n_arr2_binom = []
while element_is_zero2(n_arr1_binom, n_arr2_binom) or len(n_arr1_binom) == 0 or len(n_arr1_binom) != len(n_arr2_binom):
    n_arr1_binom = []  # Обнуляем перед каждым заполнением
    n_arr2_binom = []  # Обнуляем перед каждым заполнением
    ksi_binom, P_binom = generate_ksi_binom(n, p)
    ksi_binom.sort()
    n_arr1_binom = frequancy(ksi_binom)
    data_binom = sps.binom.rvs(n, p, size=200)
    data_binom.sort()
    n_arr2_binom = frequancy(data_binom)
print("n_arr1_binom:", n_arr1_binom)
print("n_arr2_binom:", n_arr2_binom)
x_i_binom = np.array(set(ksi_binom))
# Преобразование частот в относительные частоты
relative_frequencies1_binom = relative_frequancy(n_arr1_binom)
relative_frequencies2_binom = relative_frequancy(n_arr2_binom)
    
print("Сгенерированные значения:", ksi_binom)
print("Сгенерированные значения (scipy):", data_binom)
print(f"P_binom ={P_binom} ")
for i in range(len(P_binom)):
    s+=P_binom[i]
    S_arr_binom.append(s)
print(f"sum_P = {s}")
print(f"S_arr_binom = {S_arr_binom}")
print(f"частоты САМ биномиального распределения{n_arr1_binom}")
print(f"частоты пакетного биномиального распределения{n_arr2_binom}")
var_sum = 0
srednee = (sum(ksi_binom)/200)
print(f" среднее_биномиальное = {srednee:.6f}")
for i in ksi_binom:
    var_sum +=(i-srednee)**2
var_sum/=199
print(f" дисперсия = {var_sum}")
print(f"относительные частоты САМ биномиального распределения{relative_frequencies1_binom}")
print(f"относительные частоты пакетного биномиального распределения{relative_frequencies2_binom}")
print(x_i_binom)
print(len(relative_frequencies1_binom))
print(len(relative_frequencies2_binom))
# Подготовка данных для построения полигонов
x1 = [k for k in range(len(relative_frequencies1_binom))]
y1 = [s for s in relative_frequencies1_binom]
x2 = [k for k in range(len(relative_frequencies2_binom))]
y2 = [s for s in relative_frequencies2_binom]
N_abs_sq_div_p = [(200*(relative_frequencies1_binom[i]-P_binom[i])**2)/P_binom[i] for i in range(len(relative_frequencies1_binom))]
N_abs_sq_div_p = np.array(N_abs_sq_div_p,dtype=float)
abs_wi_p = [abs(relative_frequencies1_binom[i]-P_binom[i]) for i in range(len(relative_frequencies1_binom))]
abs_wi_p = np.array(abs_wi_p,dtype=float)
print(f"abs_wi_p = {abs_wi_p}")
print(f"N_abs_sq_div_p = {N_abs_sq_div_p}")
x_prob = np.arange(0, n + 1)
y_prob = P_binom
teor_binom_n = np.array([n_w * 200 for n_w in y_prob], dtype=float)
# Построение полигонов относительных частот
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)  # 1 row, 2 columns, first plot
plt.plot(x1, y1, marker='o', linestyle='-', label='Выборка (САМДР)')
plt.plot(x2, y2, marker='x', linestyle='--', label='Выборка (scipy)')
plt.xlabel('Значение (k)')
plt.ylabel('Относительная частота')
plt.title('Полигоны относительных частот')
plt.legend()
plt.grid(True)
plt.subplot(1, 2, 2)  # 1 row, 2 columns, second plot
plt.plot(x_prob, y_prob, marker='o', linestyle='-', color='blue', label='Теоретические вероятности')
plt.xlabel('Значение (k)')
plt.ylabel('Вероятность P(k)')
plt.title('Полигон вероятностей')
plt.legend()
plt.grid(True)
plt.tight_layout()  # Предотвращает перекрывание графиков
plt.show()


min_x1 = min(ksi_binom) 
max_x1 = max(ksi_binom) 
n_arr1_binom = pad_array(n_arr1_binom, min_x1, max_x1, n+1)
min_x2 = min(data_binom) 
max_x2 = max(data_binom) 
n_arr2_binom = pad_array(n_arr2_binom, min_x2, max_x2, n+1)
print(n_arr1_binom)
print(n_arr2_binom)
relative_frequencies1_binom = relative_frequancy(n_arr1_binom)
relative_frequencies2_binom = relative_frequancy(n_arr2_binom)
N_abs_sq_div_p = [(200*(relative_frequencies1_binom[i]-P_binom[i])**2)/P_binom[i] for i in range(len(relative_frequencies1_binom))]
N_abs_sq_div_p = np.array(N_abs_sq_div_p,dtype=float)
abs_wi_p = [abs(relative_frequencies1_binom[i]-P_binom[i]) for i in range(len(relative_frequencies1_binom))]
abs_wi_p = np.array(abs_wi_p,dtype=float)
var_sum = 0
srednee = (sum(ksi_binom)/200)
print(f" среднее_биномиальное = {srednee:.6f}")
for i in ksi_binom:
    var_sum +=(i-srednee)**2
var_sum/=199
print(f" дисперсия = {var_sum}")
for i in abs_wi_p:
    k = round(i,5)
    print(k)
for i in N_abs_sq_div_p:
    k = round(i,5)
    print(k)
print(sum(N_abs_sq_div_p))


N_abs_sq_div_p = [(200*(relative_frequencies2_binom[i]-P_binom[i])**2)/P_binom[i] for i in range(len(relative_frequencies2_binom))]
N_abs_sq_div_p = np.array(N_abs_sq_div_p,dtype=float)
abs_wi_p = [abs(relative_frequencies2_binom[i]-P_binom[i]) for i in range(len(relative_frequencies2_binom))]
abs_wi_p = np.array(abs_wi_p,dtype=float)
var_sum = 0
srednee = (sum(data_binom)/200)
print(f" среднее_биномиальное = {srednee:.6f}")
for i in data_binom:
    var_sum +=(i-srednee)**2
var_sum/=199
print(f" дисперсия = {var_sum}")
for i in abs_wi_p:
    k = round(i,5)
    print(k)
for i in N_abs_sq_div_p:
    k = round(i,5)
    print(k)
print(sum(N_abs_sq_div_p))


k = 0
wi12 = []
for a,b in zip(relative_frequencies1_binom,relative_frequencies2_binom):
    if a==0 and b==0:
        k=0
    else:
        k=(a**2 + b**2)/(a+b)
    print(k)
    wi12.append(k)
print(400*(sum(wi12)-1))


def generate_ksi_geom(p, size=200):
    ksi = []
    
    # Генерируем значения ksi
    for _ in range(size):
        k = 0
        # Генерируем случайные числа до первой удачи
        while random.random() >= p:  # Пока случайное число > p (неудача)
            k += 1
        ksi.append(k + 1)  # Добавляем 1 для учета первой удачи
    max_ksi = max(ksi)
    
    # Расчет вероятности P для всех возможных значений k
    P = np.array([p * (1 - p) ** i for i in range(max_ksi)])
    return ksi, P
s = 0
S_arr_geom = []
n_arr1_geom=[]
n_arr2_geom=[]
while len(n_arr1_geom) != len(n_arr2_geom) or len(n_arr1_geom) == 0 or element_is_zero1(n_arr1_geom) or element_is_zero1(n_arr2_geom):
    ksi_geom,P_geom = generate_ksi_geom(p)
    ksi_geom.sort()
    data_geom = sps.geom.rvs(p,size = 200)
    data_geom.sort()
    n_arr1_geom = frequancy(ksi_geom)
    n_arr2_geom = frequancy(data_geom)
    relative_frequencies1_geom = relative_frequancy(n_arr1_geom)
    relative_frequencies2_geom = relative_frequancy(n_arr2_geom)
print(f"частоты САМ геометрического распределения{n_arr1_geom}")
print(f"частоты пакетного геометрического распределения{n_arr2_geom}")
# Преобразование частот в относительные частоты
print(f"относительные частоты САМ биномгеометрического  распределения{relative_frequencies1_geom}")
print(f"относительные частоты пакетного геометрического распределения{relative_frequencies2_geom}")
for i in range(len(ksi_geom)):
    ksi_geom[i] -= 1
    data_geom[i] -= 1
print(ksi_geom)
print(data_geom)
S_arr_geom = []
s = 0
var_sum = 0
print("\nВероятности P(k) ( округленные до 5 знаков):")
for i, prob in enumerate(P_geom):
    print(f"P_geom({i}) = {round(prob, 6)}")
    s=round(s+prob,5)
    S_arr_geom.append(s)
print(f"сумма вероятностей = {round(s,6)}")
x_i_geom =list(set(ksi_geom))
print(f"x_i_geom = {len(x_i_geom)}")
print(len(relative_frequencies1_geom))
print(len(relative_frequencies2_geom))
srednee = (np.sum(data_geom)/200)
print(f" среднее_геометрического = {srednee:.6f}")
for i in data_geom:
    var_sum +=(i-srednee)**2
var_sum/=199
print(f" дисперсия = {var_sum}")
# Подготовка данных для построения полигонов
x1 = [k for k in range(len(relative_frequencies1_geom))]
y1 = [s for s in relative_frequencies1_geom]
x2 = [k for k in range(len(relative_frequencies2_geom))]
y2 = [s for s in relative_frequencies2_geom]
N_abs_sq_div_p = [(200*(relative_frequencies1_geom[i]-P_geom[i])**2)/P_geom[i] for i in range(len(relative_frequencies1_geom))]
N_abs_sq_div_p = np.array(N_abs_sq_div_p,dtype=float)
abs_wi_p = [abs(relative_frequencies1_geom[i]-P_geom[i]) for i in range(len(relative_frequencies1_geom))]
abs_wi_p = np.array(abs_wi_p,dtype=float)
print(f"abs_wi_p = {abs_wi_p}")
print(f"N_abs_sq_div_p = {N_abs_sq_div_p}")
teor_geom_n = []
x_prob = np.arange(0, len(x1))
y_prob = P_geom / np.sum(P_geom)
# Количество необходимых значений
n_w = 200
# Вычисляем teor_geom_n
teor_geom_n = n_w * y_prob
# Проверка значений в teor_geom_n:
print("Проверка значений teor_geom_n:")
for i, val in enumerate(teor_geom_n):
    print(f"Index: {i}, Value: {val}")
# Проверка суммы
print("Сумма teor_geom_n:", np.sum(teor_geom_n))
# Построение полигонов относительных частот
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)  # 1 row, 2 columns, first plot
plt.plot(x1, y1, marker='o', linestyle='-', label='Выборка (САМДР)')
plt.plot(x2, y2, marker='x', linestyle='--', label='Выборка (scipy)')
plt.xlabel('Значение (k)')
plt.ylabel('Относительная частота')
plt.title('Полигоны относительных частот')
plt.legend()
plt.grid(True)
plt.subplot(1, 2, 2)  # 1 row, 2 columns, second plot
plt.plot(x_prob, y_prob, marker='o', linestyle='-', color='green', label='Теоретические вероятности')
plt.xlabel('Значение (k)')
plt.ylabel('Вероятность P(k)')
plt.title('Полигон вероятностей')
plt.legend()
plt.grid(True)
plt.tight_layout()  # Предотвращает перекрывание графиков
plt.show()
print(len(ksi_geom))
print(len(data_geom))


var_sum = 0
srednee = (np.sum(ksi_geom)/200)
print(f" среднее_пауссоновского = {srednee:.6f}")
for i in ksi_geom:
    var_sum +=(i-srednee)**2
var_sum/=199
print(f" дисперсия = {var_sum}")
N_abs_sq_div_p = [(200*(relative_frequencies1_geom[i]-P_geom[i])**2)/P_geom[i] for i in range(len(relative_frequencies1_geom))]
N_abs_sq_div_p = np.array(N_abs_sq_div_p,dtype=float)
abs_wi_p = [abs(relative_frequencies1_geom[i]-P_geom[i]) for i in range(len(relative_frequencies1_geom))]
abs_wi_p = np.array(abs_wi_p,dtype=float)
for i in abs_wi_p:
    print(round(i,5))
for i in N_abs_sq_div_p:
    print(round(i,5))
print(np.sum(N_abs_sq_div_p))

var_sum = 0
srednee = (np.sum(data_geom)/200)
print(f" среднее_пауссоновского = {srednee:.6f}")
for i in data_geom:
    var_sum +=(i-srednee)**2
var_sum/=199
print(f" дисперсия = {var_sum}")
k = 0
summa = 0
N_abs_sq_div_p = [(200*(relative_frequencies2_geom[i]-P_geom[i])**2)/P_geom[i] for i in range(len(relative_frequencies2_geom))]
N_abs_sq_div_p = np.array(N_abs_sq_div_p,dtype=float)
abs_wi_p = [abs(relative_frequencies2_geom[i]-P_geom[i]) for i in range(len(relative_frequencies2_geom))]
abs_wi_p = np.array(abs_wi_p,dtype=float)
for i in abs_wi_p:
    print(round(i,5))
for i in N_abs_sq_div_p:
    print(round(i,5))
t=0
t = np.sum(N_abs_sq_div_p)
print(round(t,5))

for i in range(len(relative_frequencies2_geom)):
    k = (relative_frequencies1_geom[i]**2 + relative_frequencies2_geom[i]**2)/(relative_frequencies1_geom[i] + relative_frequencies2_geom[i])
    summa +=k
    print(round(k,5))
print(400*(summa-1))

def generate_ksi_poisson(lambd, size=200):
    ksi = []
    # Генерируем значения ksi
    for _ in range(size):
        k = 0
        p = math.exp(-lambd)  # Начальная вероятность P(X=0)
        alpha = random.random()
        while alpha > p: # Пока случайное число > p (неудача)
            k += 1
            alpha *= random.random()  # Вычисляем P(X=k) и сравниваем с alpha
        ksi.append(k)
    max_ksi = max(ksi)
    P = np.zeros(max_ksi + 1)
    P[0] = math.exp(-lambd)
    # Расчет вероятности P для всех возможных значений k
    for i in range(max_ksi):
        P[i+1] = P[i] * lambd / (i+1)
    return ksi, P
s = 0
S_arr_poisson = []
n_arr1_poisson=[]
n_arr2_poisson=[]
while len(n_arr1_poisson) != len(n_arr2_poisson) or len(n_arr1_poisson) == 0 :
    ksi_poisson,P_poisson = generate_ksi_poisson(lamb)
    ksi_poisson.sort()
    data_poisson = sps.poisson.rvs(lamb,size = 200)
    data_poisson.sort()
    n_arr1_poisson = frequancy(ksi_poisson)
    n_arr2_poisson = frequancy(data_poisson)
    relative_frequencies1_poisson = relative_frequancy(n_arr1_poisson)
    relative_frequencies2_poisson = relative_frequancy(n_arr2_poisson)
s=0
print(ksi_poisson)
print(data_poisson)
print("\nВероятности P(k) ( округленные до 5 знаков):")
for i, prob in enumerate(P_poisson):
    print(f"P_poisson({i}) = {round(prob, 6)}")
    s=round(s+prob,5)
    S_arr_poisson.append(s)
print(f"сумма вероятностей = {round(s,6)}")
print(f"частоты САМ пауссоновского распределения{n_arr1_poisson}")
print(f"частоты пакетного пауссоновского распределения{n_arr2_poisson}")
# Преобразование частот в относительные частоты
print(f"относительные частоты САМ пауссоновского распределения{relative_frequencies1_poisson}")
print(f"относительные частоты пакетного пауссоновского распределения{relative_frequencies2_poisson}")
x_i_poisson =list(set(ksi_poisson))
print(f"x_i_geom = {len(x_i_poisson)}")
print(len(relative_frequencies1_poisson))
print(len(relative_frequencies2_poisson))
var_sum = 0
srednee = (sum(ksi_poisson)/200)
print(f" среднее_пауссоновского = {srednee:.6f}")
for i in ksi_poisson:
    var_sum +=(i-srednee)**2
var_sum/=199
print(f" дисперсия")
# Подготовка данных для построения полигонов
x1 = [k for k in range(len(relative_frequencies1_poisson))]
y1 = [s for s in relative_frequencies1_poisson]
x2 = [k for k in range(len(relative_frequencies2_poisson))]
y2 = [s for s in relative_frequencies2_poisson]
N_abs_sq_div_p = [(200*(relative_frequencies1_poisson[i]-P_poisson[i])**2)/P_poisson[i] for i in range(len(relative_frequencies1_poisson))]
N_abs_sq_div_p = np.array(N_abs_sq_div_p,dtype=float)
abs_wi_p = [abs(relative_frequencies1_poisson[i]-P_poisson[i]) for i in range(len(relative_frequencies1_poisson))]
abs_wi_p = np.array(abs_wi_p,dtype=float)
print(f"abs_wi_p = {abs_wi_p}")
print(f"N_abs_sq_div_p = {N_abs_sq_div_p}")
teor_poisson_n = []
x_prob = np.arange(0, len(x1))
y_prob = P_poisson / np.sum(P_poisson)
# Количество необходимых значений
n_w = 200
# Вычисляем teor_geom_n
teor_poisson_n = n_w * y_prob
# Проверка значений в teor_geom_n:
print("Проверка значений teor_poisson_n:")
for i, val in enumerate(teor_poisson_n):
    print(f"Index: {i}, Value: {val}")
# Проверка суммы
print("Сумма teor_poisson_n:", np.sum(teor_poisson_n))
# Построение полигонов относительных частот
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)  # 1 row, 2 columns, first plot
plt.plot(x1, y1, marker='o', linestyle='-', label='Выборка (САМДР)')
plt.plot(x2, y2, marker='x', linestyle='--', label='Выборка (scipy)')
plt.xlabel('Значение (k)')
plt.ylabel('Относительная частота')
plt.title('Полигоны относительных частот')
plt.legend()
plt.grid(True)
plt.subplot(1, 2, 2)  # 1 row, 2 columns, second plot
plt.plot(x_prob, y_prob, marker='o', linestyle='-', color='red', label='Теоретические вероятности')
plt.xlabel('Значение (k)')
plt.ylabel('Вероятность P(k)')
plt.title('Полигон вероятностей')
plt.legend()
plt.grid(True)
plt.tight_layout()  # Предотвращает перекрывание графиков
plt.show()
print(len(ksi_poisson))
print(len(data_poisson))


var_sum = 0
srednee = (np.sum(ksi_poisson)/200)
print(f" среднее_пауссоновского = {srednee:.6f}")
for i in ksi_poisson:
    var_sum +=(i-srednee)**2
var_sum/=199
print(f" дисперсия = {var_sum}")
k = 0
summa = 0
N_abs_sq_div_p = [(200*(relative_frequencies1_poisson[i]-P_poisson[i])**2)/P_poisson[i] for i in range(len(relative_frequencies1_poisson))]
N_abs_sq_div_p = np.array(N_abs_sq_div_p,dtype=float)
abs_wi_p = [abs(relative_frequencies1_poisson[i]-P_poisson[i]) for i in range(len(relative_frequencies1_poisson))]
abs_wi_p = np.array(abs_wi_p,dtype=float)
for i in abs_wi_p:
    print(round(i,5))
for i in N_abs_sq_div_p:
    print(i)
print(sum(N_abs_sq_div_p))

var_sum = 0
srednee = (np.sum(data_poisson)/200)
print(f" среднее_пауссоновского = {srednee:.6f}")
for i in data_poisson:
    var_sum +=(i-srednee)**2
var_sum/=199
print(f" дисперсия = {var_sum}")
N_abs_sq_div_p = [(200*(relative_frequencies2_poisson[i]-P_poisson[i])**2)/P_poisson[i] for i in range(len(relative_frequencies2_poisson))]
N_abs_sq_div_p = np.array(N_abs_sq_div_p,dtype=float)
abs_wi_p = [abs(relative_frequencies2_poisson[i]-P_poisson[i]) for i in range(len(relative_frequencies2_poisson))]
abs_wi_p = np.array(abs_wi_p,dtype=float)
for i in abs_wi_p:
    print(round(i,5))
for i in N_abs_sq_div_p:
    print(i)
print(sum(N_abs_sq_div_p))

for i in range(len(relative_frequencies1_poisson)):
    k = (relative_frequencies1_poisson[i]**2 + relative_frequencies2_poisson[i]**2)/(relative_frequencies1_poisson[i] + relative_frequencies2_poisson[i])
    summa +=k
    print(k)
print(400*(summa-1))

