import os as os
import random as rn

cdir = os.getcwd()

exp_list = [1000, 5000, 10000, 25000, 50000, 100000]
n_items = len(exp_list)
file_id = 0

for k in range(0, n_items):
    file_id = file_id + 1
    iterations = exp_list[k]
    random_file = f'ran{file_id}.dat'
    for i in range(1, iterations + 1):
        print(f'iterations = {i}')
        with open(f'{cdir}/{random_file}', 'a') as random_data:
            random_value = rn.uniform(-1, 1)
            random_value = '{0:.8f}'.format(random_value)
            random_data.write(f'{random_value} \n')

