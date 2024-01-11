import numpy as np
import os
import matplotlib.pyplot as plt
import math
import re

def pBar_exact(xx, tt):
    return (tt/4.0)*xx*(1-xx)

def Sw_exact(xx, tt):
    return 0.5 + (tt/4.0)*xx*(1-xx)

def extract_float_from_filename(file_name):
    # Extract float values from file names
    return float(file_name.split('_')[2][:-4])

def custom_sort(s):
    # Extract the integer part from the string, or use 0 if there are no integers
    integer_part = int(''.join(filter(str.isdigit, s))) if any(char.isdigit() for char in s) else 0
    return (integer_part, s)

def extract_float_from_filename_t(file_name):
    # Use regular expression to extract the float part
    match = re.search(r'\d+\.\d+', file_name)
    if match:
        return float(match.group())
    else:
        return None

x = np.linspace(0.0,1.0,1024)


# PRESSURE 

file_list = os.listdir("./")

filtered_files = [filename for filename in file_list if filename.startswith('data_pBar')]
filtered_files = sorted(filtered_files, key=extract_float_from_filename) 

print(filtered_files)
times = [extract_float_from_filename_t(file_name) for file_name in filtered_files]
print(times)


for i in range(len(filtered_files)):

    file = filtered_files[i]

    x_projL2 = np.genfromtxt(file, delimiter=',', max_rows=1)
    u_projL2 = np.genfromtxt(file, delimiter=',', skip_header=1, max_rows=1)

    plt.figure()
    plt.grid(True)

    u = pBar_exact(x, times[i])
    plt.plot(x,u,c = "b",label="$u(x,t) exata$")
    plt.scatter(x_projL2, u_projL2, c="r", marker='o', label="MEF", zorder=2) # linestyle='--'
    
    plt.title(file)
    plt.xlabel("x")
    plt.ylabel("u")
    plt.ylim([-0.03,0.08])
    plt.legend()
    plt.savefig('plot'+file[:-4]+str(times[i])+'.jpeg', dpi=300, bbox_inches='tight')
    plt.close()






# SATURATION 

file_list = os.listdir("./")

filtered_files = [filename for filename in file_list if filename.startswith('data_Sw')]
filtered_files = sorted(filtered_files, key=extract_float_from_filename) 

print(filtered_files)
times = [extract_float_from_filename_t(file_name) for file_name in filtered_files]
print(times)


for i in range(len(filtered_files)):

    file = filtered_files[i]

    x_projL2 = np.genfromtxt(file, delimiter=',', max_rows=1)
    u_projL2 = np.genfromtxt(file, delimiter=',', skip_header=1, max_rows=1)

    plt.figure()
    plt.grid(True)

    u = Sw_exact(x, times[i])
    plt.plot(x,u,c = "b",label="$u(x,t) exata$")
    plt.scatter(x_projL2, u_projL2, c="r", marker='o', label="MEF", zorder=2) # linestyle='--'
    
    plt.title(file)
    plt.xlabel("x")
    plt.ylabel("u")
    plt.ylim([0.4, 0.8])
    plt.legend()
    plt.savefig('plot'+file[:-4]+str(times[i])+'.jpeg', dpi=300, bbox_inches='tight')
    plt.close()