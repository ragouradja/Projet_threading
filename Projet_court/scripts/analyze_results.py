__authors__ = "Ragousandirane RADJASANDIRANE"
__contact__ = "radja.ragou@gmail.com"
__date__ = "14/09/2021"
__version__= "1.0"

# Modules import
from multiprocessing import cpu_count
import re
import matplotlib.pyplot as plt

def extract_cpu():
    return re.compile(r"^# Nb (\w|\s+)+ : (\d)")
    
def extract_id():
    return re.compile(r"^# (Final [\w|\s]+) : (\d+)%")


def extract_zscore():
    return re.compile(r"^# (Zscore) \((\d+) \w+\): (-?\d.\d)")

def extract_hmatrix():
    return re.compile(r"^# (Last [\w|\s]+) : (\d+.\d+)")

def extract_time():
    return re.compile(r"^# (\w* time : )(\d+.\d+)s")

def extract_blosum():
    return re.compile(r"^# Blosum Matrix")

def extract_weights():
    return re.compile(r"^#( Weights \w+ ): (\d.\d) \w+ (\d.\d) \w+")
    
def extract_dssp():
    return re.compile(r"^# Secondary [\w|\s]+")



def extract_results(file_log):
    
    id = extract_id()
    zscore = extract_zscore()
    hmatrix = extract_hmatrix()
    time_align = extract_time()
    blosum = extract_blosum()
    weights = extract_weights()
    dssp = extract_dssp()
    cpu = extract_cpu()

    blosum_read = False
    dssp_read = False  
    id_list = []
    zscore_list = []
    hmatrix_list = []
    time_list = []
    weights_list = []
    cpu_list = []

    with open(file_log,"r") as results_file:
        for line in results_file:
            id_line = id.search(line)
            zscore_line = zscore.search(line)
            hmatrix_line = hmatrix.search(line)
            time_line = time_align.search(line)
            blosum_line = blosum.search(line)
            weights_line = weights.search(line)
            dssp_line = dssp.search(line)
            cpu_line = cpu.search(line)

            if id_line:
                id_list.append(int(id_line.group(2)))

            if zscore_line:
                zscore_list.append((int(zscore_line.group(2)),str(zscore_line.group(3))))

            if hmatrix_line:
                hmatrix_list.append(float(hmatrix_line.group(2)))

            if time_line:
                time_list.append(float(time_line.group(2)))

            if weights_line:
                weights_list.append((weights_line.group(2) + "_"+weights_line.group(3)))

            if cpu_line:
                cpu_list.append(int(cpu_line.group(2)))
    return id_list,zscore_list,hmatrix_list,time_list,weights_list,cpu_list

def plot_time(cpu_list,time_list):        
    plt.plot(cpu_list,time_list)
    plt.title("Execution time over number of CPU (79x79 matrix)")
    plt.xlabel("Number of CPU")
    plt.ylabel("Time (s)")    
    plt.show()

def plot_weights(weights_list,id_list,hmatrix_list):
    dico = {}
    weights_list = [("1_1")]  + weights_list
    weights_list = [("1_0")]  + weights_list


    hmatrix_list = [12158.571] + hmatrix_list
    hmatrix_list = [5396.130] + hmatrix_list


    id_list = [40] + id_list
    id_list = [36] + id_list

    identity = [f"{x}%" for x in id_list]
    print(identity)
    
    for i in range(len(weights_list)):
        if weights_list[i] not in dico:
            dico[weights_list[i]] = []
        dico[weights_list[i]] = hmatrix_list[i]
    print(weights_list)

    plt.figure(figsize = (10, 5))

    plt.bar(dico.keys(),dico.values())
    addlabels(dico.keys(),id_list )
    plt.title("Best score of High Matrix and Percent Identity by Weights for Dope and Blosum score")
    plt.xlabel("Weights Dope_Blosum")
    plt.ylabel("Best score")
    plt.show()

def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i]*50,f"{y[i]}%", ha = 'center')

if __name__ == "__main__":
    logfolder = "../log/"
    file_weights = "weight.log"
    file_time = "plt_time64.log"
    file_options = "with_options.log"
    file_zscore = "zscore.log"
    file_log = f"{logfolder}{file_zscore}"  
    id_list,zscore_list,hmatrix_list,\
    time_list,weights_list,cpu_list = extract_results(file_log)

    print(zscore_list)