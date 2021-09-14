"""Extracting informations from log file with regular expression"""

__authors__ = "Ragousandirane RADJASANDIRANE"
__contact__ = "radja.ragou@gmail.com"
__date__ = "14/09/2021"
__version__= "1.0"

# Modules import
import re
import matplotlib.pyplot as plt

def extract_cpu():
    """Regular expression for extract CPU usage"""
    return re.compile(r"^# Nb (\w|\s+)+ : (\d)")
    
def extract_id():
    """Regular expression for extract final percent identity"""
    return re.compile(r"^# (Final [\w|\s]+) : (\d+)%")

def extract_zscore():
    """Regular expression for extract zscore"""
    return re.compile(r"^# (Zscore) \((\d+) \w+\): (-?\d.\d)")

def extract_hmatrix():
    """Regular expression for extract Best score of High Matrix"""
    return re.compile(r"^# (Last [\w|\s]+) : (\d+.\d+)")

def extract_time():
    """Regular expression for extract total time of computation"""
    return re.compile(r"^# (\w* time : )(\d+.\d+)s")

def extract_weights():
    """Regular expression for extract weights used"""
    return re.compile(r"^#( Weights \w+ ): (\d.\d) \w+ (\d.\d) \w+")
    


def extract_results(file_log):
    """Extraction informations"""

    id = extract_id()
    zscore = extract_zscore()
    hmatrix = extract_hmatrix()
    time_align = extract_time()
    weights = extract_weights()
    cpu = extract_cpu()

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
            weights_line = weights.search(line)
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
    """Plotting Time by CPU usage"""

    plt.plot(cpu_list,time_list)
    plt.title("Execution time over number of CPU (79x79 matrix)")
    plt.xlabel("Number of CPU")
    plt.ylabel("Time (s)")    
    plt.show()

def plot_weights(weights_list,id_list,hmatrix_list):
    """Plotting Score and percent identity by weights used"""

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
    log_folder = "../log/analyzed/"
    file_weights = "weight.log"
    file_time = "plt_time64.log"
    file_options = "dssp.log"
    file_zscore = "zscore.log"

    # Change the file to load correct results to plot
    file_log = f"{log_folder}{file_weights}"  

    # Extracting informations from log file
    id_list,zscore_list,hmatrix_list,\
    time_list,weights_list,cpu_list = extract_results(file_log)

    # Load file_time
    #plot_time(cpu_list,time_list)

    # Load file_weights
    plot_weights(weights_list,id_list,hmatrix_list)
    