import os
import time
import re
import pandas as pd
from matplotlib import pyplot as plt


class Data:
    def __init__(self, path):
        ''' attributions:
        self.path: str
        self.acq: bool True or False
        self.ncolumn: int 1 or 2
        self.flows: [separator_equilibration, trap_equilibration/separator_loading, /trap_loading]'''
        self.flow_sep_equ = None
        self.flow_trap_equ = None
        self.flow_trap_load = None
        self.flow_sep_load = None
        self.path = path
        self.avg_flow()

    def avg_flow(self, digits=None):  # get average flows, no decimal digits limit by default
        sub_path = [d for d in os.listdir(self.path) if d.endswith('-column-separation')][0]
        log_files = [f for f in os.listdir(self.path + '/' + sub_path) if re.search(r'execution-log_\d+', f)]
        if len(log_files) > 0:
            log_path = self.path + '/' + sub_path + '/' + log_files[0]
        else:
            self.acq = False
            return
        if not os.path.exists(log_path):
            self.acq = False
            return
        with open(log_path, 'r') as f:
            if not f.read().endswith('COMPLETED\n'):  # if acquisition not done
                self.acq = False
            else:
                self.acq = True
            f.seek(0)
            try:
                self.flows = re.findall(r'average flow.+?\[uL/min\]:[\s]+?([\d\.]+)', f.read())
                if digits != None:
                    self.flows = [f'{{:.{str(digits)}f}}'.format(float(e)) for e in self.flows]
                if sub_path.endswith('Two-column-separation'):
                    self.ncolumn = 2
                    if len(self.flows) >= 1:
                        self.flow_sep_equ = self.flows[0]
                        if len(self.flows) >= 2:
                            self.flow_trap_equ = self.flows[1]
                            if len(self.flows) >= 3:
                                self.flow_trap_load = self.flows[2]
                elif sub_path.endswith('One-column-separation'):
                    self.ncolumn = 1
                    if len(self.flows) >= 1:
                        self.flow_sep_equ = self.flows[0]
                        if len(self.flows) >= 2:
                            self.flow_sep_load = self.flows[1]
            except:
                print('Extraction failed in some floder(s): {self.path}')


if __name__ == '__main__':
    # list all .d folders and sort by MS acquisition time
    ds = [d for d in os.listdir('.') if os.path.isdir(d) and d.endswith('.d')]
    ds.sort(key=lambda d: os.path.getmtime(d + '/analysis.tdf_bin') if os.path.exists(d + '/analysis.tdf_bin') else os.path.getctime(d))
    datas = [Data(d) for d in ds]

    # search average flows
    flows = [data.flows for data in datas]
    adict = {'Filename': [data.path for data in datas],
             'Separator Equilibration [uL/min]': [data.flow_sep_equ for data in datas],
             'Trap Equilibration [uL/min] (Two columns)': [data.flow_trap_equ for data in datas],
             'Trap Loading [uL/min] (Two columns)': [data.flow_trap_load for data in datas],
             'Separator Loading [uL/min] (One column)': [data.flow_sep_load for data in datas],
             'N column': [data.ncolumn for data in datas],
             'Acquisition done': [data.acq for data in datas]}
    df = pd.DataFrame(adict)
    df.set_index('Filename', inplace=True)
    t = time.strftime('%Y-%m-%d %H-%M-%S', time.localtime())
    output = 'execution-log_average_flow_summary_' + t + '.txt'
    df.to_csv(output, sep='\t')

    # plotting
    df = df.applymap(pd.to_numeric).fillna(0)
    y1 = df['Separator Equilibration [uL/min]'].tolist()
    y2 = df['Trap Equilibration [uL/min] (Two columns)'].tolist()
    y3 = df['Trap Loading [uL/min] (Two columns)'].tolist()
    y2_remove0 = [e for e in y2 if e != 0]
    y3_remove0 = [e for e in y3 if e != 0]
    x1 = range(1, len(y1) + 1)
    x2 = range(1, len(y2_remove0) + 1)
    x3 = range(1, len(y3_remove0) + 1)
    # configure figures
    plt.figure(figsize=(20, 10), dpi=80)

    # separator equilibration
    fig1 = plt.subplot(311)
    fig1.set_title('Separator Equilibration [uL/min]')
    plt.plot(x1, y1)
    for a, b in zip(x1, y1):
        plt.text(a, b, format(b, '0.2f'), ha='center', va='bottom', fontsize=10)
    plt.grid()

    # trap equilibration
    fig2 = plt.subplot(312)
    fig2.set_title('Trap Equilibration [uL/min]')
    plt.plot(x2, y2_remove0)
    for a, b in zip(x2, y2_remove0):
        plt.text(a, b, format(b, '0.2f'), ha='center', va='bottom', fontsize=10)
    plt.grid()

    # trap loading
    fig3 = plt.subplot(313)
    fig3.set_title('Trap Loading [uL/min]')
    plt.plot(x3, y3_remove0)
    for a, b in zip(x3, y3_remove0):
        plt.text(a, b, format(b, '0.2f'), ha='center', va='bottom', fontsize=10)
    plt.xlabel("File", fontsize=12)
    plt.grid()
    plt.show()
