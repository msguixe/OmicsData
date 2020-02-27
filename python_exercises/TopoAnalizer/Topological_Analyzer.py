#  ______                  __            _            __   
# /_  ______  ____  ____  / ____  ____ _(__________ _/ /   
#  / / / __ \/ __ \/ __ \/ / __ \/ __ `/ / ___/ __ `/ /    
# / / / /_/ / /_/ / /_/ / / /_/ / /_/ / / /__/ /_/ / /     
#/_/  \____/ .___/\____/_/\____/\__, /_/\___/\__,_/_/      
#         /_/__                //___/                      
#          /   |  ____  ____ _/ __  ______ ___  _____      
#         / /| | / __ \/ __ `/ / / / /_  // _ \/ ___/      
#        / ___ |/ / / / /_/ / / /_/ / / //  __/ /          
#       /_/  |_/_/ /_/\__,_/_/\__, / /___\___/_/           
#                            /____/                

# By: Laureano Tomás Daza
# Email: lauretomas@gmail.com
# Copyright 2018 Laureano Tomás Daza
# Python 3.6
# OS: Windows 10 Home version: 1709    

'''Topological Analyzer or Topo Analizer is a program to analyze mainly
the topology of proteins from Uniprot, and retrieve their motifs and
domains. The only input of the program is the Uniprot Code of the
protein, and the different outputs are:
- The aminoacid counts of each topological region
- The aminoacid analysis under 2 classification of the topological
	region(csv) and a plot(png)
- The motifs in a csv file
- The dominains in a csv file
- The sequence in fasta format'''

    
# Module Importation

'''Firstly we imported the necessary modules for the Graphical User Interface
(GUI), mainly the program is based on Tkinter but it also needs other modules 
to manage some features'''

from tkinter import *
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from PIL import ImageTk
from PIL import Image
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

'''We also imported pandas because it is used throughout the program'''

import pandas as pd

# Class Definition

class CreateToolTip(object):

    '''We defined a new class named CreateToolTip in order to place an hover
box (tooltio or pop-up) over some widget of the program to display some 
additional information '''

    def __init__(self, widget, text='widget info'):
        
        '''Here we defined the general function that runs with the program
        It reads what kind of widget are you using and bind the actions
        <Enter>(enter) and <Leave>(close) to it'''
        
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.close)
        
    def enter(self, event=None):
        
        '''Here we defined the enter function, i.e., we defined what action
        will happen when you enter with the mouse pointer in the widget.
        
        In this case the coordinates to the new window that will be displayed
        over the widget are set and in this new window we add a label with
        the text'''
        
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = tk.Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background='#C3FEFF', relief='solid', borderwidth=1,
                       font=("times", "9", "normal"))
        label.pack(ipadx=1)
        
    def close(self, event=None):
        
        '''Here we defined the complementary function to Enter, Close which
        destroys the pop-up window when we leave the widget'''
        
        if self.tw:
            self.tw.destroy()

# Function Defition

def nocode_det(*args):
    
    '''The function nocode_det detects if the inserted code is a valid
    Uniprot Code, to do that we use try/except, if the code is valid the
    label_code will be update to show the code and the protein name,
    if the code is empty it will display nothing, and if the code is invalid
    it will show an Error'''
    
    a = code.get()
    try: 
        from Bio import ExPASy
        from Bio import SwissProt
        with ExPASy.get_sprot_raw(a) as handle:
            prueba = SwissProt.read(handle)
            b = prueba.entry_name
            c = a + " - " + b
        label_code.set(c)
    except:
        if a == "":
            label_code.set(a)
        else:
            label_code.set("ERROR: invalid UNIPROT Code")

def open_window_analysis(title, df, title_save, figure, switch):
    
    '''The function open_window_analysis will open a new window when the button
    Analyze is pressed (because it used in the function topology_analisis),
    in this new window will be display a text box with the data frame from the
    analysis and a button to save the analysis.
    If the analysis is performed by any of the classifications another 2 
    buttons will be displayed, one to show a plot and the other to save it'''
    
    window = tk.Toplevel(root, background = '#eafaff')
    window.title(title)
    ttk.Button(window, text="Save Analysis",command = lambda: save_analysis(title_save,df)).grid(column=1, row=3, sticky=W)
    
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        complete_df = str(df)
   
    txt = tk.Text(window)
    txt.grid(column = 2, row = 1)
    txt.insert(tk.END,complete_df)
    txt.configure(state = "disabled")
    
    if switch == "yes":
        title_fig = title + " Stacked Barplot"
    
        ttk.Button(window, text="View Plot",command = lambda: figure_windows(title_fig,figure)).grid(column=2, row=3, sticky=W)
    
        title_save_fig = title_save + " Stacked Barplot"
    
        ttk.Button(window, text="Save Plot",command = lambda: save_figure(title_save_fig,figure)).grid(column=3, row=3, sticky=(E,W,S,N))
        
    
    
def figure_windows(title,figure):
    
    '''The function figure_windows display a new window when the button View 
    Plot is pressed. In this new window will be display only the figure, the
    figure will be in a canvas that has a scrollbar to move through the plot'''
    
    window = tk.Toplevel(root)
    window.title(title)
    fig = figure
    fig.set_size_inches(14, 14)
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.show()
    canvas_fig = canvas.get_tk_widget()
    
    canvas_fig.configure(scrollregion=(0,0,500,500))
    hbar=Scrollbar(window,orient=HORIZONTAL)
    hbar.pack(side=BOTTOM,fill=X)
    hbar.config(command=canvas_fig.xview)
    vbar=Scrollbar(window,orient=VERTICAL)
    vbar.pack(side=RIGHT,fill=Y)
    vbar.config(command=canvas_fig.yview)
    canvas_fig.config(width=300,height=300)
    canvas_fig.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
    canvas_fig.pack(side=LEFT,expand=True,fill=BOTH)
    

def save_figure(title, figure):
    
    '''The function save_figure display a window for the user to save the plot
    in the directory that he/she wants. It used when the button Save Plot is
    pressed'''
    
    filename = filedialog.asksaveasfilename(title=title,defaultextension = 'png', filetypes=[("png","*.png")])
    figure.set_size_inches(14, 14)
    figure.savefig(filename, dpi = 720)

def save_analysis(title, df):
    
    '''The function save_analysis display a window like save_figure function
    asking for the directory to save the analysis. It used when the button
    Save Analysis is pressed'''
    
    filename = filedialog.asksaveasfilename(title=title,defaultextension = 'csv', filetypes=[("csv","*.csv")])
    df.to_csv(filename)

def sequence_file(*args):
    
    '''The function sequence_file save the sequence of the protein in fasta
    format, to do so the sequence is retrieved and the other necessary
    information to make the fasta header.
    We included a try/except chunck to display an Error if the code is invalid'''
    
    a = code.get()
    try: 
        from Bio import ExPASy
        from Bio import SwissProt
        with ExPASy.get_sprot_raw(a) as handle:
            record = SwissProt.read(handle)
    except:
        if a == "":
            open_window("No Code","Please Insert an Uniprot Code", "#FFC3C3", '200x30')
        else:
            open_window("No Valid Code","Please Insert a valid Uniprot Code", "#FFC3C3", '200x30')
    
    descrip = record.description.split(";")[0]
    num = descrip.find("Full=") + 5
    descrip = descrip[num:]
    fasta_header = ">sp|" + code.get() + "|" + record.entry_name + " " + descrip + " OS=" + record.organism
    
    filename = filedialog.asksaveasfilename(defaultextension = '.fasta', filetypes=[("fasta","*.fasta")])
    TextFile = open(filename,"w")
    TextFile.write(fasta_header + '\n')
    TextFile.write(record.sequence)
    TextFile.close()

def copy_button():
    
    '''The function copy_button allows the user to copy the sequence directly
    to the clipboard to be used in other programs or whatever. We also
    included a try/except chunk to show an Error if necessary'''
    
    a = code.get()
    try: 
        from Bio import ExPASy
        from Bio import SwissProt
        with ExPASy.get_sprot_raw(a) as handle:
            record = SwissProt.read(handle)
    except:
        if a == "":
            open_window("No Code","Please Insert an Uniprot Code", "#FFC3C3", '200x30')
        else:
            open_window("No Valid Code","Please Insert a valid Uniprot Code", "#FFC3C3", '200x30')

    clip = Tk()
    clip.withdraw()
    clip.clipboard_clear()
    clip.clipboard_append(record.sequence)

def open_window(title, message, color, geometry):
    
    '''The function open_window is a generic function the open a new window
    and shows a message, it is used in other functions mainly to show Errors'''
    
    pyr = Tk()
    pyr.geometry(geometry)
    pyr.title(title)
    pyr.configure(background = color)
    lbl = Label(pyr, text=message, bg=color)
    lbl.grid(column=2, row=1)
    pyr.mainloop()

def topology_analysis(*args):
    
    '''The topology_analysis function is the main function of the program, it
    performs the aminoacyds analysis of the sequence in the different regions
    of the topology'''
    
    '''Firstly we use a try/except chunk to show an Error if necessary. If the
    code is valid we retrive from Uniprot information of the Uniprot Code'''
    
    a = code.get()
    try: 
        from Bio import ExPASy
        from Bio import SwissProt
        with ExPASy.get_sprot_raw(a) as handle:
            record = SwissProt.read(handle)
    except:
        if a == "":
            open_window("No Code","Please Insert an Uniprot Code", "#FFC3C3", '200x30')
        else:
            open_window("No Valid Code","Please Insert a valid Uniprot Code", "#FFC3C3", '200x30')

    '''From the record object we retrieve the sequence and using
    ProteinAnalysis module we analyzed this sequence and counted the
    aminoacyds of the complete sequence'''
    
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    import pandas as pd
    seq = record.sequence
    analysed_seq = ProteinAnalysis(seq)
    coun = analysed_seq.count_amino_acids()
    
    '''We also retrieve feature information and we save only the information
    regarding to the topology regions'''
    
    feat = record.features

    topo = 'INTRAMEM', 'TRANSMEM', 'TOPO_DOM'

    topology = []
    for i in feat:
        if i[0] in topo:
            topology.append(i)
    
    if topology == []:
        
        '''If there is no topology a new window showing a message will be
        displayed'''
        
        protein = str(code.get()) + " - " + record.entry_name
        open_window("No Topology", "There is no topology in the protein " + protein + 
                    "\n This protein has not any Intramembrane, Transmembrane or Topological Domain Region ","#FFC3C3", '500x40')
    else:
        
        '''If there is topology the analysis continues.
        We made the first row of the final data frame with the information
        of the complete sequence, in this case only the name and the amynoacyds
        count'''
        
        row = {}
        row["name"] = "COMPLETE_SEQ"
        row["number"] = ""
        row["end"] = ""
        row["start"] = ""

        complete_row = {**coun, **row}
    
        rows = []
        rows.append(complete_row)
        
        '''For each topology region we made a row with the name, the start and
        end position and the aminoacyds count'''
        
        for i in topology:
            row_i = {}
            row_i["number"] = 1
            row_i["end"] = i[2]
            row_i["start"] = i[1]
            row_i["name"] = i[0]
            analysed_seq_i = ProteinAnalysis(seq[(i[1]-1):(i[2]-1)])
            coun_i = analysed_seq_i.count_amino_acids()
            complete_row_i = {**coun_i, **row_i}
            rows.append(complete_row_i)
            
        '''Once the rows are made we turned them into a dataframe and order
        the columns to be displayed in an easier way'''
            
        df = pd.DataFrame.from_dict(rows)
        ordenation = list(row.keys())
            
        order = ordenation + list(coun.keys())
            
        df = df[order]
        
        '''For each topology class we made a global analysis including the
        number of times the topology class appears in the sequence'''
        
        for j in topo:
            topo_j = df[df.name == j]
            topo_j = topo_j.drop(ordenation, axis = 1)
            topo_sum = topo_j.sum(axis = 0)
            row_j = {}
            row_j["number"] = len(topo_j.index)
            row_j["end"] = ''
            row_j["start"] = ''
            row_j["name"] = 'GLOBAL_'+j
            complete_row_j = {**topo_sum, **row_j}
            rows.append(complete_row_j)

        '''The final data frame for the aminoacyd analysis is made and order
        as before, we set the name of each row with the name value'''

        final_df = pd.DataFrame.from_dict(rows)
    
        final_df.index = final_df['name']

        ordenation = list(row.keys())
    
        order = ordenation + list(coun.keys())

        final_df = final_df[order]
        
        
### CLASSIFICATION I
        
        '''If the classification I button is pressed the aminoacyd analysis
        will be performed grouping the aminoacids into 2 categories: 
        hydrophobic and hydrophilic, and the results will be displayed in a
        new window'''
        

        if (check_class1.get() == "class1") & (check_class2.get() == "noclass2"):
            
            classification1 = {'hydrophobic': ['F','W','G','A','V','L','I','M','P'],
                   'hydrophilic': ['E','D','K','H','R','Y','S','T','C','Q','N']}
            
            '''From the data frame of the aminoacyd counts we groups these
            count by the classification I and create 2 new columns, one with
            the absolute count and other with the relative count'''
            
            class_1_rows = []
            for j in final_df.iterrows():
                df_row = j[1]
                row_ji = {}
                row_ji["number"] = df_row["number"]
                row_ji["end"] = df_row["end"]
                row_ji["start"] = df_row["start"]
                row_ji["name"] = df_row["name"]
                for i in classification1:
                    class_i = df_row[classification1[i]]
                    row_ji[i] = class_i.sum(0)
                    if df_row.drop(ordenation).sum(0) == 0:
                        row_ji[i+'_rel'] = 0
                    else:
                        row_ji[i+'_rel'] = class_i.sum(0)/df_row.drop(ordenation).sum(0)
        
                class_1_rows.append(row_ji)

            final_df_class1 = pd.DataFrame.from_dict(class_1_rows)
            final_df_class1.index = final_df_class1['name']
            
            final_df_class1 = final_df_class1[ordenation + list(classification1.keys()) + 
                                  [s + '_rel' for s in list(classification1.keys())]]
            
            '''We made a stacked barplot with the information of the relative 
            frequencies for each row'''
            
            plot_data =final_df_class1[['name']+[s + '_rel' for s in list(classification1.keys())]]

            plot_image = plot_data.plot(kind="bar",stacked = True)

            figure = plot_image.get_figure()
            
            open_window_analysis("Classification I Analysis", final_df_class1, "Save Classification I Analysis", figure, "yes")

### CLASSIFICATION II
        
        elif (check_class2.get() == "class2") & (check_class1.get() == "noclass1"):
            
            '''If the classification II button is pressed the aminoacyd 
            analysis will be performed grouping the aminoacids into 5 
            categories: aromatic, alkyl, acidic, basic and neutral, and the 
            results will be displayed in a new window'''

            classification2 = {'aromatic': ['F','W'], 'alkyl': ['G','A','V','L','I','M','P'],
                   'acidic': ['E','D'], 'basic': ['K','H','R'], 
                   'neutral': ['Y','S','T','C','Q','N']}

            '''From the data frame of the aminoacyd counts we groups these
            count by the classification II and create new columns, with
            the absolute counts and with the relative counts'''

            class_2_rows = []
            for j in final_df.iterrows():
                df_row = j[1]
                row_ji = {}
                row_ji["number"] = df_row["number"]
                row_ji["end"] = df_row["end"]
                row_ji["start"] = df_row["start"]
                row_ji["name"] = df_row["name"]
                for i in classification2:
                    class_i = df_row[classification2[i]]
                    row_ji[i] = class_i.sum(0)
                    if df_row.drop(ordenation).sum(0) == 0:
                        row_ji[i+'_rel'] = 0
                    else:
                        row_ji[i+'_rel'] = class_i.sum(0)/df_row.drop(ordenation).sum(0)
        
                class_2_rows.append(row_ji)

            final_df_class2 = pd.DataFrame.from_dict(class_2_rows)

            final_df_class2.index = final_df_class2['name']

            final_df_class2 = final_df_class2[ordenation + list(classification2.keys()) + 
                                  [s + '_rel' for s in list(classification2.keys())]]
            
            '''We made a stacked barplot with the information of the relative 
            frequencies for each row'''
            
            plot_data =final_df_class2[['name']+[s + '_rel' for s in list(classification2.keys())]]

            plot_image = plot_data.plot(kind="bar",stacked = True)

            figure = plot_image.get_figure()
            
            open_window_analysis("Classification II Analysis", final_df_class2, "Save Classification II Analysis", figure, "yes")
        
#### BOTH CLASSIFICATIONS 
                    
        elif (check_class1.get() == "class1") & (check_class2.get() == "class2"):
            
            '''If both classifications buttons are pressed the aminoacyd 
            analysis will be performed grouping the aminoacids the 2 
            classification, and the results will be displayed in 2 new windows.
            Both analysis are performed as described above'''
            
            classification1 = {'hydrophobic': ['F','W','G','A','V','L','I','M','P'],
                   'hydrophilic': ['E','D','K','H','R','Y','S','T','C','Q','N']}

            class_1_rows = []
            for j in final_df.iterrows():
                df_row = j[1]
                row_1 = {}
                row_1["number"] = df_row["number"]
                row_1["end"] = df_row["end"]
                row_1["start"] = df_row["start"]
                row_1["name"] = df_row["name"]
                for i in classification1:
                    class_i = df_row[classification1[i]]
                    row_1[i] = class_i.sum(0)
                    if df_row.drop(ordenation).sum(0) == 0:
                        row_1[i+'_rel'] = 0
                    else:
                        row_1[i+'_rel'] = class_i.sum(0)/df_row.drop(ordenation).sum(0)
        
                class_1_rows.append(row_1)

            final_df_class1 = pd.DataFrame.from_dict(class_1_rows)
            final_df_class1.index = final_df_class1['name']
            
            final_df_class1 = final_df_class1[ordenation + list(classification1.keys()) + 
                                  [s + '_rel' for s in list(classification1.keys())]]
            
            plot_data_1 =final_df_class1[['name']+[s + '_rel' for s in list(classification1.keys())]]

            plot_image_1 = plot_data_1.plot(kind="bar",stacked = True)

            figure_1 = plot_image_1.get_figure()
            
            open_window_analysis("Classification I Analysis", final_df_class1, "Save Classification I Analysis", figure_1, "yes")


            classification2 = {'aromatic': ['F','W'], 'alkyl': ['G','A','V','L','I','M','P'],
                   'acidic': ['E','D'], 'basic': ['K','H','R'], 
                   'neutral': ['Y','S','T','C','Q','N']}

            class_2_rows = []
            for j in final_df.iterrows():
                df_row = j[1]
                row_2 = {}
                row_2["number"] = df_row["number"]
                row_2["end"] = df_row["end"]
                row_2["start"] = df_row["start"]
                row_2["name"] = df_row["name"]
                for i in classification2:
                    class_i = df_row[classification2[i]]
                    row_2[i] = class_i.sum(0)
                    if df_row.drop(ordenation).sum(0) == 0:
                        row_2[i+'_rel'] = 0
                    else:
                        row_2[i+'_rel'] = class_i.sum(0)/df_row.drop(ordenation).sum(0)
        
                class_2_rows.append(row_2)

            final_df_class2 = pd.DataFrame.from_dict(class_2_rows)

            final_df_class2.index = final_df_class2['name']

            final_df_class2 = final_df_class2[ordenation + list(classification2.keys()) + 
                                  [s + '_rel' for s in list(classification2.keys())]]
            
            plot_data_2 =final_df_class2[['name']+[s + '_rel' for s in list(classification2.keys())]]

            plot_image_2 = plot_data_2.plot(kind="bar",stacked = True)

            figure_2 = plot_image_2.get_figure()
            
            open_window_analysis("Classification II Analysis", final_df_class2, "Save Classification II Analysis", figure_2, "yes")
            
        elif (check_class1.get() in ("noclass1", "")) & (check_class2.get() in ("noclass2", "")):
            
            '''If no classification button is pressed a new window will be
            displayed showing the results from the aminoacids counts'''
            
            open_window_analysis("Aminoacyds Analysis", final_df, "Save No Classification Analysis", "","no")
            
def find_motifs(*args):
    
    '''The function find_motifs find the motifs of the sequences if they exists.
    It is used when the button Find Motifs is pressed'''
    
    '''Firstly we use a try/except chunk to show an Error if necessary. If the
    code is valid we retrive from Uniprot information of the Uniprot Code'''
    
    a = code.get()
    try: 
        from Bio import ExPASy
        from Bio import SwissProt
        with ExPASy.get_sprot_raw(a) as handle:
            record = SwissProt.read(handle)
    except:
        if a == "":
            open_window("No Code","Please Insert an Uniprot Code", "#FFC3C3", '200x30')
        else:
            open_window("No Valid Code","Please Insert a valid Uniprot Code", "#FFC3C3", '200x30')

    '''From the record object we retrieve the sequence and the features, and
    we save only the information regarding to motifs'''

    seq = record.sequence
    
    feat = record.features

    motifs = []
    for i in feat:
        if i[0] == "MOTIF":
            motifs.append(i)
    
    if motifs == []:
        
        '''If there is no motif a new window showing a message will be 
        displayed'''
        
        protein = str(code.get()) + " - " + record.entry_name
        open_window("No Motifs", "There is no motifs in the protein " + protein,"#FFC3C3", '330x30')
    else:
        
        '''If there are motifs the analysis continues.
        We made a row for each motif with the name/description, its position
        and the sequence. The results are displayed in a new window'''
        
        rows_m = []

        for i in motifs:
            row_i = {}
            row_i["end"] = i[2]
            row_i["start"] = i[1]
            row_i["name/description"] = i[3].split(".")[0]
            row_i["length"] = i[2]-i[1]
            row_i["sequence"] = seq[i[1]:i[2]]
            
            rows_m.append(row_i)

        final_df_motif = pd.DataFrame.from_dict(rows_m)

        order = ["name/description","length","start","end","sequence"]

        final_df_motif = final_df_motif[order]
        
        open_window_analysis("Motifs Analysis", final_df_motif, "Save Motif Analysis", "", "no")


def find_domains(*args):
    
    '''The function find_domains find the domains of the sequences if they exists.
    It is used when the button Find Domains is pressed'''
    
    '''Firstly we use a try/except chunk to show an Error if necessary. If the
    code is valid we retrive from Uniprot information of the Uniprot Code'''
    
    a = code.get()
    try: 
        from Bio import ExPASy
        from Bio import SwissProt
        with ExPASy.get_sprot_raw(a) as handle:
            record = SwissProt.read(handle)
    except:
        if a == "":
            open_window("No Code","Please Insert an Uniprot Code", "#FFC3C3", '200x30')
        else:
            open_window("No Valid Code","Please Insert a valid Uniprot Code", "#FFC3C3", '200x30')

    '''From the record object we retrieve the features, and we save only the 
    information regarding to domains'''
    
    feat = record.features

    domain = []
    for i in feat:
        if i[0] == "DOMAIN":
            domain.append(i)
    
    if domain == []:
        
        '''If there is no domain a new window showing a message will be 
        displayed'''
        
        protein = str(code.get()) + " - " + record.entry_name
        open_window("No Domain", "There is no domain in the protein " + protein,"#FFC3C3", '330x30')
    else:
        
        '''If there are domains the analysis continues.
        We made a row for each domain with the name/description, its position
        and the length. The results are displayed in a new window'''
        
        rows_d = []

        for i in domain:
            row_i = {}
            row_i["end"] = i[2]
            row_i["start"] = i[1]
            row_i["name/description"] = i[3].split(".")[0]
            row_i["length"] = i[2]-i[1]
            
            rows_d.append(row_i)

        final_df_domain = pd.DataFrame.from_dict(rows_d)

        order = ["name/description","length","start","end"]

        final_df_domain = final_df_domain[order]
        
        open_window_analysis("Domains Analysis", final_df_domain, "Save Domain Analysis", "", "no")

def help_fun():
    
    '''The function help_fun is used when the help button (Question Mark
    button) is pressed, and it displays a new window with the Intructions
    of the program'''
    
    pyr = Tk()
    pyr.title("Topological Analizer Help")
    data = 'Topological Analizer Help\n\nUsage:\n\nInsert an Uniprot Code in the entry box at the top center of the window:\n\n- if the Uniprot Code is valid it will appear below the code and\n\tthe protein name\n- if the Uniprot Code is invalid it will appear an Error\n\nAminoacyds Analysis Section:\n\nIn this section the program analyses the topological regions of the sequences,\ni.e,transmembrane, intramembrane, and topological domain regions.\n\nHere you can choose if you want to analyse the aminoacyd content of these\nregions by Classification I (Hydrophobic and Hydrophylic), Classification II\n(Aromatic, Alkyl, Acidic, Basic and Neutral), both classification or not\nclassification at all, when you click in the "Analyse" button the analysis\nwill be displayed.\n\n- If you do not choose any classification the program will display a window\nwith a data frame with the aminoacyd content by aminoacyds and you can save\nthis analysis into a csv file\n\n- If you choose one classification the program will display a window with a\ndata frame with the aminoacyd content summarize by the classification chosen.\nYou will be able to save the analysis into a csv file, display a\nstacked-barplot of the results and save the plot into a png file\n\n- If you choose both classification the program will display 2 windows, each\nof them with one analysis and all the correspoding options\n\n- If the inserted protein has not any topology the program will display an Error\n\n\nSequence Section:\n\nHere you can save the sequence in fasta format into a file or copy the complete\nsequence to the clipboard if you need to use it in other program.\n\nBesides, you can find Motifs and Domains in the sequence, if you click any of\nthese buttons a window will be display showing a data frame with the analysis\nand you can save this analysis into a csv file.\n\nIf the protein has not any Motif or Domain the program will display an Error.'    
    txt = tk.Text(pyr)
    txt.grid(column = 2, row = 1)
    txt.insert(tk.END,data)
    txt.configure(state = "disabled")
    
    pyr.mainloop()


# GUI set-up

'''Firstly we created the root windows with a title'''

root = tk.Tk()
root.title("Topological Analizer")

'''Then the main frame is created and we configure its settings'''

mainframe = tk.Frame(root,bg='#eafaff')
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)
    
mainframe.pack(expand=1, fill='both')

'''We created string variable that will store the code inserted by the user
in the code_entry created below'''

code = StringVar()

code_entry = ttk.Entry(mainframe, width=7, textvariable=code)
code_entry.grid(column=2, row=2, sticky=(W, E))

'''We used the function nocode_det with the code variable in order to detect
if the inserted code is valid or not. This result will be displayed using
the label_code string variable and the label created below'''

code.trace('w',nocode_det)

label_code = StringVar()

ttk.Label(mainframe, textvariable=label_code, background = '#eafaff').grid(column=2, row=3, sticky=(W, E))

'''We created 2 label in order to be more self explanatory, one asking the user
to insert an Uniprot Code and the other one next to where the code if valid or
the error if not will be displayed'''

ttk.Label(mainframe, text="Insert Uniprot Code", background = '#eafaff', font="Helvetica 11 bold", anchor = "center").grid(column=2, row=1, sticky=(W,E,N,S))
ttk.Label(mainframe, text="Uniprot Code", background = '#eafaff').grid(column=1, row=3, sticky=E)

'''We create the help button and used the function help_fun in it. We also
added to it a pop-up'''

'''In order to make the program in a executable file this step is a little bit
different in the script used to make the executable file, because the file
help.png here is used in my computer, but with the program that turns the
script into an executable file the png file is stored locally and the path is
a bit different. See the other version of the file '''         
          
image = Image.open(r"help.png")
photo = ImageTk.PhotoImage(image)
help_b = tk.Button(mainframe, image = photo, command = help_fun)
help_b.image = photo
help_b.grid(column=3, row=1, sticky=E)

CreateToolTip(help_b, "Help: How to use it")

## AMINOACYDS ANALYSIS SECTION

'''We created a label to set the tittle of the section'''

ttk.Label(mainframe, text="Aminoacyds Analysis", background = '#eafaff', foreground = "#0800A3",wraplength=250, font="Helvetica 11 underline", anchor = "center").grid(column=2, row=6, sticky=(W, E))

'''We create 2 string variables to save the information of the 2 check buttons 
of each classification, and label with pop-ups to set the categories of each
classification.

In the checkbuttons we set an on-value and off-value, these variables are the
ones used in the function topology_analysis to check which classification 
the user has chosen'''          
          
check_class1 = StringVar()
check_class1.set("noclass1")
check1 = Checkbutton(mainframe, text='Classification I', activebackground= '#eafaff',background = '#eafaff', state="normal", onvalue="class1", offvalue="noclass1", variable = check_class1).grid(column = 1, row = 7, sticky = (W,E))

hydrop = tk.Label(mainframe, text="     -  Hydrophobic", background = '#eafaff', wraplength=250)
hydrop.grid(column=1, row=8, sticky=(W))
CreateToolTip(hydrop, "F, W, G, A, V, L, I, M, P")

hydrof = tk.Label(mainframe, text="     -  Hydrophilic", background = '#eafaff', wraplength=250)
hydrof.grid(column=1, row=9, sticky=(W))
CreateToolTip(hydrof, "E, D, K, H, R, Y, S, T, C, Q, N")

check_class2 = StringVar()
check_class2.set("noclass2")
check2 = Checkbutton(mainframe, text='Classification II', activebackground= '#eafaff', background = '#eafaff', state="normal", onvalue="class2", offvalue="noclass2", variable = check_class2).grid(column = 3, row = 7, sticky = (W,E))

aromatic = tk.Label(mainframe, text="     -  Aromatic", background = '#eafaff', wraplength=250)
aromatic.grid(column=3, row=8, sticky=(W))
CreateToolTip(aromatic, "F, W")

alkyl = tk.Label(mainframe, text="     -  Alkyl", background = '#eafaff', wraplength=250)
alkyl.grid(column=3, row=9, sticky=(W))
CreateToolTip(alkyl, "G, A, V, L, I, M, P")

acidic = tk.Label(mainframe, text="     -  Acidic", background = '#eafaff', wraplength=250)
acidic.grid(column=3, row=10, sticky=(W))
CreateToolTip(acidic, "E, D")

basic = tk.Label(mainframe, text="     -  Basic", background = '#eafaff', wraplength=250)
basic.grid(column=3, row=11, sticky=(W))
CreateToolTip(basic, "K, H, R")

neutral = tk.Label(mainframe, text="     -  Neutral", background = '#eafaff', wraplength=250)
neutral.grid(column=3, row=12, sticky=(W))
CreateToolTip(neutral, "Y, S, T, C, Q, N")

'''Finally, we placed a button to analyze the aminoacyds'''

ttk.Button(mainframe, text="Analyze", command = topology_analysis).grid(column=2, row=13, sticky=(W,E,N,S))


## SEQUENCE SECTION

'''We created a label to place the title of the section and 4 buttons, to Find
Motifs, Save the sequence, Find Domains and Save the sequence to the clipboard'''

ttk.Label(mainframe, text="Sequence", background = '#eafaff', foreground = "#950000", wraplength=250, font="Helvetica 11 underline", anchor = "center").grid(column=2, row=15, sticky=(W,E,N,S))

ttk.Button(mainframe, text="Save Sequence", command = sequence_file).grid(column=2, row=16, sticky=(W,E,N,S))

ttk.Button(mainframe, text="Copy Sequence \n   to Clipboard", command = copy_button).grid(column=2, row=17, sticky=(W,E,N,S))

ttk.Button(mainframe, text="Find Motifs", command = find_motifs).grid(column=1, row=16, sticky=W)

ttk.Button(mainframe, text="Find Domains", command = find_domains).grid(column=3, row=16, sticky=E)

####

'''We used this for loop to add some space between the widgets in the 
different positions of the frame in order to make it looks nicer'''

for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

'''We set the pointer to be in the code_entry entry, to make it easier to the
user. And we bind to the Enter key the button Analyze, so when the user press
the Enter key the topology_analysis will be performed'''

code_entry.focus()
root.bind('<Return>', topology_analysis)

'''We set the program icon, this step is also different from the executable
file'''

root.iconbitmap(r'topo_analyzer.ico')

'''Finally we run the mainloop of the root windows to display the GUI'''

root.mainloop()
