# Libraries
## SuperVenn
from supervenn import supervenn
## Data Processing
import csv
import pathlib
## Table
import plotly.graph_objects as go
from IPython.display import display, FileLink, Markdown, HTML
%matplotlib inline
from plotnine import ggplot, aes, geom_line
import matplotlib.pyplot as plt
import os

os.chdir("c:\\Users\\vasileioubill95\\Desktop\\Pipeline\\Figures\\Figure_3\\")

# Input Settings
svenn = True
annotations = 150
svenn_file_format = ['svg']
svenn_file_name = 'supervenn'
final_svenn_file_names = [str(svenn_file_name + '.' + file_type) for file_type in svenn_file_format]

# Displaying Figures
def figure_title(label, title):
    display(HTML(f"<div style='font-size:2rem; padding;1rem 0;'><b>{label}</b>: {title}</div>"))

def figure_legend(label, title, content=""):
    display(HTML(f"<div><b>{label}</b>: <i>{title}</i>. {content} </div>"))

# Saving Figures
def save_figure(plot_name, **kwargs):
    import io
    mem = io.BytesIO()
    pyplot.savefig(mem, bbox_inches='tight')
    with open(plot_name, 'wb') as fw:
        fw.write(mem.getbuffer())
        

# Loading helper
def load_sets(*files):
    ''' Load a set of files into pairs of labeled sets
    '''
    sets = {}
    for file in map(pathlib.Path, files):
        if file.suffix == '.gmt':
            for line in map(str.strip, file.open('r')):
                line_split = line.split('\t')
                if len(line_split) < 3: continue
                term, description, *geneset = line_split
                term_description = ' '.join(filter(None, map(str.strip, [
                    file.stem if len(files) > 1 else '',
                    term,
                    description,
                ])))
                sets[term_description] = set(filter(None, map(str.strip, geneset)))
        else:
            # assumed file is newline separated genes if not a gmt
            sets[file.stem] = set(map(str.strip, file.open('r')))
    return sets



# Add the appropriate gene lists to the dictionary
gsdict = load_sets("C:/Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\6.Overlaps\\Activity\\MGH\\MGH_pos.txt",
                   "C:/Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\6.Overlaps\\Activity\\MGH\\MGH_neg.txt",
                   "C:/Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\6.Overlaps\\Activity\\Mayo\\Mayo_pos.txt",
                   "C:/Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\6.Overlaps\\Activity\\Mayo\\Mayo_neg.txt",
                   "C:/Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\6.Overlaps\\Activity\\Stanford\\Stanford_pos.txt",
                   "C:/Users\\vasileioubill95\\Desktop\\Pipeline\\Case_study_1\\6.Overlaps\\Activity\\Stanford\\Stanford_neg.txt")


# SuperVenn
if svenn:
    figure_title("Figure 2", "SuperVenn")
    supervenn(list(gsdict.values()), list(gsdict.keys()),sets_ordering= 'minimize gaps', side_plots=True,widths_minmax_ratio=0.05, min_width_for_annotation=annotations)
    for plot_name in final_svenn_file_names:
        plt.savefig(plot_name)
    figure_legend("Figure 2", "SuperVenn", "The numbers on the right represent the set sizes and the numbers on the top show how many sets the intersection is part of. The overlapping portions of the colored bars correspond to set intersections.")
