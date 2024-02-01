import numpy as np
import pandas as pd
import math 
import numpy as np
import pandas as pd
import math 
from pandas import ExcelWriter
from pandas import ExcelFile
import scipy
from scipy import stats
from scipy.stats import pearsonr
import warnings
from scipy.stats import ttest_ind
from more_itertools import unique_everseen
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from scipy.optimize import curve_fit
from lmfit.models import LinearModel, StepModel
from lmfit import Model
from matplotlib import gridspec
import re
from adjustText import adjust_text
import itertools
from pylab import *
from glob import glob
import matplotlib.style as mplstyle

from SigProfilerExtractor import sigpro as sig
# to get input from vcf files
def main_function():
	path_to_example_folder_containing_vcf_files = sig.importdata("vcf")
	data = "/gpfs/group/home/thanigan/thanigan/Sequencing/Resistant_Clones/fc-c662c8e5-11d4-4fbd-a5b3-6fa758fb5ee4/Scripps_Hanigan_WGS_DataDelivery_12samples/RP-2289/WGS/Mutect_Combined_BC/Sig_Profiler_022121_508/" # you can put the path to your folder containing the vcf !!!NOte the vcf must be in its own folder seperate from the folder that you are working in/that output will be put in.
	sig.sigProfilerExtractor("vcf", "example_output", data, reference_genome='GRCh38',minimum_signatures=1, maximum_signatures=5)
if __name__=="__main__":
	main_function()

#Wait untill the excecution is finished. The process may a couple of hours based on the size of the data.
#Check the current working directory for the "example_output" folder.
def main_function():
	path_to_example_table = sig.importdata("matrix")
	data2 = "/gpfs/group/home/thanigan/thanigan/Sequencing/Resistant_Clones/fc-c662c8e5-11d4-4fbd-a5b3-6fa758fb5ee4/Scripps_Hanigan_WGS_DataDelivery_12samples/RP-2289/WGS/Mutect_Combined_BC/Sig_Profiler_022221_508_Out/" # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	sig.sigProfilerExtractor("matrix", "example_output", data2, opportunity_genome="GRCh38", minimum_signatures=1, maximum_signatures=5)
if __name__=="__main__":
	main_function()

#output_path = "/gpfs/group/home/thanigan/thanigan/Sequencing/Resistant_Clones/fc-c662c8e5-11d4-4fbd-a5b3-6fa758fb5ee4/Scripps_Hanigan_WGS_DataDelivery_12samples/RP-2289/WGS/Mutect_Combined/SigProfiler_202021/example_output/"
#import sigProfilerPlotting as sigPlt
#sigPlt.plotSBS(data2, output_path, project, plot_type, percentage=False)