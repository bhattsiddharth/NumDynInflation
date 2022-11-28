
---
---

# Numerical simulation of inflationary dynamics for PBH formation and beyond

---
---

This code is intended to be used as per the procedure detailed in the [article](<arXiv link>) on arXiv. 

The codes given are tested using Ubuntu 20.04 and the following softwares/packages are required:  
Python 3.10+  
NumPy 1.23+  
SciPy 1.9+  
Matplotplib 3.6+  


Please make a folder named *data* in the same location as the scripts to save the simulation data in text files. 
Further, make a directory named *modes* within *data* to save simulation data of each mode of fluctuations.
Otherwise, please specify the location where you intend to save the data files wherever required in the scripts.


## Background analysis

The script *inf_dyn_background.py* concerns simulation of background dynamics and quantum fluctuations under slow-roll approximation. The values of parameters **Nt** and **v0** need to be fixed as per CMB observations. 
Once these values are fixed, you can save the data in the file *inf_bg_data.txt* within the *data* folder (or any other specified filename and location). 
The size of the file can be varied by adjusting the number of time steps for solving the ODEs.

For phase space analysis, it is recommended that you enter the initial conditions as described in the article and save the data for dynamics due to each set of initial conditions separately, labelling the text file with the value of **xi** used.
ex: *inf_bg_data_5.0.txt*
The different data files, corresponding to different sets of initial conditions can then be plotted in the same plot.


## Quantum fluctuations analysis

The script *inf_dyn_MS_full.py* concerns simulation of (background, along with ) scalar as well as tensor fluctuations. This makes use of the data from the background simulation to read values of initial conditions. There are no parameters whose values need to be fixed in this script, although the value of **v0** can be set more precisely by matching the numerical results with CMB data. In this case, please make the corresponding change in the background script as well and save background data again.
The data for each mode can be saved in files named *inf_MS_data_**Nk**.txt*, where **Nk** corresponds to different modes.


If you plan to work with multiple models of inflation, it is recommended you create separate scripts for each model and specify the same in all filenames for convenience.
ex: *quad_dyn_background.py* or *Starobinsky_MS_data_60.0.txt*


This code is meant to be pedagogical in nature. In the process of fixing the various parameters, entering suitable initial conditions, etc the user will get a clear understanding of inflationary dynamics. We are developing a much more modular and direct PYTHON package akin to a black box which would allow users to compute specific aspects and predictions of inflationary models without having to traverse through all the dynamics manually. Updates regarding the same will be posted here as and when progress is made.


Thank you for taking interest in our work!

For any queries, please contact:
bhattsid24@gmail.com
or
swagatam18@gmail.com

---
