# A 3D model of a geological structure in SE Spain

In this repository we introduce a Python application for creating a 3D model of an imbricate thrust system. The application uses the traditional geologic information to create an HTML geological map with real topography and a set of geological cross-sections with the essential structural and stratigraphic elements. We also created a model that represents a block 300 meters deep, in which all the geological elements appear.

On the basis of the high geological knowledge gained during the last three decades, the Palomeque sheets affecting the Cenozoic Malaguide Succession in the Internal Betic Zone (SE Spain) was selected to show the application. In this area, a Malaguide Cretaceous to lower Miocene Succession is deformed as an imbricate thrust system with two thrusts forming a duplex affected later by three faults with a main slipstrike kinematic. The modeled elements match well with the design of the stratigraphic intervals and the structures reported in recent scientific publications. 

## How to use
The `data` directory contains all the data needed to build the model. These data are saved in csv files that must be read during model building.

The code can be found in the repository, it can be downloaded as ZIP by clicking in the green Code button. There are three Jupyter notebooks in which the process is fully described. The order in which these notebooks should be read is as follows:
- Copernicus.ipynb. I
- 2D geological cross-sections.ipynb
- 3Dmodel_color.ipynb
- 3D Bolck.ipynb

The `figures` directory contains png image files with the cross section representations and html files with the 3d models. Note that the html files in this directory are larger than GitHug supports, so to be viewed they must be downloaded and used in a local web browser


### Download the code

The code can be found in the repository, it can be downloaded as ZIP by clicking in the green Code button. The necessary files are the four Jupyter notebook mentioned above. In these notebook we use libraries and packages like: plotlib, matplotli, pandas, numpy, sympy, etc. that are included at the begining. But we also use the files copernico.py and bezier.py from the repository https://github.com/torresjrjr/Bezier.py. This file is also hosted here.

How to run the notebook

You can run the notebooks in Jupyter-Notebook or Visual Studio Code, you also need a Python kernel installed in your computer. We recommend installing Anaconda and launching Jupyter by typing jupyter-notebook in the Anaconda Prompt.

To successfully run the notebook, you need to locate it in the same folder as the data directory and the file bezier.py. In order to do this, you may just extract the ZIP file with the whole repository. Then, launch Jupyter Notebook and select the notebook  visualizing_an_imbricate_thrust_system.ipynb. To run a cell, you can just click in the run button (next to the cell number) or click on it and press Ctrl+Enter. You're now ready to go!

Authors:

Manuel Bullejos, Manuel Martín-Martín
