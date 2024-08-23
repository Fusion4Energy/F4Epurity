# F4Epurity

## Description

F4Epurity is a python-based tool that can be used to approximate the change in activity and dose resulting from a deviation in chemical composition of a material under irradiation. 

The tool works for a point or line deviation and outputs a 3D map of the deviation in dose in *vtr* format. This dose field assumes unshielded and unscattered conditions.  

## Installation

It is recommended that the user sets up a virtual environment to install the code. The following steps detail set up of a [virtual environment](https://docs.python.org/3/tutorial/venv.html) and installation of the tool.

```bash
python3 -m venv env
```
Activate the virtual environment. When activated, the package and all of its dependancies will be installed into the virtual environment. 

```bash
source env/bin/activate
```
Clone the repository containing the code from the [GitHub repository](https://git.ccfe.ac.uk/avalenti/f4e_impurity_assessment/-/tree/master/).

```bash 
git clone git@git.ccfe.ac.uk:username/f4e_impurity_assessment.git
```

Enter the directory containing the code which can now be installed. 

```bash 
pip install .
```

Each time a user launches a new window or terminal, they need to make sure that the virtual environment is activated using the above *source* command. 

## Usage


Once installed, the tool is straightforward to use. The main tool can be ran in the command line.

```bash
f4epurity -h
```

This will show help menu in the command line which details the required and optional arguments for the tool. 

It is also possible to plot global maps of the collapsed cross section and the activity. These tools can be ran with the following commands.

```bash
f4epurity-xs -h
```
```bash
f4epurity-activity -h
```
All of the output files from running the tool are placed into a timestamped output folder. The folder will contain a *metadata.json* file which logs the input commands that were run.

## Example Usage 

Currently the tool supports specification of a point source or a line source. 

An example of how the tool would be used for a *0.05%* deviation in *Cobalt* content at the location with coordinates, *x=250, y=3124, z=10* is shown below. The path to the *VTR* file containing the input flux spectrum (ITER ex-bioshield, 5 group) is given and the assumed irradiation scenario is *SA2*. The deviation in dose is request at *1E6 seconds*.

```bash
f4epurity --element Co --delta_impurity 0.05 --input_flux ./flux_spectrum.vtr --x1 250 --y1 3124 --z1 10 --irrad_scenario SA2 --decay_time 1e6
```

With the same input parameters, but a line source that extends from *x1,y1,z1* to *x2=300, y2=3500, z2=50* is shown below.

```bash 
f4epurity --element Co --delta_impurity 0.05 --input_flux ./flux_spectrum.vtr --x1 250 --y1 3124 --z1 10 --x2 300 --y2 3500 --z2 50 --irrad_scenario SA2 --decay_time 1e6
```

Multiple source terms may also be specified in a single run of the tool. For example, the deviation in impurity content may be in multiple independent pipes. The following is an example where multiple line sources are defined. 

```bash 
f4epurity --element Co --delta_impurity 0.05 --input_flux ./flux_spectrum.vtr --x1 250 -245 --y1 3124 2341 --z1 10 0 --x2 300 412 --y2 3500 2600 --z2 50 5 --irrad_scenario SA2 --decay_time 1e6 --plot z 15
```

The output will contain the *VTR* maps of dose as before as well as a total dose map summing the contribution from the different source terms. In this example the plotting function was also used. This allows a quick method for visualising the dose map which takes two arguments, the plane, *x*, *y* or *z*, and the position. Currently only on-axis slices are supported. 

The user may also request to evaluate the change in dose at certain locations, e.g. workstations. In this case, the commands, *--workstation* and *--location* should be provided, where *workstation* is the number of the specific workstation where the dose is requested and *location* is the name of the ITER room. *'all'* can be specified for the *workstation* to output the result at all workstations in the given location. 

```bash 
f4epurity --element Co --delta_impurity 0.05 --input_flux ./flux_spectrum.vtr --x1 250 -245 --y1 3124 2341 --z1 10 0 --x2 300 412 --y2 3500 2600 --z2 50 5 --irrad_scenario SA2 --decay_time 1e6 --workstation 'all' --location 'NB cell'
```
A *csv* file is output containing the result for the deviation in dose at the requested workstations. For the list of available workstations, the user should consult the *workstations.xlsx* document in the code resources. 

**All results are ouput per gram of material**.

## Reporting Bugs

Any bugs or problems faced running the code should be raised as issues in the GitHub repository. Any proposals for new code features should also be raised as issues.

## License

The code is distributed under XXX.