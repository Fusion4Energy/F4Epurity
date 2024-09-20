[![Testing](https://github.com/Fusion4Energy/F4Epurity/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Fusion4Energy/F4Epurity/actions/workflows/ci.yml)
# F4Epurity

## Description

F4Epurity is a python-based tool that can be used to approximate the change in activity and dose resulting from a deviation in chemical composition of a material under irradiation. 

The tool works for a point or line deviation and outputs a 3D map of the deviation in dose in *vtr* format. This dose field assumes unshielded and unscattered conditions.  

## Installation

**Python versions lower than 3.10 are not supported (tested versions are 3.10 and 3.11). Both Windows and Linux OS are supported.**

It is recommended that the user sets up a virtual environment to install the code. The following steps detail set up of a [virtual environment](https://docs.python.org/3/tutorial/venv.html) and installation of the tool.

```bash
python3 -m venv env
```
Activate the virtual environment. When activated, the package and all of its dependancies will be installed into the virtual environment. 

```bash
source env/bin/activate
```
Clone the repository containing the code from the [GitHub repository](https://github.com/Fusion4Energy/F4Epurity). Note that for users without Git installed, you can instead download an archive of the repository, see [here](https://docs.github.com/en/repositories/working-with-files/using-files/downloading-source-code-archives).

```bash 
git clone git@github.com:Fusion4Energy/F4Epurity.git
```

Enter the directory containing the code which can now be installed. 

```bash 
pip install .
```

To do a developer mode install the following command can be used:

```bash
pip install -e .
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

Currently the tool supports specification of a point source or a line source. A complete list of the accepted arguments for the tool can be found [here](#Extended-description-of-admissible-options).

An example of how the tool would be used for a *0.05%* deviation in *Cobalt* content at the location with coordinates, *x=250, y=3124, z=10* is shown below. The path to the *VTR* file containing the input flux spectrum (ITER ex-bioshield, 5 group) is given. Note that the code expects the spectrum to be in the default format output by [F4Enix](https://f4enix.readthedocs.io/en/latest/index.html) e.g. "Value - ..", "ValueBin-..."). The irradiation scenario is *SA2*. The deviation in dose is requested at *1E6* seconds.

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

The output will contain a total dose map summing the contribution from the different source terms. If the ``--output_all_vtr`` flag is given then the individual *VTR* files for all source terms will also be output. In this example the plotting function was also used. This allows a quick method for visualising the dose map which takes two arguments, the plane, *x*, *y* or *z*, and the position. Currently only on-axis slices are supported. The geometry can be superimposed with this plot if a path to the folder containing the STL file of the geometry is provided. This can be specified using the ``--stl_files`` argument. 

The user may also request to evaluate the change in dose at certain locations, e.g. workstations. In this case, the commands, *--workstation* and *--location* should be provided, where *workstation* is the number of the specific workstation where the dose is requested and *location* is the name of the ITER room. *'all'* can be specified for the *workstation* to output the result at all workstations in the given location. 

```bash 
f4epurity --element Co --delta_impurity 0.05 --input_flux ./flux_spectrum.vtr --x1 250 -245 --y1 3124 2341 --z1 10 0 --x2 300 412 --y2 3500 2600 --z2 50 5 --irrad_scenario SA2 --decay_time 1e6 --workstation 'all' --location 'NB cell'
```
A *csv* file is output containing the result for the deviation in dose at the requested workstations. For the list of available workstations, the user should consult the *workstations.xlsx* document in the code resources. 

The user may have already calculated the activity of certain radionuclides and wish to simply output the resulting dose maps. This is possible using the ``--activities_file`` option. The path to a text file containing a list of the nuclides and their corresponding activities. The activity values should be specified per gram i.e. Bq/g. An example is shown below. *x1*, *y1* and *z1* are also required in this case for the location of the source.

```bash
Co060 2300000
Co060m 43000
Nb094 12000
```

Since the number of options to be provided is high, users can also provide a ``.json`` or a ``.yaml`` configuration using the ``--cfg`` option:
```bash
f4epurity --cfg path/to/config.json
```

The structure of these files is very simple and looks like the following for the .json:
```json
{
    "element": "Ta",
    "delta_impurity": 0.1,
    "input_flux": "dummy",
    "irrad_scenario": "SA2",
    "x1": [-835],
    "y1": [1994],
    "z1": [1230],
    "decay_time": 1e6
}
```
and the .yaml:
```yaml
element: Ta
delta_impurity: 0.1
input_flux: dummy
irrad_scenario: SA2
x1: [-835]
y1: [1994]
z1: [1230]
decay_time: 1e6
```

Finally, if the study requires to assess a large number of source points (or lines), e.g. hundreds of valves in a rooom where the impurity content has been changed, it may become cumbersome to list all the coordinates in either the CLI or even in the configuration file. In this case, the ``--sources_csv`` option can be used which accepts a path to a .csv file containing the locations of the sources, e.g.:
```csv
x1,y1,z1,x2,y2,z2
-835,1994,1230,-835,1994,1231
-834,1994,1230,-834,1994,1231
```
Clearly, if point sources are to be considered, [x2, y2, z2] can be omitted.
The code will throw errors if:
- the column names differ from what reported here
- not all necessary columns are present (either 3 or 6)
- the ``--sources_csv`` and ``--x1 --y2 ...`` options are tried to be provided at the same time

If multiple sources are provided, by default, it is assumed that all the components 
where the impurities are located have the same mass. In this way, the dose deviation maps can be easily sum toghether and the results are provided as mSv/h/g of component. If this assumption is not valid, masses (in g) for each component need to be provided through the ``--m`` option. This data can be also provided through the ``--sources_csv`` option adding a column ``m`` to the .csv file.


**All results are ouput per gram of material**.

## Reporting Bugs

Any bugs or problems faced running the code should be raised as [issues](https://github.com/Fusion4Energy/F4Epurity/issues) in the GitHub repository. Any proposals for new code features should also be raised as issues.

## Troubleshooting

Many issues encountered during installation can be resolved by creating a new virtual environment and reinstalling.

Some users may have some restrictions meaning that the entry points are blocked. In this case it is not possible to invoke directly ``f4epurity -h`` in the command line or run the code as a module (``python -m f4epurity -h``). In this case the tool(s) must be ran in the following way:

```bash
f4epurity ->  python -m f4epurity.main
f4epurity-xs ->  python -m f4epurity.global_effective_xs_map
f4epurity-activity ->  python -m f4epurity.global_activity_map
f4epurity-shielding -> python -m f4epurity.shielding_estimator
```

## Extended description of admissible options
This is a complete list of the available options for the tool and their description.

### Mandatory parameters

- ``--element``: the element name of the impurity to be considered by the tool (e.g. "Ta").
- ``--delta_impurity``: the percentage increase (or decrease) of impurity to be considered in the different components (sources). As an example, if the impurity in the component is increased from 0.5% to 0.6% in mass, ``delta_impurity=0.1``.
- ``--input_flux``: path to the .vtk file containing the energy binned flux to be used for activation calculations. Array name **must** start with "ValueBin-".
- ``--irrad_scenario``: irradiation scenario to be considered in the assessment. "SA2" and "DT1" are shipped with the tool. It is also possible to supply a path to file with user defined scenario.
- ``--decay_time``: decay time in seconds after which the dose deviation needs to be computed.
- ``--x1``, ``--y1``, ``--z1``: x, y, z coordinate(s) of the sources.

### Optional parameters

- `--cfg`: path to a configuration file (json or yaml) where all options can be specified instead of providing them from command line.
- ``--m``: mass(es) in grams of the component containing the impurity.
- ``--x2``, ``--y2``, ``--z2``; x, y, z coordinate(s) of the end of line source(s). They trigger x1, y1 and z1 to be intended as the start of line source(s).
- ``--sources_csv``: path to a .csv file that can help in definig multiple sources. The admissible columns are x1, y1, z1, x2, y2, z2, m.
- ``--root_output``: path to the directory where the tool outputs will be stored. If not specified, a new "output" folder is created locally where the tool is executed.
- ``--activities_file``: path to a text file containing isotopes and their pre-computed activities. This is used for instance to get a quick estimation of the dose field due to a certain distribution of ACPs. In case it is provided, "element", "delta_impurity", "input_flux", "irrad_scenario" and "decay_time" are not needed/ignored.
- ``--plot``: 2 parameters are provided here, the slice_axis and the slice_location. Output image of dose mesh at a given location is printed".
- ``--workstation``: name of the workstation for which to report the max dose e.g. 1, 3, 4 or 'all'.
- ``--location``: location of the workstation(s) e.g. Nb cell
- ``--output_all_vtr``: by default is false. If set to true, a dose map resulting from each single source is printed instead of only the global sum of all th sources. depending on the number of sources and the size of the maps this can slow down a bit the assessement. 