## Examples

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