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