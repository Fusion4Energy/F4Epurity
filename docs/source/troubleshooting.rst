############### 
Troubleshooting
###############

Many issues encountered during installation can be resolved by creating a new virtual environment and reinstalling.

Some users may have restrictions meaning that the entry points are blocked. In this case, it is not possible to invoke directly ``f4epurity -h`` in the command line or run the code as a module (``python -m f4epurity -h``). Instead, the tools must be run in the following way:

.. code-block:: bash

    f4epurity ->  python -m f4epurity.main
    f4epurity-xs ->  python -m f4epurity.global_effective_xs_map
    f4epurity-activity ->  python -m f4epurity.global_activity_map
    f4epurity-shielding -> python -m f4epurity.shielding_estimator