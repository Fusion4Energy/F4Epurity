############ 
Installation
############

**Python versions lower than 3.10 are not supported (tested versions are 3.10 and 3.11). Both Windows and Linux OS are supported.**

F4Epurity is available on PyPi. It can simply be installed using pip:

.. code-block:: bash

    pip install f4epurity

Developer installation
=======================

Clone the repository containing the code from the `GitHub repository <https://github.com/Fusion4Energy/F4Epurity>`_. Note that for users without Git installed, you can instead download an archive of the repository. See `here <https://docs.github.com/en/repositories/working-with-files/using-files/downloading-source-code-archives>`_.

.. code-block:: bash

    git clone git@github.com:Fusion4Energy/F4Epurity.git

Enter the directory containing the code, which can now be installed.

.. code-block:: bash

    pip install .

To do a developer mode install, the following command can be used:

.. code-block:: bash

    pip install -e .

Each time a user launches a new window or terminal, they need to make sure that the virtual environment is activated using the above *source* command.