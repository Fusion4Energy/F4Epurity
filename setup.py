from setuptools import setup, find_packages

setup(
    name='f4epurity',
    version='0.1.0',
    author='UKAEA - F4E',
    description='Tool to quantify the impact on the shut down dose rate of a deviation in impurity content of a material',
    install_requires=[
        'pyevtk',
        'matplotlib',
        'numpy',
        'openpyxl',
        'pyvista',
        'pandas',
    ],
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'f4epurity = f4epurity.main:main',
            'f4epurity-activity = f4epurity.global_activity_map:main',
            'f4epurity-xs = f4epurity.global_effective_xs_map:main',
        ],
    },
)
