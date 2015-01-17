# SETD8 binding

## Setup for simulation

The setup is performed by the script

```
cd setup/
python setup.py [-nmol integer]
```

The input file that describes the inhibitors contains a list of molecules. The argument *nmol* specifies which one to consider as the 0-based index of the line in the input file. Default is 0.

The receptor is taken from 4FMU.pdb.
