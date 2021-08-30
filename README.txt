This repository contains three directories, covering the following cases:
* Translocation  under fixed force
* Loop extrusion under fixed force
* Loop extrusion under fixed extension

Each directory contains the following:
* A generate.py Python script, used to produce the initial configuration of the DNA, with the SMCC attached.
* A parameters file, containing the most important system parameters.
* A parse.py Python script, used for parsing the above file.
* An input file, which should be called when executing LAMMPS.
* A run.py Python script, with an example usage.
