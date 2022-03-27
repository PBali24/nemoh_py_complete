# nemoh_py_complete
01.03.19
This is the python script to run the NEMOH Fortran exectuables
for one frequency at a time only for the moment
the script requires NEMOH to be run in the Calculation folder inside a given simulation directory with the creation
of the mesh and Results subfolders
The meshing is done via a modified script MeshTypes.py from Tim Verbreugghe
The version of the script requries a call to the MW_NEM_shared fodler that
has the shell and nemoh helper Function that are used by both this shell and the mw_nem_coupling
shell