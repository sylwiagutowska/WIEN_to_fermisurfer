# WIEN to fermisurfer


This program allows to export WIEN data to fermisurfer (*.frmsf) format in order to plot Fermi surface.

There are two modes:
1) o - orbital character - you can color your FS with orbital character (got from qtl file). It produces files FS_atom_x_orbital_x.frmsf - for each atom and each orbital
2) v - velocity - you can color your FS with fermi velocity

This program needs to  be run in a directory, where you have run your DOS calculations.
It uses the informations from *struct, *outputkgen, *energy*, *outputt* files. It additionally runs lapw2 -qtl calculations with more bands included in order to get full informations on qtl charges.
