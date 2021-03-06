AFRODITE - PR293A

A GEANT4 simulation of the AFRODITE vault at iThemba LABS, composed of the AFRODITE array of HPGe clover detectors along with various ancillary detectors

The primary objective of this simulation is to write the first standard simulation package specifically for AFRODITE that can be utilised by all students/staff that are interested. The compilation of the many different detectors (and other experimental equipment) that are available at iThemba is also of great import.

I would like to thank the numerous staff at iThemba who put up with my requests to dig through their old documents in the hopes of finding such treasures as detector drawings. Even greater are my thanks to those who have themselves delved into the depths of archives. I hope that this project for iThemba repays your effort.

////////////////////////////////////////////////////////////////////////////////////////////////////

The current implemented detectors include the following:

- Clover HPGe detectors (Eurogam type), manufactured by Canberra.
- The BGO shield complimenting the Clover HPGe detectors (Eurogam type), manufactured by Eurisys.
- LEPS HPGe detectors, manufactured by Eurisys.
- A neutron wall, composed of a total of 12 plastic scintillators (NE102).
- The "new" scattering chamber built for AFRODITE

////////////////////////////////////////////////////////////////////////////////////////////////////

A CAD model import interface called CADMesh (authored primarily by Christopher Poole) is used within this simulation for implementing complex geometries. This needs to be obtained from the following GitHub repository: https://github.com/christopherpoole/CADMesh

Alternatively, one could comment out the relevant CADMesh associated code and use only the hard-coded geometrical objects. It should be noted that an effort has been made to hard-code the geometries in GEANT4 when possible as this has computational advantages.

Adding more text here about the experiment PR239A
