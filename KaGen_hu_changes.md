# Changes and additions for KaGen_hu

This document describes changes and adaptations to [KaGen](https://github.com/sebalamm/KaGen) that
are necessary to output large graphs and their coordinates.

KaGen can produce random delaunay triangulation (rdg) and random geometric (rgg) graphs.
Some tools like Geographer, also require the coordinates of these graphs but KaGen does not stores
them. Also, it had some problems when the graphs needed to be written to a file because, although the 
graph creation is done in parallel, the graph was gathered in one PE and this PE writes the graph.
But this approach has memory problems when the graphs are too big.

For this purpose, we added functionality for writting the coordinates in a file and adapted the
way the graph is stored when using multiple PEs.

## Details

For 2d rdg and rgg graphs, the idea is to abuse the function that writes the edges `GatherPrint`
and write also the coordinates. 
Although this works for 2D, it does not work for 3D as the edge list IO takes a list
with pairs of numbers. For that we have overloaded the same function `GatherPrint` for both 
`GatherPrint( tuple<float,float>, ... )` and `GatherPrint( tuple<float,float,float>, ...)`.
A problem is that the user needs to specify the dimension at compile time using the flag
`DKAGEN_DIMENSION_2D=ON/OFF` for 2D or 3D. This way, only 2 or 3 dimensional graph can be 
created per installation.

**TODO**: create a class dedicated for the coordinates and do not abuse `GatherPrint`. This
will avoid the need to specify the supported number of dimensions at compile time.

## Examples

### Installation

Create a `build` folder and in the build folder call `cmake` as:

> cmake <br>
-DEXTERNAL_INCLUDE_DIRS=/home/Code/sparehash/include \\ <br>
-DCMAKE_INSTALL_PREFIX=../installation \\ <br>
-DKAGEN_OMIT_HEADER=OFF \\ <br>
-DKAGEN_OUTPUT_EDGES=ON \\ <br>
-DKAGEN_SINGLE_OUTPUT=ON \\<br>
-DKAGEN_DIMENSION_2D=ON \\ <br>
..

Then call 

> make -j 4 install

For 3D graph you must set `-DKAGEN_DIMENSION_2D=OFF`; this implies that 3D is ON.

### Running the generator

New parameter:

- `output_dir`: pass a folder path where both the graph file and coordinates will be stored.
- `output_format`: graph and coordinates can also be stored in binary
- `coord_output`: name of the coordinates file. If none is given, the name of the graph and coordinates
files are created automatically based in the input parameters (generator and n)

The following command uses 16 PEs to create a 2D random delaunay graph with 2^26 vertices 
and stores the graph and the coordinates files in `/scratch/usr/meshes/rdg_2d_26.bgf`
and `/scratch/usr/meshes/rdg_2d_26.bgf.xyz` respectively.

> mpirun -n 16 kagen -output_dir /scratch/usr/meshes -output_format binary  -n 26 -gen rdg_2d

The following command creates a 3D random geometric graph (remember: in this case, you
should set `-DKAGEN_DIMENSION_2D=OFF` during compilation otherwise an error message is thrown).
The `output_format` parameter is not given and the graph will be written in text using the
METIS graph format.

> mpirun -n 16 kagen -output_dir /scratch/usr/meshes -n 26 -gen rgg_3d




