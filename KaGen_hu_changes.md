# Changes and additions for KaGen_hu

This document describes changes and adaptations to [KaGen](https://github.com/sebalamm/KaGen) that
are necessary to output large graphs and their coordinates.

KaGen can produce random delaunay triangulation (rdg) and random geometric (rgg) graphs.
Some tools like Geographer, also require the coordinates of these graphs but KaGen does not stores
them. Also, it had some problems when the graphs needed to be written to a file because, although the 
graph creation is done in parallel, the graph was gathered in one PE and this PE writes the graph.
But this approach has memory problems when the graphs are too big.