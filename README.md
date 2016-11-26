# Voronoi
An implementation of Fortune's Line Sweep algorithm for Voronoi Diagram and Delaunay Triangulation generation.
![Pretty Gif](/result.gif?raw=true "Fortune's Algorithm")

## What's a Voronoi Diagram?
Wikipedia puts it best.

> In mathematics, a Voronoi diagram is a partitioning of a plane into regions based on distance to points in a specific subset of the plane. That set of points (called seeds, sites, or generators) is specified beforehand, and for each seed there is a corresponding region consisting of all points closer to that seed than to any other. These regions are called Voronoi cells. The Voronoi diagram of a set of points is dual to its Delaunay triangulation.

https://en.wikipedia.org/wiki/Voronoi_diagram

## What's a Delaunay Triangulation?
Back to Wikipedia...

> In mathematics and computational geometry, a Delaunay triangulation for a set P of points in a plane is a triangulation DT(P) such that no point in P is inside the circumcircle of any triangle in DT(P). Delaunay triangulations maximize the minimum angle of all the angles of the triangles in the triangulation; they tend to avoid sliver triangles. The triangulation is named after Boris Delaunay for his work on this topic from 1934.

https://en.wikipedia.org/wiki/Delaunay_triangulation

## Is this implementation robust?
Not yet! Here's an example of a failure case:
![Pretty Gif](/images/nonrobust.gif?raw=true "Brokenness")
