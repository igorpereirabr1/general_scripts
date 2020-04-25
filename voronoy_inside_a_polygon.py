
from shapely.geometry.polygon import Polygon,Point
from shapely.geometry import asPoint
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
import numpy as np
from collections import defaultdict

class voronoi_analysis():
  def __init__(self,polygon):
    self._polygon = polygon
    return None

  def generate_random(self,number=100):
    """Generate shapely.geometry.Point n randon objects corresponding to the
    regions of a Generate shapely.geometry.Polygon area object, based on a 
    np.random.uniform distribution.

    """
    polygon = self._polygon

    list_of_points = []
    minx, miny, maxx, maxy = polygon.bounds
    counter = 0
    while counter < number:
        pnt = Point(np.random.uniform(minx, maxx), np.random.uniform(miny, maxy))
        if polygon.contains(pnt):
            list_of_points.append(pnt)
            counter += 1
    self.list_of_points = list_of_points
    return self.list_of_points

  def coords_to_points(self,coords):
    """Convert a NumPy array of 2D coordinates `coords` to a list of shapely Point objects"""
    return list(map(asPoint, coords))

  def voronoi_polygons(self,voronoi, diameter):
    """Generate shapely.geometry.Polygon objects corresponding to the
    regions of a scipy.spatial.Voronoi object, in the order of the
    input points. The polygons for the infinite regions are large
    enough that all points within a distance 'diameter' of a Voronoi
    vertex are contained in one of the infinite polygons.

    """
    centroid = voronoi.points.mean(axis=0)

    # Mapping from (input point index, Voronoi point index) to list of
    # unit vectors in the directions of the infinite ridges starting
    # at the Voronoi point and neighbouring the input point.
    ridge_direction = defaultdict(list)
    for (p, q), rv in zip(voronoi.ridge_points, voronoi.ridge_vertices):
        u, v = sorted(rv)
        if u == -1:
            # Infinite ridge starting at ridge point with index v,
            # equidistant from input points with indexes p and q.
            t = voronoi.points[q] - voronoi.points[p] # tangent
            n = np.array([-t[1], t[0]]) / np.linalg.norm(t) # normal
            midpoint = voronoi.points[[p, q]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - centroid, n)) * n
            ridge_direction[p, v].append(direction)
            ridge_direction[q, v].append(direction)

    for i, r in enumerate(voronoi.point_region):
        region = voronoi.regions[r]
        if -1 not in region:
            # Finite region.
            yield Polygon(voronoi.vertices[region])
            continue
        # Infinite region.
        inf = region.index(-1)              # Index of vertex at infinity.
        j = region[(inf - 1) % len(region)] # Index of previous vertex.
        k = region[(inf + 1) % len(region)] # Index of next vertex.
        if j == k:
            # Region has one Voronoi vertex with two ridges.
            dir_j, dir_k = ridge_direction[i, j]
        else:
            # Region has two Voronoi vertices, each with one ridge.
            dir_j, = ridge_direction[i, j]
            dir_k, = ridge_direction[i, k]

        # Length of ridges needed for the extra edge to lie at least
        # 'diameter' away from all Voronoi vertices.
        length = 2 * diameter / np.linalg.norm(dir_j + dir_k)

        # Polygon consists of finite part plus an extra edge.
        finite_part = voronoi.vertices[region[inf + 1:] + region[:inf]]
        extra_edge = [voronoi.vertices[j] + dir_j * length,
                      voronoi.vertices[k] + dir_k * length]
        yield Polygon(np.concatenate((finite_part, extra_edge)))

  def generate_voronoy_polygons(self,number=100,plot=True):

      points = np.array([[k.x,k.y] for k in self.generate_random(number=number)])
      boundary = np.array([[k.x,k.y] for k in self.coords_to_points(self._polygon.boundary.coords.xy)])
      diameter = np.linalg.norm(boundary.ptp(axis=0))

      x, y = boundary.T

      self._voronoi_sub_polygons={}
      count = 0
      

      for p in self.voronoi_polygons(Voronoi(points), diameter):
          x, y = zip(*p.intersection(self._polygon).exterior.coords)  
          #append the sub-polygon in the dictionary
          self._voronoi_sub_polygons[count] = p.intersection(self._polygon)
          count+=1
      if plot:
          plt.figure(figsize=(8,8))
          for key,value in self._voronoi_sub_polygons.items():
              plt.plot(*value.exterior.xy,color='red')
              plt.scatter(x = value.centroid.x,y = value.centroid.y,color='blue',s=0.8)
              plt.plot(x, y, 'r-')
         #for databricks showin image
          plt.show()
          try:
              display()
          except:
              pass

      return self._voronoi_sub_polygons


#Example of usage

#coords to create a polygon
a0 = [-16.534427,	-71.606101] 
a1 = [-16.533385,	-71.607572]
a2 = [-16.532789,	-71.607054]
a3 = [-16.533399,	-71.606506]
a4 = [-16.533726,	-71.605479]

a_lats_vect = np.array([a0[0],a1[0],a2[0],a3[0],a4[0]])
a_lons_vect = np.array([a0[1],a1[1],a2[1],a3[1],a4[1]])

polyg = Polygon(np.column_stack((a_lons_vect, a_lats_vect)))

test = voronoi_analysis(polygon=polyg)

areas = test.generate_voronoy_polygons(number=150,plot=True)
