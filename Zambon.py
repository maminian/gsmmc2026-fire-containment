import numpy as np
import shapely as shp
import shapely.ops as shp_ops
import scipy.optimize as spo
import shapely.plotting as shp_plotting










def generateGeometricData(points):
  '''
  INPUT: points = a sequence of points (ex: [(0, 0), (1, 1), ...] defining a polygon
  note: do not loop back to initial point.

  OUTPUT: a shapely Polygon object representing the forest boundary
  '''
  # Define a polygon P, which is the boundary of the 'forest'.
  # cost scales with complexity of P; keep P limited in scope for now.
  # while code is in development, just take a test shape.
  V_cycle = list(V)
  V_cycle.append(V[0])
  P = shp.Polygon(V_cycle)
  dP = P.boundary # convenience

  #print(V)
  #print(V_cycle)
  #shp_plotting.plot_polygon(P)




  # generate array of barrier lines
  B = []

  for i in np.arange(0, len(V)):
    for j in np.arange(i+1, len(V)):
      line = shp.LineString((V[i], V[j]))
      if P.contains(line) and not dP.contains(shp.snap(line, dP, TOLERANCE)) and line.length >= TOLERANCE * 10:
        B.append(line)#np.append(B, line)



  #shp_plotting.plot_polygon(P)
  #shp_plotting.plot_line(shp.MultiLineString(B))
  #print(str(len(B)) + " barriers generated")

  return (P, B)






def zambon(P, B, delta, FIRESTART = shp.Point(0, 0), CONSPEED = 1, TOLERANCE = 10**(-3)):
  '''
  #INPUT:
  #P, B generated from generateGeometricData
  # FIRESTART = starting point of data, as a shp.Point() object
  #CONSPEED = how fast barriers can be built
  # TOLERANCE = largely irrelevant; some small tolerance for geometrical distance, to handle tiny numerical errors. Can ignore


  OUTPUT:
  (res, output)
  res = formal results of the optimization (not very necssary)
  output = the containment barriers to draw, in order.

  '''
  def trim(arr, ind):
    return [arr[i] for i in ind]
  # construction time vector
  p = []
  for b in B:
    p.append(b.length / CONSPEED)

  ftime = [delta[i] - p[i] for i in np.arange(len(delta))]
  #print(delta)
  # construct basic graphs:

  Gi = np.zeros((len(B), len(B)))
  Gc = np.zeros((len(B), len(B)))
  Gp = np.zeros((len(B), len(B)))


  def containsFirestart(collection):
    # input: a GeometryCollection of polygons
    # output: the (first) polygon that contains the face containing FIRESTART
    # false if none found

    for face in collection.geoms:
      if face.contains(FIRESTART):
        return face
    return False

  def noContainsFirestart(collection):
    # input: a collection of polygons
    # output: the (first) polygon NOT containing FIRESTART
    # false if none found
    for face in collection.geoms:
      if not face.contains(FIRESTART):
        return face
    return False

  # generate a record of the region protected by each barrier
  prot = [] #the ith entry is the area protected by barrier b
  for b in B:
    split1 = shp_ops.split(P, b)
    #print(len(split1.geoms))
    prot.append(noContainsFirestart(split1))

  # generate a record of the areas
  protarea = [f.area for f in prot]


  for r in np.arange(len(B)):
    for c in np.arange(r + 1, len(B)):
      l1 = B[r]
      l2 = B[c]
      # determine whether they intersect.
      # there is a tiny buffer to avoid cases where they intersect at a polygon vertex
      # and a double check just in case.
      if l1.intersects(l2):
        intsc = shp.intersection(l1, l2)
        if P.buffer(-0.01).contains(intsc):
          Gi[r][c] = 1
          Gi[c][r] = 1
      # check whether the two can possibly be mutually constructed
      if p[r] + p[c] > min(delta[r], delta[c]):
        Gc[r][c] = 1
        Gc[c][r] = 1
      # check whether one precedes the other
      face1 = prot[r]
      face2 = prot[c]

      if face1.contains(face2.buffer(-0.01)):
        Gp[r][c] = 1
      elif face2.contains(face1.buffer(-0.01)):
        Gp[c][r] = 1

  Af = shp_ops.split(P, shp.MultiLineString(B))
  protby = []

  # for each face- find which barriers protect it, and record it
  #protby[i] is a list of the INDICIES of the barriers that protect the ith face
  for i in range(len(Af.geoms)):
    arr = []
    f = Af.geoms[i]
    for j in np.arange(len(B)):
      if prot[j].contains(f.buffer(-0.01)):
        arr.append(j)
    protby.append(arr)


  # get indices for sorting the areas in descending order
  protarea_argsort = np.argsort(protarea)[::-1]
  protarea[protarea_argsort[1]]

  # RUN GREEDY HEURISTIC

  greedybans = np.array([0 for i in np.arange(len(B))])
  greedyqueue = []
  greedytime = 0
  for i in protarea_argsort:
    if greedybans[i] == 0 and delta[i] - p[i] >= greedytime:
      greedyqueue.append(B[i])
      greedytime += p[i]

      # remove any preceding or preceded barriers from consideration
      greedybans[Gp[i,:] == 1] = 1
      greedybans[Gp[:,i] == 1] = 1
  #shp_plotting.plot_polygon(P)
  #shp_plotting.plot_line(shp.MultiLineString(greedyqueue), color = "green")
  #shp_plotting.plot_points(FIRESTART, color = "red")

  greedyburned = containsFirestart(shp_ops.split(P, shp.MultiLineString(greedyqueue)))
  greedysaved = P.area - greedyburned.area
  

  #print("Area saved by Greedy Heuristic = " + str(greedysaved))


  ## Barrier Elimination Algorithm
  # practically speaking I think this does very little unless the polygon is quite complex
  # this primarily serves as a computational cost saving step
  validB = []
  for i in np.arange(len(B)):
    okBarriers = [B[j] for j in np.arange(len(B)) if Gp[j][i] == 0 and i != j]
    okBarriers.append(B[i])
    minBurned = containsFirestart(shp_ops.split(P, shp.MultiLineString(okBarriers)))
    #shp_plotting.plot_polygon(P)
    #shp_plotting.plot_line(shp.MultiLineString(okBarriers))
    if not minBurned or minBurned.area < greedyburned.area: #slightly scuffed alg: sometimes minBurned is false, which is probably a result of cutting the firestart region so finely that numerical errors or buffering cause it to not be detected. it's assumed that when this happens the burned area is very small.
      validB.append(i)

  #print(validB)

  # trim down all barrier tracking to valid barriers only
  B = trim(B, validB)
  p = trim(p, validB)
  ftime = trim(ftime, validB)
  delta = trim(delta, validB)
  Gp = Gp[validB, :]
  Gp = Gp[:,validB]
  Gc = Gc[validB, :]
  Gc = Gc[:,validB]
  Gi = Gi[validB, :]
  Gi = Gi[:,validB]
  protarea = trim(protarea, validB)


  # INTEGER PROGRAMMING


  #sco.milp
  # c = protarea


  # CONSTRAINTS 16 to 17

  def genRow(ind, n):
    x = [int(k in ind) for k in range(n)]
    return x

  A1 = [[0 for k in range(2 * len(B))]] #initialize this with a 0 row just to handle aberrant cases where no other constraints are generated
  # so that milp doesn't get mad about a empty array.


  for r in np.arange(len(B)):
    for c in np.arange(r + 1, len(B)):
      if Gi[r][c] == 1 or Gp[r][c] == 1 or Gp[c][r] == 1: # directed, so check both
        A1.append(genRow([r,c], 2 * len(B)))

  A1 = np.array(A1)
  con1 = spo.LinearConstraint(A1, ub = 1)



  # CONSTRAINT 13

  A2 =  []
  bound = []

  for j in np.arange(len(B)):
    bd = 0
    row = [0 for k in range(2 * len(B))]
    row[len(B) + j] = 1
    row[j] -= p[j]
    for i in np.arange(len(B)):
      if delta[i] <= delta[j]:
        bd -= p[i]
        row[i] -= p[i]
        row[j] -= p[i]
    bound.append(bd)
    A2.append(row)
  A2 = np.array(A2)
  #print(A2.shape)
  #print(len(bound))
  con2 = spo.LinearConstraint(A2, lb = bound)



  # CONSTRAINT 9

  A3 = []
  for j in np.arange(len(B)):
    row = [0 for k in range(2 * len(B))]
    row[len(B) + j] = 1
    row[j] = - delta[j]
    A3.append(row)

  A3 = np.array(A3)
  con3 = spo.LinearConstraint(A3, ub = 0)


  integrality_constraint = genRow(list(range(len(B))), 2*len(B))
  var_bounds = spo.Bounds([0 for i in range(2*len(B))], [np.inf if i >= len(B) else 1 for i in range(2*len(B))])

  # because MILP is a minimization algorithm; set protected areas to negative
  negprotarea = [-1 * x for x in protarea]
  cvec = negprotarea + [0 for i in range(len(B))]


  res = spo.milp(c = cvec,
          integrality = integrality_constraint,
          bounds = var_bounds,
          constraints = (con1, con2, con3)
          )
  # output whether it was successful

  print(res['message'])
  y = res['x']
  #print(y)

  sol = [B[i] for i in range(len(y)) if y[i] == 1 ]




  # get correct order of paths to draw:
  barriers = y[0:len(B)]
  c = y[len(B):len(B)*2]
  NB = sum(barriers == 1)
  sortind = np.argsort(c)[::-1]
  output = []
  for i in range(NB):
    output.append(B[sortind[i]])

  return (res, output)




  '''
  shp_plotting.plot_polygon(P)
  shp_plotting.plot_line(B[1])
  shp_plotting.plot_polygon(prot[1], color = "green")
  '''
  #print(Gp)


  
# Build Time Limits
# this will be generated via simulation in the correct run. For now we're going to let the fire expand at a constant rate and clip through walls.



# BASIC VARIABLES

TOLERANCE = 10**(-3)
CONSPEED = 20
FIRESTART = shp.Point(-4,0)

V = [
    shp.Point(1,2),
    shp.Point(3,7),
    shp.Point(5,3),
    shp.Point(6,0),
    shp.Point(4,-4),
    shp.Point(0,-2),
    shp.Point(-2,-5),
    shp.Point(-5,-1),
    shp.Point(-3,4),
    shp.Point(0, 2)
]


P, B = generateGeometricData(V)
FIRESPEED = 2
delta = [shp.distance(b, FIRESTART)/FIRESPEED for b in B]

res, out = zambon(P, B, delta, FIRESTART = FIRESTART, CONSPEED = 20)

shp_plotting.plot_polygon(P)
shp_plotting.plot_line(shp.MultiLineString(out), color = "green")

