from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
from PIL import Image, ImageStat, ImageDraw
import matplotlib.path as mplPath
import numpy as np
import itertools

def draw_circle(draw,pt,r,color):
    pt0 = (pt[0]-r, pt[1]-r)
    pt1 = (pt[0]+r, pt[1]+r)
    draw.ellipse([pt0, pt1], fill=color)

def isin_region(poly,x,y):
    bbPath = mplPath.Path(np.array(poly))
    return bbPath.contains_point([x,y])

def find_centroids(A, sz, numGen, rho):
    mX = [0]*numGen
    mY = [0]*numGen
    m  = [0]*numGen
    for y in range(sz[1]):
        for x in range(sz[0]):
            k = A[y][x]
            m[k]  +=   rho(x, y)
            mY[k] += y*rho(x, y)
            mX[k] += x*rho(x, y)
    centroids = [(mX[k]/m[k], mY[k]/m[k]) for k in range(numGen) if m[k] > 0]
    return centroids

def stipple(im, itr):
    print("Start stippling the image. (Itertation={})".format(itr))
    cutoff = ImageStat.Stat(im).mean[0]
    xsize = im.size[0]
    ysize = im.size[1]

    #初期点群の作成
    bxSz = 4
    points = []
    for x in itertools.product(range(0, xsize-int(bxSz/2), bxSz), range(0, ysize-int(bxSz/2), bxSz)):
        box = (x[0], x[1], x[0]+bxSz, x[1]+bxSz)
        region = im.crop(box)
        if (ImageStat.Stat(region).mean[0]/255 < np.random.random()):
            points.append((x[0]+int(bxSz/2), x[1]+int(bxSz/2)))
    

    # Lloyd's methodにより、Weighted Voronoi Diagramを作成（設定した繰り返し回数分回す）
    for it in range(itr):
        print("iteration = {}".format(it+1))
        im_p = Image.new("RGB",(xsize,ysize),(255, 255, 255))
        im2 = im.copy()
        vm = Voronoi(points)
        regions, vertices = voronoi_finite_polygons_2d(vm)
        matrix = [[-1 for x in range(xsize)] for y in range(ysize)]
        num = 0
        print(len(regions))
        for i,reg in enumerate(regions):
            print(i)
            polygon = vertices[reg]
            xmin = max(0,int(np.floor(polygon[:,0].min())))
            xmax = min(xsize,int(np.ceil(polygon[:,0].max())))
            ymin = max(0,int(np.floor(polygon[:,1].min())))
            ymax = min(ysize,int(np.ceil(polygon[:,1].max())))
            print(ymin,ymax,xmin,xmax)
            for y in range(ymin,ymax+1):
                for x in range(xmin,xmax+1):
                    if y >= 0 and y < ysize and x>= 0 and x< xsize:
                        if matrix[y][x] != -1:
                            continue
                        if isin_region(polygon, x,y):
                            num += 1
                            matrix[y][x] = i
        draw = ImageDraw.Draw(im_p)
        for pt in points:
            draw_circle(draw, (pt[0], pt[1]), 0.4, 0)
        im_p.save("stippled"+str(it)+".jpg")
        centroids = find_centroids(matrix, (xsize, ysize), len(points), lambda x, y : 1 - im.getpixel((x, y))/255)
        points = [(pt[0], pt[1]) for pt in centroids]
    return points

def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)