from __future__ import print_function

import pymap3d as pm
import pcl
import numpy as np
from math import cos, sin, pi, sqrt, radians

import pcl.pcl_visualization

def pointCnv(file):
    convertedPoints = []
    for line in file:
        LLA = line.split()
        # print(LLA)
        XYZ = geodetic2ecef(float(LLA[0]), float(LLA[1]), float(LLA[2]))
        # XYZ = geodetic2enu(float(LLA[0]), float(LLA[1]), float(LLA[2]))
        # XYZ = enu2camera(XYZ[0], XYZ[1], XYZ[2])
        # print(XYZ)
        convertedPoints.append([XYZ[0], XYZ[1], XYZ[2]])


    convertedPoints = np.array(convertedPoints, dtype = np.float32)

    return convertedPoints

# geodetic2ecef reference: https://stackoverflow.com/questions/8981943/lat-long-to-x-y-z-position-in-js-not-working
def geodetic2ecef(lat,lon, alt):
    cosLat = cos(lat * pi / 180.0)
    sinLat = sin(lat * pi / 180.0)
    cosLon = cos(lon * pi / 180.0)
    sinLon = sin(lon * pi / 180.0)

    rad = 6378137.0
    f = 1.0 / 298.257224
    C = 1.0 / sqrt(cosLat ** 2 + (1 - f) ** 2 * sinLat ** 2)
    S = (1.0 - f) ** 2 * C
    h = alt
    x = (rad * C + h) * cosLat * cosLon
    y = (rad * C + h) * cosLat * sinLon
    z = (rad * S + h) * sinLat
    
    return x, y, z

# ecef2enu reference: https://gist.github.com/sbarratt/a72bede917b482826192bf34f9ff5d0b
def ecef2enu(x, y, z, lat0, lon0, h0):
    a = 6378137
    b = 6356752.3142
    f = (a - b) / a
    e_sq = f * (2-f)

    lamb = radians(lat0)
    phi = radians(lon0)
    s = sin(lamb)
    N = a / sqrt(1 - e_sq * s * s)

    sin_lambda = sin(lamb)
    cos_lambda = cos(lamb)
    sin_phi = sin(phi)
    cos_phi = cos(phi)

    x0 = (h0 + N) * cos_lambda * cos_phi
    y0 = (h0 + N) * cos_lambda * sin_phi
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda

    xd = x - x0
    yd = y - y0
    zd = z - z0

    xEast = -sin_phi * xd + cos_phi * yd
    yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd

    return xEast, yNorth, zUp

def geodetic2enu(lat, lon, alt):
    # x, y, z = pm.geodetic2ecef(lat, lon, alt)
    # x, y ,z = ecef2enu(x, y, z, 45.90414414, 11.02845385, 227.5819)
    x, y, z = pm.geodetic2enu(lat, lon, alt, 45.90414414, 11.02845385, 227.5819)

    return x, y, z

def enu2camera(x, y, z):
    Qs = 0.362114
    Qx = 0.374050
    Qy = 0.592222
    Qz = 0.615007
    Rq = [[Qs*Qs+Qx*Qx-Qy*Qy-Qz*Qz, 2*Qx*Qy-2*Qs*Qz, 2*Qx*Qz+2*Qs*Qy], [2*Qx*Qy+2*Qs*Qz, Qs*Qs-Qx*Qx+Qy*Qy-Qz*Qz, 2*Qy*Qz-2*Qs*Qx], [2*Qx*Qz-2*Qs*Qy, 2*Qz*Qy+2*Qs*Qx, Qs*Qs-Qx*Qx-Qy*Qy+Qz*Qz]]
    xyz = [x, y, -z]
    p = np.matmul(xyz, Rq)
    return p

def filterOutliers(cloud, mean, sd):
    filteredCloud = cloud.make_statistical_outlier_filter()
    filteredCloud.set_mean_k(mean)
    filteredCloud.set_std_dev_mul_thresh(sd)
    filteredCloud.set_negative(False)
    filteredCloud = filteredCloud.filter()

    return filteredCloud

def filterLargeComponents(filteredCloud,index,thresh):
    remComponents = filteredCloud.make_kdtree_flann()
    indices, sqrDistances = remComponents.nearest_k_search_for_cloud(filteredCloud, index)
    distances = np.sum(sqrDistances, axis=1)
    blocks = []
    for i in range(np.shape(distances)[0]):
        if (distances[i] < float(thresh)):
            blocks.extend(indices[i])
    uniqueIndices = list(set(blocks))
    filteredCloud = filteredCloud.extract(uniqueIndices, negative=True)

    return filteredCloud

def cylinderSegmentation(filteredCloud,iter):
    seg = filteredCloud.make_segmenter_normals(ksearch=50)
    seg.set_model_type(pcl.SACMODEL_CYLINDER)
    seg.set_optimize_coefficients(True)
    seg.set_normal_distance_weight(0.1)
    seg.set_method_type(pcl.SAC_RANSAC)
    seg.set_max_iterations(iter)

    seg.set_distance_threshold(20)
    seg.set_radius_limits(0, 10)
    indices, model = seg.segment()
    filteredCloud = filteredCloud.extract(indices, negative=False)

    return filteredCloud

def planeSegmentation(filteredCloud, iter):
    fil = filteredCloud.make_passthrough_filter()
    fil.set_filter_field_name("y")
    fil.set_filter_limits(6363082.8, 6363093.0)
    cloud_filtered = fil.filter()
    print("passthrough filter:", cloud_filtered.size)

    seg = cloud_filtered.make_segmenter_normals(ksearch=50)
    seg.set_optimize_coefficients(True)
    seg.set_model_type(pcl.SACMODEL_PLANE)
    seg.set_normal_distance_weight(0.1)
    seg.set_method_type(pcl.SAC_RANSAC)
    seg.set_max_iterations(iter)
    seg.set_distance_threshold(85)
    indices, model = seg.segment()
    cloud_plane = cloud_filtered.extract(indices, negative=False)

    return cloud_plane

def meshlabSave(filteredCloud, filename):
    with open(filename, 'w') as final:
        for point in filteredCloud:
            line = "v " + str(point[0]) + " " + str(point[1]) + " "+ str(point[2])
            final.write(line)
            final.write("\n")

def readCSV():
    convertedPoints =[]
    with open('matlabData.csv', 'r') as CSV:
        for line in CSV:
            XYZ = line.split(',')
            # print(XYZ[0])
            convertedPoints.append([float(XYZ[0]), float(XYZ[1]), float(XYZ[2])])
    convertedPoints = np.array(convertedPoints, dtype = np.float32)

    return convertedPoints

def readRaw():
    with open('final_project_point_cloud.fuse', 'r') as file:
        convertedPoints = pointCnv(file)
        return convertedPoints

def main():
    # convertedPoints = readRaw()
    convertedPoints =readCSV()
    cloud = pcl.PointCloud()
    cloud.from_array(convertedPoints)
    print ("Original point cloud Data.")
    print (cloud)
    pcl.save(cloud, "OriginalCloud.pcd")
    # meshlabSave(cloud, "OriginalCloud.obj")
    # print ("output: 'OriginalCloud.pcd'")

    filteredCloud = filterOutliers(cloud, 50, 5.0)
    print ("\nFiltered without outliers cloud Data")
    print (filteredCloud)
    # meshlabSave(filteredCloud, "OutlierFilteredCloud.obj")

    filteredCloud = filterLargeComponents(filteredCloud, 1000, 5000)
    print ("\nFiltered without large components cloud Data")
    print (filteredCloud)
    # meshlabSave(filteredCloud, "LargeComponentsFilteredCloud.obj")

    filteredCloud = planeSegmentation(filteredCloud, 100)
    print ("\nCloud after plane segmentation")
    print (filteredCloud)

    # filteredCloud = cylinderSegmentation(filteredCloud, 1000)
    # print ("\nCloud after cylindrical segmentation")
    # print (filteredCloud)

    # filteredCloud = segmentation(filteredCloud, cloud_filtered)
    # print(filteredCloud)

    pcl.save(filteredCloud, "final.pcd")
    print ("output: 'final.pcd'")

    # meshlabSave(filteredCloud, 'final.obj')

if __name__ == '__main__':
    main()