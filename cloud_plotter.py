import vtk
import csv
import numpy
import math

#http://wfs-beta.geonet.org.nz/geoserver/geonet/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=geonet:quake&outputFormat=csv&cql_filter=BBOX%28origin_geom,174,-41,175,-42%29+AND+origintime%3E=%272013-07-18%27+AND+magnitude%3E3

class Points_Maker(object):
    def __init__(self, fname):
        self.fname = fname
        self.multiplyer = 1
        self.latitude = 0.001
        self.logitude = None
        self.depth = None
        self.earth_radius = 6373 #ref radius @ sea level (km)
        self.long_min = 174
        self.long_max = 175
        self.lat_min = -41
        self.lat_max = -42
        self.lat_max_distance = None
        self.long_max_distance = None
        self.points = []
        self.points_mag = []
        self.max_depth = None
        self.bounds = None
        self.get_data()

    def set_lat_and_long_max_distance(self):
        bound = (self.lat_min, self.long_max, 0)
        self.lat_max_distance = self.convert_bound_to_distance(bound)
        bound = (self.lat_max, self.long_min, 0)
        self.long_max_distance = self.convert_bound_to_distance(bound)


    def make_bounds(self):
        self.bounds = [[self.lat_min, self.long_min, 0], \
                       [self.lat_min, self.long_max, 0], \
                       [self.lat_max, self.long_min, 0], \
                       [self.lat_max, self.long_max, 0], \
                       [self.lat_min, self.long_min, self.max_depth], \
                       [self.lat_min, self.long_max, self.max_depth], \
                       [self.lat_max, self.long_min, self.max_depth], \
                       [self.lat_max, self.long_max, self.max_depth]]

    def convert_depth_to_relative_radius(self, depth):
        relative_radius = self.earth_radius - depth
        return relative_radius

    def distance_on_unit_sphere(self, lat1, long1, lat2, long2, depth):
        # Convert latitude and longitude to
        # spherical coordinates in radians.
        degrees_to_radians = math.pi/180.0

        # phi = 90 - latitude
        phi1 = (90.0 - lat1)*degrees_to_radians
        phi2 = (90.0 - lat2)*degrees_to_radians

        # theta = longitude
        theta1 = long1*degrees_to_radians
        theta2 = long2*degrees_to_radians

        # Compute spherical distance from spherical coordinates.

        # For two locations in spherical coordinates
        # (1, theta, phi) and (1, theta, phi)
        # cosine( arc length ) =
        #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
        # distance = rho * arc length

        cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) +
               math.cos(phi1)*math.cos(phi2))
        arc = math.acos( cos )
        # Remember to multiply arc by the radius of the earth
        # in your favorite set of units to get length.

        relative_radius = self.convert_depth_to_relative_radius(depth)
        distance = arc * relative_radius
        return distance

    def make_x_point(self, lat, depth):
        distance = self.distance_on_unit_sphere(self.lat_min, self.long_min, lat, self.long_min, depth)
        return distance

    def make_y_point(self, longi, depth):
        distance = self.distance_on_unit_sphere(self.lat_min, self.long_min, self.lat_min, longi, depth)
        return distance

    def get_data(self):
        reader = csv.reader(open(self.fname, "rb"))
        for row in reader:
            if "FID" in row[0]:
                pass
            else:
                self.longitude = float(row[3])
                self.latitude = float(row[4])
                self.depth = float(row[5])
                magnitude = float(row[6])
                x_point = self.make_x_point(self.latitude, self.depth)
                y_point = self.make_y_point(self.longitude, self.depth)
                point = [x_point,  y_point, self.depth]
                point_and_mag = [[x_point*self.multiplyer,  y_point*self.multiplyer, self.depth*self.multiplyer], magnitude]
                point = numpy.asarray(point)
                point_and_mag = numpy.asarray(point_and_mag)
                self.points.append(point)
                self.points_mag.append(point_and_mag)
                if self.depth > self.max_depth:
                    self.max_depth = self.depth * self.multiplyer

    def convert_bound_to_distance(self, bound):
        x, y, depth = bound
        x_point = self.make_x_point(x, depth)
        y_point = self.make_y_point(y, depth)
        bound = (x_point*self.multiplyer, y_point*self.multiplyer, depth*self.multiplyer)
        return bound

class VtkPointCloud:
    def __init__(self, pm, zMin=-0.0, zMax=156, maxNumPoints=1e6): #sets colour limits
        self.maxNumPoints = maxNumPoints
        self.vtkPolyData = vtk.vtkPolyData()
        self.clearPoints()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(self.vtkPolyData)
        mapper.SetColorModeToDefault()
        mapper.SetScalarRange(zMin, pm.max_depth)
        mapper.SetScalarVisibility(1)
        self.vtkActor = vtk.vtkActor()
        self.vtkActor.SetMapper(mapper)

    def addPoint(self, point, mag):
        if self.vtkPoints.GetNumberOfPoints() < self.maxNumPoints:
            pointId = self.vtkPoints.InsertNextPoint(point[:])
            self.vtkDepth.InsertNextValue(point[2])
            self.vtkCells.InsertNextCell(1)
            self.vtkCells.InsertCellPoint(pointId)
            self.vtkActor.GetProperty().SetPointSize(mag)
        self.vtkCells.Modified()
        self.vtkPoints.Modified()
        self.vtkDepth.Modified()
        self.vtkActor.Modified()

    def clearPoints(self):
        self.vtkPoints = vtk.vtkPoints()
        self.vtkCells = vtk.vtkCellArray()
        self.vtkDepth = vtk.vtkDoubleArray()
        self.vtkDepth.SetName('DepthArray')
        self.vtkPolyData.SetPoints(self.vtkPoints)
        self.vtkPolyData.SetVerts(self.vtkCells)
        self.vtkPolyData.GetPointData().SetScalars(self.vtkDepth)
        self.vtkPolyData.GetPointData().SetActiveScalars('DepthArray')


def main():
    pm = Points_Maker("quake.csv")
    renderer = vtk.vtkRenderer()
    pointCloud_bounds = VtkPointCloud(pm)
    for point_and_mag in pm.points_mag:
        mag = point_and_mag[1]
        point = point_and_mag[0]
        if mag > 4:
            pointCloud = VtkPointCloud(pm)
            pointCloud.addPoint(point, math.log(mag)*10)
            renderer.AddActor(pointCloud.vtkActor)
    pm.make_bounds()
    for bound in pm.bounds:
        bound = pm.convert_bound_to_distance(bound)
        pointCloud_bounds.addPoint(bound, 10)  # the number is what makes the bound markers larger
    renderer.AddActor(pointCloud_bounds.vtkActor)
    renderer.SetBackground(.2, .3, .3)
    renderer.ResetCamera()

    # apply texture
    input_image = "overlay.jpg"

    # Read the image data from a file
    reader = vtk.vtkJPEGReader()
    reader.SetFileName(input_image)

    # Create texture object
    texture = vtk.vtkTexture()
    texture.SetInput(reader.GetOutput())

    texture.SetInputConnection(reader.GetOutputPort())
    texture.InterpolateOn()
    pointCloud_overlay = VtkPointCloud(pm)
    pm.set_lat_and_long_max_distance()
    plane = vtk.vtkPlaneSource()

    print  pm.long_max_distance
    print pm.lat_max_distance
    plane.SetOrigin((0, 0, 0))
    plane.SetPoint1(pm.lat_max_distance)
    plane.SetPoint2(pm.long_max_distance)


    planeMapper = vtk.vtkPolyDataMapper()
    planeMapper.SetInputConnection(plane.GetOutputPort())
    planeActor = vtk.vtkActor()
    planeActor.SetMapper(planeMapper)
    planeActor.SetTexture(texture)

    renderer.AddActor(planeActor)
    # Render Window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # Interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Begin Interaction
    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == '__main__':
    main()
