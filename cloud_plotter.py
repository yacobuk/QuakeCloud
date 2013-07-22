import vtk
import csv
import numpy
import math

class Points_Maker(object):
    def __init__(self, fname):
        self.fname = fname
        self.points = []
        self.points_mag = []
        self.multiplyer = 100
        self.max_depth = float(155)
        self.bounds = [[0, 0, 0], [0, -self.multiplyer, 0], [self.multiplyer, 0, 0], [self.multiplyer, -self.multiplyer, 0],[0, 0, self.max_depth], [0, -self.multiplyer, self.max_depth], [self.multiplyer, 0, self.max_depth], [self.multiplyer, -self.multiplyer, self.max_depth]]
        self.get_data()

    def get_data(self):
        reader = csv.reader(open(self.fname, "rb"))
        for row in reader:
            if "FID" in row[0]:
                pass
            else:
                longitude = (float(row[3]) - 174) * self.multiplyer #de-localises value
                latitude = (float(row[4]) + 41) * self.multiplyer #de-localises value
                depth = float(row[5])
                magnitude = float(row[6])
                point = [longitude, latitude, depth]
                point_and_mag = [[longitude, latitude, depth], magnitude]
                point = numpy.asarray(point)
                point_and_mag = numpy.asarray(point_and_mag)
                self.points.append(point)
                self.points_mag.append(point_and_mag)

class VtkPointCloud:

    def __init__(self, zMin=-0.0, zMax=155, maxNumPoints=1e6): #sets colour limits
        self.maxNumPoints = maxNumPoints
        self.vtkPolyData = vtk.vtkPolyData()
        self.clearPoints()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(self.vtkPolyData)
        mapper.SetColorModeToDefault()
        mapper.SetScalarRange(zMin, zMax)
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
            #.GetProperty().SetSize(mag)
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
    pointCloud2 = VtkPointCloud()
    for bound in pm.bounds:
        pointCloud2.addPoint(bound, 10)
    for point_and_mag in pm.points_mag:
        mag = point_and_mag[1]
        point = point_and_mag[0]
        if mag > 1:
            pointCloud = VtkPointCloud()
            pointCloud.addPoint(point, math.log(mag)*1) #changes this multiplier to make the points larger
            renderer.AddActor(pointCloud.vtkActor)
    renderer.AddActor(pointCloud2.vtkActor)
    renderer.SetBackground(.2, .3, .3)  #colour
    renderer.ResetCamera()

    # Render Window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
