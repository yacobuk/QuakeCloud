__author__ = 'yacob'
import vtk
import csv
import numpy
import math
import collections
from datetime import datetime, timedelta
from itertools import groupby


class Data_Extractor(object):
    """
    collects data from a CSV, extracts the required cells, and creates an "event" list, converts lat / long
    to relative distance, applies scaling multiplyer, finds max depth
    """
    def __init__(self, session):
        self.fname = session.fname
        self.multiplyer = session.multiplier  # scalar used to increase model size
        self.event_list = []
        self.event_max_depth = -1000000  # starting height 100km high that expects to be overwritten
        self.earth_radius = 6373  # ref radius @ sea level (km)
        self.long_min = session.long_min  # defines the bounds of the sampled area
        self.long_max = session.long_max  # defines the bounds of the sampled area
        self.lat_min = session.lat_min  # defines the bounds of the sampled area
        self.lat_max = session.lat_max  # defines the bounds of the sampled area
        self.rounding_base_for_depth = 25
        self.get_data()

    def get_data(self):
        """
        main class routine, collect data, prepare, make event list, sorts list into time order
        """
        reader = csv.reader(open(self.fname, "rb"))
        for row in reader:
            if "FID" in row[0]:
                pass
            else:
                event_longitude, event_latitude = map(float, row[13].translate(None, "POINT()").strip().split())
                event_depth = float(row[5])
                event_magnitude = float(row[6])
                date, junk = row[12].split(".")
                time_of_event = datetime.strptime(date, '%Y-%m-%dT%H:%M:%S')
                event_x_point = self.make_x_point(event_latitude, event_depth)
                event_y_point = self.make_y_point(event_longitude, event_depth)
                points_mag_and_time = [event_x_point * self.multiplyer,  event_y_point * self.multiplyer,
                                       event_depth * self.multiplyer, event_magnitude, time_of_event]
                points_mag_and_time = numpy.asarray(points_mag_and_time)
                self.event_list.append(points_mag_and_time)
                if event_depth > self.event_max_depth:
                    self.event_max_depth = event_depth * self.multiplyer
                    self.event_max_depth = int(self.rounding_base_for_depth *
                                               round(float(self.event_max_depth)/self.rounding_base_for_depth)) + \
                                           self.rounding_base_for_depth

        self.event_list.sort(key=lambda point_mag_and_time: point_mag_and_time[4])

    def make_x_point(self, event_latitude, event_depth):
        """
        returns the (km) distance normalised latitude value, relative to min_lat,min.long bounds coordinates)
        """
        distance = self.distance_on_unit_sphere(self.lat_min, self.long_min, event_latitude, self.long_min, event_depth)
        return distance

    def make_y_point(self, event_longitude, event_depth):
        """
        returns the (km) distance normalised longitude value, relative to min_lat,min.long bounds coordinates)
        """
        distance = self.distance_on_unit_sphere(self.lat_min, self.long_min, self.lat_min, event_longitude, event_depth)
        return distance

    def convert_depth_to_relative_radius(self, depth):
        """
        changes the size of the perfect circle (earth) relative to the depth of the event
        """
        relative_radius = self.earth_radius - depth
        return relative_radius

    def distance_on_unit_sphere(self, lat1, long1, lat2, long2, depth):
        """
        when given reference long and lat, and event long and lat, returns the distance in Km
        """
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
        arc = math.acos(cos)
        # Remember to multiply arc by the radius of the earth
        # in your favorite set of units to get length.

        relative_radius = self.convert_depth_to_relative_radius(depth)
        distance = arc * relative_radius
        return distance


class Model_Bounds_Points_Maker(object):
    """
    creates the bounding for the model, used to populate the points limits, and overlay the map image
    """
    def __init__(self, data, session):
        self.data = data
        self.base = 25
        self.map_height = session.map_height
        self.make_bounds()
        self.make_bounds_scale()
        self.lat_max_distance = None
        self.long_max_distance = None
        self.make_bounds_distance_array()
        self.set_lat_and_long_max_distance()

    def make_bounds(self):
        """
        sets up the bounds based on the source long / lat nim / max
        """
        self.bounds = [[self.data.lat_min, self.data.long_min, self.map_height],
                       [self.data.lat_min, self.data.long_max, self.map_height],
                       [self.data.lat_max, self.data.long_min, self.map_height],
                       [self.data.lat_max, self.data.long_max, self.map_height],
                       [self.data.lat_min, self.data.long_min, self.data.event_max_depth],
                       [self.data.lat_min, self.data.long_max, self.data.event_max_depth],
                       [self.data.lat_max, self.data.long_min, self.data.event_max_depth],
                       [self.data.lat_max, self.data.long_max, self.data.event_max_depth]]

    def make_bounds_scale(self):
        self.bounds_scale = []
        for layer in xrange(0, self.data.event_max_depth, self.base):
            bound = [0, 0, layer]
            self.bounds_scale.append(bound)


    def convert_bound_to_distance(self, bound):
        """
        converts bound long / lat makers to distance vectors
        """
        x, y, depth = bound
        x_point = self.data.make_x_point(x, depth)
        y_point = self.data.make_y_point(y, depth)
        bound = (x_point * self.data.multiplyer, y_point * self.data.multiplyer, depth * self.data .multiplyer)
        return bound

    def make_bounds_distance_array(self):
        """
        refreshes the bounds array as distance vectors
        """
        bounds = []
        for bound in self.bounds:
            bound = self.convert_bound_to_distance(bound)
            bounds.append(bound)
        self.bounds = bounds

    def set_lat_and_long_max_distance(self):
        bound = (self.data.lat_min, self.data.long_max, 0)
        self.lat_max_distance = self.convert_bound_to_distance(bound)
        bound = (self.data.lat_max, self.data.long_min, 0)
        self.long_max_distance = self.convert_bound_to_distance(bound)


class Timed_Event_List_Maker(object):
    """
    converts the sequential event list into time slice grouped list
    """
    def __init__(self, event_list, session):
        self.event_list = event_list
        self.time_block_size = session.minutes_per_event_block  # size of time block  (in mins) to slice data into
        self.timed_event_dict = collections.OrderedDict()
        self.time_slicer()

    def get_key(self, d):
        """
        group by self.time_block_size
        """
        k = d + timedelta(minutes=-(d.minute % self.time_block_size))
        k = datetime(k.year, k.month, k.day, k.hour, k.minute, 0)
        return k

    def add_missing_empty_frames(self, g):
        """
        adds the missing frames into the event stream (where a time block occurs with no events)
        """
        last_key = None
        for key, items in g:
            if last_key:
                while (key-last_key).seconds > self.time_block_size*60:
                    empty_key = last_key + timedelta(minutes=self.time_block_size)
                    yield (empty_key, [])
                    last_key = empty_key
            yield (key, items)
            last_key = key

    def time_slicer(self):
        """
        slices the event list by time slots, returns a dictionary of {timeslot:[events...]}
        """

        times_only_event_list = []
        for event in self.event_list:
            times_only_event_list.append(event[4])

        g = groupby(times_only_event_list, key=self.get_key)

        counter = 0
        for key, items in self.add_missing_empty_frames(g):
            block = []
            for i, item in enumerate(items):
                block.append(self.event_list[counter])
                if item:
                    counter += 1
            self.timed_event_dict[key] = block

    def print_timed_event_list(self):
        """
        reminder of how to get to each of the timed event from the dictionary....
        """
        for key, value in self.timed_event_dict.iteritems():
            print key
            for event in value:
                print event


class VtkPointCloud:
    def __init__(self, event_max_depth, zMin=-0.0, zMax=156, maxNumPoints=1e6): #sets colour limits
        zMax = event_max_depth
        self.maxNumPoints = maxNumPoints
        self.vtkPolyData = vtk.vtkPolyData()
        self.clearPoints()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(self.vtkPolyData)
        mapper.SetColorModeToDefault()
        mapper.SetScalarRange(zMin, event_max_depth)
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


class vtkTimerCallback():
    def __init__(self, mydict, max_event_depth, ren, renderWindow, session):
        self.timer_count = 0
        self.session = session
        self.event_window_size = session.display_event_window_in_event_blocks
        self.mydict = mydict
        self.rendered_actors = collections.deque([None] * self.event_window_size, maxlen=self.event_window_size)
        self.max_event_depth = max_event_depth
        self.ren = ren
        self.renderWindow = renderWindow
        self.mag_display_limit = session.display_mag_threshold

    def execute(self, obj, event):
        pointCloud_event = VtkPointCloud(self.max_event_depth)
        current_event_mag_threshold_reached = 0
        self.ren.RemoveActor(self.rendered_actors[-1])
        try:
            self.rendered_actors[-2].GetProperty().SetOpacity(0.1)
        except:
            pass
        try:
            self.rendered_actors[-3].GetProperty().SetOpacity(0.4)
        except:
            pass
        try:
            key, event_block = self.mydict.popitem(0)
        except:
            quit()
        time_stamp = key
        print time_stamp
        if len(event_block) > 0:
            for event in event_block:
                mag = event[3]
                point = [event[0], event[1], event[2]]
                current_event_mag_threshold_reached = 0
                if mag > self.mag_display_limit:
                    current_event_mag_threshold_reached += 1
                    if self.session.display_quake_points_log:
                        pointCloud_event.addPoint(point, math.log(mag) * self.session.display_point_scalar_for_quake_points)
                    if self.session.display_quake_points_only_no_mag:
                        pointCloud_event.addPoint(point, 1)
                    if self.session.display_quake_points_linear:
                        pointCloud_event.addPoint(point, mag * self.session.display_point_scalar_for_quake_points)
                    else:
                        print "Something went wrong. Please check your quake point type settings"
                        quit()

                    #self.ren.ResetCamera()
                    self.renderWindow.Render()
        self.ren.AddActor(pointCloud_event.vtkActor)
        self.rendered_actors.appendleft(pointCloud_event.vtkActor)
        #self.ren.ResetCamera()
        self.renderWindow.Render()
        iren = obj
        iren.GetRenderWindow().Render()
        self.timer_count += 1


class Set_Render_Variables(object):
    """ 
    used to set the various controsl that allow the user to specify the supported "look", timing and map options
    """
    def __init__(self):
        # Display Vars
        self.display_mag_threshold = 1  # sets the minimum magnitude of quake to display
        self.display_event_window_in_event_blocks = 6  # sets the number of frames the event set will persist on map)
        self.minutes_per_event_block = 20  # in mins
        self.replay_frame_speed = 50  # defines the refresh rate of the replay
        self.display_point_scalar_for_quake_points = 10
        self.display_quake_points_log = False  # if True, quake points are plotted in size with magnitude in log scale
        self.display_quake_points_linear = True  # if True, quake points are plotted in size with magnitude in linear scale
        self.display_quake_points_only_no_mag = False  # if True, quake points are just plotted in size, with no magnitude data
        self.plot_boundary_markers = False  # If true, displays the corner boundary markers
        self.plot_depth_scale = True  # If true. displays the depth scale markers

        # quake_data_location
        self.fname = "quake(13).csv"

        #mapping vars
        self.overlay_image_fname = "overlay.jpg"
        self.map_height = 0  # Sets the top height of the model in km
        self.long_min = 174  # defines the minimum longitude bounds of the sampled area
        self.long_max = 175  # defines the maximum longitude bounds of the sampled area
        self.lat_min = -41  # defines the minimum latitude bounds of the sampled area
        self.lat_max = -42  # defines the maximum latitude bounds of the sampled area
        self.multiplier = 1  # used to allow the whole model to be scaled up if needed (has bugs)


def main():
    session = Set_Render_Variables()

    data = Data_Extractor(session)
    bounds = Model_Bounds_Points_Maker(data, session)
    events = Timed_Event_List_Maker(data.event_list, session)

    # sets up the VTK session
    ren = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    #renderWindow.AddRenderer(ren)
    i_ren = vtk.vtkRenderWindowInteractor()
    i_ren.SetRenderWindow(renderWindow)

    # text / legend object
    legend = vtk.vtkCornerAnnotation()
    legend.SetText(0, "          Quake Cloud 3d Visualiser\n          Depth Markers are 25Km Apart\n"
                      "          All other distances are relative \n\n\n          Copyright 2013 Jay Gattuso\n\n")
    legend.SetMaximumFontSize(15)


    # handles the overlay map image
    input_image = session.overlay_image_fname
    reader = vtk.vtkJPEGReader()
    reader.SetFileName(input_image)
    texture = vtk.vtkTexture()
    texture.SetInput(reader.GetOutput())
    texture.SetInputConnection(reader.GetOutputPort())
    texture.InterpolateOn()

    ren = vtk.vtkRenderer()
    renderWindow.AddRenderer(ren)
    plane = vtk.vtkPlaneSource()
    plane.SetOrigin((0, 0, 0))
    plane.SetPoint1(bounds.lat_max_distance)
    plane.SetPoint2(bounds.long_max_distance)

    # Create texture object
    planeMapper = vtk.vtkPolyDataMapper()
    planeMapper.SetInputConnection(plane.GetOutputPort())
    planeActor = vtk.vtkActor()
    planeActor.SetMapper(planeMapper)
    planeActor.SetTexture(texture)

    pointCloud_bounds = VtkPointCloud(data.event_max_depth)

    if session.plot_boundary_markers:
        for bound in bounds.bounds:
            pointCloud_bounds.addPoint(bound, 10)  # the number is what makes the bound markers larger

    if session.plot_depth_scale:
        for bound in bounds.bounds_scale:
            pointCloud_bounds.addPoint(bound, 10)  # the number is what makes the bound markers larger

    # add basic objects (bounds, overlay, background)
    ren.AddActor(pointCloud_bounds.vtkActor)
    ren.AddViewProp(legend)
    ren.AddActor(planeActor)
    ren.SetBackground(.2, .2, .2)
    ren.ResetCamera()

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Begin Interaction
    renderWindow.Render()

    # Initialize must be called prior to creating timer events.
    renderWindowInteractor.Initialize()

    # Sign up to receive TimerEvent
    cb = vtkTimerCallback(events.timed_event_dict, data.event_max_depth, ren, renderWindow, session)
    cb.actor = vtk.vtkActor()
    renderWindowInteractor.AddObserver('TimerEvent', cb.execute)
    timerId = renderWindowInteractor.CreateRepeatingTimer(session.replay_frame_speed)

    renderWindowInteractor.Start()

if __name__ == '__main__':
    main()
