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
        self.display_quake_points_log = session.display_quake_points_log
        self.display_quake_points_only_no_mag = session.display_quake_points_only_no_mag
        self.display_quake_points_linear = session.display_quake_points_linear
        self.display_point_scalar_for_quake_points = session.display_point_scalar_for_quake_points
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
                if self.display_quake_points_log:
                    event_magnitude = (math.log(event_magnitude) * self.display_point_scalar_for_quake_points)
                elif self.display_quake_points_only_no_mag:
                    event_magnitude = 10
                elif self.display_quake_points_linear:
                    event_magnitude = event_magnitude * self.display_point_scalar_for_quake_points
                else:
                    print "Something went wrong. Please check your \"quake point type\" settings"
                    quit()
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
        when given reference long and lat, and event long and lat, returns the distance in Km in x, y, co-ords
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
        self.bounds_x_max = None
        self.bounds_y_max = None
        self.bounds_z_max = None
        self.make_model_maxes()
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

    def make_model_maxes(self):
        """
        makes the model x, y, z, max values
        """
        self.bounds_x_max, self.bounds_y_max, self.bounds_z_max = \
            self.convert_bound_to_distance([self.data.lat_max, self.data.long_max, self.data.event_max_depth])

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

class VtkPointCloud:
    def __init__(self, event_max_depth, zMin=-0.0, zMax=156, maxNumPoints=1e6):
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
    def __init__(self, mydict, max_event_depth, rens, renderWindow, timer_text, session):
        self.timer_count = 0
        self.session = session
        self.event_window_size = session.display_event_window_in_event_blocks
        self.mydict = mydict
        self.rendered_actors = collections.deque([None] * self.event_window_size, maxlen=self.event_window_size)
        self.max_event_depth = max_event_depth
        self.rens = rens
        self.renderWindow = renderWindow
        self.mag_display_limit = session.display_mag_threshold
        self.timer_text = timer_text
        self.point_size = None
        self.set_point_size_text()

    def set_point_size_text(self):
        if self.session.display_quake_points_log:
            self.point_size = "Log"
        if self.session.display_quake_points_only_no_mag:
            self.point_size = "Location Only"
        if self.session.display_quake_points_linear:
            self.point_size = "Linear"

    def execute(self, obj, event):
        event_block = None

        key = None

        pointCloud_event = VtkPointCloud(self.max_event_depth)
        current_event_mag_threshold_reached = 0
        for i in range(4):
            self.rens[i].RemoveActor(self.rendered_actors[-1])
            try:
                self.rendered_actors[-2].GetProperty().SetOpacity(0.1)
            except:
                pass
            try:
                self.rendered_actors[-3].GetProperty().SetOpacity(0.3)
            except:
                pass
            try:
                self.rendered_actors[-4].GetProperty().SetOpacity(0.6)
            except:
                pass
            try:
                self.rendered_actors[-5].GetProperty().SetOpacity(0.75)
            except:
                pass
            try:
                key, event_block = self.mydict.popitem(0)
            except:
                quit()
            time_stamp = key
            if len(event_block) > 0:
                for event in event_block:
                    mag = event[3]
                    point = [event[0], event[1], event[2]]
                    current_event_mag_threshold_reached = 0
                    if mag > self.mag_display_limit:
                        current_event_mag_threshold_reached += 1
                        pointCloud_event.addPoint(point, mag)
            self.timer_text.SetText(1, "Date/Time: {}                    \n"
                                       "Min Quake Size Displayed: {}                  \n"
                                       "Quake Point Size Represented as: {}         \n\n\n".
                                    format(time_stamp, self.session.display_mag_threshold, self.point_size))

        #for i in range(4):
            self.rens[i].AddActor(pointCloud_event.vtkActor)
            self.rendered_actors.appendleft(pointCloud_event.vtkActor)
            self.renderWindow.Render()
            iren = obj
            iren.GetRenderWindow().Render()
        self.timer_count += 1


class VTK_Render_Parts(object):
    def __init__(self, session, event_max_depth, bounds, timed_events_dict):
        self.session = session
        self.event_max_depth = event_max_depth
        self.bounds = bounds
        self.timed_events_dict = timed_events_dict
        self.pointCloud_bounds = VtkPointCloud(self.event_max_depth)
        self.pointCloud_scale = VtkPointCloud(self.event_max_depth)
        self.actors = []
        self.renderers = []
        self.gridDimensions = 2  # makes 2 x 2 grid of viewports
        self.rendererSize = 300
        self.viewports = []
        self.camera_positions = []

        self.renderWindow = vtk.vtkRenderWindow()
        self.legend = None
        self.timer_text = None
        self.overlayActor = None
        self.pointCloud_bounds = None
        self.pointCloud_bounds = None
        self.set_camera_positions()
        self.build_model()

    def set_camera_positions(self):
        scalar = 4
        self.camera_positions = [[-self.bounds.bounds_x_max * scalar, 0, self.bounds.bounds_z_max],
                                 [0, -self.bounds.bounds_y_max * scalar/0.7, self.bounds.bounds_z_max],
                                 [self.bounds.bounds_x_max * scalar/1.7, self.bounds.bounds_y_max * scalar/1.7, self.bounds.bounds_z_max * scalar/1.7],
                                 [self.bounds.bounds_x_max * scalar/2, 0, self.bounds.bounds_z_max * scalar/2]]

        self.camera_focal_point = [[self.bounds.bounds_x_max/2, self.bounds.bounds_y_max/2, self.bounds.bounds_z_max/2],
                                 [self.bounds.bounds_x_max/2, self.bounds.bounds_y_max/2, self.bounds.bounds_z_max/2],
                                 [self.bounds.bounds_x_max/2, self.bounds.bounds_y_max/2, self.bounds.bounds_z_max/2],
                                 [self.bounds.bounds_x_max/2, self.bounds.bounds_y_max/2, 0]]

        self.camera_rolls = [-90, 180, -90, -90]
        self.camera_azimuth = [0, 0, 45, 0]
        self.camera_elevation = [15, 15, 0, 35]

    def build_model(self):
        self.make_renderers()

        for renderer in self.renderers:
            self.renderWindow.AddRenderer(renderer)
        self.make_legend_text_object()
        self.make_timed_text_object()
        self.overlayActor = self.make_texture_overly_object()
        self.pointCloud_bounds = self.make_bounds_point_cloud()
        self.pointCloud_bounds = self.make_depth_scale_pointcloud()
        self.make_render_window_populate_views()



        self.start_interactor()


    def make_renderers(self):
        for grid_id in range(4):
            if grid_id < self.gridDimensions * self.gridDimensions:
                self.renderers.append(vtk.vtkRenderer())

    def make_render_window_populate_views(self):
        self.renderWindow.SetSize(self.rendererSize * self.gridDimensions, self.rendererSize * self.gridDimensions)
        self.cameras = [vtk.vtkCamera()] * 4
        for row in range(self.gridDimensions):
            for col in range(self.gridDimensions):
                idx = row * self.gridDimensions + col
                self.viewports[:] = []
                self.renderers[idx].ResetCamera()
                self.viewports.append(float(col) * self.rendererSize / (self.gridDimensions * self.rendererSize))
                self.viewports.append(float(self.gridDimensions - (row+1)) * self.rendererSize /
                                      (self.gridDimensions * self.rendererSize))
                self.viewports.append(float(col+1) * self.rendererSize / (self.gridDimensions * self.rendererSize))
                self.viewports.append(float(self.gridDimensions - row) *
                                      self.rendererSize / (self.gridDimensions * self.rendererSize))
                if idx > (4 - 1):
                    continue
                self.renderers[idx].SetViewport(self.viewports)
                ### make camera views / view ports
                camera = vtk.vtkCamera()
                camera.SetPosition(self.camera_positions[idx])
                camera.SetRoll(self.camera_rolls[idx])
                camera.SetFocalPoint(self.camera_focal_point[idx])
                camera.Azimuth(self.camera_azimuth[idx])
                camera.Elevation(self.camera_elevation[idx])
                self.renderers[idx].SetActiveCamera(camera)
                ### add actors to model
                for actor in self.actors:
                    self.renderers[idx].AddActor(actor)
                self.renderers[idx].AddViewProp(self.legend)
                self.renderers[idx].AddViewProp(self.timer_text)
                self.renderers[idx].SetBackground(0.3, 0.3, 0.3)

    def make_legend_text_object(self):
        self.legend = vtk.vtkCornerAnnotation()
        self.legend.SetText(0,"Quake Cloud 3d Visualiser\n\
Depth Markers are 25Km Apart\n\
Maximum depth for model: {}\n\
All other distances are relative \n\n\n\
Copyright 2013 Jay Gattuso\n\n".format(self.event_max_depth))
        self.legend.SetMaximumFontSize(15)

    def make_timed_text_object(self):
        self.timer_text = vtk.vtkCornerAnnotation()
        self.timer_text.SetMaximumFontSize(20)

    def make_texture_overly_object(self):
        input_image = self.session.overlay_image_fname
        reader = vtk.vtkJPEGReader()
        reader.SetFileName(input_image)
        texture = vtk.vtkTexture()
        texture.SetInput(reader.GetOutput())
        texture.SetInputConnection(reader.GetOutputPort())
        texture.InterpolateOn()

        plane = vtk.vtkPlaneSource()
        plane.SetOrigin((0, 0, 0))
        plane.SetPoint1(self.bounds.lat_max_distance)
        plane.SetPoint2(self.bounds.long_max_distance)

        planeMapper = vtk.vtkPolyDataMapper()
        planeMapper.SetInputConnection(plane.GetOutputPort())
        self.overlayActor = vtk.vtkActor()
        self.overlayActor.SetMapper(planeMapper)
        self.overlayActor.SetTexture(texture)
        self.actors.append(self.overlayActor)

    def make_bounds_point_cloud(self):
        if self.session.plot_boundary_markers:
            for bound in self.bounds.bounds:
                self.pointCloud_bounds.addPoint(bound, 10)  # the number is what makes the bound markers larger
            self.actors.append(self.pointCloud_bounds.vtkActor)

    def make_depth_scale_pointcloud(self):
        if self.session.plot_depth_scale:
            for bound in self.bounds.bounds_scale:
                self.pointCloud_scale.addPoint(bound, 10)  # the number is what makes the bound markers larger
            self.actors.append(self.pointCloud_scale.vtkActor)

    def start_interactor(self):
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.renderWindow)
        self.renderWindow.Render()

        cb = vtkTimerCallback(self.timed_events_dict, self.event_max_depth, self.renderers, self.renderWindow, self.timer_text, self.session)
        cb.actor = vtk.vtkActor()

        #self.renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self.interactor.AddObserver('TimerEvent', cb.execute)

        self.timerId = self.interactor.CreateRepeatingTimer(self.session.replay_frame_speed)
        self.interactor.Start()

class Set_Render_Variables(object):
    """
    used to set the various controsl that allow the user to specify the supported "look", timing and map options
    """
    def __init__(self):
        # Display Vars
        self.display_mag_threshold = 2.5  # sets the minimum magnitude of quake to display
        self.display_event_window_in_event_blocks = 6  # sets the number of frames the event set will persist on map)
        self.minutes_per_event_block = 20  # in mins
        self.replay_frame_speed = 250  # defines the refresh rate of the replay in ms
        self.display_point_scalar_for_quake_points = 10
        self.display_quake_points_log = False  # if True, quake points are plotted in size with magnitude in log scale
        self.display_quake_points_linear = False  # if True, quake points are plotted in size with mag in linear scale
        self.display_quake_points_only_no_mag = True  # if True, quake points are  plotted in size, with no mag data
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
    model = VTK_Render_Parts(session, data.event_max_depth, bounds, events.timed_event_dict)

if __name__ == '__main__':
    main()
