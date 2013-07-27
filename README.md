QuakeCloud
==========

This is a 3d quake event visualiser. There are a number of supported variables in the main python code:-


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

        # Quake event data location
        self.fname = "quake_data.csv"

        # Mapping vars
        self.overlay_image_fname = "overlay.jpg"
        self.map_height = 0  # Sets the top height of the model in km
        self.long_min = 174  # defines the minimum longitude bounds of the sampled area
        self.long_max = 175  # defines the maximum longitude bounds of the sampled area
        self.lat_min = -41  # defines the minimum latitude bounds of the sampled area
        self.lat_max = -42  # defines the maximum latitude bounds of the sampled area
        self.multiplier = 1  # used to allow the whole model to be scaled up if needed (has bugs)
        

Data grabbed from geonet: http://wfs-beta.geonet.org.nz/geoserver/geonet/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=geonet:quake&outputFormat=csv&cql_filter=BBOX%28origin_geom,174,-41,175,-42%29+AND+origintime%3E=%272013-07-18%27+AND+magnitude%3E3
 
 
put the py file, the data csv file, and the image file in the same folder. 

Run the py file 

Needs vtk


