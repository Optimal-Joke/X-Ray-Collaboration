"""This module contains various data analysis methods for the Swift, Chandra, XMM, and ROSAT telescopes.
It requires the previous installation of the CIAO data analysis package, information about which can be found
here: https://cxc.harvard.edu/ciao/.
"""

import os
from collections import Counter
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from halo import halo

# data_dir should be a directory containing each object as it's own subdirectory.
# Objects should have subdirectories for each telescope's data.
data_dir = "/Users/hunterholland/Documents/Research/Laidlaw/Data/Modified"
ciao = "source ~/Documents/ciao-4.12/bin/ciao.bash"

object_dict = {}
# Creates dictionary of object names and their directories' file paths.
for entry in os.scandir(data_dir):
    if entry.is_dir() and not entry.name.startswith('.'):
        object_dict[entry.name] = entry.path
objects = list(object_dict.values())
objects.sort()


class Telescope:
    """
    """

    def __init__(self, file):
        if os.path.isfile(file):  # Assign input to filepath
            self.file_path = file
            self.file_dir, self.file_name = os.path.split(file)
        else:  # Get filepath from file input
            for directory, _subdir, files in os.walk(f"{data_dir}"):
                if file in files:
                    self.file_path = os.path.join(directory, file)
                    self.file_dir, self.file_name = directory, file
        if data_dir in self.file_path:
            # Using self.file_path, get object name from object dictionary.
            for item in object_dict.keys():
                if item in self.file_path:
                    self.obj_name = item
                    break
            # Get object path from object directory.
            try:
                obj_path = object_dict[self.obj_name]
                self.telescopes = []
                for telescope in os.scandir(obj_path):
                    if telescope.is_dir() and "." not in telescope.name:
                        # List telescopes with object data.
                        self.telescopes.append(telescope.name)
            except AttributeError:
                print("Please arrange your filesystem as described above.")
        else:
            self.obj_name = "File is not of an identifiable object."
            self.telescopes = "Object has no identifiable telescopes with data of it."

    def get_objectname(self):
        """Returns the name of the object.
        """
        return self.obj_name

    def get_filepath(self):
        """Returns the filepath of the current datafile.
        """
        return self.file_path

    def get_telescopes(self):
        """Returns the telescopes that have data for the object in question.
        """
        return self.telescopes


class Spitzer(Telescope):
    """
    """

    def __init__(self, file):
        self.telescope = "Spitzer"
        super().__init__(file)
        if not self.telescope in self.telescopes:
            init_message = f"{self.obj_name} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.obj_name} {self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.file_path}"

    def e_hist(self, min_e, max_e, nbins='auto'):
        pass


class Chandra(Telescope):
    """Takes a file path or name as input. As of now, this class can initiate data, yield 
    energy histograms, and filter files based on photon intensity (energy) and region.
    """

    def __init__(self, file):
        super().__init__(file)
        self.telescope = "Chandra"
        # Get event data from file
        self.evt_data = fits.getdata(self.file_path)
        if not self.telescope in self.telescopes:
            init_message = f"{self.obj_name} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.obj_name} {self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.file_path}"

    def e_hist(self, *args, e_list=None, e_list2=None, nbins='auto', obj=False, save=False, filename=None):
        """Makes a histogram over the specified energy range or, optionally, of up to two lists 
        of energies passed as input. Optionally saves output as PNG active file's directory.

        To specify an energy range from a minimum energy A to a maximum energy B, 
        pass the range as a tuple: (A,B).
        """
        if e_list is not None:
            # Use list of energies instead of specified ranges to make histogram
            e = e_list
        else:
            energy = self.evt_data["energy"]  # Get energy data from event data
            # Create False boolean filter list based on the energy data. This master list will be referenced and altered by each range passed as input.
            if args:
                filter_list = [False for i in range(len(energy))]
                for arg in args:
                    try:  # If a range with a max and min is passed as input.
                        min_e, max_e = arg
                    except ValueError:  # If a range with only a min is passed as input.
                        min_e, max_e = arg, float("inf")
                    except TypeError:  # If no range is passed as input.
                        min_e, max_e = 0, float("inf")
                    # New filter list with provided range.
                    r = (energy >= min_e) & (energy < max_e)
                    # Iterate through each value of the master filter list. If the corresponding index in the new filter list is True, the value at that same index in master list will change to True.
                    for i in range(len(filter_list)):
                        if r[i] == True:
                            filter_list[i] = True
            else:
                filter_list = [True for i in range(len(energy))]
            e = energy[filter_list]
        plt.hist(e, bins=nbins)
        if e_list2 is not None:
            # Overlay second histogram if data is present
            plt.hist(e_list2, bins=nbins)
        plt.xlabel("Energy [eV]")
        plt.ylabel("Count")
        if not obj:
            plt.title(f"Energy Distribution")
        elif obj:
            plt.title(f"{obj} Energy Distribution")
        if save == True:
            if obj:
                if filename is None:
                    plt.savefig(f"{self.file_dir}/ehist.png",
                                dpi=250, format="png")
                else:
                    plt.savefig(f"{self.file_dir}/{filename}",
                                dpi=250, format="png")
                print(f"Histogram saved to {self.file_dir}.")
            else:
                if filename is None:
                    print(
                        f"Histogram is not of a specified object, so you must specify a filename.\nNo file saved.")
                    return
                else:
                    plt.savefig(f"{self.file_dir}/{filename}",
                                dpi=250, format="png")
                    print(f"Histogram saved to {self.file_dir}.")
            plt.close()
        else:
            plt.show()

    def filter_energy(self, e_range, *args, newfile="filtered_energy.fits"):
        """
        """
        ranges = f"{e_range[0]},{e_range[1]}"
        for e_rng in args:
            try:
                ranges += f",{e_rng[0]}:{e_rng[1]}"
            except TypeError:
                ranges += f",{e_rng}:"

        in_filepath = self.file_path
        out_dir = os.path.dirname(self.file_path)
        command = f"{ciao}; dmcopy \"{in_filepath}[EVENTS][energy={ranges}]\" {out_dir}/{newfile}"
        print(f"shell commands executed: {command}")
        os.system(command)
        print(f"Data filtered into file {newfile} in directory {out_dir}")

    def filter_coord(self, xcenter, ycenter, shape="circle", newfile="filtered_coords.fits", **kwargs):
        """Possible shapes are circle, box, annulus, and ellipse.
        kwargs are as described at this link: https://cxc.harvard.edu/ciao/ahelp/dmregions.html
        All units are in detector coordinates.
        """
        if shape == "circle":
            radius = kwargs["radius"]
            shape_func = f"circle({xcenter},{ycenter},{radius})"
        elif shape == "box":
            width = kwargs["width"]
            height = kwargs["height"]
            angle = kwargs["angle"]
            shape_func = f"box({xcenter},{ycenter},{width}, {height}, {angle})"
        elif shape == "annulus":
            iradius = kwargs["iradius"]
            oradius = kwargs["oradius"]
            shape_func = f"annulus({xcenter},{ycenter},{iradius}, {oradius})"
        elif shape == "ellipse":
            xradius = kwargs["xradius"]
            yradius = kwargs["yradius"]
            angle = kwargs["angle"]
            shape_func = f"ellipse({xcenter},{ycenter},{xradius}, {yradius}, {angle})"

        in_filepath = self.file_path
        out_dir = os.path.dirname(self.file_path)
        command = f"{ciao}; dmcopy \"{in_filepath}[EVENTS][(x,y)={shape_func}]\" {out_dir}/{newfile}"
        print(f"shell commands executed: {command}")
        os.system(command)
        print(f"Data filtered into file {newfile} in directory {out_dir}")

    def find_halos(self):
        halo(self.file_path)

    def mosaic(self, *files):
        pass

    def find_centroids(self):
        """
        """

        # make "Data Products" directory
        try:
            os.mkdir(f"{self.file_dir}/DataProducts")
        except FileExistsError:
            pass

        # set mkpsfmap parameters
        infile1 = self.file_name
        outfile1 = "\"DataProducts/psf.fits\""
        energy = "2"
        ecf = "0.5"
        clobber = "yes"

        mkpsfmap_command = f"mkpsfmap {infile1} {outfile1} {energy} ecf={ecf} clobber={clobber}"

        # set celldetect parameters
        infile2 = self.file_name
        outfile2 = "DataProducts/celldetect.fits"
        regfile = "DataProducts/centroids.reg"
        psffile = "DataProducts/psf.fits"

        celldetect_command = f"""punlearn celldetect;
            pset celldetect regfile = {regfile};
            pset celldetect psffile = \"{psffile}\";
            pset celldetect clobber = {clobber};
            celldetect infile={infile2} outfile=\"{outfile2}\""""

        # execute shell commands
        os.system(
            f"{ciao}; cd {self.file_dir}; {mkpsfmap_command}; {celldetect_command}; exit")

        # return source xs, ys, and errors from celldetect file
        with fits.open(f"{self.file_dir}/{outfile2}") as celldetect_file:
            source_x = np.array(celldetect_file[1].data["X"])
            source_y = np.array(celldetect_file[1].data["Y"])
            source_xerr = np.array(celldetect_file[1].data["X_ERR"])
            source_yerr = np.array(celldetect_file[1].data["Y_ERR"])
            return source_x, source_y, source_xerr, source_yerr

    def add_energies(self, *files):
        pass

    def radial_profile(self, region):
        pass


class XMM(Telescope):
    """Takes a file path or name as input. As of now, this class can initiate data, yield
    energy histograms, and filter files based on photon intensity (energy) and region.
    """

    def __init__(self, file):
        super().__init__(file)
        self.telescope = "XMM"
        # Get event data from file
        self.evt_data = fits.getdata(self.file_path)
        if not self.telescope in self.telescopes:
            init_message = f"{self.obj_name} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.obj_name} {self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.file_path}"

    def e_hist(self, *args, e_list=None, e_list2=None, nbins='auto', obj=False, save=False, filename=None):
        """Makes a histogram over the specified energy range or, optionally, of up to two lists
         of energies passed as input. Optionally saves output as PNG active file's directory.

        To specify an energy range from a minimum energy A to a maximum energy B,
        pass the range as a tuple: (A,B).
        """
        if e_list is not None:
            # Use list of energies instead of specified ranges to make histogram
            e = e_list
        else:
            energy = self.evt_data["PI"]  # Get energy data from event data
            # Create False boolean filter list based on the energy data. This master list will be
            # referenced and altered by each range passed as input.
            if args:
                filter_list = [False for i in range(len(energy))]
                for arg in args:
                    try:  # If a range with a max and min is passed as input.
                        min_e, max_e = arg
                    except ValueError:  # If a range with only a min is passed as input.
                        min_e, max_e = arg, float("inf")
                    except TypeError:  # If no range is passed as input.
                        min_e, max_e = 0, float("inf")
                    # New filter list with provided range.
                    r = (energy >= min_e) & (energy < max_e)
                    # Iterate through each value of the master filter list. If the corresponding index
                    # in the new filter list is True, the value at that same index in master list will change to True.
                    for i in range(len(filter_list)):
                        if r[i] == True:
                            filter_list[i] = True
            else:
                filter_list = [True for i in range(len(energy))]
            e = energy[filter_list]
        plt.hist(e, bins=nbins)
        if e_list2 is not None:
            # Overlay second histogram if data is present
            plt.hist(e_list2, bins=nbins)
        plt.xlabel("Energy [eV]")
        plt.ylabel("Count")
        if not obj:
            plt.title(
                f"Energy Distribution")
        elif obj:
            plt.title(
                f"{obj} Energy Distribution")
        if save == True:
            if obj:
                if filename is None:
                    plt.savefig(
                        f"{self.file_dir}/ehist.png", dpi=250, format="png")
                else:
                    plt.savefig(
                        f"{self.file_dir}/{filename}", dpi=250, format="png")
                print(f"Histogram saved to {self.file_dir}.")
            else:
                if filename is None:
                    print(
                        f"Histogram is not of a specified object, so you must specify a filename.\nNo file saved.")
                    return
                else:
                    plt.savefig(
                        f"{self.file_dir}/{filename}", dpi=250, format="png")
                    print(f"Histogram saved to {self.file_dir}.")
            plt.close()
        else:
            plt.show()


class Rosat(Telescope):
    """FITS.OPEN CURRENTLY BROKEN FOR ROSAT FILES. STILL TRYING TO TROUBLESHOOT.

    Takes a file path or name as input. As of now, this class can initiate data, yield
    energy histograms, and filter files based on photon intensity (energy) and region.
    """

    def __init__(self, file):
        self.telescope = "ROSAT"
        super().__init__(file)
        if not self.telescope in self.telescopes:
            init_message = f"{self.obj_name} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.obj_name} {self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.file_path}"

    def e_hist(self, *args, e_list=None, e_list2=None, nbins='auto', obj=False, save=False, filename=None):
        """Makes a histogram over the specified energy range or, optionally, of up to two 
        lists of energies passed as input. Optionally saves output as PNG active file's directory.

        To specify an energy range from a minimum energy A to a maximum energy B, 
        pass the range as a tuple: (A,B).
        """
        if e_list is not None:
            # Use list of energies instead of specified ranges to make histogram
            e = e_list
        else:
            with fits.open(self.file_path) as hdul:
                evt_table = hdul[2]  # Get event table from file
                # Get energy data from event data
                energy = evt_table.data["PI"]
                # Create False boolean filter list based on the energy data. This master list
                # will be referenced and altered by each range passed as input.
                if args:
                    filter_list = [False for i in range(len(energy))]
                    for arg in args:
                        try:  # If a range with a max and min is passed as input.
                            min_e, max_e = arg
                        except ValueError:  # If a range with only a min is passed as input.
                            min_e, max_e = arg, float("inf")
                        except TypeError:  # If no range is passed as input.
                            min_e, max_e = 0, float("inf")
                        # New filter list with provided range.
                        r = (energy >= min_e) & (energy < max_e)
                        # Iterate through each value of the master filter list. If the corresponding index
                        # in the new filter list is True, the value at that same index in master list will change to True.
                        for i in range(len(filter_list)):
                            if r[i] == True:
                                filter_list[i] = True
                else:
                    filter_list = [True for i in range(len(energy))]
                e = energy[filter_list]
        plt.hist(e, bins=nbins)
        if e_list2 is not None:
            # Overlay second histogram if data is present
            plt.hist(e_list2, bins=nbins)
        plt.xlabel("Energy [eV]")
        plt.ylabel("Count")
        if not obj:
            plt.title(
                f"Energy Distribution")
        elif obj:
            plt.title(
                f"{obj} Energy Distribution")
        if save == True:
            if obj:
                if filename is None:
                    plt.savefig(
                        f"{self.file_dir}/ehist.png", dpi=250, format="png")
                else:
                    plt.savefig(
                        f"{self.file_dir}/{filename}", dpi=250, format="png")
                print(f"Histogram saved to {self.file_dir}.")
            else:
                if filename is None:
                    print(
                        f"Histogram is not of a specified object, so you must specify a filename.\nNo file saved.")
                    return
                else:
                    plt.savefig(
                        f"{self.file_dir}/{filename}", dpi=250, format="png")
                    print(f"Histogram saved to {self.file_dir}.")
            plt.close()
        else:
            plt.show()

    def filter_coord(self, xcenter, ycenter, shape="circle", newfile="test.fits", **kwargs):
        """Possible shapes are circle, box, annulus, and ellipse.
        kwargs are as described at this link: https://cxc.harvard.edu/ciao/ahelp/dmregions.html
        All units are in detector coordinates.
        """
        if shape == "circle":
            radius = kwargs["radius"]
            shape_func = f"circle({xcenter},{ycenter},{radius})"
        elif shape == "box":
            width = kwargs["width"]
            height = kwargs["height"]
            angle = kwargs["angle"]
            shape_func = f"box({xcenter},{ycenter},{width}, {height}, {angle})"
        elif shape == "annulus":
            iradius = kwargs["iradius"]
            oradius = kwargs["oradius"]
            shape_func = f"annulus({xcenter},{ycenter},{iradius}, {oradius})"
        elif shape == "ellipse":
            xradius = kwargs["xradius"]
            yradius = kwargs["yradius"]
            angle = kwargs["angle"]
            shape_func = f"ellipse({xcenter},{ycenter},{xradius}, {yradius}, {angle})"

        in_filepath = self.file_path
        out_dir = os.path.dirname(self.file_path)
        os.system(
            f"{ciao}; dmcopy \"{in_filepath}[STDEVT][(x,y)={shape_func}]\" {out_dir}/{newfile}")


class Swift(Telescope):
    """Takes a file path or name as input. As of now, this class can initiate data, yield
    energy histograms, and filter files based on photon intensity (energy) and region.
    """

    def __init__(self, file):
        super().__init__(file)
        self.telescope = "Swift"
        # Get event data from file
        self.evt_data = fits.getdata(self.file_path)
        if not self.telescope in self.telescopes:
            init_message = f"{self.obj_name} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.obj_name} {self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.file_path}"

    def e_hist(self, *args, e_list=None, e_list2=None, nbins='auto', obj=False, save=False, filename=None):
        """Makes a histogram over the specified energy range or, optionally, of up to two 
        lists of energies passed as input. Optionally saves output as PNG active file's directory.

        To specify an energy range from a minimum energy A to a maximum energy B, pass the range as a tuple: (A,B).
        """
        if e_list is not None:
            # Use list of energies instead of specified ranges to make histogram
            e = e_list
        else:
            energy = self.evt_data["PI"]  # Get energy data from event data
            # Create False boolean filter list based on the energy data. This master list
            # will be referenced and altered by each range passed as input.
            if args:
                filter_list = [False for i in range(len(energy))]
                for arg in args:
                    try:  # If a range with a max and min is passed as input.
                        min_e, max_e = arg
                    except ValueError:  # If a range with only a min is passed as input.
                        min_e, max_e = arg, float("inf")
                    except TypeError:  # If no range is passed as input.
                        min_e, max_e = 0, float("inf")
                    # New filter list with provided range.
                    r = (energy >= min_e) & (energy < max_e)
                    # Iterate through each value of the master filter list. If the corresponding
                    # index in the new filter list is True, the value at that same index in master list will change to True.
                    for i in range(len(filter_list)):
                        if r[i] == True:
                            filter_list[i] = True
            else:
                filter_list = [True for i in range(len(energy))]
            e = energy[filter_list]
        plt.hist(e, bins=nbins)
        if e_list2 is not None:
            # Overlay second histogram if data is present
            plt.hist(e_list2, bins=nbins)
        plt.xlabel("Energy [eV]")
        plt.ylabel("Count")
        if not obj:
            plt.title(
                f"Energy Distribution")
        elif obj:
            plt.title(
                f"{obj} Energy Distribution")
        if save == True:
            if obj:
                if filename is None:
                    plt.savefig(
                        f"{self.file_dir}/ehist.png", dpi=250, format="png")
                else:
                    plt.savefig(
                        f"{self.file_dir}/{filename}", dpi=250, format="png")
                print(f"Histogram saved to {self.file_dir}.")
            else:
                if filename is None:
                    print(
                        f"Histogram is not of a specified object, so you must specify a filename.\nNo file saved.")
                    return
                else:
                    plt.savefig(
                        f"{self.file_dir}/{filename}", dpi=250, format="png")
                    print(f"Histogram saved to {self.file_dir}.")
            plt.close()
        else:
            plt.show()

    def filter_energy(self, e_range: tuple, *args: tuple, newfile="filtered_energy.fits"):
        """Energy ranges should be input as tuples. For example, a range from 100-400ev would be
        represented as (100,400), and ranges from 100-400eV, 500-900ev, and 1000-2000eV would be
        represented as (100,400), (500,900), (1000,2000).
        """
        ranges = f"{e_range[0]}:{e_range[1]}"
        for e_rng in args:
            try:
                ranges += f",{e_rng[0]}:{e_rng[1]}"
            except TypeError:
                ranges += f",{e_rng}:"

        in_filepath = self.file_path
        out_dir = os.path.dirname(self.file_path)
        command = f"{ciao}; dmcopy \"{in_filepath}[EVENTS][PI={ranges}]\" {out_dir}/{newfile}"
        print(f"shell commands executed: {command}")
        os.system(command)
        print(f"Data filtered into file {newfile} in directory {out_dir}")

    def filter_coord(self, xcenter, ycenter, shape="circle", newfile="test.fits", **kwargs):
        """Possible shapes are circle, box, annulus, and ellipse.
        kwargs are as described at this link: https://cxc.harvard.edu/ciao/ahelp/dmregions.html
        All units are in detector coordinates.
        """
        if shape == "circle":
            radius = kwargs["radius"]
            shape_func = f"circle({xcenter},{ycenter},{radius})"
        elif shape == "box":
            width = kwargs["width"]
            height = kwargs["height"]
            angle = kwargs["angle"]
            shape_func = f"box({xcenter},{ycenter},{width}, {height}, {angle})"
        elif shape == "annulus":
            iradius = kwargs["iradius"]
            oradius = kwargs["oradius"]
            shape_func = f"annulus({xcenter},{ycenter},{iradius}, {oradius})"
        elif shape == "ellipse":
            xradius = kwargs["xradius"]
            yradius = kwargs["yradius"]
            angle = kwargs["angle"]
            shape_func = f"ellipse({xcenter},{ycenter},{xradius}, {yradius}, {angle})"

        in_filepath = self.file_path
        out_dir = os.path.dirname(self.file_path)
        os.system(
            f"{ciao}; dmcopy \"{in_filepath}[EVENTS][(x,y)={shape_func}]\" {out_dir}/{newfile}")

    def find_centroids(self):
        """This (crude) algorithm iterates through the pixels of an image and returns the
        RA and Dec coordinates of the maxima of the image as two 1D numpy arrays.
        """
        with fits.open(f"{self.file_path}") as file:
            # load in RA/Dec coordinates
            xs = np.array(file[1].data["X"])
            ys = np.array(file[1].data["Y"])

            # determine x and y bins
            xbinedges = np.histogram_bin_edges(xs, bins='fd')
            ybinedges = np.histogram_bin_edges(ys, bins='fd')

            # iterate through the bins on each axis
            image_maxes = []
            for j in range(len(ybinedges)-1):
                for i in range(len(xbinedges)-1):
                    cell_xmin, cell_xmax = xbinedges[i], xbinedges[i+1]
                    cell_ymin, cell_ymax = ybinedges[j], ybinedges[j+1]

                    # isolate the events within each bin
                    cell_constraints = np.where((xs < cell_xmax) & (
                        xs >= cell_xmin) & (ys < cell_ymax) & (ys >= cell_ymin))

                    cell_xs = xs[cell_constraints]
                    cell_ys = ys[cell_constraints]
                    coord_pairs = list(zip(cell_xs, cell_ys))

                    # find the maximum point in the bin
                    coord_counts = Counter(coord_pairs).most_common(1)
                    for result in coord_counts:
                        image_maxes.append(result)

            # sort the identified maxima by count and return their x(RA) and y(Dec) coordinates
            image_maxes = sorted(image_maxes, key=lambda x: x[1], reverse=True)
            print(len(image_maxes))
            centroid_xs = np.array([item[0][0] for item in image_maxes])
            centroid_ys = np.array([item[0][1] for item in image_maxes])
            return centroid_xs, centroid_ys


test = Swift("/Users/hunterholland/Documents/Research/Laidlaw/Data/Modified/L1517/Swift/xrt/event/sw00034249004xpcw3po_cl.evt.gz")
test.find_centroids()
