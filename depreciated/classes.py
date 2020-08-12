import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# data_dir should be a directory containing each object as it's own subdirectory. Objects should have subdirectories for each telescope's data.
data_dir = "/Users/hunterholland/Documents/Research/Laidlaw/Data/Modified"

object_dict = {}
# Creates dictionary of object names and their directories' file paths.
for entry in os.scandir(data_dir):
    if entry.is_dir() and not entry.name.startswith('.'):
        object_dict[entry.name] = entry.path
objects = list(object_dict.values())
objects.sort()


class Telescope:
    def __init__(self, file):
        if os.path.isfile(file):  # Assign input to filepath
            self.path = file
            self.file_dir, self.filename = os.path.split(file)
        else:  # Get filepath from file input
            for directory, _subdir, files in os.walk(f"{data_dir}"):
                if file in files:
                    self.path = os.path.join(directory, file)
                    self.file_dir = directory
                    self.filename = file
        # Using self.path, get object name from object dictionary.
        for item in object_dict.keys():
            if item in self.path:
                self.objectname = item
        # Get object path from object directory.
        obj_path = object_dict[self.objectname]
        self.telescopes = []
        for telescope in os.scandir(obj_path):
            if telescope.is_dir() and "." not in telescope.name:
                # List telescopes with object data.
                self.telescopes.append(telescope.name)

    def get_objectname(self):
        """Returns the name of the object.
        """
        return self.objectname

    def get_filepath(self):
        """Returns the filepath of the current datafile.
        """
        return self.path

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
            init_message = f"{self.objectname} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.path}"

    def e_hist(self, min_e, max_e, nbins='auto'):
        pass


class Chandra(Telescope):
    """Takes a file path or name as input. As of now, this class can initiate data, yield energy histograms, and filter files based on photon intensity (energy) and position.
    """

    def __init__(self, file):
        self.telescope = "Chandra"
        super().__init__(file)
        if not self.telescope in self.telescopes:
            init_message = f"{self.objectname} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.path}"

    def e_hist(self, e_range=None, e_list=None, e_list2=None, nbins='auto', object=False, save=False, filename=None):
        """Makes a histogram over the specified energy range or, optionally, of up to two lists of energies passed as input. Optionally saves output as PNG active file's directory.

        To specify an energy range from a minimum energy A to a maximum energy B, use list notation: [A, B].
        """
        if e_list is not None:
            # Use list of energies instead of specified ranges to make histogram
            e_band = e_list
        else:
            try:  # If a range with a max and min is passed as input
                min_e, max_e = e_range
            except ValueError:  # If a range with only a min is passed as input
                min_e, max_e = e_range, float("inf")
            except TypeError:  # If no range is passed as input
                min_e, max_e = 0, float("inf")
            evt_data = fits.getdata(self.path)  # Get data from event file
            energy = evt_data["energy"]  # Extract energy data
            min_thresh = energy >= min_e  # Establish min filter
            max_thresh = energy < max_e  # Establish max filter
            # Filter energy for everything between specified min and max values
            e_band = energy[min_thresh & max_thresh]
        plt.hist(e_band, bins=nbins)
        if e_list2 is not None:
            # Overlay second histogram if data is present
            plt.hist(e_list2, bins=nbins)
        plt.xlabel("Energy (eV)")
        plt.ylabel("Count")
        if not object:
            plt.title(f"Energy Distribution")
        elif object:
            plt.title(f"{object} Energy Distribution")
        if save == True:
            if object:
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

    def e_mask(self, *args, newfile=False, filename=None):
        """If newfile=False, method returns a list of masked energies that fall within the given ranges. Each range should be specifed as a list of length 2. For example, a range of energies from 600-1000eV would be denoted [600, 1000]. The final range may be of length 1 — for example, [1000]. This defaults to a range [1000, infinity].

        If newfile=True, a new file will be created that only contains rows with energies within the specified ranges.
        """
        evt_data = fits.getdata(self.path)  # Get event data from file
        energy = evt_data["energy"]  # Get energy data from event data
        # Create False boolean filter list based on the energy data. This master list will be referenced and altered by each range passed as input.
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
            # Master filter list (filter_list) is now an ordered list of bools. For each value in the energy data, if it is within the specified energy ranges, the corresponding index in filter_list is now True.
        if newfile == False:
            energies = energy[filter_list]
            return energies  # Return list of filtered energies.
        else:
            with fits.open(self.path) as hdul:
                evt_table = hdul[1]
                evt_data = evt_table.data
                # Use energy filter to filter events.
                evt_data = evt_data[filter_list]
                if filename is None:
                    filename = f"{self.filename}"
                try:
                    # Save new file to directory of original file.
                    hdul.writeto(f"{self.file_dir}/{filename}")
                    print(f"Data filtered into file {filename}")
                except NameError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")
                except OSError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")

    def coord_mask(self, center_RA, center_Dec, shape="box", size_RA=None, size_Dec=None, radius=None, radius_RA=None, radius_Dec=None, rotation=None, newfile=False, filename=None):
        """Given a set of values describing the shape of a region, this method returns the list of energies within that region, which can then be used for analysis (filter by energy, make a histogram, etc). If newfile=True, a new FITS file will be created from the original, containing only the event data for everything within the region.

        shapes: box, circle, ellipse (only box has functionality, currently)

        This method was designed with reference to the application "SAOImage ds9" and its ability to overlay regions on a given image. Once a region is created and selected, go to Region->Get Information to view its size parameters. Change the units to "detector," and those numbers may be used as input for this method.
        """
        evt_data = fits.getdata(self.path)  # Get event data from file
        ra = evt_data["x"]  # Get RA coordinates from event data
        dec = evt_data["y"]  # Get Dec coordinates from event data
        # Create True boolean filter list the length of the coordinate data. This will be referenced and altered by each successive coordinate filter (ra_range and dec_range).
        filter_list = [True for i in range(len(ra))]
        if shape == "box":
            radius_RA, radius_Dec = size_RA/2, size_Dec/2
            # Get min and max RA and Dec values based on the given parameters and shape geometry.
            min_RA, max_RA = center_RA-radius_RA, center_RA+radius_RA
            min_Dec, max_Dec = center_Dec-radius_Dec, center_Dec+radius_Dec
            # Define bool lists the length of the coordinate data —— one for each coordinate. Values are True if they lie in the region and False if otherwise.
            ra_range = (ra >= min_RA) & (ra <= max_RA)
            dec_range = (dec >= min_Dec) & (dec <= max_Dec)
        elif shape == "circle":
            pass
        elif shape == "ellipse":
            pass
        else:
            print("Shape must be either \"circle\", \"ellipse\", or \"box\".")
            return
        coord_ranges = [ra_range, dec_range]
        # For each coordinate filter list, iterate through each value of the master filter list. This results in filter_list being a list of bools that reads True for any event that falls within the specified region.
        for coordinate in coord_ranges:
            for i in range(len(filter_list)):
                if coordinate[i] == False:
                    filter_list[i] = False
        if newfile == False:
            with fits.open(self.path) as hdul:
                evt_table = hdul[1]
                evt_table.data = evt_table.data[filter_list]
                # Energies present in the specified region.
                energies_in_region = evt_table.data["energy"]
            return energies_in_region
        else:
            with fits.open(self.path) as hdul:
                evt_table = hdul[1]
                evt_table.data = evt_table.data[filter_list]
                if filename is None:
                    filename = f"{self.filename}"
                try:
                    hdul.writeto(f"{self.file_dir}/{filename}")
                    print(f"Data filtered into file {filename}")
                except NameError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")
                except OSError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")


class XMM(Telescope):
    """Takes a file path or name as input. As of now, this class can initiate data, yield energy histograms, and filter files based on photon intensity (energy) and position.
    """

    def __init__(self, file):
        self.telescope = "XMM"
        super().__init__(file)
        if not self.telescope in self.telescopes:
            init_message = f"{self.objectname} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.path}"

    def e_hist(self, e_range=None, e_list=None, e_list2=None, nbins='auto', object=False, save=False, filename=None):
        """Makes a histogram over the specified energy range or, optionally, of up to two lists of energies passed as input. Optionally saves output as PNG active file's directory.

        To specify an energy range from a minimum energy A to a maximum energy B, use list notation: [A, B].
        """
        if e_list is not None:
            # Use list of energies instead of specified ranges to make histogram
            e_band = e_list
        else:
            try:  # If a range with a max and min is passed as input
                min_e, max_e = e_range
            except ValueError:  # If a range with only a min is passed as input
                min_e, max_e = e_range, float("inf")
            except TypeError:  # If no range is passed as input
                min_e, max_e = 0, float("inf")
            evt_data = fits.getdata(self.path)  # Get data from event file
            energy = evt_data["PI"]  # Extract energy data
            min_thresh = energy >= min_e  # Establish min filter
            max_thresh = energy < max_e  # Establish max filter
            # Filter energy for everything between specified min and max values
            e_band = energy[min_thresh & max_thresh]
        plt.hist(e_band, bins=nbins)
        if e_list2 is not None:
            # Overlay second histogram if data is present
            plt.hist(e_list2, bins=nbins)
        plt.xlabel("Energy (eV)")
        plt.ylabel("Count")
        if not object:
            plt.title(
                f"Energy Distribution")
        elif object:
            plt.title(
                f"{object} Energy Distribution")
        if save == True:
            if object:
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

    def e_mask(self, *args, newfile=False, filename=None):
        """If newfile=False, method returns a list of masked energies that fall within the given ranges. Each range should be specifed as a list of length 2. For example, a range of energies from 600-1000eV would be denoted [600, 1000]. The final range may be of length 1 — for example, [1000]. This defaults to a range [1000, infinity].

        If newfile=True, a new file will be created that only contains rows with energies within the specified ranges.
        """
        evt_data = fits.getdata(self.path)  # Get event data from file
        energy = evt_data["PI"]  # Get energy data from event data
        # Create False boolean filter list based on the energy data. This master list will be referenced and altered by each range passed as input.
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
            # Master filter list (filter_list) is now an ordered list of bools. For each value in the energy data, if it is within the specified energy ranges, the corresponding index in filter_list is now True.
        if newfile == False:
            energies = energy[filter_list]
            return energies  # Return list of filtered energies.
        else:
            with fits.open(self.path) as hdul:
                evt_data = hdul[1].data
                # Use energy filter to filter events.
                evt_data = evt_data[filter_list]
                if filename is None:
                    filename = f"{self.filename}"
                try:
                    # Save new file to directory of original file.
                    hdul.writeto(f"{self.file_dir}/{filename}")
                    print(f"Data filtered into file {filename}")
                except NameError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")
                except OSError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")

    def coord_mask(self, center_RA, center_Dec, shape="box", size_RA=None, size_Dec=None, radius=None, radius_RA=None, radius_Dec=None, rotation=None, newfile=False, filename=None):
        """Given a set of values describing the shape of a region, this method returns the list of energies within that region, which can then be used for analysis (filter by energy, make a histogram, etc). If newfile=True, a new FITS file will be created from the original, containing only the event data for everything within the region.

        shapes: box, circle, ellipse (only box has functionality, currently)

        This method was designed with reference to the application "SAOImage ds9" and its ability to overlay regions on a given image. Once a region is created and selected, go to Region->Get Information to view its size parameters. Change the units to "detector," and those numbers may be used as input for this method.
        """
        evt_data = fits.getdata(self.path)  # Get event data from file
        ra = evt_data["X"]  # Get RA coordinates from event data
        dec = evt_data["Y"]  # Get Dec coordinates from event data
        # Create True boolean filter list the length of the coordinate data. This will be referenced and altered by each successive coordinate filter (ra_range and dec_range).
        filter_list = [True for i in range(len(ra))]
        if shape == "box":
            radius_RA, radius_Dec = size_RA/2, size_Dec/2
            # Get min and max RA and Dec values based on the given parameters and shape geometry.
            min_RA, max_RA = center_RA-radius_RA, center_RA+radius_RA
            min_Dec, max_Dec = center_Dec-radius_Dec, center_Dec+radius_Dec
            # Define bool lists the length of the coordinate data —— one for each coordinate. Values are True if they lie in the region and False if otherwise.
            ra_range = (ra >= min_RA) & (ra <= max_RA)
            dec_range = (dec >= min_Dec) & (dec <= max_Dec)
        elif shape == "circle":
            pass
        elif shape == "ellipse":
            pass
        else:
            print("Shape must be either \"circle\", \"ellipse\", or \"box\".")
            return
        coord_ranges = [ra_range, dec_range]
        # For each coordinate filter list, iterate through each value of the master filter list. This results in filter_list being a list of bools that reads True for any event that falls within the specified region.
        for coordinate in coord_ranges:
            for i in range(len(filter_list)):
                if coordinate[i] == False:
                    filter_list[i] = False
        if newfile == False:
            with fits.open(self.path) as hdul:
                evt_table = hdul[1]
                evt_table.data = evt_table.data[filter_list]
                # Energies present in the specified region.
                energies_in_region = evt_table.data["PI"]
            return energies_in_region
        else:
            with fits.open(self.path) as hdul:
                evt_table = hdul[1]
                evt_table.data = evt_table.data[filter_list]
                if filename is None:
                    filename = f"{self.filename}"
                try:
                    hdul.writeto(f"{self.file_dir}/{filename}")
                    print(f"Data filtered into file {filename}")
                except NameError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")
                except OSError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")


class Rosat(Telescope):
    """FITS.OPEN CURRENTLY BROKEN FOR ROSAT FILES. STILL TRYING TO TROUBLESHOOT.

    Takes a file path or name as input. As of now, this class can initiate data, yield energy histograms, and filter files based on photon intensity (energy) and position.
    """

    def __init__(self, file):
        self.telescope = "ROSAT"
        super().__init__(file)
        if not self.telescope in self.telescopes:
            init_message = f"{self.objectname} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.path}"

    def e_hist(self, e_range=None, e_list=None, e_list2=None, nbins='auto', object=False, save=False, filename=None):
        """Makes a histogram over the specified energy range or, optionally, of up to two lists of energies passed as input. Optionally saves output as PNG active file's directory.

        To specify an energy range from a minimum energy A to a maximum energy B, use list notation: [A, B].
        """
        if e_list is not None:
            # Use list of energies instead of specified ranges to make histogram
            e_band = e_list
        else:
            try:  # If a range with a max and min is passed as input
                min_e, max_e = e_range
            except ValueError:  # If a range with only a min is passed as input
                min_e, max_e = e_range, float("inf")
            except TypeError:  # If no range is passed as input
                min_e, max_e = 0, float("inf")
            # Open data file
            hdul = fits.open(self.path, ignore_missing_end=True)
            evt_data = hdul[2].data  # Get event data from file
            energy = evt_data["PI"]  # Get energy data from event data
            min_thresh = energy >= min_e  # Establish min filter
            max_thresh = energy < max_e  # Establish max filter
            # Filter energy for everything between specified min and max values
            e_band = energy[min_thresh & max_thresh]
        plt.hist(e_band, bins=nbins)
        if e_list2 is not None:
            # Overlay second histogram if data is present
            plt.hist(e_list2, bins=nbins)
        plt.xlabel("Energy (eV)")
        plt.ylabel("Count")
        if not object:
            plt.title(
                f"Energy Distribution")
        elif object:
            plt.title(
                f"{object} Energy Distribution")
        if save == True:
            if object:
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

    def e_mask(self, *args, newfile=False, filename=None):
        """If newfile=False, method returns a list of masked energies that fall within the given ranges. Each range should be specifed as a list of length 2. For example, a range of energies from 600-1000eV would be denoted [600, 1000]. The final range may be of length 1 — for example, [1000]. This defaults to a range [1000, infinity].

        If newfile=True, a new file will be created that only contains rows with energies within the specified ranges.
        """
        hdul = fits.open(self.path)  # Open data file
        evt_data = hdul[2].data  # Get event data from file
        energy = evt_data["PI"]  # Get energy data from event data
        # Create False boolean filter list based on the energy data. This master list will be referenced and altered by each range passed as input.
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
            # Master filter list (filter_list) is now an ordered list of bools. For each value in the energy data, if it is within the specified energy ranges, the corresponding index in filter_list is now True.
        if newfile == False:
            energies = energy[filter_list]
            return energies  # Return list of filtered energies.
        else:
            with fits.open(self.path) as hdul:
                evt_data = hdul[2].data
                # Use energy filter to filter events.
                evt_data = evt_data[filter_list]
                if filename is None:
                    filename = f"{self.filename}"
                try:
                    # Save new file to directory of original file.
                    hdul.writeto(f"{self.file_dir}/{filename}")
                    print(f"Data filtered into file {filename}")
                except NameError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")
                except OSError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")

    def coord_mask(self, center_RA, center_Dec, shape="box", size_RA=None, size_Dec=None, radius=None, radius_RA=None, radius_Dec=None, rotation=None, newfile=False, filename=None):
        """Given a set of values describing the shape of a region, this method returns the list of energies within that region, which can then be used for analysis (filter by energy, make a histogram, etc). If newfile=True, a new FITS file will be created from the original, containing only the event data for everything within the region.

        shapes: box, circle, ellipse (only box has functionality, currently)

        This method was designed with reference to the application "SAOImage ds9" and its ability to overlay regions on a given image. Once a region is created and selected, go to Region->Get Information to view its size parameters. Change the units to "detector," and those numbers may be used as input for this method.
        """
        hdul = fits.open(self.path)  # Open data file
        evt_data = hdul[2].data  # Get event data from file
        ra = evt_data["X"]  # Get RA coordinates from event data
        dec = evt_data["Y"]  # Get Dec coordinates from event data
        # Create True boolean filter list the length of the coordinate data. This will be referenced and altered by each successive coordinate filter (ra_range and dec_range).
        filter_list = [True for i in range(len(ra))]
        if shape == "box":
            radius_RA, radius_Dec = size_RA/2, size_Dec/2
            # Get min and max RA and Dec values based on the given parameters and shape geometry.
            min_RA, max_RA = center_RA-radius_RA, center_RA+radius_RA
            min_Dec, max_Dec = center_Dec-radius_Dec, center_Dec+radius_Dec
            # Define bool lists the length of the coordinate data —— one for each coordinate. Values are True if they lie in the region and False if otherwise.
            ra_range = (ra >= min_RA) & (ra <= max_RA)
            dec_range = (dec >= min_Dec) & (dec <= max_Dec)
        elif shape == "circle":
            pass
        elif shape == "ellipse":
            pass
        else:
            print("Shape must be either \"circle\", \"ellipse\", or \"box\".")
            return
        coord_ranges = [ra_range, dec_range]
        # For each coordinate filter list, iterate through each value of the master filter list. This results in filter_list being a list of bools that reads True for any event that falls within the specified region.
        for coordinate in coord_ranges:
            for i in range(len(filter_list)):
                if coordinate[i] == False:
                    filter_list[i] = False
        if newfile == False:
            with fits.open(self.path) as hdul:
                evt_table = hdul[2]
                evt_table.data = evt_table.data[filter_list]
                # Energies present in the specified region.
                energies_in_region = evt_table.data["PI"]
            return energies_in_region
        else:
            with fits.open(self.path) as hdul:
                evt_table = hdul[2]
                evt_table.data = evt_table.data[filter_list]
                if filename is None:
                    filename = f"{self.filename}"
                try:
                    hdul.writeto(f"{self.file_dir}/{filename}")
                    print(f"Data filtered into file {filename}")
                except NameError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")
                except OSError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")


class Swift(Telescope):
    """Takes a file path or name as input. As of now, this class can initiate data, yield energy histograms, and filter files based on photon intensity (energy) and position.
    """

    def __init__(self, file):
        self.telescope = "Swift"
        super().__init__(file)
        if not self.telescope in self.telescopes:
            init_message = f"{self.objectname} has no data from {self.telescope} telescope."
        else:
            init_message = f"{self.telescope} data initiated."
        print(init_message)

    def __repr__(self):
        return f"{self.telescope} object from path {self.path}"

    def e_hist(self, e_range=None, e_list=None, e_list2=None, nbins='auto', object=False, save=False, filename=None):
        """Makes a histogram over the specified energy range or, optionally, of up to two lists of energies passed as input. Optionally saves output as PNG active file's directory.

        To specify an energy range from a minimum energy A to a maximum energy B, use list notation: [A, B].
        """
        if e_list is not None:
            # Use list of energies instead of specified ranges to make histogram
            e_band = e_list
        else:
            try:  # If a range with a max and min is passed as input
                min_e, max_e = e_range
            except ValueError:  # If a range with only a min is passed as input
                min_e, max_e = e_range, float("inf")
            except TypeError:  # If no range is passed as input
                min_e, max_e = 0, float("inf")
            evt_data = fits.getdata(self.path)  # Get data from event file
            energy = evt_data["PI"]  # Extract energy data
            min_thresh = energy >= min_e  # Establish min filter
            max_thresh = energy < max_e  # Establish max filter
            # Filter energy for everything between specified min and max values
            e_band = energy[min_thresh & max_thresh]
        plt.hist(e_band, bins=nbins)
        if e_list2 is not None:
            # Overlay second histogram if data is present
            plt.hist(e_list2, bins=nbins)
        plt.xlabel("Energy (eV)")
        plt.ylabel("Count")
        if not object:
            plt.title(
                f"Energy Distribution")
        elif object:
            plt.title(
                f"{object} Energy Distribution")
        if save == True:
            if object:
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

    def e_mask(self, *args, newfile=False, filename=None):
        """If newfile=False, method returns a list of masked energies that fall within the given ranges. Each range should be specifed as a list of length 2. For example, a range of energies from 600-1000eV would be denoted [600, 1000]. The final range may be of length 1 — for example, [1000]. This defaults to a range [1000, infinity].

        If newfile=True, a new file will be created that only contains rows with energies within the specified ranges.
        """
        evt_data = fits.getdata(self.path)  # Get event data from file
        energy = evt_data["PI"]  # Get energy data from event data
        # Create False boolean filter list based on the energy data. This master list will be referenced and altered by each range passed as input.
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
            # Master filter list (filter_list) is now an ordered list of bools. For each value in the energy data, if it is within the specified energy ranges, the corresponding index in filter_list is now True.
        if newfile == False:
            energies = energy[filter_list]
            return energies  # Return list of filtered energies.
        else:
            with fits.open(self.path) as hdul:
                evt_data = hdul[1].data
                # Use energy filter to filter events.
                evt_data = evt_data[filter_list]
                if filename is None:
                    filename = f"{self.filename}"
                try:
                    # Save new file to directory of original file.
                    hdul.writeto(f"{self.file_dir}/{filename}")
                    print(f"Data filtered into file {filename}")
                except NameError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")
                except OSError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")

    def coord_mask(self, center_RA, center_Dec, shape="box", size_RA=None, size_Dec=None, radius=None, radius_RA=None, radius_Dec=None, rotation=None, newfile=False, filename=None):
        """Given a set of values describing the shape of a region, this method returns the list of energies within that region, which can then be used for analysis (filter by energy, make a histogram, etc). If newfile=True, a new FITS file will be created from the original, containing only the event data for everything within the region.

        shapes: box, circle, ellipse (only box has functionality, currently)

        This method was designed with reference to the application "SAOImage ds9" and its ability to overlay regions on a given image. Once a region is created and selected, go to Region->Get Information to view its size parameters. Change the units to "detector," and those numbers may be used as input for this method.
        """
        evt_data = fits.getdata(self.path)  # Get event data from file
        ra = evt_data["x"]  # Get RA coordinates from event data
        dec = evt_data["y"]  # Get Dec coordinates from event data
        # Create True boolean filter list the length of the coordinate data. This will be referenced and altered by each successive coordinate filter (ra_range and dec_range).
        filter_list = [True for i in range(len(ra))]
        if shape == "box":
            radius_RA, radius_Dec = size_RA/2, size_Dec/2
            # Get min and max RA and Dec values based on the given parameters and shape geometry.
            min_RA, max_RA = center_RA-radius_RA, center_RA+radius_RA
            min_Dec, max_Dec = center_Dec-radius_Dec, center_Dec+radius_Dec
            # Define bool lists the length of the coordinate data —— one for each coordinate. Values are True if they lie in the region and False if otherwise.
            ra_range = (ra >= min_RA) & (ra <= max_RA)
            dec_range = (dec >= min_Dec) & (dec <= max_Dec)
        elif shape == "circle":
            pass
        elif shape == "ellipse":
            pass
        else:
            print("Shape must be either \"circle\", \"ellipse\", or \"box\".")
            return
        coord_ranges = [ra_range, dec_range]
        # For each coordinate filter list, iterate through each value of the master filter list. This results in filter_list being a list of bools that reads True for any event that falls within the specified region.
        for coordinate in coord_ranges:
            for i in range(len(filter_list)):
                if coordinate[i] == False:
                    filter_list[i] = False
        if newfile == False:
            with fits.open(self.path) as hdul:
                evt_table = hdul[1]
                evt_table.data = evt_table.data[filter_list]
                # Energies present in the specified region.
                energies_in_region = evt_table.data["PI"]
            return energies_in_region
        else:
            with fits.open(self.path) as hdul:
                evt_table = hdul[1]
                evt_table.data = evt_table.data[filter_list]
                if filename is None:
                    filename = f"{self.filename}"
                try:
                    hdul.writeto(f"{self.file_dir}/{filename}")
                    print(f"Data filtered into file {filename}")
                except NameError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")
                except OSError:
                    print(
                        f"""A file with the name {filename} already exists in {self.file_dir}\nPlease use the \"filename\" parameter to give it another name""")
