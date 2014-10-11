import numpy, netCDF4, argparse

def distribution_csv_to_nc(input_filename, output_filename):
    data = numpy.loadtxt(input_filename, delimiter=",")

    # throw away the first column (crystal index)
    data = data[:,1:]

    nc = netCDF4.Dataset(output_filename, "w")
    nc.createDimension("distribution_index", 1)
    nc.createDimension("crystal_index", data.shape[0])
    nc.createDimension("parameter_index", data.shape[1])

    var = nc.createVariable("distributions", "f8",
                            ("distribution_index",
                             "crystal_index",
                             "parameter_index"))
    labels = ["crystal axis x", "crystal axis y", "crystal axis z", "size",
              "dislocation density", "time of last recrystallization",
              "size at last recrystallization"]
    var.parameter_labels = ", ".join(labels)

    # save this as the first distribution
    var[0,:] = data

    nc.close()

parser = argparse.ArgumentParser()

parser.add_argument("INPUT", nargs=1)
parser.add_argument("OUTPUT", nargs=1)

options = parser.parse_args()

distribution_csv_to_nc(options.INPUT[0], options.OUTPUT[0])
