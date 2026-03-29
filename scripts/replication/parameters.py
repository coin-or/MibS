# Script to replicate experiments of paper 
# "Improving Directions in Mixed Integer Bilevel Linear Optimization"
# Battista F. & Ted K. Ralphs
# Last edited 2026

import os
from argparse import ArgumentDefaultsHelpFormatter
from jsonargparse import ArgumentParser

def _abs_path(path):
    if path is None:
        return None
    
    return os.path.abspath(os.path.expanduser(path))


def _parseCliPairs(argument, cliValues):
    """
    Convert alternating key-value pairs CLI values into a dictionary.

    Example input:
        ["KEY1", "VALUE1", "KEY2", "VALUE2"]

    Example output:
        {
            "KEY1": "VALUE1",
            "KEY2": "VALUE2",
        }
    """
    if cliValues is None:
        return {}

    if len(cliValues) % 2 != 0:
        raise ValueError(f"Incorrect number of argument for {argument}. KEY VALUE pairs required.")

    resultDict = {}
    for i in range(0, len(cliValues), 2):
        name = cliValues[i]
        value = cliValues[i + 1]
        resultDict[name] = _abs_path(value)

    return resultDict


class Parameters:
    """
    Shared CLI parser and parameter container.

    The class defines the common command-line arguments used by the
    replication scripts, parses them, and stores the resulting values in
    `self.parameters`.

    This class performs only syntactic validation of the arguments.
    Semantic validation (e.g., checking that required arguments are present, 
    or that file paths exist) should be performed separately after parsing 
    by the user.
    """

    def __init__(self):
        self.parser = ArgumentParser(
            description="Command-line arguments for the replication scripts.",
            formatter_class=ArgumentDefaultsHelpFormatter,
        )
        self.parameters = dict()

        self.addArguments()

    def __getitem__(self, key):
        """
        Return a parsed parameter by name.

        This allows dictionary-style access such as `params["outputDir"]`
        after `parse()` has been called.
        """
        return self.parameters[key]

    def addArguments(self):
        # run.py parameters
        self.parser.add_argument(
            "--binariesPath",
            nargs="+",
            default=None,
            help="Executable binaries to use for the experiments, provided as VERSION PATH pairs.",
        )

        self.parser.add_argument(
            "--instanceDirs",
            nargs="+",
            default=None,
            help="Instance-set directories to process, provided as DATASET DIRECTORY pairs.",
        )

        self.parser.add_argument(
            "--outputDir",
            default=os.path.join(".", "results"),
            help="Directory where experiment outputs, summaries, and generated files will be written.",
        )

        self.parser.add_argument(
            "--writeParams",
            action="store_true",
            help="Write parameter files to <outputDir>/parameters instead of passing parameters on the \
                    command line.",
        )

        self.parser.add_argument(
            "--testName",
            default="mibs",
            help="Job name to use when submitting experiments through qsub.",
        )

        self.parser.add_argument(
            "--pbsFile",
            default=None,
            help="Path to the PBS submission script to use with qsub. If omitted, experiments are run \
                    locally.",
        )

        # make_plots.py parameters

        self.parser.add_argument(
            "--dataSets",
            nargs="+",
            default=None,
            help="Dataset names to include when generating plots and summaries.",
        )

        self.parser.add_argument(
            "--aggregateDatasets",
            action="store_true",
            help="Aggregate the selected datasets into a single combined summary and set of plots.",
        )

        self.parser.add_argument(
            "--fileCsvIn",
            default=None,
            help="Path to an existing summary CSV file to read instead of parsing raw output files \
                    in <outputDir>.",
        )

        self.parser.add_argument(
            "--makePlotsImprovingDirectionsPaper",
            action="store_true",
            help="If set, use the specific plotting configuration for the 'Improving Directions in \
                    Mixed Integer Bilevel Linear Optimization' paper.",
        )


    def parse(self):
        """
        Parse and syntactically validate the command-line arguments and store the results.

        Raises:
            ValueError: If required arguments are missing, or
                if any key-value pair argument has an invalid format.
        """
        args = self.parser.parse_args()

        try:
            self.parameters["binariesPath"] = _parseCliPairs("--binariesPath", args.binariesPath)
            self.parameters["instanceDirs"] = _parseCliPairs("--instanceDirs", args.instanceDirs)
        except ValueError as exc:
            self.parser.error(str(exc))

        # Used in run.py, but not yet in make_plots.py 
        self.parameters['versions'] = list(self.parameters['binariesPath'].keys())

        for binName, path in self.parameters['binariesPath'].items():
            self.parameters['binariesPath'][binName] = _abs_path(path)

        for dataset, path in self.parameters['instanceDirs'].items():
            self.parameters['instanceDirs'][dataset] = _abs_path(path)

        self.parameters['outputDir'] = _abs_path(args.outputDir)
        self.parameters['writeParams'] = args.writeParams
        self.parameters['testName'] = args.testName
        self.parameters['pbsFile'] = _abs_path(args.pbsFile)
        self.parameters["dataSets"] = args.dataSets
        self.parameters['aggregateDatasets'] = args.aggregateDatasets
        self.parameters['fileCsvIn'] = _abs_path(args.fileCsvIn)
        self.parameters['makePlotsImprovingDirectionsPaper'] = args.makePlotsImprovingDirectionsPaper

        # TODO: read mibsParamsInputs and commonParams from CLI as well?
        # Each configuration could be a file with parameters, and we can have
        # a CLI argument for the directory where these files are stored.
        # Then we can read them in and pass to the runExperiments function.

if __name__ == "__main__":
    # Example usage
    params = Parameters()
    params.parse()

    print("Parsed parameters:")
    for key, value in params.parameters.items():
        print(f"{key}: {value}")