from pathlib import Path

class Source:
    """Source "abstract" class

    Not truly an abstract class, because otherwise all methods would need to be defined when inheriting
    and some sources do not make use of all of these methods
    """
    def __init__(self, sdir=Path('.'), **kwargs) -> None:
        self.sdir = sdir
        for k, v in kwargs.items():
            setattr(self, k, v)

    def check_files(self):
        """If processed files are present, return 0. If raw files are present, return 1. Otherwise return 2.
        """
        pass

    def get_dir(self):
        """Get directory of Source

        Returns:
            Path: Source directory
        """
        return self.sdir

    def download(self):
        """Download necessary files
        """
        pass

    def parse_drugs(self):
        """Parse drug file(s) into json
        """
        pass

    def parse_targets(self):
        """Parse target file(s) into json
        """
        pass

    def parse_associations(self):
        """Parse association file(s) into json
        """
        pass

    def load_drugs(self):
        """Load drugs from json
        """
        pass

    def load_targets(self):
        """Load targets from json
        """
        pass

    def load_associations(self):
        """Load associations from json
        """
        pass

    def run(self):
        """Execute necessary methods to fully pre-process Source data
        """
        pass