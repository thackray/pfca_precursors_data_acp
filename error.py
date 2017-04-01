"""Exceptions for use in chembox."""

class ChemboxError(Exception):
    """Base class for chembox exceptions."""
    pass

class ParameterError(ChemboxError):
    """Exception raised when model parameter is undefined or invalid."""
    pass

class YamlBooleanError(ChemboxError):
    """Exception raised when yaml erroneously reads something as a boolean."""
    def __init__(self, string):
        self.spec = string
    def __str__(self):
        return repr('yaml is reading something as a boolean where it '+\
                    'shouldn\'t.  '+self.spec)

class ExternalForcingError(ChemboxError):
    """Exception raised when there is no external forcing."""
    pass
