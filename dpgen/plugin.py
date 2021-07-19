"""Base of plugin systems."""


class Plugin:
    """A class to register plugins.

    Examples
    --------
    >>> Plugin = Register()
    >>> @Plugin.register("xx")
        def xxx():
            pass
    >>> print(Plugin.plugins['xx'])
    """
    def __init__(self):
        self.plugins = {}

    def register(self, key):
        """Register a plugin.
        
        Parameter
        ---------
        key: str
            Key of the plugin.
        """
        def decorator(object):
            self.plugins[key] = object
            return object
        return decorator
    
    def get_plugin(self, key):
        return self.plugins[key]

    def __add__(self, other):
        self.plugins.update(other.plugins)
        return self
