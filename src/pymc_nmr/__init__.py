#the initial __init__.py as a "holding" file

from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)
