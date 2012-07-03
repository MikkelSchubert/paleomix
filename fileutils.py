import os
import uuid
import threading
import functools


class GlobalCache:
    _enabled = 0
    _cache   = {}
    _hits    = 0
    _misses  = 0

    def __init__(self):
        pass

    def __enter__(self):
        assert threading.current_thread().name == "MainThread"
        GlobalCache._enabled += 1


    def __exit__(self, _exc_type, _exc_value, _traceback):
        assert threading.current_thread().name == "MainThread"
        assert GlobalCache._enabled

        GlobalCache._enabled -= 1
        if not GlobalCache._enabled:
            GlobalCache._cache = {}


    @classmethod
    def cache_enabled(cls, obj):
        @functools.wraps(obj)
        def memoizer(*args, **kwargs):
            if cls._enabled:
                assert threading.current_thread().name == "MainThread"

                key = (obj, tuple(args), frozenset(kwargs.iteritems()))
                if key not in GlobalCache._cache:
                    cls._misses += 1
                    result = GlobalCache._cache[key] = obj(*args, **kwargs)
                    return result
                else:
                    cls._hits += 1
                    return GlobalCache._cache[key]
            else:
                return obj(*args, **kwargs)

        return memoizer


@GlobalCache.cache_enabled
def file_is_empty(filename):
    if not os.path.exists(filename):
        return True

    length = os.stat(filename).st_size
    if length > 20:
        return False

    with open(filename, "rb") as fobj:
        magic = fobj.read(2)
        if not magic:
            return True
        elif magic == "BZ":
            return (length <= 14) # BZip2
        elif magic == "\x1f\x8b":
            return (length <= 20) # GZip
       
        return False


def valid_file(filename):
    return not file_is_empty(filename)


def add_postfix(filename, postfix):
    filename, ext = os.path.splitext(filename)
    return filename + postfix + ext


def swap_ext(filename, ext):
    filename, _ = os.path.splitext(filename)
    if not ext.startswith("."):
        ext = "." + ext

    return filename + ext    


def create_temp_dir(root):
    while True:
        uuid4 = str(uuid.uuid4())
        path = os.path.join(root, uuid4)
    
        if not os.path.exists(path):
            os.makedirs(path, mode = 0700)
            return path


def missing_files(filenames):
    """Given a list of filenames, returns a list of those that
    does not exist. Note that this function does not differentiate
    between files and folders."""
    result = []
    for filename in filenames:
        if not valid_file(filename):
            result.append(filename)
            
    return result


@GlobalCache.cache_enabled
def modified_after(younger, older):
    """Returns true any of the files expected to be 'younger' have
    been modified after any of the files expected to be 'older'."""
    def get_mtimes(filenames):
        for filename in filenames:
            yield os.path.getmtime(os.path.realpath(filename))

    return max(get_mtimes(younger)) >= min(get_mtimes(older))


@GlobalCache.cache_enabled
def is_executable(filename):
    return os.path.isfile(filename) and os.access(filename, os.X_OK)


def executable_exists(filename):
    if os.path.dirname(filename):
        return is_executable(filename)

    for path in os.environ["PATH"].split(os.pathsep):
        if is_executable(os.path.join(path, filename)):
            return True

    return False
