import os
import os.path
import uuid
import threading
import functools


class GlobalCache:
    _enabled = 0
    _cache   = {}
    _hits    = 0
    _misses  = 0

    @classmethod
    def __enter__(cls):
        assert threading.current_thread().name == "MainThread"
        cls._enabled += 1

    @classmethod
    def __exit__(cls, _exc_type, _exc_value, _traceback):
        assert threading.current_thread().name == "MainThread"
        assert cls._enabled
        cls._enabled -= 1

        if not cls._enabled:
            cls._cache = {}

    @classmethod
    def cache_enabled(cls, obj):
        cache = obj.cache = {}

        @functools.wraps(obj)
        def memoizer(*args, **kwargs):
            if cls._enabled:
                assert threading.current_thread().name == "MainThread"

                key = (tuple(args), frozenset(kwargs.iteritems()))
                if key not in cache:
                    cls._misses += 1
                    result = cache[key] = obj(*args, **kwargs)
                    return result
                else:
                    cls._hits += 1
                    return cache[key]
            else:
                return obj(*args, **kwargs)

        return memoizer


@GlobalCache.cache_enabled
def file_is_empty(filename):
    if not os.path.exists(filename):
        return True

    with open(filename, "rb") as fobj:
        magic = fobj.read(2)
        if not magic:
            return True

        fobj.seek(0, os.SEEK_END)
        length = fobj.tell()
        
        if magic == "BZ":
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
