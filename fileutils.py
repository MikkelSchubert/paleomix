import os
import os.path
import uuid


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


def is_executable(filename):
    return os.path.isfile(filename) and os.access(filename, os.X_OK)


def executable_exists(filename):
    if os.path.dirname(filename):
        print "dirname"
        return is_executable(filename)

    for path in os.environ["PATH"].split(os.pathsep):
        if is_executable(os.path.join(path, filename)):
            return True

    return False
