import os
import uuid
import shutil



def add_postfix(filename, postfix):
    filename, ext = os.path.splitext(filename)
    return filename + postfix + ext


def swap_ext(filename, ext):
    filename, _ = os.path.splitext(filename)
    if not ext.startswith("."):
        ext = "." + ext

    return filename + ext    


def reroot_path(root, filename):
    """Returns the basename of filename joined to root."""
    return os.path.join(root, os.path.basename(filename))


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
        if not os.path.exists(filename):
            result.append(filename)
            
    return result


def modified_after(younger, older):
    """Returns true any of the files expected to be 'younger' have
    been modified after any of the files expected to be 'older'."""
    def get_mtimes(filenames):
        for filename in filenames:
            yield os.path.getmtime(os.path.realpath(filename))

    return max(get_mtimes(younger)) >= min(get_mtimes(older))


def is_executable(filename):
    return os.path.isfile(filename) and os.access(filename, os.X_OK)


def executable_exists(filename):
    if os.path.dirname(filename):
        return is_executable(filename)

    for path in os.environ["PATH"].split(os.pathsep):
        if is_executable(os.path.join(path, filename)):
            return True

    return False


def missing_executables(filenames):
    result = []
    for filename in filenames:
        if not executable_exists(filename):
            result.append(filename)
            
    return result


def move_file(source, destination):
    directory = os.path.dirname(destination)
    if not os.path.exists(directory):
        os.makedirs(directory)

    shutil.move(source, destination)


def copy_file(source, destination):
    directory = os.path.dirname(destination)
    if not os.path.exists(directory):
        os.makedirs(directory)

    shutil.copy(source, destination)
