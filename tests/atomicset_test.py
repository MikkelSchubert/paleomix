import nose

from flexmock import flexmock

from pypeline.atomiccmd import AtomicCmd, CmdError
from pypeline.atomicset import ParallelCmds, SequentialCmds


################################################################################
################################################################################
## Functions and properties with same expected behavior for each set

_COMMON_RETURNS = [1, 2, 3]
_COMMON_FUNCTIONS = [
    ("join",         [],      _COMMON_RETURNS),
    ("commit",       ["TMP"], None) ]
_COMMON_PROPERTIES = [
    ("executables",  _COMMON_RETURNS),
    ("input_files",  _COMMON_RETURNS),
    ("output_files", _COMMON_RETURNS) ]
    

def test_sets__common_functions():
    def build_func_mocks(func_name, arguments, return_value):
        for return_value in _COMMON_RETURNS:
            cmd_mock = flexmock(AtomicCmd(["ls"]))
            cmd_mock.should_receive(func_name).with_args(*arguments).and_return([return_value]).once
            yield cmd_mock

    def test_function(cmdclass, func_name, arguments, return_value):
        commandset = cmdclass(build_func_mocks(func_name, arguments, return_value))
        cmd_prop   = getattr(commandset, func_name)(*arguments)              
        assert cmd_prop == return_value

    for cmdclass in (ParallelCmds, SequentialCmds):
        for (func_name, arguments, return_value) in _COMMON_FUNCTIONS:
            yield test_function, cmdclass, func_name, arguments, return_value


def test_sets__common_properties():
    def build_prop_mocks(prop_name, return_value):
        for return_value in _COMMON_RETURNS:
            properties = dict([(prop_name, [return_value])])
            yield flexmock(AtomicCmd(["ls"]), **properties)

    def test_property(cmdclass, prop_name, return_value):
        commandset = cmdclass(build_prop_mocks(prop_name, return_value))
        cmd_prop   = getattr(commandset, prop_name)
        assert cmd_prop == return_value
        
    for cmdclass in (ParallelCmds, SequentialCmds):
        for (prop_name, return_values) in _COMMON_PROPERTIES:
            yield test_property, cmdclass, prop_name, return_values


def test_sets__str__():
    def test_function(cls):
        mocks = []
        for ii in range(3):
            cmd_mock = flexmock(AtomicCmd(["ls"]))
            cmd_mock.should_receive("__str__").with_args().and_return(str(10 ** ii))
            mocks.append(cmd_mock)
        commands = cls(mocks)
        description = str(commands)

        assert description == "[1, 10, 100]"

    yield test_function, ParallelCmds
    yield test_function, SequentialCmds



################################################################################
################################################################################
## Parallel commands

def test_parallel_commands__atomiccmds():
    mocks = []
    for _ in range(3):
        cmd_mock = flexmock(AtomicCmd(["ls"]))
        cmd_mock.should_receive('run').with_args("xTMPx").once
        mocks.append(cmd_mock)

    cmds = ParallelCmds(mocks)
    cmds.run("xTMPx")

@nose.tools.raises(CmdError)
def test_parallel_commands__reject_sequential():
    command = AtomicCmd(["ls"])
    seqcmd  = SequentialCmds([command])
    ParallelCmds([seqcmd])

def test_parallel_commands__accept_parallel():
    command = AtomicCmd(["ls"])
    parcmd  = ParallelCmds([command])
    ParallelCmds([parcmd])

@nose.tools.raises(CmdError)
def test_parallel_commands__reject_noncommand():
    ParallelCmds([object()])

@nose.tools.raises(CmdError)
def test_parallel_commands__reject_empty_commandset():
    ParallelCmds([])



################################################################################
################################################################################
## Sequential commands

def test_sequential_commands__atomiccmds():
    mocks = []
    for _ in range(3):
        cmd_mock = flexmock(AtomicCmd(["ls"]))
        cmd_mock.should_receive('run').with_args("xTMPx").ordered.once
        cmd_mock.should_receive('join').with_args().and_return([0]).ordered.once
        mocks.append(cmd_mock)
    
    cmds = SequentialCmds(mocks)
    cmds.run("xTMPx")

def test_sequential_commands__abort_on_error():
    cmd_mock_1 = flexmock(AtomicCmd(["ls"]))
    cmd_mock_1.should_receive('run').with_args("xTMPx").ordered
    cmd_mock_1.should_receive('join').with_args().and_return([0])
    cmd_mock_2 = flexmock(AtomicCmd(["ls"]))
    cmd_mock_2.should_receive('run').with_args("xTMPx").ordered
    cmd_mock_2.should_receive('join').with_args().and_return([1])
    cmd_mock_3 = flexmock(AtomicCmd(["ls"]))
    cmd_mock_3.should_receive('run').with_args("xTMPx").never
    cmd_mock_3.should_receive('join').with_args().and_return([None])
   
    cmds = SequentialCmds([cmd_mock_1, cmd_mock_2, cmd_mock_3])
    cmds.run("xTMPx")

    assert cmds.join() == [0, 1, None]


def test_sequential_commands__accept_parallel():
    command = AtomicCmd(["ls"])
    parcmd  = ParallelCmds([command])
    SequentialCmds([parcmd])

def test_sequential_commands__accept_sequential():
    command = AtomicCmd(["ls"])
    seqcmd  = SequentialCmds([command])
    SequentialCmds([seqcmd])

@nose.tools.raises(CmdError)
def test_sequential_commands__reject_noncommand():
    SequentialCmds([object()])


@nose.tools.raises(CmdError)
def test_sequential_commands__reject_empty_commandset():
    SequentialCmds([])
    
