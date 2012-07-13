#!/usr/bin/python

class TaskError(RuntimeError):
    pass
    


class TaskGraph:
    class Node:
        DONE, RUNNING, RUNABLE, QUEUED, OUTDATED, ERROR = range(6)

        def __init__(self, task):
            self.task          = task
            self._state        = None
            self._fixed_state  = None
            self.subnodes      = frozenset()
            self.dependencies  = frozenset()


        @property
        def state(self):
            if self._fixed_state is not None:
                return self._fixed_state

            assert self._state is not None
            return self._state


        def __str__(self):
            return str(self.task)


    def __init__(self, tasks):
        self._graph = None
        self._tasks = {}

        def collapse(subtasks):
            for task in subtasks:
                if task not in self._tasks:
                    self._tasks[task] = TaskGraph.Node(task)
                    collapse(task.subnodes)
                    collapse(task.dependencies)
        collapse(tasks)
                    
        # TODO: Check that all are 'Node's


    def set_task_state(self, node, state):
        if state not in (None, TaskGraph.Node.RUNNING, TaskGraph.Node.ERROR):
            raise ValueError("Cannot set states other than RUNNING and ERROR, or cleared (None).")
            
        self._tasks[node.task]._fixed_state = state
        self._graph = None


    def __iter__(self):
        """Returns a graph of tasks."""
        if self._graph is None:
            self._graph = self._build_graph(self._tasks)
        return iter(self._graph)


    def iterflat(self):
        if self._graph is None:
            self._graph = self._build_graph(self._tasks)
        return iter(self._tasks.values())


    @classmethod
    def _build_graph(cls, tasks):
        for (task, node) in tasks.iteritems():
            node._state       = None
            node.subnodes     = frozenset(tasks[subnode] for subnode in task.subnodes)
            node.dependencies = frozenset(tasks[dependency] for dependency in task.dependencies)

        top_nodes = cls._get_top_nodes(tasks)
        for node in top_nodes:
            cls._update_states(node)

        return top_nodes


    @classmethod
    def _get_top_nodes(cls, tasks):
        """Returns a sequence of nodes that are not depended upon by any other nodes."""

        dependencies = set()
        for (_, node) in tasks.iteritems():
            dependencies.update(node.subnodes)
            dependencies.update(node.dependencies)
        
        return set(tasks.itervalues()) - dependencies
        
 
    @classmethod
    def _update_states(cls, node):
        if node._state is not None:
            # Possibly return a fixed state
            return node.state

        # Update sub-tasks, before checking for fixed states
        state = node.DONE
        for subnode in node.subnodes:
            state = max(state, cls._update_states(subnode))
        for subnode in node.dependencies:
            state = max(state, cls._update_states(subnode))

        if state == node.DONE:
            if not node.task.is_done or node.task.is_outdated:
                state = node.RUNABLE
        elif state in (node.RUNNING, node.RUNABLE, node.QUEUED):
            if node.task.is_done:
                state = node.OUTDATED
            else:
                state = node.QUEUED
        node._state = state

        # Possibly return a fixed state
        return node.state
