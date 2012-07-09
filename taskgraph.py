

class TaskGraph:
    DONE, RUNNING, RUNABLE, QUEUED, OUTDATED, ERROR = range(6)


    def __init__(self, tasks):
        self._tasks        = list(tasks)
        self._states       = None
        self._fixed_states = {}


    def mark_task(self, task, state):
        self._fixed_states[task] = state
        self._states = None


    def clear_task_state(self, task):
        del self._fixed_states[task]
        self._states = None


    def get_state(self, task):
        return self._get_states()[task]


    def get_runable_tasks(self, max_tasks):
        max_tasks = max(0, max_tasks - self.running_tasks())

        for (task, state) in self._get_states().items():
            if not max_tasks:
                break
            elif state == self.RUNABLE:
                yield task
                max_tasks -= 1


    def any_runable_left(self):
        for state in self._get_states().itervalues():
            if state in (self.RUNNING, self.RUNABLE):
                return True

        return False

    
    def running_tasks(self):
        return self._get_states().values().count(self.RUNNING)


    def _get_states(self):
        if self._states is None:
            self._states = self._build_state_map(self._tasks, self._fixed_states)

        return self._states


    @classmethod
    def _build_state_map(cls, tasks, fixed_states):
        states = {}
        for task in tasks:
            cls._update_state_map(task, states, fixed_states)

        return states


    @classmethod
    def _update_state_map(cls, task, states, fixed_states):
        if task in states:
            return states[task]

        # Update sub-tasks, before checking for fixed states
        state = cls.DONE
        for subtask in task.subnodes:
            state = max(state, cls._update_state_map(subtask, states, fixed_states))

        if task in fixed_states:
            state = fixed_states[task]
        elif state == cls.DONE:
            if not task.is_done or task.is_outdated:
                state = cls.RUNABLE
        elif state in (cls.RUNNING, cls.RUNABLE, cls.QUEUED):
            if task.is_done:
                state = cls.OUTDATED
            else:
                state = cls.QUEUED

        states[task] = state

        return state
