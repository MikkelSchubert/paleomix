import time
import multiprocessing

import ui


def _call_run(node):
    """Wrapper function, required in order to call Node.run()
    in subprocesses, since it is not possible to pickle 
    bound functions (e.g. self.run)"""
    return node.run()
        


class Pypeline:
    RUNNING, FINISHED = "Running", "Finished"

    def __init__(self, config):
        self._config = config
        self._nodes = []


    def add_node(self, node):
        self._nodes.append(node)


    def run(self, max_running = 4):
        running, last_running = {}, None
        pool = multiprocessing.Pool(max_running)
        states = self._check_states(self._nodes, {})

        self.print_nodes(self._nodes, states)
        while any((states.get(node) != Pypeline.FINISHED) for node in self._nodes):
            if not self._check_running_nodes(states, running):
                break

            while (len(running) < max_running):
                node = self._get_runnable_node(self._nodes, states)
                if not node:
                    break

                states[node] = Pypeline.RUNNING
                running[node] = pool.apply_async(_call_run, args = (node,))
                
            if running != last_running:
                self.print_nodes(self._nodes, states)
                last_running = dict(running)
                
            time.sleep(1)

        pool.close()
        pool.join()

        self._check_running_nodes(states, running)

        ui.print_msg("Done ...")


    @classmethod
    def _check_running_nodes(cls, states, running):
        errors = False
        for (node, proc) in running.items():
            if not proc.ready():
                continue

            states[node] = Pypeline.FINISHED
            running.pop(node)
                    
            try:
                error = "Run() returned false."
                result = proc.get()
            except Exception, error:
                result = False

            if not result:
                errors = True
                print "%s: Error occurred running command (terminating gracefully): %s" \
                    % (node, error)

        return not errors
 

    @classmethod
    def _check_states(cls, nodes, states):
        for node in nodes:
            if node not in states:
                if node.output_exists():
                    states[node] = cls.FINISHED
            
            cls._check_states(node.dependencies, states)

        return states

    @classmethod
    def _get_runnable_node(cls, nodes, states):
        for node in nodes:
            if states.get(node):
                continue
            elif all((states.get(dep) == cls.FINISHED) for dep in node.dependencies):
                return node
            elif any((states.get(node) is None) for dep in node.dependencies):
                node = cls._get_runnable_node(node.dependencies, states)
                if node:
                    return node
        return None


    @classmethod
    def print_nodes(cls, nodes, states):
        print
        ui.print_msg("Pipeline:")
        cls.print_sub_nodes(nodes, states, "   ")
        

    @classmethod
    def print_sub_nodes(cls, nodes, states, prefix = "", live = True):
        for node in nodes:
            is_live = live and not node.output_exists()
            if not is_live:
                print_func = ui.print_disabled
            elif states.get(node) == Pypeline.RUNNING:
                print_func = ui.print_info
            else:
                print_func = ui.print_msg
            
            print_func(prefix + "+ " + str(node))
            current_prefix = prefix + ("  " if (node == nodes[-1]) else "|  ")

            if node.dependencies:
                cls.print_sub_nodes(node.dependencies, states, current_prefix + "   ", is_live)
            else:
                print_func(current_prefix)
