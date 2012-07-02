import time
import multiprocessing

import ui


def _call_run(node, config):
    """Wrapper function, required in order to call Node.run()
    in subprocesses, since it is not possible to pickle 
    bound functions (e.g. self.run)"""
    return node.run(config)
        


class Pypeline:
    RUNNING, FINISHED = "Running", "Finished"

    def __init__(self, config):
        self._top_nodes = []
        self._config = config


    def add_node(self, node):
        self._top_nodes.append(node)


    def run(self, max_running = 4, dry_run = False):
        running, last_running = {}, None

        remaining_nodes = self._get_unfinished_nodes(self._top_nodes)
        if dry_run or not remaining_nodes:
            self.print_nodes(self._top_nodes, running, remaining_nodes)
            ui.print_msg("Done ...")
            return            

        pool = multiprocessing.Pool(max_running)
        while remaining_nodes:
            if not self._check_running_nodes(running):
                break

            while (len(running) < max_running):
                node = self._get_runnable_node(remaining_nodes)
                if not node:
                    break

                running[node] = pool.apply_async(_call_run, args = (node, self._config))
                remaining_nodes.remove(node)
                
            if running != last_running:
                self.print_nodes(self._top_nodes, running, remaining_nodes)
                last_running = dict(running)
                
            time.sleep(1)

        pool.close()
        pool.join()

        self._check_running_nodes(running)

        ui.print_msg("Done ...")


    @classmethod
    def _check_running_nodes(cls, running):
        errors = False
        for (node, proc) in running.items():
            if not proc.ready():
                continue

            running.pop(node)
                   
            try:
                error = "\tRun() returned false."
                result = proc.get()
            except Exception, error:
                result = False

            if not result:
                errors = True
                ui.print_err("%s: Error occurred running command (terminating gracefully):\n%s" \
                    % (node, error))

        return not errors
 
    
    @classmethod
    def _get_unfinished_nodes(cls, nodes):
        unfinished = set()
        for node in nodes:
            if node.output_status() != node.EXISTS:
                unfinished.add(node)
            unfinished.update(cls._get_unfinished_nodes(node.dependencies))
        return unfinished


    @classmethod
    def _get_runnable_node(cls, nodes):
        """Returns a node for which all dependencies are met, or None if 
        no such nodes exist."""

        for node in nodes:
            if node.output_status() == node.EXISTS:
                continue
            elif all((dep.output_status() == node.EXISTS) for dep in node.dependencies):
                return node

        return None


    @classmethod
    def print_nodes(cls, top_nodes, running, remaining):
        print
        ui.print_msg("Pipelin (%i nodes running, %i left):" \
                         % (len(running), len(remaining)))
        cls.print_sub_nodes(list(top_nodes), running, "   ")
        

    @classmethod
    def print_sub_nodes(cls, nodes, running, prefix = ""):
        for node in nodes:
            print_func = ui.print_info
            if node not in running:
                status = node.output_status()
                if status == node.EXISTS:
                    print_func = ui.print_disabled
                elif status == node.OUTDATED:
                    print_func = ui.print_warn
                else:
                    print_func = ui.print_msg
            
            print_func(prefix + "+ " + str(node))
            current_prefix = prefix + ("  " if (node == nodes[-1]) else "|  ")

            if node.dependencies:
                cls.print_sub_nodes(node.dependencies, running, current_prefix + "   ")
            else:
                print_func(current_prefix)
