import time
import signal
import multiprocessing

import ui
import taskgraph
import fileutils




class Pypeline:
    def __init__(self, config):
        self._top_nodes = []
        self._config = config


    def add_node(self, node):
        self._top_nodes.append(node)


    def run(self, max_running = 4, dry_run = False):
        running, last_running = {}, None
        pool = multiprocessing.Pool(max_running, _init_worker)
        nodes = taskgraph.TaskGraph(self._top_nodes)
        
        if dry_run:
            ui.print_node_tree(self._top_nodes, nodes)
            ui.print_msg("Dry run done ...")
            return 0
    
        try:
            while self._check_running_nodes(running, nodes):
                if running == last_running:
                    time.sleep(1)
                    continue
                elif not nodes.any_runable_left():
                    break

                for task in nodes.get_runable_tasks(max_running):
                    running[task] = pool.apply_async(_call_run, args = (task, self._config))
                    nodes.mark_task(task, nodes.RUNNING)
                last_running = dict(running)
                
                ui.print_node_tree(self._top_nodes, nodes)
        except KeyboardInterrupt:
            ui.print_err("Keyboard interrupt detected, terminating gracefully ...")
            pool.close()
            pool.join()
            return False

        pool.close()
        pool.join()

        if not self._check_running_nodes(running, nodes):
            ui.print_err("Errors were detected ...")
            return False

        ui.print_msg("Done ...")
        return True


    @classmethod
    def _check_running_nodes(cls, running, graph):
        errors = False
        for (node, proc) in running.items():
            if not proc.ready():
                continue

            running.pop(node)
            graph.clear_task_state(node)
                   
            try:
                error = "\tRun() returned false."
                result = proc.get()
            except Exception, error:
                result = False

            if not result:
                errors = True
                ui.print_err("%s: Error occurred running command (terminating gracefully):\n%s\n" \
                    % (node, "\n".join(("\t" + line) for line in str(error).strip().split("\n"))))

        return not errors
 



def _init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def _call_run(node, config):
    """Wrapper function, required in order to call Node.run()
    in subprocesses, since it is not possible to pickle 
    bound functions (e.g. self.run)"""
    return node.run(config)


def _collect_nodes(nodes):
    result = set()
    for node in nodes:
        result.add(node)
        result.update(_collect_nodes(node.subnodes))

    return result


def _get_runable_node(nodes, running):
    """Returns a node for which all subnodes are done, or None if 
    no such nodes exist."""

    for node in nodes:
        if (node in running) or node.is_done:
            continue
        elif node.is_runable:
            return node
        else:
            subnode = _get_runable_node(node.subnodes, running)
            if subnode:
                return subnode

    return None


def _any_nodes_left(nodes, running):
    for node in nodes:
        if (node not in running) and not node.is_done:
            return True
    return False
