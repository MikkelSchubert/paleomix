import time
import multiprocessing

import ui




class Pypeline:
    def __init__(self, config):
        self._top_nodes = []
        self._config = config


    def add_node(self, node):
        self._top_nodes.append(node)


    def run(self, max_running = 4, dry_run = False):
        nodes = []
        if not dry_run:
            nodes = _collect_nodes(self._top_nodes)

        running, last_running = {}, None
        pool = multiprocessing.Pool(max_running)
        while True:
            if not self._check_running_nodes(running):
                break

            if running != last_running:
                while (len(running) < max_running):
                    node = _get_runable_node(nodes, running)
                    if not node:
                        break

                    running[node] = pool.apply_async(_call_run, args = (node, self._config))
                
                ui.print_node_tree(self._top_nodes, running)
                last_running = dict(running)

                if not _any_nodes_left(nodes, running):
                    break
                
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
