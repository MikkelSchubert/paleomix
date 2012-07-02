import time
import multiprocessing

import ui




class Pypeline:
    def __init__(self, config):
        self._top_nodes = []
        self._config = config


    def add_node(self, node):
        self._top_nodes.append(node)


    def run(self, max_running = 4, shallow_run = True, dry_run = False):
        nodes = []
        if not dry_run:
            nodes = set(self._top_nodes) 
            if not shallow_run:
                nodes = _collect_nodes(self._top_nodes)

        running, last_running = {}, None
        pool = multiprocessing.Pool(min(len(nodes), max_running))
        while True:
            if not self._check_running_nodes(running):
                break

            if running != last_running:
                while (len(running) < max_running):
                    node = _get_runable_node(nodes, running)
                    if not node:
                        break

                    running[node] = pool.apply_async(_call_run, args = (node, self._config))
                
                self.print_nodes(self._top_nodes, running)
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
 
    

    @classmethod
    def print_nodes(cls, top_nodes, running):
        print
        ui.print_msg("Pipeline (%i nodes running):" % (len(running),))
        cls.print_sub_nodes(list(top_nodes), running, "   ")
        

    @classmethod
    def print_sub_nodes(cls, nodes, running, prefix = ""):
        for node in nodes:
            print_func = ui.print_info
            if node not in running:
                if node.is_done:
                    print_func = ui.print_disabled
                elif node.is_outdated:
                    print_func = ui.print_warn
                else:
                    print_func = ui.print_msg
            
            print_func(prefix + "+ " + str(node))
            current_prefix = prefix + ("  " if (node == nodes[-1]) else "|  ")

            if node.subnodes:
                cls.print_sub_nodes(node.subnodes, running, current_prefix + "   ")
            else:
                print_func(current_prefix)





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
