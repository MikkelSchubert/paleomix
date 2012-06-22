import ui


class Pypeline:
    RUNNING, FINISHED = "Running", "Finished"

    def __init__(self, config):
        self._config = config
        self._nodes = []


    def add_node(self, node):
        self._nodes.append(node)


    def run(self, max_threads = 1):
        states = self._check_states(self._nodes, {})

        self.print_nodes(self._nodes)
        while any((states.get(node) != Pypeline.FINISHED) for node in self._nodes):
            node = self._get_runnable_node(self._nodes, states)

            states[node] = Pypeline.RUNNING
            if not node.run():
                print "Error occurred, terminating gracefully."
                break
            states[node] = Pypeline.FINISHED

            self.print_nodes(self._nodes)


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
    def print_nodes(cls, nodes):
        print
        ui.print_msg("Pipeline:")
        cls.print_sub_nodes(nodes, "   ")
        

    @classmethod
    def print_sub_nodes(cls, nodes, prefix = "", live = True):
        for node in nodes:
            is_live = live and not node.output_exists()
            print_func = (ui.print_msg if is_live else ui.print_disabled)
            
            print_func(prefix + "+ " + str(node))
            current_prefix = prefix + ("  " if (node == nodes[-1]) else "|  ")

            if node.dependencies:
                cls.print_sub_nodes(node.dependencies, current_prefix + "   ", is_live)
            else:
                print_func(current_prefix)

        
