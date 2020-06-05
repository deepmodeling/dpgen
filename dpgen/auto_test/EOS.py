from Property import Property

class EOS (Property) :
    def __init__ (self,
                  parameter) :
        """
        Constructor

        Parameters
        ----------
        parameters : dict
                A dict that specifies the interaction.
        """
        pass
        
    def make_confs(self, 
                   work_path,
                   path_to_equi):
        """
        Make configurations needed to compute the property. 
        The tasks directory will be named as work_path/task.xxxxxx

        Parameters
        ----------
        work_path : str
                The path where the task are located 
        path_to_equi : str
                The path of the task that equilibrated the configuration.
        """
        pass
        
    def _cmpt(self, 
              output_file,
              all_tasks,
              all_res):
        """
        Compute the property.

        Parameters
        ----------
        output_file:
                The file to output the property
        all_tasks : list of str
                The list of directories of tasks
        all_res : list of str
                The list of results
        """
        pass
