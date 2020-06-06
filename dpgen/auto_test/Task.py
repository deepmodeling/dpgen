from abc import ABC,abstractmethod

class Task(ABC):
    @abstractmethod
    def __init__ (self, 
                  paramters, 
                  path_to_poscar) :
        """
        Constructor

        Parameters
        ----------
        parameters : dict
                A dict that specifies the interaction.
        path_to_poscar : str
                The path to POSCAR. Indicating in which system the task will be initialized.
        """
        pass

    @abstractmethod
    def make_potential_files(self, 
                             output_dir):
        """
        Prepare potential files for a computational task.
        For example, the VASP prepares POTCAR. 
        DeePMD prepares frozen model(s).

        Parameters
        ----------
        output_dir : str
                The directory storing the potential files.
        """
        pass

    @abstractmethod
    def make_input_file(self,
                        output_dir, 
                        task_type, 
                        task_param):
        """
        Prepare input files for a computational task
        For example, the VASP prepares INCAR.
        LAMMPS (including DeePMD, MEAM...) prepares in.lammps.
        IMPORTANT: Task parameters should be stored in output_dir/task.json

        Parameters
        ----------
        output_dir : str
                The directory storing the input files.
        task_type : str
                Can be
                - "relaxation:": structure relaxation
                - "static": static computation calculates the energy, force... of a strcture
        task_parame: dict
                The parameters of the task.
                For example the VASP interaction can be provided with
                { "ediff": 1e-6, "ediffg": 1e-5 }        
        """
        pass

    @abstractmethod
    def compute (self):
        """
        Compute output of the task. 
        IMPORTANT: The output configuration should be converted and stored in a CONTCAR file.

        Returns
        -------
        result_dict: dict
                A dict that storing the result. For example:
                { "energy": xxx, "force": [xxx] }
        """
        pass

    @property
    @staticmethod
    @abstractmethod
    def forward_files(self):
        """
        Return forward files.
        """
        pass

    @property
    @staticmethod
    @abstractmethod
    def forward_common_files(self):
        """
        Return forward common files.
        """
        pass

    @property
    @staticmethod
    @abstractmethod
    def backward_files(self):
        """
        Return backward files.
        """
        pass

