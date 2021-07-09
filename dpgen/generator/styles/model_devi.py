from abc import ABC, abstractmethod
from typing import List, TYPE_CHECKING, Tuple

if TYPE_CHECKING:
    import numpy as np
    import dpdata

class ModelDeviEngien(ABC):
    """This is the base class of the model deviation engien."""
    def __init__(self, jdata: dict, mdata: dict):
        self.jdata = jdata
        self.mdata = mdata

    @abstractmethod
    def make_input(self, directory: str, system: dpdata.System, models: List[str]):
        """Make a simulation input from a initial system.
        
        Parameters
        ----------
        system: dpdata.System
            initial system to run simulations
        models: list[str]
            list of models
        """
        raise NotImplementedError("Not implemented")

    @abstractmethod
    def get_running_parameters(self) -> Tuple[str, List[str]]:
        """Get running parameters.

        Returns
        -------
        command: str
            the running command
        forward_files: list of str
            forward files
        backward_files: list of str
            backward files
        common_files: list of str
            common files
        """
        raise NotImplementedError("Not implemented")


    @abstractmethod
    def extract_trajectory(self, directory) -> 'Trajectory':
        """Obtain a Trajectory object from a directory.
        
        Returns
        -------
        trajectory: Trajectory
            the Trajectory object
        """
        return Trajectory(directory)


class Trajectory(ABC):
    """A trajectory containing lots of frames or clusters.
    
    Parameters
    ----------
    directory: str
        The directory of a simulation
    """
    def __init__(self, directory: str) -> None:
        self.directory = directory
        self.status = 0
    
    @abstractmethod
    def get_model_deviations(self) -> np.ndarray:
        """Get model deviations of all frames or clusters.

        Returns
        -------
        model_deviations: np.ndarray[float], 1D or 2D
            Model deviation of all frames or all clusters. 
            - frame: 1D, the maximum model deviation of all atoms
            - clusters: 2D, the atomic model deviation of centeral atom
        
        Notes
        -----
        The index of model deviation should be as the same as frames.
        """
        # for performance, returns a numpy array instead of a list
        raise NotImplementedError("Not implemented")
    
    @abstractmethod
    def get_frames(self, idx: List[Tuple[int]]) -> List["Frame"]:
        """Get list of frames from idx.

        Parameters
        ----------
        idx: list[tuple[int]]
            List of indexes.
        
        Returns
        -------
        frames: list[Frame]
            List of frames.
        """
        return [Frame(self, ii) for ii in idx]
    
    def open_trajectory(self):
        """Open the trajectory.
        
        Notes
        -----
        For example, open the file stream, initialize the parameters,
        and others.
        """
        self.status = 1

    def close_trajectory(self):
        """Close the trajectory. Call this method when all frames have
        been read.
        """
        self.status = 0
    
    def __lt__(self, other: "Trajectory") -> bool:
        return self.directory < other.directory


class Frame(ABC):
    """A frame or cluster of trajectory, which should have a model deviation.
    
    Parameters
    ----------
    trajectory: Trajectory
        The trajectory that this frame belongs to.
    idx: tuple[int]
        The index of the frame or cluster in the trajectory. Usually a frame needs
        1 index and a cluster needs 2.
    """
    def __init__(self, trajectory: Trajectory, idx: Tuple[int]) -> None:
        self.trajectory = trajectory
        self.idx = idx

    @abstractmethod
    def read_frame(self) -> dpdata.System:
        """Read the frame.
        
        Returns
        -------
        system: dpdata.System
            returns a system that contain the frame.
        """
        if self.trajectory.status == 0:
            self.trajectory.open_trajectory()
        raise NotImplementedError("Not implemented")

    def __lt__(self, other: "Frame") -> bool:
        # to only open a file at the same time
        if self.trajectory is other.trajectory:
            return self.idx < other.idx
        return self.trajectory < other.trajectory