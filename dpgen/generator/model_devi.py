from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Iterator, List, TYPE_CHECKING, Tuple, Iterator, Dict
from dpgen.plugin import Plugin

if TYPE_CHECKING:
    import numpy as np
    import dpdata

class ModelDeviEngien(ABC):
    """This is the base class of the model deviation engien."""

    __ModelDeviPlugin=Plugin()

    @staticmethod
    def register(key):
        return ModelDeviEngien.__ModelDeviPlugin.register(key)

    @staticmethod
    def engiens():
        return ModelDeviEngien.__ModelDeviPlugin.plugins
    
    @staticmethod
    def get_engien(key):
        try:
            return ModelDeviEngien.engiens()[key]
        except IndexError as e:
            raise RuntimeError("Unsupported engien!") from e

    def __init__(self, jdata: dict, mdata: dict):
        self.jdata = jdata
        self.mdata = mdata

    @abstractmethod
    def make_input(self, iter_idx:int, sys_idx: int, directory: Iterator[str], conf_name: str, models: List[str]):
        """Make a simulation input from a directory.

        The length of system should be only 1.
        
        Parameters
        ----------
        iter_idx: int
            index of iteration
        sys_idx: int
            index of system
        directory: Iterator
            generate the working directory to run simulations
        conf_name: str
            file name of initial system
        models: list[str]
            list of models
        """
        raise NotImplementedError("Not implemented")

    @abstractmethod
    def get_running_parameters(self, work_path: str) -> Dict[str]:
        """Get running parameters.

        Parameters
        ----------
        work_path: str
            the work path

        Returns
        -------
        parameters: dict[str]
            a dict which should contain the following keys
                command: str
                    the running command
                forward_files: list of str
                    forward files
                backward_files: list of str
                    backward files
                common_files: list of str
                    common files
        """
        return {
            "command": None,
            "forward_files": None,
            "backward_files": None,
            "common_files": None,
        }


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
    engien: ModelDeviEngien
        Engien
    directory: str
        The directory of a simulation
    """
    def __init__(self, engien: ModelDeviEngien, directory: str) -> None:
        self.engien = engien
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
    def get_frame(self, idx: int) -> "Frame":
        """Get a frame from idx.

        Parameters
        ----------
        idx: int
            Index.
        
        Returns
        -------
        frames: Frame
            Frames.
        
        Notes
        -----
        A new trajectory should implement thier own Frame, as the current
        Frame do nothing in `read_frame` method.
        """
        return Frame(self, idx)
    
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