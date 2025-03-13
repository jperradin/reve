from tqdm import tqdm

from .config.settings import Settings
from .io.reader.reader_factory import ReaderFactory
from .core.system import System
from .utils import *

def main(settings: Settings):
    """
    Main function to test the package.
    """

    if settings.verbose:
        print(settings)
    reader = ReaderFactory().get_reader(settings.file_location)
    reader.set_verbose(settings.verbose)

    system = System(reader, settings)

    if settings.range_of_frames[1] == -1:
        total = system.get_num_frames()
    else:
        total = settings.range_of_frames[1] - settings.range_of_frames[0]

    progress_bar_kwargs = {
        "disable": not settings.verbose,
        "leave": True,
        "ncols": 10000
    }


    progress_bar = tqdm(enumerate(system.iter_frames()), desc="Parsing frames", unit="frame", **progress_bar_kwargs)
    for i, frame in progress_bar:
        if settings.lattice.apply_custom_lattice:
            frame.set_lattice(settings.lattice.custom_lattice)
        wpositions = frame.get_wrapped_positions()
        j = 22023
        k = 1500
        l = 22398
        p1 = wpositions[j]
        p2 = wpositions[k]
        p3 = wpositions[l]
        pbc_angle = geometry.calculate_pbc_angle(p1, p2, p3, frame.get_lattice())
        dir_angle = geometry.calculate_direct_angle(p1, p2, p3)

        atom1 = frame.get_atom_by_id(j)
        atom2 = frame.get_atom_by_id(k)
        atom3 = frame.get_atom_by_id(l)
        
        print(f"Frame {i}: angle {atom1.symbol} - {atom2.symbol} - {atom3.symbol} is {pbc_angle:.2f}")
        print(f"Frame {i}: angle {atom1.symbol} - {atom2.symbol} - {atom3.symbol} is {dir_angle:.2f}")

