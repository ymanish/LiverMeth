"""
Pipeline for nucleosome free energy landscape calculation.

Submodules
----------
sequence_io     : FASTA / BED loading, nucleosome-occupied region extraction
methylation     : Bismark .cov loading, MN encoding, terminal reversion
sliding_window  : Sliding window generation and batching
parallel        : Pool-based parallel free energy calculation
"""

from .sequence_io import (
    load_fasta,
    load_atac_peaks,
    nucleosome_occupied_segments,
    atac_peak_segments,
    write_nucleosome_fasta,
    write_atac_fasta,
)
from .methylation import (
    find_cpg_positions,
    load_methylated_positions,
    apply_methylation_to_sequence,
)
from .sliding_window import (
    Subsequence,
    sliding_window_sequence,
    nuc_window_generator,
    batcher,
)
from .parallel import (
    run_pool_energy,
    results_to_dataframe,
    fill_from_reference,
    fill_from_unmethylated,
    filter_unchanged_windows,
    copy_energies_from_reference,
)
from .states import (
    total_states_index,
    total_open_states,
    get_states,
)
from .sequence_design import (
    position_map,
    widom_backbone,
    inject_cpgs,
)
