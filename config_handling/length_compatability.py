from typing import List
from config_handling.formatting import InclRange, AdapterPair

PRIMER3_MAX_LEN = 60


def potential_length_overflow(adapters: List[AdapterPair],
                              hetero_len: int, binding_len: InclRange) -> bool:
    """Returns whether a primer produced by this program might exceed 60bp in
    length."""
    adapter_lens = []
    for adapter_pair in adapters:
        adapter_lens.append(len(adapter_pair.forward))
        adapter_lens.append(len(adapter_pair.reverse))

    max_adapter_len = max(adapter_lens)
    max_len = max_adapter_len + hetero_len + binding_len.stop

    if max_len > PRIMER3_MAX_LEN:
        return max_len
    else:
        return False
