# === General Defaults ===
NUM_SPACERS = 4
# The default maximum length of any spacer created
MAX_SPACER_LENGTH = 12
ABSOLUTE_MAX_SPACER_LENGTH = 50
# The default size of the heterogeneity region to proceed the binding region of
# the primer
NUM_HETERO = 12
ABSOLUTE_MAX_NUM_HETERO = 50

# === Random sampling defaults ===
# The number of forward and reverse primers to compare against each other.
NUM_PAIRINGS_TO_COMP = 35
# The number of primers in the initial sets of potential primers.
INITIAL_PRIMER_SET_SIZE = 1000

# === SpacerAlignmentGen Specific defaults ===
# Begin default criteria weighting.
GET_SMALLEST_TOTAL_LEN_DEFAULT = 1
# We make this 4, the corresponding method is likely to return values a fourth
# of the size of get_smallest_total_len_list().
GET_SMALLEST_OF_ANY_SPACER_DEFAULT = 4

