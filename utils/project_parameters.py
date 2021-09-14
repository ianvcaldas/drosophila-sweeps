"""Parameters of the current project and analysis."""

from .popgen_summary_statistics import (
    NumberOfSNPs,
    Pi,
    NumberOfHaplotypes,
    HaplotypeHomozygosity,
    H12,
    H2byH1,
    TajimasD,
)

default_summary_statistics = [
    NumberOfSNPs(),
    Pi(fraction=True),
    NumberOfHaplotypes(),
    HaplotypeHomozygosity(),
    H12(),
    H2byH1(),
    TajimasD(),
]

summary_statistic_order = ["pi", "num_snps", "num_haps", "H1", "H12", "H2/H1", "taj_D"]

locus_size = 1_000_000
data_dimension = 21
smallest_window = 1000

default_simulation_parameters = {
    # Biology
    "locus-size": locus_size,
    "sample-size": 205,
    # SLiM options
    "last-generation": 5000,
    "simplification-interval": 500,
    "restarts-limit": 100,
    # Feature calculation and normalizing
    "data-dimension": data_dimension,
    "smallest-window": smallest_window,
}
