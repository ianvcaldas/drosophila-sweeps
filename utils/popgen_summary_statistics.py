"""Workers calculate a certain set of features from one simulated dataset.

They take as input the binary variants matrix and a list of genome positions.
As output, they provide one value per feature.

Very important: every single feature here assumes that the data contains ONLY
segregating sites, no monomorphic sites. Besides, all loci have two alleles,
coded as 0 and 1. Any violation of these assumptions will produce unexpected
results.
"""

from collections import Counter, defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.spatial.distance import pdist


class Worker:
    def __init__(self):
        # Used for setting parameters
        pass

    def run(self, positions, data, window_sizes):
        """Perform the task.
        Should return a dict where the keys are each feature and the values are
        the results of calculating that feature.
        """
        pass

    def __str__(self):
        return self.name


class WindowMeasureBase(Worker):
    """Parent base class for different measures across windows.
    Each output is the statistic calculated in a subwindow of specified size,
    emanating from a specified position.
    """

    def run_bp(self, positions, data, center_pos, window_sizes, last_position=1):
        self.prepare(positions, data)
        results = {ws: {} for ws in window_sizes}
        for ws in window_sizes:
            for center in center_pos[ws]:
                # This allows workers to check the window size they're working
                # with without passing it as an explicit argument to
                # self.calculate().
                # We are measuring the EFFECTIVE window size, that is, the size
                # of the window that remains within the locus boundaries.
                self.curr_eff_window_size = min(
                    ws,  # window fully within locus
                    0.5 * ws + center,  # touches start of locus
                    0.5 * ws + (last_position - center),  # touches end of locus
                )
                subpositions, subwindow = self._make_window(positions, data, center, ws)
                results[ws][center] = self.calculate(subpositions, subwindow)
        return results

    def run_snp(self, positions, data, window_size, window_step):
        self.prepare(positions, data)
        results = defaultdict(dict)
        for i in range(0, len(data), window_step):
            subpositions = positions[i : (i + window_size)]
            subwindow = data[i : (i + window_size), :]
            center = subpositions[len(subpositions) // 2]
            # There's no way for a SNP window to extend beyond the locus
            ws = subpositions[-1] - subpositions[0]
            self.curr_eff_window_size = ws
            results[ws][center] = self.calculate(subpositions, subwindow)
        return results

    def prepare(self, positions, data):
        """Do any preparations necessary before things are run on multiple
        subwindows and with multiple central positions. Can be used to set up
        calculations so they don't have to be repeated every window size/center
        position.
        """
        pass

    def _make_window(self, positions, data, center, window_size):
        """Create subwindow emanating from a specified position.
        Window sizes are given in relative locus size, from 0 to 1.
        """
        upper = center + window_size / 2
        lower = center - window_size / 2
        indices = np.where((positions > lower) & (positions < upper))[0]
        subpositions = positions[indices]
        subwindow = data[indices, :]
        return subpositions, subwindow

    def calculate(self, positions, data):
        """Takes in haplotypes as a numpy array with [position, samples].
        Positions are a list of physical genomic locations.
        This function calculates a single feature.
        """
        pass


class HaploCounter(WindowMeasureBase):
    """For statistics that deal with the number and frequency of different
    haplotypes.
    """

    def get_haplotype_counts(self, haps):
        haps_strings = [hap.tostring() for hap in haps.T]
        counts = Counter(haps_strings)
        return counts


class NumberOfSNPs(WindowMeasureBase):
    """The number of variant sites in the window."""

    def __init__(self):
        self.name = "num_snps"

    def calculate(self, positions, data):
        return data.shape[0]


class Pi(WindowMeasureBase):
    """Average nucleotide diversity.

    Defined as the probability that, if you pick a random pair of haplotypes
    from the sample, a random site will be different between them.

    Defined as the average number of pairwise differences between samples. This
    is with option fraction=False.

    If fraction=True, we divide the value of Pi by the window size. If the
    window size is given in bp, this produces an intuitive interpretation: the
    probability that, if you pick two sequences from the population, any given
    site differs between them. If the window size is given in chromosome
    fraction, you need to divide the resulting value by the entire locus size
    to recover the intuitive interpretation.

    Cite: Nei and Li (1979). Mathematical model for studying genetic variation
    in terms of restriction endonucleases. Proc Natl Acad Sci USA 76 (10),
    pp. 5269–5273.
    """

    def __init__(self, fraction=False, with_replacement=True):
        self.name = "pi"
        self.fraction = fraction
        self.with_replacement = with_replacement

    def calculate(self, positions, data):
        if self.with_replacement:
            return self._calc_with_replacement(positions, data)
        else:
            return self._calc_without_replacement(positions, data)

    def _calc_with_replacement(self, positions, haps):
        num_seg_sites = haps.shape[0]
        if num_seg_sites == 0:
            return 0
        diff_sites = pdist(haps.transpose(), metric="hamming") * num_seg_sites
        num_comp = len(diff_sites)
        pi = np.sum(diff_sites) / (num_comp)
        if np.isnan(pi):
            pi = 0
        # Only perform fraction if window size is in base pairs, not if it's a fraction
        # of the genome.
        if self.fraction and self.curr_eff_window_size > 1:
            pi = pi / self.curr_eff_window_size
        return pi

    def _calc_without_replacement(self, positions, window):
        raise NotImplementedError


class NumberOfHaplotypes(HaploCounter):
    """Counts how many different haplotypes are present in the window.

    Expected to increase with window size and sample size; this makes it hard
    to compare.
    """

    def __init__(self):
        self.name = "num_haps"

    def calculate(self, positions, data):
        counts = self.get_haplotype_counts(data)
        return len(counts)


class WattersonsTheta(WindowMeasureBase):
    """Measures the number of segregating sites. Since we expect the number
    of segregating sites to increase with the sample size, Watterson's theta
    has a correction for sample size.

    Cite: Watterson, G. (1975). On the number of segregating sites in genetical
    models without recombination. Theoretical Population Biology 7 (2), pp.
    256–276.
    """

    def __init__(self):
        self.name = "theta_w"
        self.correction = None

    def calculate(self, positions, data):
        num_seg_sites = data.shape[0]
        num_samples = data.shape[1]
        if num_seg_sites == 0:
            return 0
        if self.correction is None:
            self.correction = np.sum([1 / k for k in range(1, num_samples)])
        theta_hat = num_seg_sites / self.correction
        return theta_hat


class HaplotypeHomozygosity(HaploCounter):
    """Also known as H1.

    The probability that two haplotypes, taken randomly from the population,
    are identical.
    """

    def __init__(self):
        self.name = "H1"

    def calculate(self, positions, data):
        num_samples = data.shape[1]
        if data.shape[0] == 0:  # 0 segregating sites
            return 1
        counts = self.get_haplotype_counts(data)
        return np.sum([(count / num_samples) ** 2 for count in counts.values()])


class H12(HaploCounter):
    """The haplotype homozygosity, only the two most common haplotypes are
    treated as one. Useful for looking for genetic sweeps.

    If there is only one fixed haplotype, H12 = haplotype homozygosity = 1.

    Cite: Garud, N. R. et al. (2015). Recent Selective Sweeps in North American
    Drosophila melanogaster Show Signatures of Soft Sweeps. PLOS Genetics 11
    (2), pp. 1–32.
    """

    def __init__(self):
        self.name = "H12"

    def calculate(self, positions, data):
        num_samples = data.shape[1]
        if data.shape[0] <= 1:  # 0 or 1 segregating sites
            return 1
        counts = self.get_haplotype_counts(data)
        top_two_freqs = [count[1] / num_samples for count in counts.most_common(2)]
        h1 = np.sum([(count / num_samples) ** 2 for count in counts.values()])
        p1 = top_two_freqs[0]
        try:
            p2 = top_two_freqs[1]
        except IndexError:
            # only one haplotype, fixed
            p2 = 0
        return h1 + 2 * p1 * p2


class H2byH1(HaploCounter):
    """The haplotype homozygosity, but ignoring the most common haplotype. In
    conjunction with H12, can be used to distinguish between soft and hard
    selective sweeps.

    Cite: Garud, N. R. et al. (2015). Recent Selective Sweeps in North American
    Drosophila melanogaster Show Signatures of Soft Sweeps. PLOS Genetics 11
    (2), pp. 1–32.
    """

    def __init__(self):
        self.name = "H2/H1"

    def calculate(self, positions, data):
        num_samples = data.shape[1]
        if data.shape[0] == 0:  # 0 segregating sites
            return np.nan
        counts = self.get_haplotype_counts(data)
        try:
            top_count = counts.most_common(1)[0][1]
        except IndexError:
            print(data, positions)
        top_freq = top_count / num_samples
        h1 = np.sum([(count / num_samples) ** 2 for count in counts.values()])
        h2 = h1 - top_freq ** 2
        return h2 / h1


class TajimasD(WindowMeasureBase):
    """A statistic to measure signals of selection.

    Cite: Tajima, F. (1989). Statistical method for testing the neutral
    mutation hypothesis by DNA poly- morphism. Genetics 123 (3), pp. 585–595.
    """

    def __init__(self):
        self.name = "taj_D"
        self.pi_calc = Pi(fraction=False)
        self.a1 = None
        self.e1 = None
        self.e2 = None

    def _prepare_values(self, n):
        """Get auxiliary values in the computation of D.

        The notation here is identical to Tajima's original paper.
        n is the number of samples.
        """
        a1 = np.sum([1 / k for k in range(1, n)])
        a2 = np.sum([1 / (k ** 2) for k in range(1, n)])
        b1 = (n + 1) / (3 * (n - 1))
        b2 = (2 * (n ** 2 + n + 3)) / (9 * n * (n - 1))
        c1 = b1 - 1 / a1
        c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 ** 2)
        e1 = c1 / a1
        e2 = c2 / (a1 ** 2 + a2)
        self.a1 = a1
        self.e1 = e1
        self.e2 = e2

    def calculate(self, positions, data):
        num_samples = data.shape[1]
        if num_samples <= 1:
            return np.nan
        if self.a1 is None:
            self._prepare_values(num_samples)
        S = data.shape[0]  # number of segregating sites
        if S == 0:
            return np.nan
        a1 = self.a1
        e1 = self.e1
        e2 = self.e2
        big_pi = self.pi_calc.calculate(positions, data)
        var_D = np.sqrt(e1 * S + e2 * S * (S - 1))
        D = (big_pi - S / a1) / var_D
        return D


class LinkageDisequilibrium(WindowMeasureBase):
    """
    Calculates measures of linkage disequilibrium between pairs of loci and
    summarizes them across pairs.
    The pairs of loci for which r^2 is calculated is restricted to pairs whose
    loci are at a certain minimum distance, defined as a fraction of the
    current window size.
    """

    def __init__(self, min_dist=0):
        """
        min_dist: the fraction of the current window size that serves as the
        minimum distance between loci. For instance, if min_dist = 1/2, then
        only pairs of loci at least window_size/2 apart will have r2 calculated
        for them.
        """
        self.min_dist = min_dist

    def prepare(self, positions, data):
        self._calculated_pairs = dict()

    def _ld_measure(self, locus_A, locus_B, kind="r2"):
        """
        Calculate D or r^2 between the two loci A and B.
        """
        assert kind in ["r2", "D"]
        a1 = np.sum(locus_A) / len(locus_A)  # Frequency of allele 1 in locus A
        b1 = np.sum(locus_B) / len(locus_B)  # Frequency of allele 1 in locus B
        # This probably is the speed bottleneck:
        count_a1b1 = sum(hap == (1, 1) for hap in zip(locus_A, locus_B))
        a1b1 = count_a1b1 / len(
            locus_A
        )  # Frequency of haplotype 11 between the two loci.
        D = a1b1 - a1 * b1
        if kind == "D":
            return D
        r2 = (D ** 2) / (a1 * (1 - a1) * b1 * (1 - b1))
        return r2

    def calculate(self, positions, data):
        num_sites = data.shape[0]
        min_abs_dist = self.curr_eff_window_size * self.min_dist
        if num_sites == 0:
            return 1  # complete non-random assortment
        pairs = [
            pair
            for pair in combinations(range(num_sites), 2)
            if abs(positions[pair[0]] - positions[pair[1]]) > min_abs_dist
        ]
        sum_r_squared = 0
        for pair in pairs:
            pos = (positions[pair[0]], positions[pair[1]])
            if pos in self._calculated_pairs:
                sum_r_squared += self._calculated_pairs[pos]
            else:
                sum_r_squared += self._ld_measure(data[pair[0], :], data[pair[1], :])
                self._calculated_pairs[pos] = sum_r_squared
        return sum_r_squared / len(pairs)


class Zns(LinkageDisequilibrium):
    """
    Calculate the Zns statistic of Kelly (1997). It is just the average r2
    among every pair of loci in the genomic window.

    Cite:
    Kelly, J. K. (1997). A Test of Neutrality Based on Interlocus Associations.
    Genetics 146 (3), pp. 1197– 1206. url:
    http://www.genetics.org/content/146/3/1197.
    """

    def __init__(self):
        self.min_dist = 0


class SiteFrequencySpectrum(Worker):
    """Calculates the SFS of the data. Corresponds to a graph where the x-axis is
    the frequency of the alternate allele and the y-axis is the number of positions
    in the data with that frequency.
    """

    def __init__(self, folded=False, relative_freq=False):
        """If folded=True, return the folded SFS.
        If relative_freq=True, make the x-axis the proportion of sites instead
        of the number of sites.
        """
        self.name = "sfs"
        self.folded = folded
        self.relative_freq = relative_freq

    def run(self, positions, data, window_sizes):
        num_samples = data.shape[1]
        freqs = np.bincount(data.sum(axis=1, dtype=int), minlength=num_samples + 1)
        freqs[0] = 0  # constant sites are not implemented
        if self.folded:
            rev_freqs = np.flip(freqs, axis=0)
            freqs = freqs + rev_freqs
            midpoint = int(np.ceil(len(freqs) / 2))
            freqs = freqs[:midpoint]
        if self.relative_freq:
            num_sites = data.shape[0]
            freqs = freqs / num_sites
        digits = len(str(len(freqs)))
        return {f"sfs_{str(i).zfill(digits)}": freqs[i] for i in range(len(freqs))}


class FourGameteTest(WindowMeasureBase):
    """
    Does the four gamete test on the outermost SNPs of a window. Under an
    infinite sites model, between any pair of loci that have all four gametes
    --- 00, 01, 10, and 11 --- there must have occurred at least one
    recombination event. This statistic just returns 1 if all gametes are
    present and 0 if they aren't.
    """

    def __init__(self):
        self.name = "four_gamete_test"

    def calculate(self, positions, data):
        if data.shape[0] <= 1:  # 0 or 1 segregating sites
            return 0
        locus_A = data[0, :]
        locus_B = data[-1, :]
        num_gametes = len(Counter(zip(locus_A, locus_B)))
        return int(num_gametes == 4)


class OutermostLD(LinkageDisequilibrium):
    """
    Returns measure of LD (D or r2) between the two outermost SNPs within a
    window.
    """

    def __init__(self, kind):
        assert kind in ["D", "r2"]
        self.name = f"Outermost {kind}"
        self.kind = kind

    def calculate(self, positions, data):
        if data.shape[0] <= 1:  # 0 or 1 segregating sites
            return np.nan
        locus_A = data[0, :]
        locus_B = data[-1, :]
        return self._ld_measure(locus_A, locus_B, self.kind)


class HaplotypeLDBound(HaploCounter):
    """
    Uses a local estimator of the minimum number of recombination events in a
    window.

    Cite:
    Myers, S. R. e R. C. Griffiths (2003). Bounds on the Minimum Number of
    Recombination Events in a Sample History. Genetics 163 (1), pp. 375–394.
    url: http://www.genetics.org/content/163/1/375.
    """

    def __init__(self):
        self.name = "haplotype_ld_bound"

    def calculate(self, positions, data):
        H = len(self.get_haplotype_counts(data))
        S = data.shape[0]
        if S <= 1:
            return np.nan
        return H - S - 1


class IBDMoment(WindowMeasureBase):
    """Calculates the mean, variance, skew, kurtosis and other standardized
    moments of the distribution of IBD segments within the provided window.

    This implementation assumes that the number of samples in the data is
    bigger than 2, and an even number. This should be true for any true
    application."""

    def __init__(self, moment=1):
        names = {1: "mean", 2: "var", 3: "skew", 4: "kurtosis"}
        try:
            this_name = names[moment]
        except KeyError:
            this_name = f"std_moment{moment}"
        self.name = f"IBD-{this_name}"
        self.moment = moment

    def calculate(self, positions, data):
        num_samples = data.shape[1]
        ibd_distr = []
        for i in range(0, num_samples, 2):
            j = i + 1
            # If there's an odd sample size, skip the last sample that doesn't
            # have a pair
            if j == num_samples:
                continue
            heteroz_positions = positions[data[:, i] != data[:, j]]
            heteroz_tracts = np.diff(heteroz_positions)
            ibd_distr.extend(heteroz_tracts)
        # If there are no SNPs in window for pair of genomes, the tract is
        # as big as the window. It could be bigger still, but we have no
        # information outside the window.
        if ibd_distr == []:
            ibd_distr = [self.curr_eff_window_size, self.curr_eff_window_size]
        if len(ibd_distr) == 1:
            ibd_distr = [ibd_distr[0], ibd_distr[0]]
        desc = stats.describe(ibd_distr)
        if self.moment == 1:
            return desc.mean
        elif self.moment == 2:
            return desc.variance
        elif self.moment == 3:
            return desc.skewness
        elif self.moment == 4:
            return desc.kurtosis
        else:
            # If moment is weird, return the standardized moment.
            return stats.moment(ibd_distr, moment=self.moment) / (
                np.std(ibd_distr) ** self.moment
            )
