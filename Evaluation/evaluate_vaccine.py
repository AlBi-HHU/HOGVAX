import gzip
import re
import numpy as np
import pandas as pd
import read_peptides, binarize_entries


class EvaluationMetrices(object):
    def __init__(self, peptides, binding_affinities, allele_freq):
        self.peptides = peptides
        self.binding_affinities = binding_affinities
        self.allele_freq = allele_freq

    def pred_at_least_one(self, population):
        """
        Calculate the probability that any person of the given population presents at least one peptide
        (based on Gifford EvalVax-Unlinked)
        """
        self.population = population
        loci_prob = 1
        for locus in self.binding_affinities.columns.levels[0]:
            loci_prob *= (1 - self.locus_prob(locus))
        loci_prob = 1 - loci_prob
        # print('Total prob', loci_prob)
        return loci_prob

    def locus_prob(self, locus):
        """ based on Gifford """
        # Equation 5
        alleles_ba = self.binding_affinities[locus]
        alleles_prob = 1 - np.prod(1-alleles_ba, axis=0)
        # Equation 6
        diploid_ba = 1 - np.outer((1-alleles_prob), (1-alleles_prob))
        np.fill_diagonal(diploid_ba, alleles_prob)
        diploid_ba = pd.DataFrame(diploid_ba, index=alleles_prob.index, columns=alleles_prob.index)
        # Equation 7
        diploid_af = np.outer(self.allele_freq[locus].loc[self.population], self.allele_freq[locus].loc[self.population])
        diploid_af = pd.DataFrame(diploid_af, index=self.allele_freq[locus].columns, columns=self.allele_freq[locus].columns)
        # sum up the multiplied values for each allele pair (done by np.sum automatically)
        prob = np.sum(np.sum(diploid_af * diploid_ba))
        # print(prob)
        return prob

    def get_all_genotypes(self, genotypes):
        # generate each genotype only once
        indices = genotypes.index
        for i in range(len(indices)):
            hap1 = genotypes.index[i]
            for j in range(0, len(indices)):#range(i, len(indices)):
                hap2 = genotypes.index[j]
                yield hap1, hap2

    def pred_genotype_hits(self, population, single_allele_ba, max_interest, mhc_type):
        """
        Calculate the probability that an individual of a specific population has exactly n peptide hla hits
        (based on Gifford EvalVax-Robust)
        """
        # Equation 1
        genotype_af = np.outer(self.allele_freq.loc[population], self.allele_freq.loc[population])
        genotype_af = pd.DataFrame(genotype_af, index=self.allele_freq.columns, columns=self.allele_freq.columns)
        # compute number of peptides binding to each allele e(a)
        hits_per_allele = np.sum(single_allele_ba)
        # Equation 2
        hits_per_genotype = {}
        for hap1, hap2 in self.get_all_genotypes(genotype_af):
            unique_alleles = set(hap1).union(set(hap2))
            count, no_predictions = self.get_genotype_counts(unique_alleles, hits_per_allele, mhc_type)
            if count in hits_per_genotype:
                hits_per_genotype[count].append((hap1, hap2))
            else:
                hits_per_genotype[count] = [(hap1, hap2)]
        # Equation3 count exact hits
        m = int(max(hits_per_genotype))
        probs = {'count': [], 'prob': []}
        for i in range(m + 1):
            if i in hits_per_genotype:
                frq = np.sum(genotype_af[h1][h2] for h1, h2 in hits_per_genotype[i])
                probs['count'].append(i)
                probs['prob'].append(frq)
            else:
                probs['count'].append(i)
                probs['prob'].append(0)
        # calculate expected value
        expected_value = np.sum(np.array(probs['count']) * np.array(probs['prob']))
        # Equation 4 for k in [1, max_interest] filter data frame for all counts >= k
        cur_df = pd.DataFrame.from_dict(probs)
        prob_n = {'n': [], 'prob': []}
        for n in range(1, max_interest + 1):
            prob_n['n'].append(n)
            # have to use converse probability
            prob_n['prob'].append(1-sum(cur_df[cur_df['count'] < n]['prob']))
        return probs, prob_n, expected_value, no_predictions

    def get_genotype_counts(self, alleles, hits_per_allele, mhc_type):
        count = 0
        no_ba_prediction = []
        for allele in alleles:
            if mhc_type == 'mhc1':
                reg = re.search(r'(HLA-[A-Z]+)(\d)', allele)
                locus = reg.group(1)
            else:
                reg = re.search(r'HLA-DP|HLA-DQ|DRB1', allele)
                locus = reg.group(0)
            if (locus, allele) in hits_per_allele:
                count += hits_per_allele[locus, allele]
            else:
                no_ba_prediction.append((locus, allele))
        return int(count), no_ba_prediction

    def number_covered_alleles(self, population, locus):
        """
        Calculate the number of covered alleles of the population and with this metric the population coverage;
        Track information about alleles for which we don't have prediction data from netmhc.
        """
        # binding affinities only contain peptides from vaccine
        alleles_ba = self.binding_affinities[locus]
        alleles_fa = self.allele_freq[locus]

        # count hits per allele
        tmp = self.binding_affinities[locus].sum().to_frame().reset_index()
        tmp = tmp.rename(columns={'genotype':'count', 0:'hits'})
        allele_hit_counts = tmp.groupby('hits').count().reset_index()
        fractions = allele_hit_counts['count'] / sum(allele_hit_counts['count'])
        allele_hit_counts = allele_hit_counts.assign(frac=fractions)

        # check for all alleles with (entry >= 0.638) / non-zero entry -> sum up corresponding frequencies
        hit_allele_freq = {}
        filter = (alleles_ba != 0).any()
        covered_alleles = filter.index[filter]
        covered_fraction = 0
        for allele in covered_alleles:
            if allele in alleles_fa:
                covered_fraction += alleles_fa[allele][population]
                hit = tmp[tmp['count'] == allele]['hits'].values[0]
                if hit in hit_allele_freq:
                    hit_allele_freq[hit] += alleles_fa[allele][population]
                else:
                    hit_allele_freq[hit] = alleles_fa[allele][population]

        no_predictions = []
        fraction_no_predictions = 0
        for allele in alleles_fa:
            if allele == 'unknown':
                no_predictions.append(allele)
                fraction_no_predictions += alleles_fa[allele][population]
            if allele not in alleles_ba:
                no_predictions.append(allele)
                fraction_no_predictions += alleles_fa[allele][population]

        total = sum(alleles_fa.loc[population])
        return covered_fraction, no_predictions, fraction_no_predictions, total, allele_hit_counts, hit_allele_freq
