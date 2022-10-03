import altair_saver as saver
import altair as alt
import pandas as pd


class EvaluationPlots(object):
    def __init__(self, data, outdir):
        self.colors = ['#9AD6E5', '#66A6D9', '#BBB465', '#FFC12C', '#DA0034', '#3A7A9B', '#FFAD00']
        self.data = data
        self.outdir = outdir

    def ilp_opitvax_bar_plot(self):
        min_ = round(min(self.data['Population Coverage']), 2) - 0.01
        bar = alt.Chart(self.data).mark_bar().encode(
            x='Method',
            y=alt.Y('Population Coverage:Q', scale=alt.Scale(domain=[min_,1.0]), axis=alt.Axis(titleColor='#3a7a9b')),
            color=alt.Color('Method',
                            scale=alt.Scale(
                                domain=self.data['Method'].tolist(),
                                range=['#9ad6e5', '#3a7a9b', '#9ad6e5', '#3a7a9b']),
                            legend=None)
        )

        crosses = alt.Chart(self.data).mark_circle(size=60, color='#FFAD00').encode(
            x='Method',
            y=alt.Y('Number of Peptides:Q', axis=alt.Axis(titleColor='#FFAD00'))
        )

        layered = alt.layer(bar, crosses).resolve_scale(
            y='independent'
        ).configure_axis(
            labelFontSize=20,
            titleFontSize=20
        )

        layered.save(self.outdir + 'unlinked_pop_coverage_bar.html')

    def covered_locus_pie(self, loci, mhc_type, method):
        charts = []
        for locus in loci:
            df = pd.DataFrame.from_dict(self.data[locus])
            pie = alt.Chart(df).encode(
                theta=alt.Theta('values:Q', stack=True),
                color=alt.Color('allele_category:N', scale=alt.Scale(range=self.colors))
            ).properties(
                title='Coverage of alleles at ' + locus
            )

            c1 = pie.mark_arc(stroke='#fff')
            c2 = pie.mark_text(outerRadius=168).encode(
                text='values:Q'
            )

            chart = (c1 + c2)
            charts.append(chart)

        layered = alt.hconcat(*[chart for chart in charts])
        layered.configure_view(
            strokeWidth=0
        ).configure_axis(
            labelFontSize=20,
            titleFontSize=20
        ).save(self.outdir + 'pies_covered_alleles_' + mhc_type + '_' + method + '.html')

    def hit_covered_locus_pie(self, loci, mhc_type, method):
        charts = []
        for locus in loci:
            df = pd.DataFrame.from_dict(self.data[locus])
            order = sorted(df.allele_category[:-2], key=float) + ['No-BA-prediction', 'Uncovered']
            pie = alt.Chart(df).transform_calculate(
                order=f"-indexof({order}, datum.Origin)"
            ).encode(
                theta=alt.Theta('values:Q', stack=True),
                color=alt.Color('allele_category:N', scale=alt.Scale(scheme='category20b')),
                order='order:Q'
            ).properties(
                title='Coverage of alleles at ' + locus
            )

            chart = pie.mark_arc(stroke='#fff')
            charts.append(chart)

        layered = alt.hconcat(*[chart for chart in charts]).resolve_scale(color='independent')
        layered.configure_view(
            strokeWidth=0
        ).configure_axis(
            labelFontSize=20,
            titleFontSize=20
        ).save(self.outdir + 'pies_hit_covered_alleles_' + mhc_type + '_' + method + '.html')

    def stacked_hit_covered_alleles(self, mhc_key, method_key):
        order = ['No-BA-prediction', 'Uncovered'] + sorted([x for x in self.data['allele_category'] if x!='No-BA-prediction' and x!='Uncovered'], key=float, reverse=True)
        chart = alt.Chart(self.data).mark_bar().encode(
            x='method',
            y='sum(values)',
            color=alt.Color('allele_category',
                            scale=alt.Scale(scheme='viridis'),
                            sort=order,
                            legend=alt.Legend(columns=3, symbolLimit=0)),
            order='order:Q'
        )
        chart.show()
        chart.save(self.outdir + 'stacked_bar_hit_alleles_' + mhc_key + '_' + method_key + '.html')

    def hit_pies(self, loci, mhc_type, method):
        if mhc_type == 'mhc1':
            bins = int(max([self.data[loc]['hits'].max() for loc in loci])) + 1
        else:
            bins = 10
        for locus in loci:
            df = self.data[locus]
            pie = alt.Chart(df).encode(
                theta=alt.Theta('count', stack=True),
                color=alt.Color('hits', scale=alt.Scale(range=self.colors), bin=alt.Bin(maxbins=bins))
            ).properties(
                title='Coverage of alleles at ' + locus
            )

            # c1 = pie.mark_arc(stroke='#fff')
            c2 = pie.mark_text(outerRadius=168).encode(
                text='hits'
            )

            # chart = (c1 + c2)
            chart = pie.mark_arc() + c2
            chart.configure_view(
                strokeWidth=0
            ).configure_axis(
                labelFontSize=20,
                titleFontSize=20
            ).save(self.outdir + 'allele_hit_pie_' + locus + '_' + mhc_type + '_' + method + '.html')

    def hit_histogram(self, mhc, method):
        charts = []
        for population, data in self.data.items():
            df = pd.DataFrame.from_dict(data[0])
            hist = alt.Chart(df).mark_bar(color='#9AD6E5').encode(
                x=alt.X('count:Q',
                        title='Number of peptide-HLA hits',
                        bin=alt.Bin(maxbins=30)),
                y=alt.Y('prob:Q', title='Frequency')
            ).properties(
                title=population
            )

            exp = pd.DataFrame(data[1])
            line = alt.Chart(exp).mark_rule(size=4, color='#FFAD00').encode(
                x='expected:Q'
            )

            text = line.mark_text(
                align='left',
                baseline='middle',
                dx=7,
                color='#FFAD00'
            ).encode(
                text=alt.Text('expected', format=',.3f')
            )

            chart = (hist + line + text)
            charts.append(chart)

        layered = alt.hconcat(*[chart for chart in charts])
        layered.configure_view(
            strokeWidth=0
        ).configure_axis(
            labelFontSize=20,
            titleFontSize=20
        ).save(self.outdir + 'hit_histogram_' + mhc + '_' + method + '.html')

    def population_hit_bar_plot(self, mhc, method):
        sorting = sorted(list(set([p for p in self.data['population'] if not p == 'Average']))) + ['Average']
        bar = alt.Chart(self.data).mark_bar().encode(
            x=alt.X('population', axis=alt.Axis(title=None, labels=False, ticks=False), sort=sorting),
            y=alt.Y('prob:Q', title='Population coverage'),
            color=alt.Color('population', scale=alt.Scale(range=self.colors)),
            column=alt.Column('n', title='Minimum number of peptide-Hla hits')
        ).configure_axis(
            grid=False,
            labelFontSize=20,
            titleFontSize=20
        ).configure_view(
            strokeOpacity=0
        )

        bar.save(self.outdir + 'min_hit_number_bar_' + mhc + '_' + method + '.html')

    def composition_plot(self, mhc, method):
        bar = alt.Chart(self.data).mark_bar().encode(
            x=alt.X('Protein:N', axis=alt.Axis(labelAngle=0)),
            y='Count',
            color=alt.Color('Protein', scale=alt.Scale(range=self.colors))
        ).properties(
            width=230
        ).configure_axis(
            labelFontSize=20,
            titleFontSize=20
        )
        bar.save(self.outdir + 'composition_bar_' + mhc + '_' + method + '.html')