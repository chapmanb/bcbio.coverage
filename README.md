# bcbio.coverage

Summarize coverage of high throughput sequencing experiments, emphasizing
approaches that identify poorly covered genes or pathways of biological
interest. Current variant calling and prioritization approaches identify
potentially deleterious mutations, but do not provide an easy way to summarize
regions with low or no coverage where we may have missed variants due to
coverage issues.

The inputs are:

- BigWig file of coverage or BAM file from which we can calculate coverage.
- BED file or list of gene identifiers of genes of interest.

and it outputs:

- Ranked list of gene regions with low or no coverage, prioritized by: coding
  region coverage and total bases missed.

Longer term goals include handling non-coding region prioritization.

## Resources

- Average coverage data for exomes from the [NHLBI Exome Sequencing Project (ESP)][esp]
  converted into BigWig format: [ESP6500SI-V2-coverage.bw][esp-bw].

[esp]: http://evs.gs.washington.edu/EVS/
[esp-bw]: https://s3.amazonaws.com/biodata/coverage/ESP6500SI-V2-coverage.bw

## License

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html
