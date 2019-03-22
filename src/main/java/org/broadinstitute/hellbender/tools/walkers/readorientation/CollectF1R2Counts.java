package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * At each genomic locus, count the number of F1R2/F2R1 alt reads.
 * {@link LearnReadOrientationModel} uses the tsv output of this tool
 *
 * <h3>Usage Example</h3>
 *
 * gatk CollectF1R2Counts \
 *   -R GRCh38.fasta \
 *   -I tumor.bam \
 *   -O tumor-artifact-prior-table.tsv \
 *   -alt-table tumor-alt.tsv \
 *   -ref-hist tumor-ref.metrics \
 *   -alt-hist tumor-alt.metrics
 */

@CommandLineProgramProperties(
        summary = "Collect F1R2 read counts for the Mutect2 orientation bias mixture model filter",
        oneLineSummary = "Collect F1R2 read counts for the Mutect2 orientation bias mixture model filter",
        programGroup = CoverageAnalysisProgramGroup.class
)

public class CollectF1R2Counts extends LocusWalker {
    @ArgumentCollection
    protected CollectF1R2CountsArgumentCollection CF1R2Args = new CollectF1R2CountsArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Base path for output files")
    private String outputBaseName;

    private F1R2CountsCollector f1R2CountsCollector;

    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public void onTraversalStart() {
        f1R2CountsCollector = new F1R2CountsCollector(CF1R2Args, getHeaderForReads(), outputBaseName);
    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        f1R2CountsCollector.process(alignmentContext, referenceContext);
    }

    @Override
    public Object onTraversalSuccess() {
        f1R2CountsCollector.writeHistograms();
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        f1R2CountsCollector.close();
    }
}
