package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.primitives.Ints;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.Histogram;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class F1R2CountsCollector {

    public static final String ALT_TABLE_EXTENSION = ".alt_table";
    public static final String ALT_HIST_EXTENSION = ".alt_histogram";
    public static final String REF_HIST_EXTENSION = ".ref_histogram";

    private static final Logger logger = LogManager.getLogger(F1R2CountsCollector.class);

    private final Set<String> samples;

    private final SAMFileHeader header;

    private final CollectF1R2CountsArgumentCollection CF1R2Args;

    // For each sample and for reference context, count ref sites in a histogram keyed by depth
    private final Map<String, Map<String, Histogram<Integer>>> refSiteHistograms;

    // For each sample store the total depths of alt sites with alt depth = 1 separately to save memory
    private Map<String, DepthOneHistograms> depthOneAltHistograms;

    // alt table writer for each sample
    private Map<String, AltSiteRecord.AltSiteRecordTableWriter> altTableWriters;

    private final String outputBaseName;

    public F1R2CountsCollector(final CollectF1R2CountsArgumentCollection CF1R2Args, final SAMFileHeader header, final String outputBaseName) {
        this.CF1R2Args = CF1R2Args;
        samples = ReadUtils.getSamplesFromHeader(header);
        this.header = header;
        this.outputBaseName = outputBaseName;

        refSiteHistograms = new HashMap<>(samples.size());
        depthOneAltHistograms = new HashMap<>(samples.size());
        altTableWriters = new HashMap<>(samples.size());

        for (final String sample : samples) {
            final Map<String, Histogram<Integer>> refSitesHistogramsForThisSample = new HashMap<>(F1R2FilterConstants.ALL_KMERS.size());

            // Initialize for each reference the histogram of the counts of reference sites by depth
            F1R2FilterConstants.ALL_KMERS.forEach(context -> {
                Histogram<Integer> emptyRefHistogram = F1R2FilterUtils.createRefHistogram(context, CF1R2Args.maxDepth);
                refSitesHistogramsForThisSample.put(context, emptyRefHistogram);
            });

            refSiteHistograms.put(sample, refSitesHistogramsForThisSample);

            depthOneAltHistograms.put(sample, new DepthOneHistograms(CF1R2Args.maxDepth));

            // Intentionally not use try-with-resources so that the writer stays open outside of the try block
            final String altTablePath = outputBaseName + IOUtils.urlEncode(sample) + ALT_TABLE_EXTENSION;
            try {
                altTableWriters.put(sample, new AltSiteRecord.AltSiteRecordTableWriter(IOUtils.fileToPath(new File(altTablePath)), sample));
            } catch (IOException e) {
                throw new UserException(String.format("Encountered an IO exception creating a writer for %s", altTablePath), e);
            }
        }
    }

    public void process(final AlignmentContext alignmentContext, final ReferenceContext referenceContext) {
        final int position = referenceContext.getInterval().getStart();
        final String refContext = referenceContext.getKmerAround(position, F1R2FilterConstants.REF_CONTEXT_PADDING);
        if (refContext == null) {
            return;
        }
        final Nucleotide refBase = F1R2FilterUtils.getMiddleBase(refContext);

        if (refContext.contains("N") || refContext.length() != F1R2FilterConstants.REFERENCE_CONTEXT_SIZE) {
            return;
        }

        if (refContext == null) {
            logger.warn(String.format("Skipped a site with null reference at interval %s, k-mer = %s",
                    referenceContext.getInterval().toString(), refContext));
            return;
        }

        for (final Map.Entry<String, AlignmentContext> entry : alignmentContext.splitContextBySampleName(header).entrySet()) {
            final String sample = entry.getKey();
            final ReadPileup pileup = entry.getValue().getBasePileup().makeFilteredPileup(pe -> pe.getQual() > CF1R2Args.minBaseQuality);
            final int[] baseCounts = pileup.getBaseCounts();
            final int depth = (int) MathUtils.sum(baseCounts);

            if (!isPileupGood(pileup)) {
                return;
            }

            // Make a copy of base counts and update the counts of ref to -1. Now the maxElementIndex of the array gives us
            // the alt base.
            final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
            baseCountsCopy[refBase.ordinal()] = -1;
            final int altBaseIndex = MathUtils.maxElementIndex(baseCountsCopy);
            final boolean referenceSite = baseCounts[altBaseIndex] == 0;

            // If the site is ref, we simply update the coverage histogram
            if (referenceSite) {
                refSiteHistograms.get(sample).get(refContext).increment(Math.min(depth, CF1R2Args.maxDepth));
                return;
            }

            // If we got here, we have an alt site with a single alt base
            final Nucleotide altBase = Nucleotide.decode(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

            final int refCount = baseCounts[refBase.ordinal()];
            final int altCount = baseCounts[altBaseIndex];
            Utils.validate(altCount > 0, "We must have a nonzero alt read but got " + altCount);

            final int refF1R2 = pileup.getNumberOfElements(pe -> Nucleotide.decode(pe.getBase()) == refBase && ReadUtils.isF1R2(pe.getRead()));
            final int altF1R2 = pileup.getNumberOfElements(pe -> Nucleotide.decode(pe.getBase()) == altBase && ReadUtils.isF1R2(pe.getRead()));

            if (altCount == 1) {
                final ReadOrientation type = altF1R2 == 1 ? ReadOrientation.F1R2 : ReadOrientation.F2R1;
                depthOneAltHistograms.get(sample).increment(refContext, altBase, type, depth);
                return;
            }

            try {
                altTableWriters.get(sample).writeRecord(new AltSiteRecord(refContext, refCount, altCount, refF1R2, altF1R2, altBase));
            } catch (IOException e) {
                throw new UserException("Encountered an IO Exception writing to the alt data table", e);
            }
        }
    }

    public void writeHistograms() {
        for (final String sample : samples) {
            final MetricsFile<?, Integer> refMetricsFile = new MetricsFile<>();
            refMetricsFile.addHeader(new StringHeader(sample));
            refSiteHistograms.get(sample).values().forEach(refMetricsFile::addHistogram);
            refMetricsFile.write(new File(outputBaseName + IOUtils.urlEncode(sample) + REF_HIST_EXTENSION));

            final MetricsFile<?, Integer> altMetricsFile = new MetricsFile<>();
            altMetricsFile.addHeader(new StringHeader(sample));
            depthOneAltHistograms.get(sample).getHistograms().forEach(altMetricsFile::addHistogram);
            altMetricsFile.write(new File(outputBaseName + IOUtils.urlEncode(sample) + ALT_HIST_EXTENSION));
        }
    }

    public void close() {
        if (altTableWriters != null) {
            for (final AltSiteRecord.AltSiteRecordTableWriter writer : altTableWriters.values()) {
                if (writer != null) {
                    try {
                        writer.close();
                    } catch (IOException e) {
                        throw new UserException("Encountered an IO exception while closing the alt table writer", e);
                    }
                }
            }
        }
    }

    /**
     * Use a series of heuristics to detect a bad pileup.
     */
    private boolean isPileupGood(final ReadPileup pileup){
        final int[] baseCounts = pileup.getBaseCounts();
        final int depth = (int) MathUtils.sum(baseCounts);

        List<Integer> mappingQualities = Ints.asList(pileup.getMappingQuals());

        // If more than 1% of the reads is indel then consider this site an indel
        final int indelThreshold = depth/100;
        boolean isIndel = pileup.getNumberOfElements(pe -> pe.isDeletion() || pe.isAfterInsertion() || pe.isBeforeDeletionStart()) > indelThreshold;

        // If depth (the sum of base counts) is 0 but the pileup is non-empty, that means all the reads
        // have deleted bases at this particular locus
        isIndel = isIndel || depth == 0 && pileup.size() > 0;

        return depth > 0 && ! isIndel && MathUtils.median(mappingQualities) >= CF1R2Args.minMedianMapQual;
    }
}
