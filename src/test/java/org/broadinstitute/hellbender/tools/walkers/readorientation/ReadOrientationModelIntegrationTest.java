package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class ReadOrientationModelIntegrationTest extends CommandLineProgramTest {
    public static final String testDir = toolsTestDir + "read_orientation_filter/";
    public static final String hapmapBamSnippet = testDir + "hapmap-20-plex-chr-20-21-read-orientation.bam";
    public static final String intervalList = testDir + "hapmap-20-plex-chr-20-21-read-orientation.interval_list";

    @DataProvider(name = "scatterCounts")
    public Object[][] getScatterCounts(){
        return new Object[][]{{1}, {5}};
    }
    /**
     * Test that the tool sites of orientation bias that are manually picked out.
     * Also tests scattering CollectF1R2Counts
     */
    @Test(dataProvider = "scatterCounts")
    public void testOnRealBam(final int scatterCount) throws IOException {
        final File scatteredDir = createTempDir("scattered");

        // Step 1: SplitIntervals
        final File intervalDir = createTempDir("intervals");
        runCommandLine(
            Arrays.asList(
                    "-R", b37_reference_20_21,
                    "-L", intervalList,
                    "-O", intervalDir.getAbsolutePath(),
                    "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount)),
            SplitIntervals.class.getSimpleName());

        // Step 2: CollectF1R2Counts
        final File[] intervals = intervalDir.listFiles();
        for (int i = 0; i < intervals.length; i++){
            runCommandLine(Arrays.asList(
                    "-R", b37_reference_20_21,
                    "-I", hapmapBamSnippet,
                    "-L", intervals[i].getAbsolutePath(),
                    "-O", scatteredDir.getAbsolutePath() + "/scatter_" + i),
                    CollectF1R2Counts.class.getSimpleName());

            final int iFinal = i;
            final File refHist = Arrays.stream(scatteredDir.listFiles())
                    .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.REF_HIST_EXTENSION))
                    .filter(file -> file.getAbsolutePath().contains("scatter_" + iFinal ))
                    .findFirst().get();
            // Ensure that we print every bin, even when the count is 0
            final int lineCount = (int) Files.lines(Paths.get(refHist.getAbsolutePath())).filter(l -> l.matches("^[0-9].+")).count();
            Assert.assertEquals(lineCount, F1R2FilterConstants.DEFAULT_MAX_DEPTH);
        }

        // Step 3: LearnReadOrientationModel
        final File priorTable = createTempFile("prior", ".tsv");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, priorTable.getAbsolutePath());

        final List<File> refMetricsFiles = Arrays.stream(scatteredDir.listFiles())
                .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.REF_HIST_EXTENSION)).collect(Collectors.toList());

        final List<File> altMetricsFiles = Arrays.stream(scatteredDir.listFiles())
                .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.ALT_HIST_EXTENSION)).collect(Collectors.toList());

        final List<File> altTableFiles = Arrays.stream(scatteredDir.listFiles())
                .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.ALT_TABLE_EXTENSION)).collect(Collectors.toList());

        refMetricsFiles.forEach(f ->
                args.addArgument(CollectF1R2CountsArgumentCollection.REF_SITE_METRICS_LONG_NAME, f.getAbsolutePath()));
        altMetricsFiles.forEach(f ->
                args.addArgument(CollectF1R2CountsArgumentCollection.ALT_DEPTH1_HISTOGRAM_LONG_NAME, f.getAbsolutePath()));
        altTableFiles.forEach(f ->
                args.addArgument(CollectF1R2CountsArgumentCollection.ALT_DATA_TABLE_LONG_NAME, f.getAbsolutePath()));
        runCommandLine(args.getArgsList(), LearnReadOrientationModel.class.getSimpleName());

        final ArtifactPriorCollection artifactPriorCollection = ArtifactPriorCollection.readArtifactPriors(priorTable);

        // Step 4: Mutect 2
        final File unfilteredVcf = GATKBaseTest.createTempFile("unfiltered", ".vcf");
        final File filteredVcf = GATKBaseTest.createTempFile("filtered", ".vcf");
        final File bamout = GATKBaseTest.createTempFile("SM-CEMAH", ".bam");

        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                    "-I", hapmapBamSnippet,
                    "-R", b37_reference_20_21,
                    "-O", unfilteredVcf.getAbsolutePath(),
                    "-bamout", bamout.getAbsolutePath()),
                Mutect2.class.getSimpleName()));

        new Main().instanceMain(makeCommandLineArgs(
                Arrays.asList(
                        "-V", unfilteredVcf.getAbsolutePath(),
                        "-R", b37_reference_20_21,
                        "--" + M2FiltersArgumentCollection.ARTIFACT_PRIOR_TABLE_NAME, priorTable.getAbsolutePath(),
                        "-O", filteredVcf.getAbsolutePath()),
                FilterMutectCalls.class.getSimpleName()));


        // These artifacts have been verified manually
        // The pair is of type (Position, Expected Source of Prior Probability)
        // Prior for e.g. TGA->A F2R1 should come from TCA->T F1R2
        final List<Triple<Integer, ReadOrientation, ArtifactState>> knownArtifacts = Arrays.asList(
                new ImmutableTriple<>(23421079, ReadOrientation.F2R1, ArtifactState.F2R1_G), // CAC->G F2R1
                new ImmutableTriple<>(34144749, ReadOrientation.F2R1, ArtifactState.F2R1_A), // TGA->A F2R1, ref context is not canonical
                new ImmutableTriple<>(62165528, ReadOrientation.F1R2, ArtifactState.F1R2_A)); // CGC->A F1R2

        for (final Triple<Integer, ReadOrientation, ArtifactState> artifact : knownArtifacts) {
            final int position = artifact.getLeft();
            final ReadOrientation expectedReadOrientaiton = artifact.getMiddle();
            final ArtifactState expectedSourceOfPrior = artifact.getRight();

            Optional<VariantContext> variant = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                    .filter(vc -> vc.getStart() == position).findFirst();
            Assert.assertTrue(variant.isPresent());

            // Check that the expected filters were applied
            Assert.assertTrue(variant.get().getFilters().contains(GATKVCFConstants.READ_ORIENTATION_ARTIFACT_FILTER_NAME));
        }
    }
}