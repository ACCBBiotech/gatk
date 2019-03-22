package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class LearnReadOrientationModelUnitTest extends CommandLineProgramTest {
    @Test
    public void testCombineHistograms(){
        final File unscatteredDir = createTempDir("unscattered");
        final File scatteredDir = createTempDir("scattered");

        // Step 0: Run CollectF1R2Counts without scatter-gather
        runCommandLine(Arrays.asList(
                "-R", b37_reference_20_21,
                "-I", ReadOrientationModelIntegrationTest.hapmapBamSnippet,
                "-L", ReadOrientationModelIntegrationTest.intervalList,
                "-O", unscatteredDir.getAbsolutePath() + "/"),
                CollectF1R2Counts.class.getSimpleName());


        // Step 1: SplitIntervals
        final File intervalDir = createTempDir("intervals");
        final int scatterCount = 5;
        runCommandLine(Arrays.asList(
            "-R", b37_reference_20_21,
            "-L", ReadOrientationModelIntegrationTest.intervalList,
            "-O", intervalDir.getAbsolutePath(),
            "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, Integer.toString(scatterCount)),
            SplitIntervals.class.getSimpleName());

        // Step 2: CollectF1R2Counts
        final File[] intervals = intervalDir.listFiles();
        for (int i = 0; i < intervals.length; i++){
            runCommandLine(Arrays.asList(
                    "-R", b37_reference_20_21,
                    "-I", ReadOrientationModelIntegrationTest.hapmapBamSnippet,
                    "-L", intervals[i].getAbsolutePath(),
                    "-O", scatteredDir.getAbsolutePath() + "/scatter_" + i),
                    CollectF1R2Counts.class.getSimpleName());
        }

        final List<File> scatteredRefMetricsFiles = Arrays.stream(scatteredDir.listFiles())
                .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.REF_HIST_EXTENSION)).collect(Collectors.toList());

        final List<File> scatteredAltMetricsFiles = Arrays.stream(scatteredDir.listFiles())
                .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.ALT_HIST_EXTENSION)).collect(Collectors.toList());

        final List<File> scatteredAltTableFiles = Arrays.stream(scatteredDir.listFiles())
                .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.ALT_TABLE_EXTENSION)).collect(Collectors.toList());

        final List<Histogram<Integer>> ref = LearnReadOrientationModel.sumHistogramsFromFiles(scatteredRefMetricsFiles, true);
        final List<Histogram<Integer>> alt = LearnReadOrientationModel.sumHistogramsFromFiles(scatteredAltMetricsFiles, false);
        final List<AltSiteRecord> altSites = LearnReadOrientationModel.gatherAltSiteRecords(scatteredAltTableFiles).getRight();

        final File refHistUnscattered = Arrays.stream(unscatteredDir.listFiles())
                .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.REF_HIST_EXTENSION)).findFirst().get();

        final File altHistUnscattered = Arrays.stream(unscatteredDir.listFiles())
                .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.ALT_HIST_EXTENSION)).findFirst().get();

        final File altTableUnscattered = Arrays.stream(unscatteredDir.listFiles())
                .filter(file -> file.getAbsolutePath().endsWith(F1R2CountsCollector.ALT_TABLE_EXTENSION)).findFirst().get();

        final List<Histogram<Integer>> refTruth = LearnReadOrientationModel.readMetricsFile(refHistUnscattered).getAllHistograms();
        final List<Histogram<Integer>> altTruth = LearnReadOrientationModel.readMetricsFile(altHistUnscattered).getAllHistograms();
        final List<AltSiteRecord> altSitesTruth = AltSiteRecord.readAltSiteRecords(altTableUnscattered.toPath()).getRight();


        for (Histogram<Integer> truth : refTruth){
            final Histogram<Integer> eval = ref.stream().filter(h -> h.getValueLabel().equals(truth.getValueLabel())).findAny().get();
            Assert.assertEquals(eval.getSumOfValues(), truth.getSumOfValues());
        }

        for (Histogram<Integer> truth : altTruth){
            final Histogram<Integer> eval = alt.stream().filter(h -> h.getValueLabel().equals(truth.getValueLabel())).findAny().get();
            Assert.assertEquals(eval.getSum(), truth.getSum());
            Assert.assertEquals(eval.getSumOfValues(), truth.getSumOfValues());
        }

        Assert.assertEquals(altSites.size(), altSitesTruth.size());

    }
}
