package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.barclay.argparser.Argument;

public class CollectF1R2CountsArgumentCollection {
    public static final String ALT_DATA_TABLE_LONG_NAME = "alt-table";
    public static final String ALT_DEPTH1_HISTOGRAM_LONG_NAME = "alt-hist";
    public static final String REF_SITE_METRICS_LONG_NAME = "ref-hist";
    public static final String MIN_MEDIAN_MQ_LONG_NAME = "median-mq";
    public static final String MIN_BASE_QUALITY_LONG_NAME = "min-bq";
    public static final String MAX_DEPTH_LONG_NAME = "max-depth";

    @Argument(fullName = MIN_MEDIAN_MQ_LONG_NAME, doc = "skip sites with median mapping quality below this value", optional = true)
    public int minMedianMapQual = 30;

    @Argument(fullName = MIN_BASE_QUALITY_LONG_NAME, doc = "exclude bases below this quality from pileup", optional = true)
    public int minBaseQuality = 20;

    @Argument(fullName = MAX_DEPTH_LONG_NAME, doc = "sites with depth higher than this value will be grouped", optional = true)
    public int maxDepth = F1R2FilterConstants.DEFAULT_MAX_DEPTH;
}
