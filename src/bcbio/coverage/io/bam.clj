(ns bcbio.coverage.io.bam
  "Reading from sequence reads in BAM format, leveraging Picard."
  (:import [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency]
           [net.sf.picard.sam BuildBamIndex]
           [org.broadinstitute.sting.commandline Tags]
           [ org.broadinstitute.sting.gatk.arguments ValidationExclusion]
           [org.broadinstitute.sting.gatk.datasources.reads
            LocusShard SAMDataSource SAMReaderID]
           [org.broadinstitute.sting.gatk.datasources.providers
            CoveredLocusView LocusShardDataProvider]
           [org.broadinstitute.sting.gatk.executive WindowMaker]
           [org.broadinstitute.sting.gatk.filters DuplicateReadFilter
            FailsVendorQualityCheckFilter NotPrimaryAlignmentFilter
            MalformedReadFilter UnmappedReadFilter]
           [org.broadinstitute.sting.gatk.resourcemanagement ThreadAllocation]
           [org.broadinstitute.sting.utils GenomeLocParser]
           )
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [bcbio.align.gref :as gref]
            [bcbio.run.itx :as itx]))

(defn index-bam
  "Generate BAM index, skipping if already present."
  [in-bam]
  (let [index-choices [(str in-bam ".bai") (string/replace in-bam ".bam" ".bai")]]
    (if-let [index-file (->> index-choices
                             (remove itx/needs-run?)
                             (filter #(itx/up-to-date? % in-bam))
                             first)]
      index-file
      (do
        (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
        (itx/with-tx-file [tx-index-file (first index-choices)]
          (BuildBamIndex/createIndex (SAMFileReader. (io/file in-bam)) (io/file tx-index-file)))
        (first index-choices)))))

(defn get-bam-source
  "Retrieve Picard SAMFileReader for input BAM, with querying."
  [bam-file]
  (let [index-file (index-bam bam-file)]
    (SAMFileReader. (io/file bam-file) (io/file index-file))))

(defn get-sample-names
  [bam-file]
  (with-open [rdr (SAMFileReader. (io/file bam-file))]
    (->> rdr
         .getFileHeader
         .getReadGroups
         (map #(.getSample %)))))

(defn- get-gatk-bam-sources
  "Create required source references for GATK mapping to BAM files."
  [bam-files loc-parser]
  (let [file-ids (map #(SAMReaderID. (io/file %) (Tags.)) bam-files)
        file-spans (zipmap file-ids (map #(-> % io/file SAMFileReader. .getFilePointerSpanningReads)
                                         bam-files))
        filters [(DuplicateReadFilter.) (FailsVendorQualityCheckFilter.)
                 (NotPrimaryAlignmentFilter.) (UnmappedReadFilter.)]
        data-source (SAMDataSource. file-ids (ThreadAllocation.) nil loc-parser
                                    false SAMFileReader$ValidationStringency/LENIENT
                                    nil nil (ValidationExclusion.) filters false)]
    [data-source file-spans]))

(defprotocol GATKIteration
  (get-align-contexts [this]))

(defrecord GATKLocationIterator [view shard window-maker shard-dp]
  GATKIteration
  (get-align-contexts [_]
    (iterator-seq view))
  java.io.Closeable
  (close [_]
    (.close shard)
    (.close window-maker)
    (.close shard-dp)
    (.close view)))

(defn prep-bam-region-iter
  "Create an iterator returning GATK AlignmentContexts over the provided coordinates.
   Navigates the depths of GATK classes to provide a simple way to get pileups at
   each position in a set of BAM files."
  [bam-files ref-file coord]
  (let [sample-names (mapcat get-sample-names bam-files)
        loc-parser (GenomeLocParser. (gref/get-seq-dict ref-file))
        loc (.createGenomeLoc loc-parser (:chr coord) (:start coord) (:end coord))
        [data-source file-spans] (get-gatk-bam-sources bam-files loc-parser)
        shard (LocusShard. loc-parser data-source [loc] file-spans)
        window-maker (WindowMaker. shard loc-parser (.seek data-source shard) [loc] sample-names)
        shard-dp (LocusShardDataProvider. shard (.getReadsInfo data-source) loc-parser loc
                                          (.next window-maker) nil nil)
        view (CoveredLocusView. shard-dp)]
    (GATKLocationIterator. view shard window-maker shard-dp)))
