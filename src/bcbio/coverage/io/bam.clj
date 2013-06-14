(ns bcbio.coverage.io.bam
  "Reading from sequence reads in BAM format, leveraging Picard."
  (:import [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency]
           [net.sf.picard.sam BuildBamIndex]
           [org.broadinstitute.sting.commandline Tags]
           [ org.broadinstitute.sting.gatk.arguments ValidationExclusion]
           [org.broadinstitute.sting.gatk.datasources.reads
            LocusShardBalancer SAMDataSource SAMReaderID]
           [org.broadinstitute.sting.gatk.datasources.providers
            CoveredLocusView LocusShardDataProvider]
           [org.broadinstitute.sting.gatk.downsampling
            DownsamplingMethod DownsampleType]
           [org.broadinstitute.sting.gatk.executive WindowMaker]
           [org.broadinstitute.sting.gatk.filters DuplicateReadFilter
            FailsVendorQualityCheckFilter NotPrimaryAlignmentFilter
            MalformedReadFilter UnmappedReadFilter]
           [org.broadinstitute.sting.gatk.resourcemanagement ThreadAllocation]
           [org.broadinstitute.sting.utils GenomeLocParser GenomeLocSortedSet])
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

(def ^{:private true} gatk-filters
  {:duplicate (DuplicateReadFilter.)
   :vendor-quality (FailsVendorQualityCheckFilter.)
   :not-primary (NotPrimaryAlignmentFilter.)
   :unmapped (UnmappedReadFilter.)})

(defn- get-gatk-bam-source
  "Create required source references for GATK mapping to BAM files."
  [bam-files loc-parser filters downsample]
  (let [file-ids (map #(SAMReaderID. (io/file %) (Tags.)) bam-files)
        gatk-filters (map #(get gatk-filters %) filters)
        data-source (SAMDataSource. file-ids (ThreadAllocation.) nil loc-parser
                                    false SAMFileReader$ValidationStringency/LENIENT
                                    nil
                                    (when downsample (DownsamplingMethod. DownsampleType/BY_SAMPLE
                                                                          (int downsample) nil))
                                    (ValidationExclusion.) gatk-filters false)]
    data-source))

(defn- get-gatk-shards
  "Retrieve an iterator over GATK shard region of the genomes for specified locations."
  [data-source locs loc-parser]
  (-> data-source
      (.createShardIteratorOverIntervals (GenomeLocSortedSet. loc-parser locs)
                                         (LocusShardBalancer.))
      .iterator
      iterator-seq))

(defn- get-gatk-views
  "Retrieve GATK LocusView objects for a given shard region of the genome."
  [shard data-source loc-parser sample-names]
  (let [windows (iterator-seq (WindowMaker. shard loc-parser (.seek data-source shard)
                                            (.getGenomeLocs shard) sample-names))]
    (map #(CoveredLocusView.
           (LocusShardDataProvider. shard (.getSourceInfo %) loc-parser
                                    (.getLocus %) % nil nil))
         windows)))

(defprotocol GATKIteration
  (get-align-contexts [this]))

(defrecord GATKLocationIterator [views]
  GATKIteration
  (get-align-contexts [_]
    (mapcat #(iterator-seq %) views))
  java.io.Closeable
  (close [_]
    (doseq [view views]
      (.close view))))

(defn prep-bam-region-iter
  "create an iterator returning GATK AlignmentContexts over the provided locus coordinates.
   Navigates the depths of GATK classes to provide a simple way to get pileups at
   each position in a set of BAM files."
  [bam-files ref-file coords & {:keys [filters downsample]
                               :or {filters [:duplicate :vendor-quality :not-primary :unmapped]}}]
  (let [sample-names (mapcat get-sample-names bam-files)
        loc-parser (GenomeLocParser. (gref/get-seq-dict ref-file))
        locs (map #(.createGenomeLoc loc-parser (:chr %) (:start %) (:end %)) coords)
        data-source (get-gatk-bam-source bam-files loc-parser filters downsample)
        shards (get-gatk-shards data-source locs loc-parser)
        views (mapcat #(get-gatk-views % data-source loc-parser sample-names) shards)]
    (GATKLocationIterator. views)))
