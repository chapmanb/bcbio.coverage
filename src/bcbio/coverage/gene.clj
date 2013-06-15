(ns bcbio.coverage.gene
  "Identify gene regions with low and no coverage blocks."
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clojure.tools.cli :refer [cli]]
            [incanter.stats :as istat]
            [me.raynes.fs :as fs]
            [org.bioclojure.bio.ensembl.core :as ens]
            [bcbio.coverage.io.bam :as bam]
            [bcbio.coverage.io.bed :as bed]
            [bcbio.coverage.io.bigwig :as bigwig]
            [bcbio.run.itx :as itx]
            [bcbio.run.parallel :refer [rmap]]))

;; ## Regions from gene names

(defn- transcript->translated-coords
  "Retrieve chromosomal coordinates of translated protein coding sequence in a transcript.
   Works around issue with Ensembl/jEnsembl where translated coordinates
   map into introns, creation large apparent exons that include introns."
  [t]
  (let [chr (-> t .getChromosomeMapping .getTarget .getChromosomeName)
        strand (-> t .getChromosomeMapping .getTargetCoordinates .getStrandInt)
        trl (.getCanonicalTranslation t)]
    (remove nil?
            (map (fn [m]
                   (let [raws (-> m .getSourceCoordinates .getStart)
                         rawe (-> m .getSourceCoordinates .getEnd)
                         c1 (.getChromosomePositionFromBASE trl raws)
                         c2 (.getChromosomePositionFromBASE trl rawe)
                         [s e] (if (neg? strand) [c2 c1] [c1 c2])]
                     ; Sanity check our remapping matches expected
                     ; size, remove regions that do not
                     (when (= (- e s) (Math/abs (- raws rawe)))
                       {:start s :end e :chr chr :strand strand})))
                 (.getTranslationMappings trl)))))

(defn- gene-name-matches
  "Annotate gene match with information on all possible names for a gene.
   Tries to catch synonyms while preferring the primary name for multiple matches."
  [g name-orig]
  (let [name (string/upper-case name-orig)
        display (string/upper-case (.getDisplayName g))
        aliases (set (map string/upper-case (.getAllSynonyms g)))]
    {:g g
     :count (cond
             (= name display) 2
             (contains? aliases name) 1
             :else 0)}))

(defn- fetch-coding-coords
  "Retrieve a flattened set of coordinates for coding regions attached to a gene."
  [gene-name params]
  (if-let [gene (->> (ens/gene-name->genes (get params :species "human") gene-name)
                     (map #(gene-name-matches % gene-name))
                     (filter #(pos? (:count %)))
                     (sort-by :count >)
                     (map :g)
                     first)]
    (->> gene
         ens/gene-transcripts
         (filter #(= "protein_coding" (.getBiotype %)))
         (mapcat transcript->translated-coords)
         bed/merge-intervals)
    (throw (Exception. (format "Did not find Gene information for %s in Ensembl" gene-name)))))

(defn- genes->coding-bed
  "Convert a file of gene names into a BED file with coding coordinates."
  [gene-file params]
  (let [bed-file (str (itx/file-root gene-file) "-exons.bed")]
    (when (itx/needs-run? bed-file)
      (itx/with-tx-file [tx-bed-file bed-file]
        (with-open [rdr (io/reader gene-file)
                    wtr (io/writer tx-bed-file)]
          (ens/with-registry (ens/registry :ensembldb)
            (doseq [gene-name (map #(-> % (string/split #"\s+") first string/trimr) (line-seq rdr))]
              (doseq [coord (fetch-coding-coords gene-name params)]
                (.write wtr (format "%s\t%s\t%s\t%s\n"
                                    (:chr coord) (dec (:start coord)) (:end coord) gene-name))))))))
    bed-file))

(defn get-coord-bed
  "Retrieve a BED file of coordinates, handling gene name or coordinate input."
  [in-file params]
  (if (bed/is-bed? in-file)
    in-file
    (genes->coding-bed in-file params)))

;; ## Coverage retrieval

(defn- ext->ftype
  "Extract file type to process from input file extension"
  [in-file]
  (-> in-file fs/extension string/lower-case (subs 1) keyword))

(defmulti get-coverage*
  (fn [ftype & args]
    ftype))

(defn- counts->cov
  [coord counts]
  (map (fn [i] {:i i :n (get counts i 0)}) (range (:start coord) (:end coord))))

(defmethod get-coverage* :bam
  ^{:doc "Retrieve average coverage information over a genomic region using GATK pileups."}
  [_ in-files coord ref-file params]
  ;; Perform initial setup of BAM files, including indexing
  (doseq [source (rmap bam/get-bam-source in-files
                       (get params :cores 1) 1)]
    (.close source))
  (with-open [iter (bam/prep-bam-region-iter in-files ref-file [coord]
                                             :downsample (:downsample params))]
    (let [bam-counts (reduce (fn [coll x]
                               (assoc coll (dec (.getPosition x))
                                      (/ (.size x) (count in-files))))
                             {} (bam/get-align-contexts iter))]
      (counts->cov coord bam-counts))))

(defmethod get-coverage* :bw
  ^{:doc "Retrieve coverage information from input BigWig file"}
  [_ in-files coord _ _]
  (when (> (count in-files) 1)
    (throw (Exception. "Do not currently handle multiple BigWig inputs.")))
  (with-open [source (bigwig/get-source (first in-files))]
    (let [bw-counts (reduce (fn [coll x]
                              (if-let [v (.getWigValue x)]
                                (assoc coll (.getStartBase x) v)
                                coll))
                            {} (iterator-seq (.getBigWigIterator source (:chr coord) (:start coord)
                                                                 (:chr coord) (inc (:end coord)) false)))]
      (counts->cov coord bw-counts))))

(defn get-coverage
  "Retrieve coverage in a region, handling multiple input files
  and heterogeneous file types."
  [input-info coord ref-file params]
  (let [in-files (if (instance? java.lang.String input-info) [input-info] input-info)
        ftypes (set (map ext->ftype in-files))]
    (if (> (count ftypes) 1)
      (throw (Exception. "Cannot retrieve from multiple heterogeneous input sources."))
      (get-coverage* (first ftypes) in-files coord ref-file params))))

(defn split-into-blocks
  "Split a sequence of positions into blocks of consecutive items within n bases."
  [n in-xs]
  (loop [xs (sort in-xs)
         cur-block []
         final []]
    (cond
     (empty? xs)
     (conj final cur-block)
     (or (empty? cur-block)
         (<= (- (first xs) (last cur-block)) n))
     (recur (rest xs) (conj cur-block (first xs)) final)
     :else
     (recur (rest xs) [(first xs)] (conj final cur-block)))))

(defn- identify-nocoverage-blocks
  "Extract blocks of nocoverage bases from raw positions"
  [low-is params]
  (->> low-is
       (split-into-blocks (get-in params [:block :distance] 100.0))
       (remove empty?)
       (map (juxt first last))
       (filter (fn [[s e]] (> (- e s) (get-in params [:block :min] 10.0))))))


(defn coverage-report-region
  "Calculate stats for problematic coverage in a chromosome region."
  [input-files coord ref-file params]
  (let [cov (->> (get-coverage input-files coord ref-file params)
                 (map #(assoc % :low (< (:n %) (get params :coverage 10.0)))))]
    {:count (count (filter :low cov))
     :coord coord
     :blocks (identify-nocoverage-blocks (->> cov (filter :low) (map :i)) params)
     :coverages (map :n (filter #(not (:low %)) cov))
     :size (- (:end coord) (:start coord))}))

(defn coverage-report-bygene
  "Calculate coverage reports for all passed regions, grouped by gene."
  [input-files coords ref-file params]
  (let [all-regions (rmap #(coverage-report-region input-files % ref-file params)
                          coords (get params :cores 1) (get params :chunk-size 1))]
    (for [[name regions] (group-by #(get-in % [:coord :name]) all-regions)]
      {:name name
       :coords (->> regions
                    (map :coord)
                    (map #(select-keys % [:chr :start :end])))
       :avg-reads (format "%.1f" (istat/mean (mapcat :coverages regions)))
       :blocks (mapcat :blocks regions)
       :size (apply + (map :size regions))
       :percent-nocoverage (format "%.1f" (* 100.0 (/ (apply + (map :count regions))
                                                      (apply + (map :size regions)))))})))

(defn problem-coverage
  "Identify problematic coverage for all supplied gene regions.
   params is a map of attributes to influence what we identify
    :coverage - Minimum coverage required for a position to be considered covered.
    :block - Map of parameters for calculating blocks of nocoverage
      :min - Minimum size of a block to report
      :distance - Allowed distance between nocoverage bases to be considered in a block
    :organism - Organism Ensembl name"
  [coverage-input gene-file ref-file params]
  (let [gene-coord-file (get-coord-bed gene-file params)]
    (with-open [rdr (io/reader gene-coord-file)]
      (doseq [report (coverage-report-bygene coverage-input (bed/get-iterator rdr)
                                             ref-file params)]
        (println report)))))

(defn -main [& args]
  (let [[options args banner] (cli args)]
    (when (not= (count args) 3)
      (println "Usage: gene <BAM or BigWig coverage file> <BED file of gene regions>" "<Reference FASTA file>")
      (System/exit 1))
    (apply problem-coverage (concat args [{:coverage 10
                                           :block {:min 50 :distance 5}}]))))
