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
            [bcbio.run.itx :as itx]))

;; ## Regions from gene names

(defn- fetch-coding-coords
  "Retrieve a flattened set of coordinates for coding regions attached to a gene."
  [gene-name params]
  (let [species (get params :species "human")
        gene (first (ens/gene-name->genes species gene-name))
        coords (mapcat ens/transcript->exon-coords (ens/gene-transcripts gene))]
    (bed/merge-intervals coords)))

(defn- genes->coding-bed
  "Convert a file of gene names into a BED file with coding coordinates."
  [gene-file params]
  (let [bed-file (str (itx/file-root gene-file) "-exons.bed")]
    (when (itx/needs-run? bed-file)
      (with-open [rdr (io/reader gene-file)
                  wtr (io/writer bed-file)]
        (ens/with-registry (ens/registry :ensembldb)
          (doseq [gene-name (map #(-> % (string/split #"\s+") first string/trimr) (line-seq rdr))]
            (doseq [coord (fetch-coding-coords gene-name params)]
              (.write wtr (format "%s\t%s\t%s\t%s\n"
                                  (:chr coord) (dec (:start coord)) (:end coord) gene-name)))))))
    bed-file))

(defn get-coord-bed
  "Retrieve a BED file of coordinates, handling gene name or coordinate input."
  [in-file params]
  (if (bed/is-bed? in-file)
    in-file
    (genes->coding-bed in-file params)))

;; ## Coverage retrieval

(defmulti get-coverage*
  "Retrieve coverage, dispatching on type of input file."
  (fn [ftype & args] ftype))

(defmethod get-coverage* :bw
  ^{:doc "Retrieve coverage at a position from a BigWig file of by-position coverage"}
  [_ source contig pos]
  (let [bw-item (first (iterator-seq
                        (.getBigWigIterator source contig pos contig (inc pos) false)))]
    (if bw-item
      (.getWigValue bw-item)
      0)))

(defmethod get-coverage* :bam
  ^{:doc "Retrieve position coverage directly from a BAM file"}
  [_ source contig pos]
  (with-open [iter (.queryOverlapping source contig (inc pos) (inc pos))]
    (count (iterator-seq iter))))

(defprotocol RetrieveCoverage
  (get-coverage [this contig pos]))

(defrecord CoverageRetriever [ftype source]
  RetrieveCoverage
  (get-coverage [_ contig pos]
    (get-coverage* ftype source contig pos))
  java.io.Closeable
  (close [_]
    (.close source)))

(defrecord PopCoverageRetriever [ftype sources]
  RetrieveCoverage
  (get-coverage [_ contig pos]
    (istat/median (map #(get-coverage* ftype % contig pos))))
  java.io.Closeable
  (close [_]
    (doseq [x sources]
      (.close x))))

(defmulti get-coverage-retriever* class)

(defn- ext->ftype
  "Extract file type to process from input file extension"
  [in-file]
  (-> in-file fs/extension string/lower-case (subs 1) keyword))

(defn- get-source
  [in-file ftype]
  (case ftype
    :bam (bam/get-bam-source in-file)
    :bw (bigwig/get-source in-file)))

(defmethod get-coverage-retriever* java.lang.String
  ^{:doc "Get a retriever that calculates coverage per position, handling alternative file types."}
  [in-file]
  (let [ftype (ext->ftype in-file)
        source (get-source in-file ftype)]
    (CoverageRetriever. ftype source)))

(defmethod get-coverage-retriever* clojure.lang.Seqable
  ^{:doc "Prepare a retriever calculating the average over multiple input files."}
  [in-files]
  (let [ftypes (set (map ext->ftype in-files))]
    (if (> (count (ftypes)) 1)
      (throw (Exception. "Cannot retrieve from multiple heterogeneous input sources."))
      (PopCoverageRetriever. (first ftypes) (map #(get-source % (first ftypes)) in-files)))))

(defn get-coverage-retriever
  "Retrieve coverage, handling multiple input files and different coverage types."
  [file-info]
  (let [test-file (if (and (not (instance? java.lang.String file-info))
                           (= 1 (count file-info)))
                    (first file-info)
                    file-info)]
    (get-coverage-retriever* test-file)))

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
  [cov-by-pos params]
  (->> (keys cov-by-pos)
       (split-into-blocks (get-in params [:block :distance] 100.0))
       (remove empty?)
       (map (juxt first last))
       (filter (fn [[s e]] (> (- e s) (get-in params [:block :min] 10.0))))))

(defn region-problem-coverage
  "Calculate stats for problematic coverage in a chromosome region."
  [retriever contig start end params]
  (let [cov (reduce (fn [coll i]
                      (let [n (get-coverage retriever contig i)]
                        (if (< n (:coverage params 10.0))
                          (assoc-in coll [:low i] n)
                          (assoc-in coll [:high i] n))))
                    {:high {} :low {}} (range start end))]
    {:count (count (:low cov))
     :blocks (identify-nocoverage-blocks (:low cov) params)
     :coverages (concat (vals (:high cov)) (vals (:low cov)))
     :size (- end start)}))

(defn gene-problem-coverage
  "Calculate stats for problematic coverage across a set of gene regions."
  [retriever coords params]
  (let [regions (map #(region-problem-coverage retriever (:chr %) (:start %) (:end %) params)
                     coords)]
    {:name (-> coords first :name)
     :coords (map #(select-keys % [:chr :start :end]) coords)
     :avg-reads (istat/mean (mapcat :coverages regions))
     :blocks (mapcat :blocks regions)
     :size (apply + (map :size regions))
     :percent-nocoverage (* 100.0 (/ (apply + (map :count regions))
                                     (apply + (map :size regions))))}))

(defn problem-coverage
  "Identify problematic coverage for all supplied gene regions.
   params is a map of attributes to influence what we identify
    :coverage - Minimum coverage required for a position to be considered covered.
    :block - Map of parameters for calculating blocks of nocoverage
      :min - Minimum size of a block to report
      :distance - Allowed distance between nocoverage bases to be considered in a block
    :organism - Organism Ensembl name"
  [coverage-input gene-file params]
  (let [retriever (get-coverage-retriever coverage-input)
        gene-coord-file (get-coord-bed gene-file params)]
    (with-open [rdr (io/reader gene-coord-file)]
      (doseq [coords (map second (group-by :name (bed/get-iterator rdr)))]
        (println (gene-problem-coverage retriever coords params))))))

(defn -main [& args]
  (let [[options args banner] (cli args)]
    (when (not= (count args) 2)
      (println "Usage: gene <BAM or BigWig coverage file> <BED file of gene regions>")
      (System/exit 1))
    (apply problem-coverage (concat args [{:coverage 10
                                           :block {:min 50 :distance 5}}]))))
