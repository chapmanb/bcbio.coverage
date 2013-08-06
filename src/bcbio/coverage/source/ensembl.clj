(ns bcbio.coverage.source.ensembl
  "Retrieve gene and exon coordinate information from Ensembl using API."
  (:require [clojure.string :as string]
            [clojure.java.io :as io]
            [org.bioclojure.bio.ensembl.core :as ens]
            [bcbio.coverage.io.bed :as bed]
            [bcbio.run.itx :as itx]
            [bcbio.run.parallel :refer [rmap]]
            ))

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
    {:gene-name gene-name
     :coords (->> gene
                  ens/gene-transcripts
                  (filter #(or (not= "canonical" (:transcripts params))
                               (.isCanonical %)))
                  (filter #(= "protein_coding" (.getBiotype %)))
                  (mapcat transcript->translated-coords)
                  bed/merge-intervals)}
    (throw (Exception. (format "Did not find Gene information for %s in Ensembl" gene-name)))))

(defn- gene-str->coding-coords
  "Identify coding coordinates for gene information from line of text."
  [line params]
  (ens/with-registry (ens/registry :ensembldb)
    (-> line
        (string/split #"\s+")
        first string/trimr
        (fetch-coding-coords params))))

(defn- genes->coding-bed
  "Convert a file of gene names into a BED file with coding coordinates."
  [gene-file params]
  (let [approach-str (if-let [x (:transcripts params)] (str "-" x) "")
        bed-file (str (itx/file-root gene-file) "-exons" approach-str ".bed")]
    (when (itx/needs-run? bed-file)
      (itx/with-tx-file [tx-bed-file bed-file]
        (with-open [rdr (io/reader gene-file)
                    wtr (io/writer tx-bed-file)]
          (doseq [{:keys [gene-name coords]} (rmap #(gene-str->coding-coords % params) (line-seq rdr)
                                                   (get params :cores 1) (get params :chunk-size 1))]
            (doseq [coord coords]
              (.write wtr (format "%s\t%s\t%s\t%s\n"
                                  (:chr coord) (dec (:start coord)) (:end coord) gene-name)))))))
    bed-file))

(defn get-coord-bed
  "Retrieve a BED file of coordinates, handling gene name or coordinate input."
  [in-file params]
  (if (bed/is-bed? in-file)
    in-file
    (genes->coding-bed in-file params)))
