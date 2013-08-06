(ns bcbio.coverage.source.ensembl
  "Retrieve gene and exon coordinate information from Ensembl using API."
  (:require [clojure.string :as string]
            [clojure.java.io :as io]
            [org.bioclojure.bio.ensembl.core :as ens]
            [bcbio.coverage.io.bed :as bed]
            [bcbio.run.itx :as itx]
            [bcbio.run.parallel :refer [rmap]]
            [me.raynes.fs :as fs])
  (:import [uk.ac.roslin.ensembl.model.core Species]))

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

(defn- name->gene
  "Retrieve a single gene from a gene-name, finding best representation for multiple."
  [gene-name params]
  (->> (ens/gene-name->genes (get params :species "human") gene-name)
       (map #(gene-name-matches % gene-name))
       (filter #(pos? (:count %)))
       (sort-by :count >)
       (map :g)
       first))

(defn- fetch-coding-coords
  "Retrieve a flattened set of coordinates for coding regions attached to a gene."
  [gene-name params]
  (if-let [gene (if (.startsWith gene-name "ENS")
                  (ens/gene (:species params) gene-name)
                  (name->gene gene-name params))]
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
  [lines params]
  (ens/with-registry (ens/registry :ensembldb)
    (vec (map
          #(-> %
               (string/split #"\s+")
               first string/trimr
               (fetch-coding-coords params))
          lines))))

(defn- genes->coding-bed
  "Convert a file of gene names into a BED file with coding coordinates."
  [gene-file params]
  (let [approach-str (if-let [x (:transcripts params)] (str "-" x) "")
        bed-file (str (itx/file-root gene-file) "-exons" approach-str ".bed")
        cores (get params :cores 1)
        work-size 20]
    (when (itx/needs-run? bed-file)
      (itx/with-tx-file [tx-bed-file bed-file]
        (with-open [rdr (io/reader gene-file)
                    wtr (io/writer tx-bed-file)]
          (doseq [line-block (partition-all (* cores work-size) (line-seq rdr))]
            (doseq [coord-block (rmap #(gene-str->coding-coords % params)
                                      (partition-all work-size line-block)
                                      (get params :cores 1) (get params :chunk-size 1))]
              (doseq [{:keys [gene-name coords]} coord-block]
                (doseq [coord coords]
                  (.write wtr (format "%s\t%s\t%s\t%s\n"
                                      (:chr coord) (dec (:start coord)) (:end coord) gene-name)))))))))
    bed-file))

(defn get-coord-bed
  "Retrieve a BED file of coordinates, handling gene name or coordinate input."
  [in-file params]
  (if (bed/is-bed? in-file)
    in-file
    (genes->coding-bed in-file params)))

;; ## Gene names for organisms

(defn- chr->exons
  "Retrieve transcript exons on a specific chromosome. Enables parallelization."
  [chr-name base-out-file params]
  (let [out-file (str (io/file (fs/parent base-out-file)
                               (str (fs/name base-out-file) "-parts")
                               (format "%s-%s.txt" (fs/name base-out-file) chr-name)))]
    (when (itx/needs-run? out-file)
      (when-not (fs/exists? (fs/parent out-file))
        (fs/mkdirs (fs/parent out-file)))
      (itx/with-tx-file [tx-out-file out-file]
        (with-open [wtr (io/writer tx-out-file)]
          (ens/with-registry (ens/registry :ensembldb)
            (let [species (ens/species (:species params))
                  chr (get (.getChromosomes species) chr-name)]
              (println (.getChromosomeName chr) (.getLength chr))
              (doseq [gene (.getGenesOnRegion chr (int 1) (int (.getLength chr)))]
                (when (first (.getHomologies gene))
                  (.write wtr (str (.getStableID gene) "\n")))))))))
    out-file))

(defn coding-genes
  "Retrieve list of evolutionarily conserved coding genes to use for
   full whole genome assessment of coverage."
  [out-file params]
  (when (itx/needs-run? out-file)
    (ens/with-registry (ens/registry :ensembldb)
      (let [chr-files (rmap #(chr->exons % out-file params)
                            (ens/list-chromosomes (:species params))
                            (get params :cores 1) (get params :chunk-size 1))]
        (itx/with-tx-file [tx-out-file out-file]
          (with-open [wtr (io/writer tx-out-file)]
            (doseq [chr-file chr-files]
              (with-open [rdr (io/reader chr-file)]
                (doseq [line (line-seq rdr)]
                  (.write wtr (str line "\n")))))))))))

(defn species-exon-coords
  "Retrieve the exon coordinates for conserved transcript regions in an organism.
   params is a map with: species, transcripts (specify canonical, if desired),
   cores (for parallelization)"
  [out-file params]
  (-> out-file
      (coding-genes params)
      (get-coord-bed params)))
