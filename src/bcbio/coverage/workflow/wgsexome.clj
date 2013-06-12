(ns bcbio.coverage.workflow.wgsexome
  "Compare whole genome sequencing data to exome data, assessing coverage and variant calls."
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [clojure.pprint :refer [pprint]]
            [clojure.tools.cli :refer [cli]]
            [clj-yaml.core :as yaml]
            [incanter.charts :as icharts]
            [incanter.core :as icore]
            [bcbio.coverage.gene :as gene]
            [bcbio.coverage.io.bed :as bed]
            [bcbio.coverage.source.esp :as esp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.variantcontext :as gvc]))

;; ## Reports

(defn- variant-coverage
  "Calculate coverage of a defined set of variant coordinates."
  [cov-get vrn-coords params]
  (reduce (fn [coll coord]
            (let [n (gene/get-coverage cov-get (:chr coord) (:start coord))
                  kw (if (< n (:coverage params 10.0)) :uncovered :covered)]
              (assoc coll kw (inc (get coll kw)))))
          {:uncovered 0 :covered 0} vrn-coords))

(defn- coverage-report
  "Provide summary of coverage for a set of coordinates.
   Calculates:
   - Summary metrics of coverage in gene.
   - High frequency variant calls in comparison without coverage in current."
  [cov-get coords vrn-coords params]
  (-> (gene/gene-problem-coverage cov-get coords params)
      (assoc :variant-compare (variant-coverage cov-get vrn-coords params))
      (dissoc :coords)))

(defmulti sample-variants
  "Retrieve variants from different input types for a group of samples."
  (fn [config _ _ _]
    (keyword (get config :type "vcf"))))

(defmethod sample-variants :esp
  ^{:doc "Retrieval from summarized ESP variant files."}
  [config region ref-file params]
  (esp/variants-in-region (:variant config) ref-file region (get params :maf-thresh 0.0)))

(defn- vc-has-sample?
  "Does the current variant context have a called variant in the provided samples."
  [samples vc af-thresh]
  (letfn [(non-ref-alleles [g]
            (filter #(and (.isCalled %) (.isNonReference %)) (:alleles g)))]
    (let [sample-gs (filter #(contains? samples (:sample-name %)) (:genotypes vc))
          ac (count (mapcat non-ref-alleles sample-gs))
          total (count (mapcat :alleles sample-gs))
          af (/ ac total)]
      (> af af-thresh))))

(defmethod sample-variants :vcf
  ^{:doc "Retrieval from standard VCF per sample variant files."}
  [config region ref-file params]
  (let [samples (set (map :name (:samples config)))]
    (with-open [vcf-get (gvc/get-vcf-retriever ref-file (:variant config))]
      (vec (filter #(vc-has-sample? samples % (get params :af-thresh 0.0))
                   (gvc/variants-in-region vcf-get region))))))

(defn- variant-report
  "Identify variants, optionally subset by other existing variants.
   Want to retrieve:
    - Total counts in current.
    - High frequency in other comparison not seen in current."
  [config coords vrn-coords ref-file params]
  (let [vcs (mapcat #(sample-variants config % ref-file
                                      (-> params (dissoc :maf-thresh) (dissoc :af-thresh)))
                    coords)]
    {:total (count vcs)
     :unique (count (remove #(contains? vrn-coords (select-keys % [:chr :start])) vcs))}))

; ## Framework

(defn- high-freq-variant-coords
  "Retrieve coordinates of high frequency variants in the given exon regions."
  [config coords ref-file params]
  (->> coords
       (mapcat #(sample-variants config % ref-file params))
       (map #(select-keys % [:chr :start]))
       set))

(defn- compare-region
  "Compare coverage based metrics and variations in a set of gene exons."
  [wgs-config exome-config coords ref-file params]
  (with-open [wgs-c-get (gene/get-coverage-retriever (map :coverage (:samples wgs-config)) params)
              exome-c-get (gene/get-coverage-retriever (map :coverage (:samples exome-config)) params)]
    (let [wgs-vrns (high-freq-variant-coords wgs-config coords ref-file params)
          exome-vrns (high-freq-variant-coords exome-config coords ref-file params)]
      {:wgs
       {:coverage (coverage-report wgs-c-get coords exome-vrns params)
        :variant (variant-report wgs-config coords exome-vrns ref-file params)}
       :exome
       {:coverage (coverage-report exome-c-get coords wgs-vrns params)
        :variant (variant-report exome-config coords wgs-vrns ref-file params)}})))

(defn- do-compare
  "Compare a WGS and exome experiment in specific gene regions of interest."
  [wgs-config exome-config region-file ref-file params]
  (let [gene-coord-file (gene/get-coord-bed region-file params)]
    (with-open [rdr (io/reader gene-coord-file)]
      (->> (bed/get-iterator rdr)
           (group-by :name)
           (map second)
           (map #(compare-region wgs-config exome-config % ref-file params))))))

(defn- cmp->csv
  "Convert coverage comparison information into flattened CSV output."
  [cmp]
  (map #(get-in cmp %)
       [[:wgs :coverage :name]
        [:wgs :coverage :size]
        [:wgs :coverage :avg-reads]
        [:exome :coverage :avg-reads]
        [:wgs :coverage :percent-nocoverage]
        [:exome :coverage :percent-nocoverage]
        [:wgs :coverage :variant-compare :covered]
        [:exome :coverage :variant-compare :covered]
        [:wgs :coverage :variant-compare :uncovered]
        [:exome :coverage :variant-compare :uncovered]
        [:wgs :variant :unique]
        [:exome :variant :unique]]))

(defn- write-to-csv
  "Flatten comparisons between WGS and exome into CSV output, one line per gene."
  [cmps out-file]
  (let [header ["gene", "size" "wgs_avg_reads" "exome_avg_reads" "wgs_pct_nocov"
                "exome_pct_nocov", "wgs_highaf_var_shared" "exome_highaf_var_shared"
                "wgs_highaf_var_unique" "exome_highaf_var_unique"
                "wgs_var_unique" "exome_var_unique"]]
    (with-open [wtr (io/writer out-file)]
      (csv/write-csv wtr (cons header (map cmp->csv cmps))))))

(defn- cmps->dataset
  "Convert comparison information into an Incanter dataset with top differences labeled."
  [cmps params]
  (let [[top-inputs main-inputs]
        (->> cmps
             (map (fn [x]
                    {:gene (get-in x [:wgs :coverage :name])
                     :wgs (get-in x [:wgs :coverage :percent-nocoverage])
                     :exome (get-in x [:exome :coverage :percent-nocoverage])}))
             (sort-by #(Math/abs (- (:wgs %) (:exome %))))
             (split-at (get params :top-diffs 10)))]
    {:ds (icore/to-dataset (concat
                        (map #(assoc % :top true) top-inputs)
                        (map #(assoc % :top false) main-inputs)))
     :top top-inputs}))

(defn- plot-top-differences
  "Identify variants with top coverage differences between two samples.
   Plot a scatter plot of coverage with top changed genes labeled."
  [cmps csv-file params]
  (let [out-file (str (itx/file-root csv-file) ".png")
        cmp-ready (cmps->dataset cmps params)
        plot (icharts/scatter-plot :wgs :exome :data (:ds cmp-ready) :group-by :top
                                   :x-label "WGS" :y-label "exome")]
    (doseq [x (:top cmp-ready)]
      (icharts/add-text plot (:wgs x) (- (:exome x) 3) (:gene x)))
    (icore/save plot out-file)
    out-file))

(defn compare-from-config
  "Perform comparison reading input information from YAML configuration file."
  ([config-file out-file cores]
     (let [config (-> config-file slurp yaml/parse-string
                      (assoc-in [:params :cores] cores))
           [wgs-config exome-config] (take 2 (:experiments config))
           cmps (do-compare wgs-config exome-config (:regions config) (:ref-file config)
                            (:params config))]
       (write-to-csv cmps out-file)
       (plot-top-differences cmps out-file (:params config)))
     out-file)
  ([config-file out-file]
     (compare-from-config config-file out-file 1)))

(defn -main [& args]
  (let [[options args banner]
        (cli args
             ["-c" "--cores" "Number of cores to use" :default 1
              :parse-fn #(Integer. %)])]
    (when (not= (count args) 2)
      (println banner)
      (println "Usage: wgsexome <config file> <output CSV file>")
      (System/exit 1))
    (compare-from-config (first args) (second args) (:cores options))))
