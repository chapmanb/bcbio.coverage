(ns bcbio.coverage.workflow.wgsexome
  "Compare whole genome sequencing data to exome data, assessing coverage and variant calls."
  (:require [clojure.java.io :as io]
            [clj-yaml.core :as yaml]
            [bcbio.coverage.gene :as gene]
            [bcbio.coverage.io.bed :as bed]
            [bcbio.coverage.source.esp :as esp]
            [bcbio.variation.variantcontext :as gvc]))

;; ## Reports

(defn- coverage-report
  "Provide summary of coverage for a set of coordinates.
   Calculates:
   - Summary metrics of coverage in gene.
   - High frequency variant calls in comparison without coverage in current."
  [cov-get coords params]
  (gene/gene-problem-coverage cov-get coords params))

(defmulti sample-variants
  "Retrieve variants from different input types for a group of samples."
  (fn [config _ _ _]
    (keyword (get config :type "vcf"))))

(defmethod sample-variants :esp
  ^{:doc "Retrieval from summarized ESP variant files."}
  [config region ref-file params]
  (esp/variants-in-region (:variant config) ref-file region (get :maf-thresh params 0.0)))

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
      (vec (filter #(vc-has-sample? samples % 0.0)
                   (gvc/variants-in-region vcf-get region))))))

(defn- variant-report
  "Identify variants, optionally subset by other existing variants.
   Want to retrieve:
    - Total counts in current.
    - High frequency in other comparison not seen in current."
  [config coords ref-file params]
  (let [vcs (mapcat #(sample-variants config % ref-file params) coords)]
    (println (count vcs))))

; ## Framework

(defn- compare-region
  "Compare coverage based metrics and variations in a set of gene exons."
  [wgs-config exome-config coords ref-file params]
  (let [wgs-c-get (gene/get-coverage-retriever (map :coverage (:samples wgs-config)))
        exome-c-get (gene/get-coverage-retriever (map :coverage (:samples exome-config)))]
    (coverage-report wgs-c-get coords params)
    (coverage-report exome-c-get coords params)
    (variant-report wgs-config coords ref-file params)
    (variant-report exome-config coords ref-file params)))

(defn- do-compare
  "Compare a WGS and exome experiment in specific gene regions of interest."
  [wgs-config exome-config region-file ref-file params]
  (let [gene-coord-file (gene/get-coord-bed region-file params)]
    (with-open [rdr (io/reader gene-coord-file)]
      (doseq [coords (map second (group-by :name (bed/get-iterator rdr)))]
        (compare-region wgs-config exome-config coords ref-file params)))))

(defn compare-from-config
  "Perform comparison reading input information from YAML configuration file."
  [config-file]
  (let [config (-> config-file slurp yaml/parse-string)
        [wgs-config exome-config] (take 2 (:experiments config))]
    (do-compare wgs-config exome-config (:regions config) (:ref-file config)
                (:params config))))
