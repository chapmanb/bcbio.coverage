(ns bcbio.align.remap
  "Remap chromosome names between UCSC and GRC/Ensembl/UCSC conventions"
  (:require [clojure.set :refer [map-invert]]
            [clojure.string :as string]
            [bcbio.align.gref :as gref]))

(def ^{:doc "Cached maps of UCSC to Ensembl coordinate conversion"
       :private true}
  base-ucsc-maps
  {:GRCm38
   {"chrM" "MT", "chrUn_GL456359" "GL456359.1", "chrUn_GL456239"
   "GL456239.1", "chr7_GL456219_random" "GL456219.1",
   "chr1_GL456213_random" "GL456213.1", "chr4_GL456216_random"
   "GL456216.1", "chr1" "1", "chr1_GL456212_random" "GL456212.1",
   "chr2" "2", "chr1_GL456211_random" "GL456211.1",
   "chr1_GL456221_random" "GL456221.1", "chr1_GL456210_random"
   "GL456210.1", "chr3" "3", "chr4" "4", "chr5" "5", "chr6" "6",
   "chr5_GL456354_random" "GL456354.1", "chr7" "7", "chr8" "8", "chrX"
   "X", "chr9" "9", "chrY" "Y", "chr5_JH584299_random" "JH584299.1",
   "chrY_JH584303_random" "JH584303.1", "chrX_GL456233_random"
   "GL456233.1", "chr10" "10", "chr4_GL456350_random" "GL456350.1",
   "chr5_JH584298_random" "JH584298.1", "chrY_JH584302_random"
   "JH584302.1", "chr11" "11", "chr5_JH584297_random" "JH584297.1",
   "chrY_JH584301_random" "JH584301.1", "chr12" "12",
   "chrY_JH584300_random" "JH584300.1", "chr4_JH584295_random"
   "JH584295.1", "chr5_JH584296_random" "JH584296.1", "chr13" "13",
   "chr14" "14", "chr4_JH584294_random" "JH584294.1",
   "chr4_JH584293_random" "JH584293.1", "chr15" "15",
   "chr4_JH584292_random" "JH584292.1", "chr16" "16", "chrUn_GL456390"
   "GL456390.1", "chr17" "17", "chr18" "18", "chrUn_GL456370"
   "GL456370.1", "chrUn_GL456381" "GL456381.1", "chrUn_GL456392"
   "GL456392.1", "chr19" "19", "chrUn_GL456360" "GL456360.1",
   "chrUn_GL456382" "GL456382.1", "chrUn_GL456393" "GL456393.1",
   "chrUn_GL456372" "GL456372.1", "chrUn_GL456383" "GL456383.1",
   "chrUn_GL456394" "GL456394.1", "chrUn_GL456385" "GL456385.1",
   "chrUn_GL456396" "GL456396.1", "chrUn_GL456387" "GL456387.1",
   "chrUn_GL456366" "GL456366.1", "chrUn_JH584304" "JH584304.1",
   "chrUn_GL456367" "GL456367.1", "chrUn_GL456378" "GL456378.1",
   "chrUn_GL456389" "GL456389.1", "chrUn_GL456368" "GL456368.1",
   "chrUn_GL456379" "GL456379.1"}
   :GRCh37
  {"chrM" "MT" "chrMT" "MT" "chrUn_gl000211" "GL000211", "chrUn_gl000222" "GL000222",
   "chrUn_gl000233" "GL000233", "chrUn_gl000244" "GL000244", "chrUn_gl000212" "GL000212",
   "chrUn_gl000223" "GL000223", "chrUn_gl000234" "GL000234", "chrUn_gl000245" "GL000245",
   "chrUn_gl000213" "GL000213", "chrUn_gl000224" "GL000224", "chrUn_gl000235" "GL000235",
   "chrUn_gl000246" "GL000246", "chr6_mcf_hap5" "HSCHR6_MHC_MCF", "chrUn_gl000214" "GL000214",
   "chrUn_gl000225" "GL000225", "chrUn_gl000236" "GL000236", "chrUn_gl000247" "GL000247",
   "chr1" "1", "chr6_cox_hap2" "HSCHR6_MHC_COX", "chrUn_gl000215" "GL000215",
   "chrUn_gl000226" "GL000226", "chrUn_gl000237" "GL000237", "chrUn_gl000248" "GL000248",
   "chr2" "2", "chrUn_gl000216" "GL000216", "chrUn_gl000227" "GL000227",
   "chrUn_gl000238" "GL000238", "chrUn_gl000249" "GL000249", "chr3" "3",
   "chrUn_gl000217" "GL000217", "chrUn_gl000228" "GL000228", "chrUn_gl000239" "GL000239",
   "chr9_gl000201_random" "GL000201", "chr4" "4", "chr11_gl000202_random" "GL000202",
   "chrUn_gl000218" "GL000218", "chrUn_gl000229" "GL000229", "chr9_gl000200_random" "GL000200",
   "chr19_gl000209_random" "GL000209", "chr5" "5", "chrUn_gl000219" "GL000219",
   "chr1_gl000192_random" "GL000192", "chr18_gl000207_random" "GL000207", "chr6" "6",
   "chr21_gl000210_random" "GL000210", "chr17_gl000206_random" "GL000206",
   "chr9_gl000199_random" "GL000199", "chr1_gl000191_random" "GL000191",
   "chr4_gl000194_random" "GL000194", "chr19_gl000208_random" "GL000208",
   "chr17_gl000205_random" "GL000205", "chr7" "7", "chr9_gl000198_random" "GL000198",
   "chr8_gl000197_random" "GL000197", "chr4_gl000193_random" "GL000193",
   "chr17_gl000204_random" "GL000204", "chr8" "8", "chrX" "X", "chr8_gl000196_random" "GL000196",
   "chr7_gl000195_random" "GL000195", "chr20" "20", "chr9" "9", "chrY" "Y",
   "chr17_gl000203_random" "GL000203", "chr10" "10", "chr21" "21", "chr6_dbb_hap3" "HSCHR6_MHC_DBB",
   "chr11" "11", "chr22" "22", "chr6_ssto_hap7" "HSCHR6_MHC_SSTO", "chr17_ctg5_hap1" "HSCHR17_1",
   "chr12" "12", "chr13" "13", "chr14" "14", "chr15" "15", "chr16" "16",
   "chr6_mann_hap4" "HSCHR6_MHC_MANN", "chr17" "17", "chr18" "18", "chr19" "19",
   "chr6_qbl_hap6" "HSCHR6_MHC_QBL", "chr6_apd_hap1" "HSCHR6_MHC_APD",
   "chrUn_gl000240" "GL000240", "chrUn_gl000230" "GL000230", "chrUn_gl000241" "GL000241",
   "chr4_ctg9_hap1" "HSCHR4_1", "chrUn_gl000220" "GL000220", "chrUn_gl000231" "GL000231",
   "chrUn_gl000242" "GL000242", "chrUn_gl000221" "GL000221", "chrUn_gl000232" "GL000232",
   "chrUn_gl000243" "GL000243"}})

(def ucsc-maps
  (-> base-ucsc-maps
      (assoc :hg19 (map-invert (:GRCh37 base-ucsc-maps)))
      (assoc :mm10 (map-invert (:GRCm38 base-ucsc-maps)))))

(defn- fix-non-version-names
  "Convert any non-versioned names into the representative version in ref-dict."
  [base-map ref-dict]
  (letfn [(find-best-match [x check]
            (first (filter #(.startsWith % x) check)))]
    (reduce (fn [coll [k v]]
              (assoc coll k
                     (if (contains? ref-dict v)
                       v
                       (find-best-match v (keys ref-dict)))))
            {} base-map)))

(defn- add-alt-keys
  "Add alternative key variations:
    - underscore to dash in UCSC names
    - chr added to all Ensembl names instead of UCSC names"
  [base-map modtype]
  {:pre [(= modtype :underscore)]}
  (reduce (fn [coll [k v]]
            (-> coll
                (assoc k v)
                (assoc (string/replace k "_" "-") v)
                (assoc (str "chr" v) v)))
          {} base-map))

(defn prep-rename-map
  "Fix Ensembl/UCSC name mappings to handle common problem cases."
  [map-key ref-file]
  (-> (get ucsc-maps map-key)
      (fix-non-version-names (gref/get-seq-name-map ref-file))
      (add-alt-keys :underscore)))

(defn- chrs-from-fasta-file
  "Retrieve a list of all chromosome names from a reference FASTA file."
  [ref-file]
  (map #(.getSequenceName %) (-> ref-file gref/get-seq-dict .getSequences)))

(defn chrom-name-remap
  "Remap chromosome names from UCSC to Ensembl names"
  [map-key ref-file orig-ref-file]
  (let [rename-map (prep-rename-map map-key ref-file)
        ref-chrs (set (chrs-from-fasta-file ref-file))
        vcf-chrs (when (and orig-ref-file (not= orig-ref-file ref-file))
                   (chrs-from-fasta-file orig-ref-file))]
    (reduce (fn [coll x]
              (let [remap-x (get coll x)]
                (if (and remap-x (contains? ref-chrs remap-x))
                  (assoc coll x remap-x)
                  (assoc coll x x))))
            rename-map vcf-chrs)))

(defn match-to-ref
  "Match coordinate contigs to a reference file, ensuring UCSC/Ensembl compatibility."
  [coords build ref-file]
  {:pre [(not (nil? build))]}
  (let [rename-map (prep-rename-map (keyword build) ref-file)]
    (letfn [(is-patch? [coord]
              (.endsWith (:chr coord) "_PATCH"))
            (fix-coord [coord]
              (cond
               (is-patch? coord) nil
               (contains? rename-map (:chr coord)) (assoc coord :chr (get rename-map (:chr coord)))
               :else coord))]
      (remove nil? (map fix-coord coords)))))
