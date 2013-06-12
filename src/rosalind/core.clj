(ns rosalind.core
  (:require [clojure.string :as str]))


;; Introduction to the Bioinformatics Armory

(defn freqs [strand]
  (for [sym [\A \C \G \T]]
    ((frequencies strand) sym)))


;; Introduction to Protein Databases

(defn processes-from-db [protein]
  (->> (str "http://www.uniprot.org/uniprot/" protein ".txt")
       slurp
       (re-seq #"(?m)DR\s+GO;\s+GO:\d+; P:(.+?)\;")
       (map second)))

(doseq [s (processes-from-db "B1IMN7")]
  (println s))


;; Transcribing DNA into RNA

(defn transcribe [string]
  (str/replace string "T" "U"))


;; Complementing a Strand of DNA

(defn revcomp [string]
  (let [rmap {\A \T
              \T \A
              \C \G
              \G \C}]
    (apply str (map rmap (reverse string)))))


;; Computing GC Content

(defn max-gc-content [fasta]
  (letfn [(gc-ratio-exact [s]
            (/ (->> s
                    (filter #{\G \C})
                    count)
               (count s)))
          
          (gc-content-percent [s]
            (* 100 (float (gc-ratio-exact s))))
          
          (get-gc-content [s]
            (let [lines (str/split s #"\n")
                  name (first lines)
                  string (apply str (apply concat (rest lines)))]
              [name (gc-content-percent string)]))]
    
    (let [record (fn [s])]
      (apply (partial max-key second)
             (map get-gc-content
                  (rest (str/split fasta #">")))))))


;; Counting Point Mutations

(defn hamming-cp-mutations [a b]
  (count (filter false? (map = a b))))


;;; Translating RNA into Protein

(defmacro deftable [tname & rest]
  `(def ~tname (apply hash-map '(~@rest))))

(deftable proteins
  UUU F      CUU L      AUU I      GUU V
  UUC F      CUC L      AUC I      GUC V
  UUA L      CUA L      AUA I      GUA V
  UUG L      CUG L      AUG M      GUG V
  UCU S      CCU P      ACU T      GCU A
  UCC S      CCC P      ACC T      GCC A
  UCA S      CCA P      ACA T      GCA A
  UCG S      CCG P      ACG T      GCG A
  UAU Y      CAU H      AAU N      GAU D
  UAC Y      CAC H      AAC N      GAC D
  UAA Stop   CAA Q      AAA K      GAA E
  UAG Stop   CAG Q      AAG K      GAG E
  UGU C      CGU R      AGU S      GGU G
  UGC C      CGC R      AGC S      GGC G
  UGA Stop   CGA R      AGA R      GGA G
  UGG W      CGG R      AGG R      GGG G) 

(defn to-protein [s]
  (->> s
       (partition 3)
       (map (partial apply str))
       (map symbol)
       (map proteins)
       (take-while #(not= % 'Stop))
       (apply str)))


;; Finding a motif in DNA

(defn motif-indices [target motif]
  (let [m (seq motif)
        parts (partition (count m) 1 target)]
    (for [i (range (count parts)) :when (= (nth parts i) m)] (inc i))))

(motif-indices "GATATATGCATATACTT" "ATAT") ;=> 2 4 10


;; Rabbits and Recurrence Relations

(defn rabbits [extra-months k]
  (loop [months (- extra-months 2)
         pairs-last-month 1
         pairs-this-month 1]
    (if (= months 0)
      pairs-this-month
      (recur (dec months)
             pairs-this-month
             (+' pairs-this-month (*' pairs-last-month k))))))

(rabbits 5 3) ;=> 19
