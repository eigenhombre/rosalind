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
            (* 100 (float (gc-content s))))
          
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
