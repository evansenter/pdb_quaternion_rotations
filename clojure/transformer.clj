; To use, open a REPL and enter the following lines (from the directory where this file resides)
; 
; (load-file "transformer.clj")
; (in-ns 'spatial-transformations)
; (def data1 (loadData "./../json/2EZM_backbone.json"))
; (def data2 (loadData "./../json/2EZM_backbone_1496_1497_20_degree_rotation.json"))
; (def data3 (loadData "./../json/2EZM_backbone_561_562_20_degree_rotation.json"))

(ns spatial-transformations
  (:use [clojure.contrib.json :only (read-json write-json)])
  (:import (java.io FileReader PrintWriter)))

; Import and I/O code.

(defn updateKeys [hash keys convertor]
  (reduce (fn [memo key] (assoc memo key (convertor (memo key)))) hash keys))

(defn parseHash [pdbHash] 
  (-> pdbHash (updateKeys [:x :y :z] #(Float/parseFloat %)) (updateKeys [:atom_id] #(Integer/parseInt %))))

(defn parseJson [fileReader] 
  (map parseHash (read-json fileReader)))
  
(defn writeJson [object]
  (write-json object (new PrintWriter *out*)))
  
(defn loadData [filePath]
  (filter #(re-find #"^(C|CA|N)$" (% :atom_name)) (parseJson (new FileReader filePath))))
  
(defn pointFromHash [pdbHash]
  (map #(% pdbHash) [:x :y :z]))
  
(defn quadPoints [array]
  (let [quads (partition 4 1 array)]
    (map #(map pointFromHash %) quads)))
    
(defn dataAfterAtomId [array atomId]
  (drop-while #(< (% :atom_id) atomId) array))

; Spatial helper functions.

(defn conjugate [array]
  (cons (first array) (map (partial * -1) (rest array))))

(defn magnitude [array] 
  (Math/sqrt (apply + (map #(Math/pow % 2) array))))

(defn dot [& arrays] 
  (apply + (apply map * arrays)))

(defn cross [[a b c] [d e f]] 
  [(- (* b f) (* c e)) (- (* c d) (* a f)) (- (* a e) (* b d))])

(defn normalize [array]
  (map #(/ % (magnitude array)) array))
  
(defn centerPoint [array]
  (map #(/ % (count array)) (reduce #(map + %1 %2) array)))

(defn toRadians [degrees]
  (* degrees (/ Math/PI 180)))

(defn toDegrees [radians]
  (* radians (/ 180 Math/PI)))

(defn vectorFromPoints [& arrays] 
  (apply map - (reverse arrays)))

(defn normalFromPoints [a b c]
  (cross (vectorFromPoints a b) (vectorFromPoints b c)))
  
(defn distance [& arrays]
  (->> arrays (apply vectorFromPoints) magnitude))
  
; Dihedral angle calculations.

(defn rotationQuaternion [angle array]
  (cons (Math/cos (/ (toRadians angle) 2)) (map (partial * (Math/sin (/ (toRadians angle) 2))) (normalize array))))

(defn quaternionProduct [[a1 & v1] [a2 & v2]]
  (cons 
    (- (* a1 a2) (dot v1 v2)) 
    (map + (map (partial * a1) v2) (map (partial * a2) v1) (cross v1 v2))))

(defn vectorAngle [array1 array2] 
  (toDegrees (Math/acos (/ (dot array1 array2) (* (magnitude array1) (magnitude array2))))))

(defn dihedralAngle [[a b c d]]
  (* 
    (vectorAngle (normalFromPoints a b c) (normalFromPoints b c d)) 
    (if (> (dot (normalFromPoints a b c) (vectorFromPoints c d)) 0) 1 -1)))
    
(defn dihedralAnglesFromHash [array]
  (map dihedralAngle (quadPoints array)))
  
(defn dihedralAnglesFromPoints [array]
  (map dihedralAngle (partition 4 1 array)))

; Translations and rotations.

(defn translate [array transform] 
  (map + array transform))

(defn rotatePoint [angle axis1 axis2 point]
  (let [originTranslation (translate point (map (partial * -1) axis1))]
    (let [rotationArray (rotationQuaternion angle (vectorFromPoints axis1 axis2))]
      (let [rotatedPoint (quaternionProduct (quaternionProduct rotationArray (cons 0 originTranslation)) (conjugate rotationArray))]
        (translate (rest rotatedPoint) axis1)))))
        
(defn rotateAfter [array angle atomId]
  (let [[axis1 axis2] (map pointFromHash (take 2 (dataAfterAtomId array atomId)))]
    (->> (dataAfterAtomId array atomId) (map pointFromHash) (map (partial rotatePoint angle axis1 axis2)))))
    
(defn rotatedCoords [array angle atomId]
  (let [rotatedSection (rotateAfter array angle atomId)]
    (concat (map pointFromHash (reverse (drop (count rotatedSection) (reverse array)))) rotatedSection)))

(defn jsonFormatRotation [array angle atomId]
  (writeJson (zipmap (take (count (dataAfterAtomId atomId)) (iterate inc atomId)) (rotateAfter array angle atomId))))

(defn rmsd [& arrays]
  (let [distances (apply map distance arrays)]
    (/ (magnitude distances) (Math/sqrt (count distances)))))

; Energy functions.

(defn distantPairs [array]
  (for [atom1 array atom2 array :when 
    (and
      (< (atom1 :atom_id) (atom2 :atom_id)) 
      (not-any? #(= (atom2 :atom_id) (% :atom_id)) (take 4 (drop-while #(< (% :atom_id) (atom1 :atom_id)) array))))]
    [atom1 atom2]))
  
(defn energyMap [array]
  (map (fn [[pdbHash1 pdbHash2]]
    (let [
      atom1 (pointFromHash pdbHash1) 
      atom2 (pointFromHash pdbHash2) 
      atomicDistance (distance atom1 atom2)
      atomicEnergy (if (< atomicDistance 3.4) (/ 100 (Math/pow atomicDistance 2)) 0)
      ] { :atom1 pdbHash1 :atom2 pdbHash2 :atomicDistance atomicDistance :atomicEnergy atomicEnergy })) (distantPairs array)))
      
(defn energy [array]
  (apply + (map :atomicEnergy (energyMap array))))
  
(defn energyOver [array threshold]
  (let [highEnergyPairs (filter #(> (% :atomicEnergy) threshold) (energyMap array))]
    (map #(updateKeys % [:atom1 :atom2] :atom_id) highEnergyPairs)))