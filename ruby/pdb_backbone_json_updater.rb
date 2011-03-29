# Takes two arguments, the first is the path to the JSON transformation file and the second is to the PDB file to update.

require "json"

PDB_COLUMNS = {
  record_name:       1..6,
  atom_id:           7..11,
  atom_name:         13..16,
  alt_location:      17..17,
  residue_name:      18..20,
  chain_id:          22..22,
  residue_id:        23..26,
  residue_insertion: 27..27,
  x:                 31..38,
  y:                 39..46,
  z:                 47..54,
  occupancy:         55..60,
  temp_factor:       61..66,
  segment_id:        73..76,
  element:           77..78,
  charge:            79..80
}

json_file, pdb_file = ARGV
transformation      = JSON.parse(File.read(json_file))

print File.read(pdb_file).split(/\n/).map(&:strip).map { |line|
  if line =~ /^ATOM/ && key = transformation.keys.find { |json_key| json_key == " #{line}"[PDB_COLUMNS[:atom_id]].strip }
    " #{line}".tap do |updated_line|
      [:x, :y, :z].zip(transformation[key]).each do |axis, value|
        updated_line[PDB_COLUMNS[axis]] = (("%.3f" % value).reverse + (" " * PDB_COLUMNS[axis].to_a.size)).reverse[-PDB_COLUMNS[axis].to_a.size..-1]
      end
    end.strip
  else
    line
  end
}.join("\n")