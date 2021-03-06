# Takes a location to a PDB file, and outputs that PDB file with *just* the N, CA and C atoms included.

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

print File.read(ARGV.first).split(/\n/).map(&:strip).reject { |line| 
  line =~ /^ATOM/ && " #{line}"[PDB_COLUMNS[:atom_name]].strip !~ /^(C|N|CA)$/
}.map(&:strip).join("\n")