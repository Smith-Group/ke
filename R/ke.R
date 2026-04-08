#' Convert between one and three letter residue names 
#'
#' @return Named character vector mapping one-letter and three-letter amino-acid
#'   residue codes.
#'
#' @docType data
#' @rdname residue_names
#' @examples
#' aa1_to_aa3["W"]
#' aa1_to_aa3[c("W", "F", "Y")]
#'
#' @export
aa1_to_aa3 <- c("TRP", "PHE", "TYR", "MET", "LEU", "ILE", "VAL", "ALA", "GLY", "SER", "THR", "ARG", "LYS", "HIS", "ASN", "GLN", "ASP", "GLU", "PRO", "CYS")

#' @docType data
#' @rdname residue_names
#' @examples
#' aa3_to_aa1["TRP"]
#' aa3_to_aa1[c("TRP", "PHE", "TYR")]
#'
#' @export
aa3_to_aa1 <- c("W", "F", "Y", "M", "L", "I", "V", "A", "G", "S", "T", "R", "K", "H", "N", "Q", "D", "E", "P", "C")

names(aa1_to_aa3) <- aa3_to_aa1
names(aa3_to_aa1) <- aa1_to_aa3

pdb_atomres <-    c(" H   ALA", " HA  ALA", " HB1 ALA", " HB2 ALA", " HB3 ALA", " H   ARG", " HA  ARG", " HB2 ARG", " HB3 ARG", " HD2 ARG", " HD3 ARG", " HE  ARG", " HG2 ARG", " HG3 ARG", "HH11 ARG", "HH12 ARG", "HH21 ARG", "HH22 ARG", " H   ASN", " HA  ASN", " HB2 ASN", " HB3 ASN", "HD21 ASN", "HD22 ASN", " H   ASP", " HA  ASP", " HB2 ASP", " HB3 ASP", " H   CYS", " HA  CYS", " HB2 CYS", " HB3 CYS", " HG  CYS", " H   GLN", " HA  GLN", " HB2 GLN", " HB3 GLN", " HG2 GLN", " HG3 GLN", "HE21 GLN", "HE22 GLN", " H   GLU", " HA  GLU", " HB2 GLU", " HB3 GLU", " HG2 GLU", " HG3 GLU", " H   GLY", " HA2 GLY", " HA3 GLY", " H   HIS", " HA  HIS", " HB2 HIS", " HB3 HIS", " HD1 HIS", " HD2 HIS", " HE1 HIS", " HE2 HIS", " H   ILE", " HA  ILE", " HB  ILE", " CD1 ILE", "HD11 ILE", "HD12 ILE", "HD13 ILE", "HG12 ILE", "HG13 ILE", "HG21 ILE", "HG22 ILE", "HG23 ILE", " H   LEU", " HA  LEU", " HB2 LEU", " HB3 LEU", " HG  LEU", "HD11 LEU", "HD12 LEU", "HD13 LEU", "HD21 LEU", "HD22 LEU", "HD23 LEU", " H   LYS", " HA  LYS", " HB2 LYS", " HB3 LYS", " HD2 LYS", " HD3 LYS", " HE2 LYS", " HE3 LYS", " HG2 LYS", " HG3 LYS", " HZ1 LYS", " HZ2 LYS", " HZ3 LYS", " H1  MET", " H2  MET", " H3  MET", " HA  MET", " HB2 MET", " HB3 MET", " HE1 MET", " HE2 MET", " HE3 MET", " HG2 MET", " HG3 MET", " H   PHE", " HA  PHE", " HB2 PHE", " HB3 PHE", " HD1 PHE", " HD2 PHE", " HE1 PHE", " HE2 PHE", " HZ  PHE", " HA  PRO", " HB2 PRO", " HB3 PRO", " HD2 PRO", " HD3 PRO", " HG2 PRO", " HG3 PRO", " H   SER", " HA  SER", " HB2 SER", " HB3 SER", " HG  SER", " H   THR", " HA  THR", " HB  THR", " HG1 THR", "HG21 THR", "HG22 THR", "HG23 THR", " H   TRP", " HA  TRP", " HB2 TRP", " HB3 TRP", " HD1 TRP", " HE1 TRP", " HE3 TRP", " HZ2 TRP", " HZ3 TRP", " HH2 TRP", " H   TYR", " HA  TYR", " HB2 TYR", " HB3 TYR", " HD1 TYR", " HD2 TYR", " HE1 TYR", " HE2 TYR", " HH  TYR", " H   VAL", " HA  VAL", " HB  VAL", "HG11 VAL", "HG12 VAL", "HG13 VAL", "HG21 VAL", "HG22 VAL", "HG23 VAL")
amber_atomres <-  c(" H   ALA", " HA  ALA", " HB1 ALA", " HB2 ALA", " HB3 ALA", " H   ARG", " HA  ARG", " HB1 ARG", " HB2 ARG", " HD1 ARG", " HD2 ARG", " HE  ARG", " HG1 ARG", " HG2 ARG", "1HH1 ARG", "2HH1 ARG", "1HH2 ARG", "2HH2 ARG", " H   ASN", " HA  ASN", " HB1 ASN", " HB2 ASN", "1HD2 ASN", "2HD2 ASN", " H   ASP", " HA  ASP", " HB1 ASP", " HB2 ASP", " H   CYS", " HA  CYS", " HB1 CYS", " HB2 CYS", " HG  CYS", " H   GLN", " HA  GLN", " HB1 GLN", " HB2 GLN", " HG1 GLN", " HG2 GLN", "1HE2 GLN", "2HE2 GLN", " H   GLU", " HA  GLU", " HB1 GLU", " HB2 GLU", " HG1 GLU", " HG2 GLU", " H   GLY", " HA1 GLY", " HA2 GLY", " H   HIS", " HA  HIS", " HB1 HIS", " HB2 HIS", " HD1 HIS", " HD2 HIS", " HE1 HIS", " HE2 HIS", " H   ILE", " HA  ILE", " HB  ILE", " CD  ILE", " HD1 ILE", " HD2 ILE", " HD3 ILE", "1HG1 ILE", "2HG1 ILE", "1HG2 ILE", "2HG2 ILE", "3HG2 ILE", " H   LEU", " HA  LEU", " HB1 LEU", " HB2 LEU", " HG  LEU", "1HD1 LEU", "2HD1 LEU", "3HD1 LEU", "1HD2 LEU", "2HD2 LEU", "3HD2 LEU", " H   LYS", " HA  LYS", " HB1 LYS", " HB2 LYS", " HD1 LYS", " HD2 LYS", " HE1 LYS", " HE2 LYS", " HG1 LYS", " HG2 LYS", " HZ1 LYS", " HZ2 LYS", " HZ3 LYS", " H1  MET", " H2  MET", " H3  MET", " HA  MET", " HB1 MET", " HB2 MET", " HE1 MET", " HE2 MET", " HE3 MET", " HG1 MET", " HG2 MET", " H   PHE", " HA  PHE", " HB1 PHE", " HB2 PHE", " HD1 PHE", " HD2 PHE", " HE1 PHE", " HE2 PHE", " HZ  PHE", " HA  PRO", " HB1 PRO", " HB2 PRO", " HD1 PRO", " HD2 PRO", " HG1 PRO", " HG2 PRO", " H   SER", " HA  SER", " HB1 SER", " HB2 SER", " HG  SER", " H   THR", " HA  THR", " HB  THR", " HG1 THR", "1HG2 THR", "2HG2 THR", "3HG2 THR", " H   TRP", " HA  TRP", " HB1 TRP", " HB2 TRP", " HD1 TRP", " HE1 TRP", " HE3 TRP", " HZ2 TRP", " HZ3 TRP", " HH2 TRP", " H   TYR", " HA  TYR", " HB1 TYR", " HB2 TYR", " HD1 TYR", " HD2 TYR", " HE1 TYR", " HE2 TYR", " HH  TYR", " H   VAL", " HA  VAL", " HB  VAL", "1HG1 VAL", "2HG1 VAL", "3HG1 VAL", "1HG2 VAL", "2HG2 VAL", "3HG2 VAL")
charmm_atomres <- c(" HN  ALA", " HA  ALA", " HB1 ALA", " HB2 ALA", " HB3 ALA", " HN  ARG", " HA  ARG", " HB1 ARG", " HB2 ARG", " HD1 ARG", " HD2 ARG", " HE  ARG", " HG1 ARG", " HG2 ARG", "HH11 ARG", "HH12 ARG", "HH21 ARG", "HH22 ARG", " HN  ASN", " HA  ASN", " HB1 ASN", " HB2 ASN", "HD21 ASN", "HD22 ASN", " HN  ASP", " HA  ASP", " HB1 ASP", " HB2 ASP", " HN  CYS", " HA  CYS", " HB1 CYS", " HB2 CYS", " HG1 CYS", " HN  GLN", " HA  GLN", " HB1 GLN", " HB2 GLN", " HG1 GLN", " HG2 GLN", "HE21 GLN", "HE22 GLN", " HN  GLU", " HA  GLU", " HB1 GLU", " HB2 GLU", " HG1 GLU", " HG2 GLU", " HN  GLY", " HA1 GLY", " HA2 GLY", " HN  HSP", " HA  HSP", " HB1 HSP", " HB2 HSP", " HD1 HSP", " HD2 HSP", " HE1 HSP", " HE2 HSP", " HN  ILE", " HA  ILE", " HB  ILE", " CD  ILE", " HD1 ILE", " HD2 ILE", " HD3 ILE", "HG11 ILE", "HG12 ILE", "HG21 ILE", "HG22 ILE", "HG23 ILE", " HN  LEU", " HA  LEU", " HB1 LEU", " HB2 LEU", " HG  LEU", "HD11 LEU", "HD12 LEU", "HD13 LEU", "HD21 LEU", "HD22 LEU", "HD23 LEU", " HN  LYS", " HA  LYS", " HB1 LYS", " HB2 LYS", " HD1 LYS", " HD2 LYS", " HE1 LYS", " HE2 LYS", " HG1 LYS", " HG2 LYS", " HZ1 LYS", " HZ2 LYS", " HZ3 LYS", " HT1 MET", " HT2 MET", " HT3 MET", " HA  MET", " HB1 MET", " HB2 MET", " HE1 MET", " HE2 MET", " HE3 MET", " HG1 MET", " HG2 MET", " HN  PHE", " HA  PHE", " HB1 PHE", " HB2 PHE", " HD1 PHE", " HD2 PHE", " HE1 PHE", " HE2 PHE", " HZ  PHE", " HA  PRO", " HB1 PRO", " HB2 PRO", " HD1 PRO", " HD2 PRO", " HG1 PRO", " HG2 PRO", " HN  SER", " HA  SER", " HB1 SER", " HB2 SER", " HG1 SER", " HN  THR", " HA  THR", " HB  THR", " HG1 THR", "HG21 THR", "HG22 THR", "HG23 THR", " HN  TRP", " HA  TRP", " HB1 TRP", " HB2 TRP", " HD1 TRP", " HE1 TRP", " HE3 TRP", " HZ2 TRP", " HZ3 TRP", " HH2 TRP", " HN  TYR", " HA  TYR", " HB1 TYR", " HB2 TYR", " HD1 TYR", " HD2 TYR", " HE1 TYR", " HE2 TYR", " HH  TYR", " HN  VAL", " HA  VAL", " HB  VAL", "HG11 VAL", "HG12 VAL", "HG13 VAL", "HG21 VAL", "HG22 VAL", "HG23 VAL")

#' Convert between proton atom/residue name pairs
#'
#' @return Named character vector mapping atom/residue identifiers between naming
#'   conventions.
#'
#' @rdname atom_residue_names
#' @export
amber_to_pdb <- pdb_atomres
names(amber_to_pdb)	<- amber_atomres

#' @rdname atom_residue_names
#' @export
pdb_to_amber <- amber_atomres
names(pdb_to_amber) <- pdb_atomres

#' @rdname atom_residue_names
#' @export
charmm_to_pdb <- pdb_atomres
names(charmm_to_pdb) <- charmm_atomres

#' @rdname atom_residue_names
#' @export
pdb_to_charmm <- charmm_atomres
names(pdb_to_charmm) <- pdb_atomres

#' Convert PDB formatted atom lines into a matrix
#'
#' Generate a matrix containing coordinates in a set of PDB atom lines
#'
#' The column names in the resulting matrix correspond to columns 13 through 27 of the 
#' ATOM lines. The definitions of the PDB columns and corresponding colnames columns are
#' as follows:
#' \tabular{lllll}{
#'   PDB \tab colnames \tab Data Type \tab Field \tab Definition \cr
#'   13-16 \tab 1-4 \tab Atom \tab name \tab Atom name \cr
#'   17 \tab 5 \tab Character \tab altLoc \tab Alternate location indicator \cr
#'   18-20 \tab 6-8 \tab Residue name \tab resName \tab Residue name \cr
#'   22 \tab 10 \tab Character \tab chainID \tab Chain identifier \cr
#'   23-26 \tab 11-4 \tab Integer \tab resSeq \tab Residue sequence number \cr
#'   27 \tab 15 \tab AChar \tab iCode \tab Code for insertion of residues
#' }
#'
#' See \url{http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM}
#' for more information.
#'
#' @param atom_lines character vector of ATOM lines
#'
#' @return 3xN matrix with a column for every atom
#'
#' @export
atom_lines_to_coord <- function(atom_lines) {

	coord <- rbind(
		as.numeric(substr(atom_lines, 31, 38)),
		as.numeric(substr(atom_lines, 39, 46)),
		as.numeric(substr(atom_lines, 47, 54))
	)
	
	colnames(coord) <- substr(atom_lines, 13, 27)
	
	coord
}

#' Convert a matrix of coordinates into PDB formatted atom lines
#'
#' Make PDB ATOM lines from a 3xN matrix of coordinates
#'
#' The names of the columns are expected to correspond to columns 13 to 27 of the PDB
#' atom records.
#'
#' @param coord 3xN matrix with a column for every atom
#'
#' @return a character vector with the PDB ATOM lines
#'
#' @export
coord_to_atom_lines <- function(coord) {

	sprintf("ATOM  %5i %s   %8.3f%8.3f%8.3f  1.00  0.00", seq_len(ncol(coord)), colnames(coord), coord[1,], coord[2,], coord[3,])
}

#' Get or set four-character atom names
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB ATOM records
#'
#' @return Getter: character vector of four-character atom names. Setter:
#'   modified `coord` with updated atom names.
#'
#' @rdname coord_atomnames
#' @export
coord_atomnames <- function(coord) {

	if (is.array(coord) && is.character(dimnames(coord)[[2]])) {
		substr(dimnames(coord)[[2]], 1, 4)
	} else if (is.character(coord)) {
		substr(coord, 1, 4)
	} else {
		stop("coord not named array or character")
	}
}

#' @param value character vector with four characters per element
#'
#' @rdname coord_atomnames
#' @export
`coord_atomnames<-` <- function(coord, value) {

	stopifnot(nchar(value) == 4)

	if (is.array(coord) && is.character(dimnames(coord)[[2]])) {
		substr(dimnames(coord)[[2]], 1, 4) <- value
	} else if (is.character(coord)) {
		substr(coord, 1, 4) <- value
	} else {
		stop("coord not named array or character")
	}
	
	coord
}

#' Get or set three-character residue names
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB res records
#'
#' @return Getter: character vector of three-letter residue names. Setter:
#'   modified `coord` with updated residue names.
#'
#' @rdname coord_resnames
#' @export
coord_resnames <- function(coord) {

	substr(dimnames(coord)[[2]], 6, 8)
}

#' @param value character vector with three characters per element
#'
#' @rdname coord_resnames
#' @export
`coord_resnames<-` <- function(coord, value) {

	stopifnot(nchar(value) == 4)

	substr(dimnames(coord)[[2]], 6, 8) <- value
	
	coord
}

#' Get or set atom/residue name pairs
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB ATOM records
#'
#' @return Getter: character vector of atom/residue identifiers. Setter:
#'   modified `coord` with updated atom/residue identifiers.
#'
#' @rdname coord_atomresnames
#' @export
coord_atomresnames <- function(coord) {

	atomresnames <- substr(dimnames(coord)[[2]], 1, 8)
	substr(atomresnames, 5, 5) <- " "
	
	atomresnames
}

#' @param value character vector with eight characters per element (1-4: atom name, 6-8: residue name)
#'
#' @rdname coord_atomresnames
#' @export
`coord_atomresnames<-` <- function(coord, value) {

	stopifnot(nchar(value) == 4)

	substr(dimnames(coord)[[2]], 1, 4) <- substr(value, 1, 4)
	substr(dimnames(coord)[[2]], 6, 8) <- substr(value, 6, 8)
	
	coord
}

#' Get or set residue sequence number
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB ATOM records
#'
#' @return Getter: integer vector of residue sequence numbers. Setter:
#'   modified `coord` with updated residue sequence numbers.
#'
#' @rdname coord_resseq
#' @export
coord_resseq <- function(coord) {

	if (is.array(coord) && is.character(dimnames(coord)[[2]])) {
		as.integer(substr(dimnames(coord)[[2]], 11, 14))
	} else if (is.character(coord)) {
		as.integer(substr(coord, 11, 14))
	} else {
		stop("coord not named array or character")
	}
}

#' @param value integer or character vector with no more than four digits
#'
#' @rdname coord_resseq
#' @export
`coord_resseq<-` <- function(coord, value) {

	value <- as.character(value)
	stopifnot(nchar(value) <= 4)

	if (is.array(coord) && is.character(dimnames(coord)[[2]])) {
		substr(dimnames(coord)[[2]], 11, 14) <- sprintf("%-4s", value)
	} else if (is.character(coord)) {
		substr(coord, 11, 14) <- sprintf("%-4s", value)
	} else {
		stop("coord not named array or character")
	}
	
	coord
}

#' Determine whether an atom is rapidly exchangeable at neutral pH
#'
#' @param coord matrix or array with second dimension having columns 13 to 27 of the PDB ATOM records
#' @param amber logical indicating whether atom names are from AMBER force field
#'
#' @return Logical vector with one value per atom identifier in `coord`,
#'   indicating whether the proton is treated as rapidly exchangeable.
#'
#' @export
coord_rapidly_exchangeable <- function(coord, amber=FALSE) {

	rapidly_exchangeable_re <- c(
		"^ HH  TYR", # tyrosine hydroxyl
		"^ HZ[123] LYS", # lysine amino group
		"^ HE  ARG", # arginine epsilon NH
		"^HH[12][12] ARG", # arginine eta NH2
		"^ HG  SER", # serine hydroxyl
		"^ HG1 THR", # threonine hydroxyl
		"^ HD1 HIS", # histidine delta amide hydrogen
		"^ HE2 HIS", # histidine epsilon amide hydrogen
		"^ H[123]  ...", # N-terminal amino group
		"^ HG  CYS" # cysteine hydrogen
	)
	
	if (amber) {
		rapidly_exchangeable_re <- c(
			"^ HH  TYR", # tyrosine hydroxyl
			"^ HZ[123] LYS", # lysine amino group
			"^ HE  ARG", # arginine epsilon NH
			"^[12]HH[12] ARG", # arginine eta NH2
			"^ HG  SER", # serine hydroxyl
			"^ HG1 THR", # threonine hydroxyl
			"^ HD1 HIS", # histidine delta amide hydrogen
			"^ HE2 HIS", # histidine epsilon amide hydrogen
			"^ H[123]  ..." # N-terminal amino group
		)
	}
	
	match_mat <- sapply(rapidly_exchangeable_re, grepl, dimnames(coord)[[2]])
	
	apply(match_mat, 1, any)
}

#' Read an ensemble of PDB files
#'
#' @param pdb_files character vector paths to PDB format file(s)
#' @param model_idx vector with indices of models to read
#' @param proton_only logical indicating whether to only read hydrogen ATOM records
#'
#' @return Numeric coordinate array. For a single input file the dimensions are
#'   `(xyz, atoms, models)`; for multiple input files an additional trailing file
#'   dimension is retained.
#'
#' @export
read_ensemble <- function(pdb_files, model_idx=NULL, proton_only=FALSE) {

	coord_list <- lapply(pdb_files, function(pdb_file) {
	
		pdb_lines <- readLines(pdb_file)
		model_line_idx <- grep("^MODEL", pdb_lines)
		
		if (length(model_line_idx)) {
			start_line_idx <- model_line_idx+1
			end_line_idx <- c(model_line_idx[-1]-1, length(pdb_lines))
			model_idx_list <- lapply(seq_along(start_line_idx), function(i) seq(start_line_idx[i], end_line_idx[i]))
		} else {
			model_idx_list <- list(seq_along(pdb_lines))
		}
		
		if (!is.null(model_idx)) {
			model_idx_list <- model_idx_list[model_idx]
		}
		
		coord_model_list <- lapply(model_idx_list, function(idx) {
			atom_lines <- grep("^ATOM  ..... ....[ A]", pdb_lines[idx], value=TRUE)
			coord <- atom_lines_to_coord(atom_lines)
			substr(colnames(coord), 5, 5) <- " "
			
			coord
		})
		
		coord_model_list
	})
	
	if (is.null(model_idx)) {
		model_idx <- seq_along(coord_list[[1]])
	}
	
	coord_list <- unlist(coord_list, recursive=FALSE)
	
	unique_atoms <- colnames(coord_list[[1]])
	
	for (i in seq_along(coord_list)[-1]) {
		unique_atoms <- intersect(unique_atoms, colnames(coord_list[[i]]))
	}
	
	#unique_atoms <- intersect(unique_atoms, c(exchangeable_atomids(unique_atoms), n_atomids))
	
	if (proton_only) {
		unique_atoms <- grep("^(H|.H)", unique_atoms, value=TRUE)
	}
	
	for (i in seq_along(coord_list)) {
		coord_list[[i]] <- coord_list[[i]][,unique_atoms]
	}
	
	coord_array <- simplify2array(coord_list)
	
	if (dim(coord_array)[3] == length(pdb_files)) {
	
		dimnames(coord_array)[[3]] <- pdb_files
	
	} else {
	
		coord_array <- array(
			coord_array,
			dim=c(dim(coord_array)[1:2], dim(coord_array)[3]/length(pdb_files), length(pdb_files)),
			dimnames=c(dimnames(coord_array)[1:2], list(paste0("M", model_idx), pdb_files))
		)
		
		if (dim(coord_array)[4] == 1) {
			
			new_dim <- dim(coord_array)[-4]
			new_dimnames <- dimnames(coord_array)[-4]
			dim(coord_array) <- new_dim
			dimnames(coord_array) <- new_dimnames
		}
	}
	
	coord_array
}

#' Select groups of protons that are typically detectable with unique chemical shifts
#'
#' @param coord_mat 3xN matrix with a column for every atom
#' @param alpha_group logical indicating whether to group aliphatic alpha protons
#' @param beta_group logical indicating whether to group aliphatic beta protons
#' @param gamma_group logical indicating whether to group aliphatic gamma protons
#' @param delta_group logical indicating whether to group aliphatic delta protons
#' @param epsilon_group logical indicating whether to group aliphatic epsilon protons
#' @param methyl_group logical indicating whether to group methyl protons
#' @param aromatic_group logical indicating whether to group phenylalanine/tyrosine delta/epsilon protons
#' @param amine_group logical indicating whether to group asparagine/glutamine amine protons
#' @param exclude_exchangeable logical indicating whether rapidly exchangeable protons should be excluded
#' @param include character vector of exchangeable protons to be explicitly included
#' @param amber logical indicating if names are from Amber force field
#'
#' @return Vector with selected protons. The names give the group the proton belongs to. The name will be the heavy atom name if all of the protons bound to that heavy atom are in the same group. Otherwise, the name will be the proton name.
#'
#' @export
coord_proton_groups <- function(coord_mat, alpha_group=FALSE, beta_group=FALSE, gamma_group=FALSE, delta_group=FALSE, epsilon_group=FALSE, amine_group=FALSE, methyl_group=FALSE, aromatic_group=TRUE, exclude_exchangeable=TRUE, include=character(), amber=FALSE) {

	heavy_idx <- substr(colnames(coord_mat), 2, 2) %in% c("C","N","O","S")
	
	atom_dist <- as.matrix(stats::dist(t(coord_mat)))
	
	group_names <- colnames(coord_mat)[heavy_idx][apply(atom_dist[!heavy_idx,heavy_idx], 1, which.min)]
	names(group_names) <- colnames(coord_mat)[!heavy_idx]
	
	if (exclude_exchangeable) {
		exclude_idx <- which(coord_rapidly_exchangeable(coord_mat[,!heavy_idx], amber))
		exclude_idx <- setdiff(exclude_idx, which(names(group_names) %in% include))
		if (length(exclude_idx)) {
			group_names <- group_names[-exclude_idx]
		}
	}
	
	if (aromatic_group) {
		group_names <- sub("( CD| CE)[12] (PHE|TYR)", "\\1* \\2", group_names)
	}
	
	if (methyl_group) {
		group_names <- sub("( CD)[12] (LEU)", "\\1* \\2", group_names)
		group_names <- sub("( CG)[12] (VAL)", "\\1* \\2", group_names)
	}
	
	groups_with_nonequivalent_protons <- character()
	
	if (!alpha_group) {
		groups_with_nonequivalent_protons <- c(" CA  GLY")
	}
	
	if (!beta_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " CB  ARG", " CB  ASN", " CB  ASP", " CB  CYS", " CB  GLN", " CB  GLU", " CB  HIS", " CB  LEU", " CB  LYS", " CB  MET", " CB  PHE", " CB  PRO", " CB  SER", " CB  TRP", " CB  TYR")
	}
	
	if (!gamma_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " CG  ARG", " CG  GLN", " CG  GLU", " CG1 ILE", " CG  LYS", " CG  MET", " CG  PRO")
	}
	
	if (!delta_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " CD  ARG", " CD  LYS", " CD  PRO")
	}
	
	if (!epsilon_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " CE  LYS")
	}
	
	if (!amine_group) {
		groups_with_nonequivalent_protons <- c(groups_with_nonequivalent_protons, " ND2 ASN", " NE2 GLN")
	}
	
	nonequivalent_idx <- substr(group_names, 1, 8) %in% groups_with_nonequivalent_protons
	group_names[nonequivalent_idx] <- names(group_names)[nonequivalent_idx]
	
	group_names
}

#' Calculate bond distances between a set of atoms
#'
#' @param coord_mat 3xN matrix with a column for every atom
#' @param depth maximum bond distance to calculate
#'
#' @return NxN matrix with bond distances between atoms (or NA if not calculated)
#'
#' @export
bond_separation <- function(coord_mat, depth=4) {

	bond_sep_mat <- matrix(NA_integer_, nrow=ncol(coord_mat), ncol=ncol(coord_mat), dimnames=list(colnames(coord_mat), colnames(coord_mat)))
	diag(bond_sep_mat) <- 0

	heavy_idx <- substr(colnames(coord_mat), 2, 2) %in% c("C","N","O","S")
	
	atom_dist <- as.matrix(stats::dist(t(coord_mat)))

	heavy_names <- colnames(coord_mat)[heavy_idx][apply(atom_dist[!heavy_idx,heavy_idx], 1, which.min)]
	proton_names <- colnames(coord_mat)[!heavy_idx]

	bond_sep_mat[cbind(heavy_names,proton_names)] <- bond_sep_mat[cbind(proton_names,heavy_names)] <- 1

	heavy_bond_mat <- which(atom_dist[heavy_idx,heavy_idx] > 0 & atom_dist[heavy_idx,heavy_idx] < 2, arr.ind=TRUE)
	heavy_bond_mat <- matrix(colnames(coord_mat)[heavy_idx][heavy_bond_mat], ncol=2)

	bond_sep_mat[heavy_bond_mat] <- 1

	for (i in seq_len(depth-1)) {
		for (j in seq_len(nrow(bond_sep_mat))) {
			bond_idx <- which(bond_sep_mat[j,] == i)
			onebond_idx <- unique(which(bond_sep_mat[,bond_idx,drop=FALSE] == 1, arr.ind=TRUE)[,1])
			nbond_idx <- onebond_idx[is.na(bond_sep_mat[j,onebond_idx])]
			bond_sep_mat[j,nbond_idx] <- i+1
		}
	}

	bond_sep_mat
}

#' Collapse first two dimensions of an array
#'
#' @param arr two or three dimensional array with the first two dimensions having the same length
#' @param idx_list list of character or integer indices to accumulate using func
#' @param func function to apply
#'
#' @return Matrix or array with the first two dimensions collapsed according to
#'   `idx_list`, preserving any trailing dimensions.
#'
#' @export
collapse_array <- function(arr, idx_list, func=`+`) {

	idx_length <- sapply(idx_list, length)
	
	source_idx <- lapply(seq_len(max(idx_length)), function(idx) sapply(idx_list, function(x) x[idx]))
	dest_idx <- lapply(seq_along(source_idx), function(idx) which(!is.na(source_idx[[idx]])))
	source_idx <- lapply(seq_along(source_idx), function(idx) source_idx[[idx]][dest_idx[[idx]]])
	
	stopifnot(length(dim(arr)) %in% c(2L, 3L))
	
	if (length(dim(arr)) == 2) {
		
		arr_rows_collapsed <- arr[source_idx[[1]],,drop=FALSE]
		for (i in seq_along(source_idx)[-1]) {
			arr_rows_collapsed[dest_idx[[i]],] <- func(arr_rows_collapsed[dest_idx[[i]],],arr[source_idx[[i]],])
		}
	
		arr_rows_cols_collapsed <- arr_rows_collapsed[,source_idx[[1]],drop=FALSE]
		for (i in seq_along(source_idx)[-1]) {
			arr_rows_cols_collapsed[,dest_idx[[i]]] <- func(arr_rows_cols_collapsed[,dest_idx[[i]]],arr_rows_collapsed[,source_idx[[i]]])
		}
	
	} else if (length(dim(arr)) == 3) {
		
		arr_rows_collapsed <- arr[source_idx[[1]],,,drop=FALSE]
		for (i in seq_along(source_idx)[-1]) {
			arr_rows_collapsed[dest_idx[[i]],,] <- func(arr_rows_collapsed[dest_idx[[i]],,],arr[source_idx[[i]],,])
		}
	
		arr_rows_cols_collapsed <- arr_rows_collapsed[,source_idx[[1]],,drop=FALSE]
		for (i in seq_along(source_idx)[-1]) {
			arr_rows_cols_collapsed[,dest_idx[[i]],] <- func(arr_rows_cols_collapsed[,dest_idx[[i]],],arr_rows_collapsed[,source_idx[[i]],])
		}
	}
	
	dimnames(arr_rows_cols_collapsed)[1:2] <- list(names(idx_list), names(idx_list))
	
	arr_rows_cols_collapsed
}

#' Check analytical derivatives via finite difference approximation
#'
#' @param func function to evaluate
#' @param value first argument to pass to `func`
#' @param dv delta used to increment elements of `value`
#' @param vdims vector of dimensions that must be incremented separately
#' @param gdims vector of dimensions in gradient corresponding to `vdims`
#' @param ... additional arguments passed to `func`
#'
#' @return difference between `gradient` attribute returned by `func` and finite difference approximation
#'
#' @export
deriv_check <- function(func, value, dv, vdims, gdims, ...) {

	result <- func(value, ..., gradient=TRUE)
	
	# gradient calculated analytically to check
	gradient <- attr(result, "gradient")
	attr(result, "gradient") <- NULL
	
	# gradient calculated using finite difference
	gradient_fd <- gradient
	gradient_fd[] <- NA

	stopifnot(dim(value)[vdims] == dim(gradient)[gdims])

	idx_mat <- as.matrix(expand.grid(lapply(vdims, function(i) seq_len(dim(value)[i]))))
	v_list <- lapply(seq_along(dim(value)), function(i) seq_len(dim(value)[i]))
	g_list <- lapply(seq_along(dim(gradient)), function(i) seq_len(dim(gradient)[i]))

	for (i in seq_len(nrow(idx_mat))) {
	
		v_list[vdims] <- as.list(idx_mat[i,])
		v_idx <- as.matrix(expand.grid(v_list))
		g_list[gdims] <- as.list(idx_mat[i,])
		g_idx <- as.matrix(expand.grid(g_list))
		
		new_value <- value
		new_value[v_idx] <- new_value[v_idx] + dv
		
		new_result <- func(new_value, ..., gradient=FALSE)
		
		gradient_fd[g_idx] <- (new_result-result)/dv
	}
	
	gradient-gradient_fd
}

#' Calculate internuclear vectors from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`)
#'
#' @return 3D array (pairs, models, xyz) with internuclear vectors. Atom pair names follow the format `resSeq:Atom-resSeq:Atom`.
#;
#' @export
coord_array_to_r_array <- function(coord_array, atom_pairs) {

	r_array <- aperm(coord_array[atom_pairs[,2],,,drop=FALSE]-coord_array[atom_pairs[,1],,,drop=FALSE], c(1,3,2))
	
	if (!is.null(dimnames(coord_array))) {
		atom_names <- sub(" *([^ ]+) *([A-Z]{3}) (.) +([^ ]+) *", "\\4:\\1", dimnames(coord_array)[[1]])
		names(atom_names) <- dimnames(coord_array)[[1]]
		dimnames(r_array)[[1]] <- paste(atom_names[atom_pairs[,1]], atom_names[atom_pairs[,2]], sep="-")
		dimnames(r_array)[[3]] <- c("x", "y", "z")
	}
	
	r_array
}

#' Back-propagate energy derivative from r array to coordinates
#'
#' @param d_energy_d_coord_array 3D array (atoms, xyz, models) to accumulate derivatives into
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`
#' @param d_energy_d_r_array 3D array (pairs, models, xyz)
#'
#' @return 3D array (pairs, xyz, models) with d_energy_d_coord_array
coord_array_to_r_array_backprop <- function(d_energy_d_coord_array, atom_pairs, d_energy_d_r_array) {

	#print(str(d_energy_d_coord_array))
	#print(str(d_energy_d_r_array))
	# propagate the internuclear vector derivatives back onto the atomic coordinates
	for (i in seq_len(dim(d_energy_d_r_array)[1])) {
		
		pair_grad <- t(d_energy_d_r_array[i,,])
		# first component is subtracted so derivative is too
		d_energy_d_coord_array[atom_pairs[i,1],,] <- d_energy_d_coord_array[atom_pairs[i,1],,] - pair_grad
		# second component is added so derivative is too
		d_energy_d_coord_array[atom_pairs[i,2],,] <- d_energy_d_coord_array[atom_pairs[i,2],,] + pair_grad
	}
	
	d_energy_d_coord_array
}

#' Calculate dipole-dipole interaction tensors from internuclear vectors
#'
#' @details
#' Given an internuclear vector with components (`x`, `y`, `z`) and distance (`r`), the vector form (`d`) of the dipole-dipole interaction tensor is defined:
#'
#' \eqn{d = 1/r^5 [ z^2 - 1/2 (x^2+y^2), \sqrt{3}/2 (x^2-y^2), \sqrt{3} xz, \sqrt{3} yz, \sqrt{3} xy) ]}
#'
#' Note that the coefficients of the last four elements are different than what was given in Smith 2020 (\doi{10.1007/s10858-019-00288-8}). Those published coefficients do not account for the 3/2 factor in the cartesian tensor.
#'
#' @param r_array 3D array (pairs, models, xyz) with internuclear vectors
#' @param dist logical indicating whether to return distance-dependent form
#' @param unit logical indicating whether to return unit (distance-independent) form
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return 3D array (pairs, models, tensor elements) with interaction tensors
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 4D array (pairs, models, tensor elements, xyz).
#'
#' @export
r_array_to_d_array <- function(r_array, dist=TRUE, unit=FALSE, gradient=FALSE) {

	stopifnot(dist || unit)

	input_dim <- dim(r_array)
	input_dimnames <- dimnames(r_array)
	
	# convert to matrix if necessary
	if (length(dim(r_array)) == 3) {
		dim(r_array) <- c(prod(dim(r_array)[1:2]), dim(r_array)[3])
	}

	x <- r_array[,1]
	y <- r_array[,2]
	z <- r_array[,3]
	sr3x <- 1.73205080756888*x #sqrt(3)
	sr3y <- 1.73205080756888*y #sqrt(3)
	xsq <- x*x
	ysq <- y*y
	zsq <- z*z
	xsq_ysq <- xsq+ysq
	rsq <- xsq_ysq+zsq
	r <- sqrt(rsq)
	db <- cbind(
		zsq-0.5*(xsq_ysq),
		0.866025403784439*(xsq-ysq), # sqrt(3)/2
		sr3x*z,
		sr3y*z,
		sr3x*y
	)
	
	if (dist) {
		r5 <- rsq*rsq*r
		da_dist <- 1/r5
		value_dist <- da_dist*db
		colnames(value_dist) <- c("d1", "d2", "d3", "d4", "d5")
		
		# convert back to array if necessary
		if (length(input_dim) == 3) {
			dim(value_dist) <- c(input_dim[1:2], 5)
			dimnames(value_dist) <- c(input_dimnames[1:2], list(c("d1", "d2", "d3", "d4", "d5")))
		}
	}
	
	if (unit) {
		da_unit <- 1/rsq
		value_unit <- da_unit*db
		colnames(value_unit) <- c("d1", "d2", "d3", "d4", "d5")
		
		# convert back to array if necessary
		if (length(input_dim) == 3) {
			dim(value_unit) <- c(input_dim[1:2], 5)
			dimnames(value_unit) <- c(input_dimnames[1:2], list(c("d1", "d2", "d3", "d4", "d5")))
		}
	}
	
	if (gradient) {
	
		# it is numerically more efficient to calculate the derivative with the product rule
		sr3z <- 1.73205080756888*z
		zero <- numeric(nrow(r_array))
		ddbdxyz <- array(c(
			-x, sr3x, sr3z, zero, sr3y,
			-y, -sr3y, zero, sr3z, sr3x,
			2*z, zero, sr3x, sr3y, zero
		), dim=c(nrow(r_array), 5, 3))
		
		# calculate derivative using product rule: db*ddadxyz + da*ddbdxyz
		# da could be premultiplied during construction of ddbdxyz
		
		if (dist) {
			r7inv <- 1/(r5*rsq)
			#ddadxyz <- c(-5*x/r7, -5*y/r7, -5*z/r7)
			ddadxyz_dist <- -5*r_array*r7inv
		
			grad_dist <- as.vector(db)*as.vector(ddadxyz_dist[,c(1L,1L,1L,1L,1L,2L,2L,2L,2L,2L,3L,3L,3L,3L,3L)]) + da_dist*ddbdxyz
		
			# convert to array if necessary
			if (length(input_dim) == 3) {
				dim(grad_dist) <- c(input_dim[1:2], 5, input_dim[3])
				dimnames(grad_dist) <- c(input_dimnames[1:2], list(c("d1", "d2", "d3", "d4", "d5")), input_dimnames[3])
			}
			
			attr(value_dist, "gradient") <- grad_dist
		}
		
		if (unit) {
			r4inv <- 1/(rsq*rsq)
			#ddadxyz <- c(-5*x/r7, -5*y/r7, -5*z/r7)
			ddadxyz_unit <- -2*r_array*r4inv
		
			grad_unit <- as.vector(db)*as.vector(ddadxyz_unit[,c(1L,1L,1L,1L,1L,2L,2L,2L,2L,2L,3L,3L,3L,3L,3L)]) + da_unit*ddbdxyz
		
			# convert to array if necessary
			if (length(input_dim) == 3) {
				dim(grad_unit) <- c(input_dim[1:2], 5, input_dim[3])
				dimnames(grad_unit) <- c(input_dimnames[1:2], list(c("d1", "d2", "d3", "d4", "d5")), input_dimnames[3])
			}
			
			attr(value_unit, "gradient") <- grad_unit
		}
	}
	
	if (dist) {
		if (unit) {
			list(dist=value_dist, unit=value_unit)
		} else {
			value_dist
		}
	} else {
		value_unit
	}
}

#' Back-propagate energy derivative from d array to r array
#'
#' @param d_d_array_d_r_array 4D array (pairs, models, tensor elements, xyz) from
#'    `gradient` attribute of `r_array_to_d_array()`
#' @param d_energy_d_d_array 3D array (pairs, models, tensor elements)
#'
#' @return 3D array (pairs, models, xyz) with d_energy_d_r_array
r_array_to_d_array_backprop <- function(d_d_array_d_r_array, d_energy_d_d_array) {

	# calculate de/dr = dd/dr * de/dd for each xyz component of the internuclear vectors
	# coerce d_energy_d_d_array to a vector so it gets repeated over dd/dr
	d_energy_d_r_array_all <- d_d_array_d_r_array*as.vector(d_energy_d_d_array)
	# sum the individual interaction tensor component derivatives associated with x, y, and z
	apply(d_energy_d_r_array_all, c(1,2,4), sum)
}

#' Calculate group norm squared from dipole-dipole interaction tensors
#'
#' @param d_array 3D array (pairs, models, tensor elements) with interaction tensors
#' @param grouping list of integer vectors giving groupings of models to average interaction tensors
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return vector with norm squared for each atom pair. 
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array (pairs, models, tensor elements).
#'
#' @export
d_array_to_g <- function(d_array, grouping, gradient=FALSE) {

	value <- numeric(dim(d_array)[1])
	
	if (is.integer(grouping)) {
		
		# convert integer vector grouping to list grouping if necessary
		grouping <- lapply(lapply(seq_len(max(grouping)), "==", grouping), which)
	}
	
	if (gradient) {
	
		grad <- d_array
		grad[] <- 0
		attr(grad, "gradient") <- NULL
	}
	
	for (i in seq_along(grouping)) {
		
		# average dipole interaction tensor
		d_matrix <- 0
		
		# sum the dipole interaction tensors within each group
		for (j in seq_along(grouping[[i]])) {
			d_matrix <- d_matrix + d_array[,grouping[[i]][[j]],,drop=FALSE]
		}
		
		if (gradient) {
		
			# D[((dself + dother)/grouplength)^2*grouplength/totallength, {dself}]
			d_matrix_grad <- d_matrix*(2/dim(d_array)[2]/length(grouping[[i]]))
			# the gradients are the same for all members of the group because they contribute equally
			for (j in seq_along(grouping[[i]])) {
				grad[,grouping[[i]][[j]],] <- d_matrix_grad
			}
		}
		
		# calculate the mean by dividing by the length
		d_matrix <- d_matrix*(1/length(grouping[[i]]))
		
		# calculate self dot product (norm squared) and accumulate group's contribution to mean g
		value <- value + rowSums(d_matrix^2)*(length(grouping[[i]])/dim(d_array)[2])
	}
	
	if (gradient) {
	
		attr(value, "gradient") <- grad
	}
	
	value
}

#' Calculate group norm squared from dipole-dipole interaction tensors
#'
#' @param d_array 3D array (pairs, models, tensor elements) with interaction tensors
#' @param grouping_mat integer matrix (groups, models) giving groupings of models to average interaction tensors
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return vector with norm squared for each atom pair. 
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array (pairs, models, tensor elements).
#'
#' @export
d_array_to_g_matrix <- function(d_array, grouping_mat, gradient=FALSE) {

	g_list <- apply(grouping_mat, 1, d_array_to_g, d_array=d_array, gradient=gradient, simplify=FALSE)
	
	value <- simplify2array(g_list)
	
	if (!is.matrix(value)) {
		value <- matrix(value, nrow=dim(d_array)[1])
	}
	
	if (gradient) {
	
		attr(value, "gradient") <- simplify2array(lapply(g_list, attr, "gradient"))
	}
	
	value
}

#' Back-propagate energy derivative from g matrix to d array
#'
#' @param d_g_matrix_d_d_array 4D array (pairs, models, tensor components, groups) fom
#'    `gradient` attribute of `d_array_to_g_matrix()`
#' @param d_energy_d_g_matrix matrix (pairs, groups)
#'
#' @return 3D array (pairs, models, tensor elements) with d_energy_d_d_array
d_array_to_g_matrix_backprop <- function(d_g_matrix_d_d_array, d_energy_d_g_matrix) {

	# calculate de/dd = dg/dd * de/dg for all individual interaction tensor components
	d_energy_d_d_array <- 0
	# sum the contributions from the different norm squared values
	for (i in seq_len(dim(d_g_matrix_d_d_array)[4])) {
		d_energy_d_d_array <- d_energy_d_d_array + d_g_matrix_d_d_array[,,,i]*d_energy_d_g_matrix[,i]
	}
	
	d_energy_d_d_array
}

#' Calculate matrix of a values from a matrix of g values
#'
#' @param g_matrix matrix of g values with columns associated with different groupings
#' @param a_coef coefficients used for calculating a values
#'
#' @return Numeric matrix of `a` values with one row per atom pair and one
#'   column per eigenvalue combination in `a_coef`.
g_matrix_to_a_matrix <- function(g_matrix, a_coef) {

	a_matrix <- matrix(0, nrow=nrow(g_matrix), ncol=ncol(a_coef), dimnames=list(rownames(g_matrix), colnames(a_coef)))
	
	for (i in seq_len(ncol(a_coef))) {
		for (j in which(a_coef[,i] != 0)) {
			a_matrix[,i] <- a_matrix[,i] + a_coef[j,i]*g_matrix[,j]
		}
	}
	
	a_matrix
}

#' Back-propagate energy derivative from a matrix to g matrix
#'
#' @param a_coef coefficients used for calculating a values
#' @param d_energy_d_a_matrix matrix (pairs, eigenvalues)
#'
#' @return matrix (pairs, groupings) with d_energy_d_g_matrix
g_matrix_to_a_matrix_backprop <- function(a_coef, d_energy_d_a_matrix) {

	d_energy_d_g_matrix <- matrix(0, nrow=nrow(d_energy_d_a_matrix), ncol=nrow(a_coef), dimnames=list(rownames(d_energy_d_a_matrix), NULL))

	for (i in seq_len(nrow(a_coef))) {
		for (j in which(a_coef[i,] != 0)) {
			d_energy_d_g_matrix[,i] <- d_energy_d_g_matrix[,i] + d_energy_d_a_matrix[,j]*a_coef[i,j]
		}
	}
	
	d_energy_d_g_matrix
}

#' Calculate sigma from a matrix of a values
#'
#' @param a_matrix matrix of a values with columns associated with eigenvalues
#' @param lambda_prime_vec eigenvalues augmented with tumbling rate
#' @param proton_mhz spectrometer proton field strength in MHz
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return Numeric vector of sigma rates, one per row of `a_matrix`. If
#'   `gradient = TRUE`, the `"gradient"` attribute is a matrix of derivatives
#'   with respect to `a_matrix`.
a_matrix_to_sigma <- function(a_matrix, lambda_prime_vec, proton_mhz, gradient=FALSE) {

	# K^2
	# ((1.2566370614e-6 kilogram meters / (ampere^2 second^2))/(4*pi)*(1.054571726e-34 meter^2 kilograms / second)*(267.513e6 radian s^-1 T^-1)^2)^2
	# 5.69549944e-49 rad m^6 / s^2
	K_sq <- 5.69549944e-49
	omega0 <- proton_mhz*1e6*2*pi

	a_matrix <- t(a_matrix)
	
	J0 <- -2*colSums(a_matrix/lambda_prime_vec)
	#Jomega <- -2*colSums(a_matrix*lambda_prime_vec/(lambda_prime_vec^2 + (omega0)^2))
	J2omega <- -2*colSums(a_matrix*lambda_prime_vec/(lambda_prime_vec^2 + (2*omega0)^2))
	
	value <- 0.1*K_sq*(3*J2omega - 0.5*J0)
	
	if (gradient) {
	
		d_J0_d_a_matrix <- -2/lambda_prime_vec
		d_J2omega_d_a_matrix <- -2*lambda_prime_vec/(lambda_prime_vec^2 + (2*omega0)^2)
		
		d_sigma_d_a_matrix <- 0.1*K_sq*(3*d_J2omega_d_a_matrix - 0.5*d_J0_d_a_matrix)
		
		attr(value, "gradient") <- matrix(d_sigma_d_a_matrix, nrow=ncol(a_matrix), ncol=nrow(a_matrix), byrow=TRUE)
	}
	
	value
}

#' Back-propagate energy derivative from sigma to a matrix
#'
#' @param d_sigma_d_a_matrix matrix (pairs, eigenvalues) fom `gradient` attribute of
#'    `a_matrix_to_sigma()`
#' @param d_energy_d_sigma vector (pairs)
#'
#' @return matrix (pairs, eigenvalues) with d_energy_d_a_matrix
a_matrix_to_sigma_backprop <- function(d_sigma_d_a_matrix, d_energy_d_sigma) {

	d_sigma_d_a_matrix*d_energy_d_sigma
}

#' Calculate overall tumbling modes from diffusion tensor principal values
#'
#' This helper constructs the overall tumbling factor for use in
#'    `a_matrix_to_relax()` from the rank-2 rotational diffusion generator in
#'    the tesseral basis
#'    \eqn{(Y_{2,0}, Y_{2,2}^{c}, Y_{2,1}^{c}, Y_{2,1}^{s}, Y_{2,2}^{s})}.
#'
#' Let \eqn{(D_x, D_y, D_z)} denote the principal values of the rotational
#' diffusion tensor in the diffusion-frame principal-axis system. In this
#' basis, the rank-2 diffusion generator is
#' \deqn{
#' \mathbf{L}^{(2)} =
#' \begin{pmatrix}
#' 3(D_x + D_y) & \sqrt{3}(D_x - D_y) & 0 & 0 & 0 \\
#' \sqrt{3}(D_x - D_y) & D_x + D_y + 4D_z & 0 & 0 & 0 \\
#' 0 & 0 & D_x + 4D_y + D_z & 0 & 0 \\
#' 0 & 0 & 0 & 4D_x + D_y + D_z & 0 \\
#' 0 & 0 & 0 & 0 & D_x + D_y + 4D_z
#' \end{pmatrix}.
#' }
#' The returned `lambda_overall_vec` contains the negative decay rates
#' associated with the distinct eigenvalues of that matrix, so that they can be
#' combined with the negative internal-motion eigenvalues used elsewhere in the
#' package.
#' The current implementation assumes that `dxyz_vec` is already ordered in the
#' same diffusion-frame convention as the generator above. In particular, the
#' axial-symmetry shortcut is only triggered for \eqn{D_x = D_y}, so the unique
#' diffusion axis should correspond to the `z` axis of the supplied diffusion
#' frame.
#'
#' For a matrix input `(pairs, 5)`, each row is interpreted either as a single
#' unit rank-2 tensor or as an averaged rank-2 unit tensor in the diffusion
#' frame. Its squared norm
#' \deqn{S^2 = \lVert \mathbf d^{\mathrm u} \rVert^2}
#' is used as a measure of residual second-rank orientational order.
#' `S^2 = 1` corresponds to no orientational averaging, while `S^2 = 0`
#' corresponds to complete isotropic averaging of the rank-2 interaction.
#' Matrix cross-correlation is supported only when every row of both inputs has
#' \eqn{S^2 = 1} within `tol`, so that the rows can be interpreted as explicit
#' unit-tensor states rather than averaged tensors.
#'
#' When `S^2 > s2_min`, the row is normalized to a unit rank-2 direction and
#' projected onto the overall tumbling eigenmodes. For autocorrelation, if
#' \eqn{\hat{\mathbf d}} is the normalized row and \eqn{\mathbf v_i} are the
#' orthonormal eigenvectors of \eqn{\mathbf{L}^{(2)}}, the directional
#' overall-mode weights are
#' \deqn{a_{pi}^{\mathrm{dir}} = (\hat{\mathbf d}_p \cdot \mathbf v_i)^2.}
#' Because the eigenvectors are orthonormal, these directional weights sum to
#' one for each row:
#' \deqn{\sum_i a_{pi}^{\mathrm{dir}} = 1.}
#' For matrix cross-correlation with explicit unit-tensor rows
#' \eqn{\hat{\mathbf d}^A_p} and \eqn{\hat{\mathbf d}^B_p}, the corresponding
#' directional amplitudes are
#' \deqn{
#' a_{pi}^{\mathrm{dir}} =
#' (\hat{\mathbf d}^A_p \cdot \mathbf v_i)
#' (\hat{\mathbf d}^B_p \cdot \mathbf v_i).
#' }
#'
#' For axially symmetric diffusion with \eqn{D_x = D_y}, the generator is
#' already diagonal in three grouped mode classes. In the autocorrelation case,
#' the directional weights reduce to
#' \deqn{a_{p1}^{\mathrm{dir}} = \hat d_{p1}^2,}
#' \deqn{a_{p2}^{\mathrm{dir}} = \hat d_{p2}^2 + \hat d_{p5}^2,}
#' \deqn{a_{p3}^{\mathrm{dir}} = \hat d_{p3}^2 + \hat d_{p4}^2,}
#' while for cross-correlation they become
#' \deqn{a_{p1}^{\mathrm{dir}} = \hat d^A_{p1}\hat d^B_{p1},}
#' \deqn{a_{p2}^{\mathrm{dir}} = \hat d^A_{p2}\hat d^B_{p2} + \hat d^A_{p5}\hat d^B_{p5},}
#' \deqn{a_{p3}^{\mathrm{dir}} = \hat d^A_{p3}\hat d^B_{p3} + \hat d^A_{p4}\hat d^B_{p4}.}
#'
#' For fully anisotropic diffusion, the coupled
#' \eqn{(Y_{2,0}, Y_{2,2}^{c})} block is diagonalized analytically. Writing
#' \deqn{
#' A = 3(D_x + D_y), \qquad
#' B = \sqrt{3}(D_x - D_y), \qquad
#' C = D_x + D_y + 4D_z,
#' }
#' and
#' \deqn{
#' \theta = \frac{1}{2}\operatorname{atan2}(2B, A - C),
#' }
#' the first two directional weights are
#' \deqn{
#' a_{p1}^{\mathrm{dir}} =
#' (\cos\theta\,\hat d_{p1} + \sin\theta\,\hat d_{p2})^2,
#' }
#' \deqn{
#' a_{p2}^{\mathrm{dir}} =
#' (-\sin\theta\,\hat d_{p1} + \cos\theta\,\hat d_{p2})^2,
#' }
#' with
#' \deqn{a_{p3}^{\mathrm{dir}} = \hat d_{p3}^2, \quad
#' a_{p4}^{\mathrm{dir}} = \hat d_{p4}^2, \quad
#' a_{p5}^{\mathrm{dir}} = \hat d_{p5}^2.}
#' In the matrix cross-correlation case, the corresponding formulas replace
#' these squares by products of the projected \eqn{A} and \eqn{B} rows.
#'
#' The current implementation regularizes the poorly defined limit
#'    \eqn{S^2 \to 0} by blending the directional weights with symmetry-based
#'    limiting weights using `S^2` itself as the blend coefficient. For each
#'    row,
#' \deqn{
#' \mathbf a_p^{\mathrm{overall}} =
#' S_p^2 \, \mathbf a_p^{\mathrm{dir}} +
#' (1 - S_p^2)\, \mathbf a^{\mathrm{lim}}.
#' }
#' The limiting weights are chosen to respect the symmetry of the corresponding
#' diffusion model:
#' \itemize{
#'   \item Isotropic diffusion: one mode with weight \eqn{(1)}.
#'   \item Axially symmetric diffusion: three grouped weights
#'      \eqn{(1/5, 2/5, 2/5)}.
#'   \item Fully anisotropic diffusion: five equal weights
#'      \eqn{(1/5, 1/5, 1/5, 1/5, 1/5)}.
#' }
#' These fallback values were chosen so that when the residual rank-2
#' direction becomes negligible, no artificial orientation is introduced by the
#' regularization. In the isotropic case there is only one overall mode, so its
#' weight must be one. In the fully anisotropic case, loss of directional
#' information implies that no one eigenmode should be preferred over any
#' other, so the five weights are taken to be equal. In the axially symmetric
#' case, the corresponding grouped weights are obtained by summing the equal
#' fully anisotropic weights over the degenerate mode pairs, giving
#' \eqn{(1/5, 2/5, 2/5)}. This makes the fallback consistent with the symmetry
#' of the diffusion operator and with the normalization
#' \eqn{\sum_i a_{pi}^{\mathrm{overall}} = 1} used for the autocorrelation
#' implementation here.
#' Thus the returned rows sum to one in the autocorrelation case, and
#' `s2_min` controls the point below which directional information is deemed too
#' weak to normalize reliably.
#'
#' For a 3D array input `(pairs, models, 5)`, the overall-mode amplitudes are
#' averaged directly over the ensemble of unit tensors. This path supports both
#' auto- and cross-correlation. For each pair \eqn{p}, model \eqn{m}, and
#' overall mode \eqn{i}, let \eqn{\mathbf d^A_{pm}} and
#' \eqn{\mathbf d^B_{pm}} denote the unit rank-2 tensors for the two
#' interactions in that model. The returned amplitudes are then
#' \deqn{
#' a_{pi}^{\mathrm{overall}} =
#' \frac{1}{M}\sum_{m=1}^M
#' (\mathbf d^A_{pm} \cdot \mathbf v_i)
#' (\mathbf d^B_{pm} \cdot \mathbf v_i),
#' }
#' where \eqn{\mathbf v_i} are the rank-2 diffusion eigenvectors and
#' \eqn{M} is the number of ensemble members. For autocorrelation,
#' \eqn{\mathbf d^A = \mathbf d^B}, so this reduces to the average of squared
#' projections. For cross-correlation, corresponding rows and models of
#' `dunit_a_array` and `dunit_b_array` are paired directly before averaging.
#' This avoids the normalization ambiguity that arises for near-zero averaged
#' rank-2 tensors.
#'
#' @param dxyz_vec numeric vector with diffusion tensor principal values
#'    `c(Dx, Dy, Dz)`
#' @param dunit_a_array array of unit dipole-dipole tensors for vector A in the
#'    diffusion frame. A matrix `(pairs, 5)` uses the averaged-`dunit`
#'    regularized model, while a 3D array `(pairs, models, 5)` averages
#'    overall-mode weights directly over the ensemble of unit vectors.
#' @param dunit_b_array optional array of unit dipole-dipole tensors for vector
#'    B in the diffusion frame. If `NULL`, the autocorrelation case is assumed
#'    and `dunit_a_array` is used for both inputs.
#' @param s2_min numeric threshold below which anisotropic directional weights
#'    are replaced by their limiting fallback values. Defaults to `1e-4` as a
#'    conservative choice that is also intended to be suitable for later
#'    single-precision implementations.
#' @param tol numeric tolerance used to collapse degenerate overall tumbling
#'    modes
#'
#' @return List with elements `a_overall_matrix` and `lambda_overall_vec`.
#'    Equal overall decay rates are collapsed so the number of columns is the
#'    minimum needed to represent isotropic, axially symmetric, or fully
#'    anisotropic diffusion. Matrix input supports autocorrelation and also
#'    supports cross-correlation when all rows of both inputs satisfy
#'    \eqn{S^2 \approx 1}; 3D array input supports both auto- and
#'    cross-correlation by direct ensemble averaging.
#'
#' @export
dxyz_dunit_to_overall_modes <- function(dxyz_vec, dunit_a_array, dunit_b_array = NULL, s2_min = 1e-4, tol = sqrt(.Machine$double.eps)) {
	is_auto <- is.null(dunit_b_array)
	if (is_auto) {
		dunit_b_array <- dunit_a_array
	}

	stopifnot(
		is.numeric(dxyz_vec),
		length(dxyz_vec) == 3,
		is.numeric(s2_min),
		length(s2_min) == 1,
		s2_min >= 0,
		is.numeric(tol),
		length(tol) == 1,
		tol >= 0
	)

	dx <- dxyz_vec[1]
	dy <- dxyz_vec[2]
	dz <- dxyz_vec[3]
	sr3 <- 1.73205080756888

	if (is.matrix(dunit_a_array)) {
		if (!is_auto) {
			stopifnot(
				is.matrix(dunit_b_array),
				identical(dim(dunit_a_array), dim(dunit_b_array))
			)
			s2_a_vec <- rowSums(dunit_a_array^2)
			s2_b_vec <- rowSums(dunit_b_array^2)
			s2_scale <- pmax(1, abs(s2_a_vec), abs(s2_b_vec))
			if (any(abs(s2_a_vec - 1) > tol * s2_scale) ||
				any(abs(s2_b_vec - 1) > tol * s2_scale)) {
				stop("Matrix cross-correlation is only supported when every row of both inputs has S^2 = 1 within `tol`")
			}
		}
		return(.dxyz_dunit_matrix_to_overall_modes(dxyz_vec, dunit_a_array, dunit_b_array, s2_min = s2_min, tol = tol, regularize = is_auto))
	}
	if (is.array(dunit_a_array) && length(dim(dunit_a_array)) == 3 &&
		is.array(dunit_b_array) && length(dim(dunit_b_array)) == 3) {
		return(.dxyz_dunit_array_to_overall_modes(dxyz_vec, dunit_a_array, dunit_b_array, tol = tol))
	}

	stop("`dunit_a_array` must be either a matrix `(pairs, 5)` or a 3D array `(pairs, models, 5)`")
}

.dxyz_dunit_matrix_to_overall_modes <- function(dxyz_vec, dunit_a_matrix, dunit_b_matrix = dunit_a_matrix, s2_min = 1e-4, tol = sqrt(.Machine$double.eps), regularize = TRUE) {
	stopifnot(
		is.matrix(dunit_a_matrix),
		ncol(dunit_a_matrix) == 5,
		is.matrix(dunit_b_matrix),
		identical(dim(dunit_a_matrix), dim(dunit_b_matrix))
	)

	dx <- dxyz_vec[1]
	dy <- dxyz_vec[2]
	dz <- dxyz_vec[3]
	sr3 <- 1.73205080756888

	# Clip S^2 to [0, 1] to guard against tiny numerical excursions outside the
	# physically meaningful range for averaged rank-2 unit tensors.
	s2_a_vec <- pmin(1, pmax(0, rowSums(dunit_a_matrix^2)))
	s2_b_vec <- pmin(1, pmax(0, rowSums(dunit_b_matrix^2)))
	s2_vec <- s2_a_vec
	dunit_a_norm_vec <- sqrt(s2_a_vec)
	dunit_b_norm_vec <- sqrt(s2_b_vec)
	dunit_a_unit_matrix <- dunit_a_matrix
	valid_a_idx <- dunit_a_norm_vec > 0
	if (any(valid_a_idx)) {
		dunit_a_unit_matrix[valid_a_idx, ] <- dunit_a_unit_matrix[valid_a_idx, , drop = FALSE] / dunit_a_norm_vec[valid_a_idx]
	}
	dunit_b_unit_matrix <- dunit_b_matrix
	valid_b_idx <- dunit_b_norm_vec > 0
	if (any(valid_b_idx)) {
		dunit_b_unit_matrix[valid_b_idx, ] <- dunit_b_unit_matrix[valid_b_idx, , drop = FALSE] / dunit_b_norm_vec[valid_b_idx]
	}
	blend_vec <- if (regularize) s2_vec else rep(1, length(s2_vec))

	# Isotropic
	if (abs(dx - dy) <= tol * max(1, abs(dx), abs(dy)) &&
		abs(dx - dz) <= tol * max(1, abs(dx), abs(dz))) {
		a_overall_matrix <- matrix(
			rowSums(dunit_a_unit_matrix * dunit_b_unit_matrix),
			nrow = nrow(dunit_a_matrix),
			ncol = 1,
			dimnames = list(NULL, "overall1")
		)
		lambda_overall_vec <- c(overall1 = -6 * dx)
		return(list(
			a_overall_matrix = a_overall_matrix,
			lambda_overall_vec = lambda_overall_vec
		))
	}

	# Axially symmetric
	if (abs(dx - dy) <= tol * max(1, abs(dx), abs(dy))) {
		a_overall_matrix_dir <- cbind(
			dunit_a_unit_matrix[, 1] * dunit_b_unit_matrix[, 1],
			dunit_a_unit_matrix[, 2] * dunit_b_unit_matrix[, 2] +
				dunit_a_unit_matrix[, 5] * dunit_b_unit_matrix[, 5],
			dunit_a_unit_matrix[, 3] * dunit_b_unit_matrix[, 3] +
				dunit_a_unit_matrix[, 4] * dunit_b_unit_matrix[, 4]
		)
		a_overall_matrix <- blend_vec * a_overall_matrix_dir +
			(1 - blend_vec) * matrix(rep(c(1 / 5, 2 / 5, 2 / 5), each = nrow(dunit_a_matrix)), nrow = nrow(dunit_a_matrix))
		colnames(a_overall_matrix) <- paste0("overall", 1:3)
		lambda_overall_vec <- c(
			overall1 = -3 * (dx + dy),
			overall2 = -(2 * dx + 4 * dz),
			overall3 = -(5 * dx + dz)
		)
		return(list(
			a_overall_matrix = a_overall_matrix,
			lambda_overall_vec = lambda_overall_vec
		))
	}

	# Fully anisotropic
	block_a <- 3 * (dx + dy)
	block_b <- sr3 * (dx - dy)
	block_c <- dx + dy + 4 * dz
	block_trace_half <- 0.5 * (block_a + block_c)
	block_diff_half <- 0.5 * (block_a - block_c)
	block_radius <- sqrt(block_diff_half^2 + block_b^2)

	lambda_overall_vec <- -c(
		block_trace_half + block_radius,
		block_trace_half - block_radius,
		dx + 4 * dy + dz,
		4 * dx + dy + dz,
		dx + dy + 4 * dz
	)

	theta <- 0.5 * atan2(2 * block_b, block_a - block_c)
	cos_theta <- cos(theta)
	sin_theta <- sin(theta)

	a_overall_matrix_dir <- cbind(
		(cos_theta * dunit_a_unit_matrix[, 1] + sin_theta * dunit_a_unit_matrix[, 2]) *
			(cos_theta * dunit_b_unit_matrix[, 1] + sin_theta * dunit_b_unit_matrix[, 2]),
		(-sin_theta * dunit_a_unit_matrix[, 1] + cos_theta * dunit_a_unit_matrix[, 2]) *
			(-sin_theta * dunit_b_unit_matrix[, 1] + cos_theta * dunit_b_unit_matrix[, 2]),
		dunit_a_unit_matrix[, 3] * dunit_b_unit_matrix[, 3],
		dunit_a_unit_matrix[, 4] * dunit_b_unit_matrix[, 4],
		dunit_a_unit_matrix[, 5] * dunit_b_unit_matrix[, 5]
	)
	a_overall_matrix <- blend_vec * a_overall_matrix_dir +
		(1 - blend_vec) * matrix(1 / 5, nrow = nrow(dunit_a_matrix), ncol = length(lambda_overall_vec))
	colnames(a_overall_matrix) <- paste0("overall", 1:5)
	names(lambda_overall_vec) <- colnames(a_overall_matrix)

	list(
		a_overall_matrix = a_overall_matrix,
		lambda_overall_vec = lambda_overall_vec
	)
}

.dxyz_dunit_array_to_overall_modes <- function(dxyz_vec, dunit_a_array, dunit_b_array, tol = sqrt(.Machine$double.eps)) {
	stopifnot(
		is.array(dunit_a_array),
		length(dim(dunit_a_array)) == 3,
		dim(dunit_a_array)[3] == 5,
		is.array(dunit_b_array),
		length(dim(dunit_b_array)) == 3,
		identical(dim(dunit_a_array), dim(dunit_b_array))
	)

	n_pairs <- dim(dunit_a_array)[1]
	n_models <- dim(dunit_a_array)[2]
	dx <- dxyz_vec[1]
	dy <- dxyz_vec[2]
	dz <- dxyz_vec[3]
	sr3 <- 1.73205080756888

	a1 <- matrix(dunit_a_array[, , 1], nrow = n_pairs, ncol = n_models)
	a2 <- matrix(dunit_a_array[, , 2], nrow = n_pairs, ncol = n_models)
	a3 <- matrix(dunit_a_array[, , 3], nrow = n_pairs, ncol = n_models)
	a4 <- matrix(dunit_a_array[, , 4], nrow = n_pairs, ncol = n_models)
	a5 <- matrix(dunit_a_array[, , 5], nrow = n_pairs, ncol = n_models)
	b1 <- matrix(dunit_b_array[, , 1], nrow = n_pairs, ncol = n_models)
	b2 <- matrix(dunit_b_array[, , 2], nrow = n_pairs, ncol = n_models)
	b3 <- matrix(dunit_b_array[, , 3], nrow = n_pairs, ncol = n_models)
	b4 <- matrix(dunit_b_array[, , 4], nrow = n_pairs, ncol = n_models)
	b5 <- matrix(dunit_b_array[, , 5], nrow = n_pairs, ncol = n_models)

	# Isotropic
	if (abs(dx - dy) <= tol * max(1, abs(dx), abs(dy)) &&
		abs(dx - dz) <= tol * max(1, abs(dx), abs(dz))) {
		a_overall_matrix <- matrix(
			rowMeans(a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4 + a5 * b5),
			nrow = n_pairs,
			ncol = 1,
			dimnames = list(NULL, "overall1")
		)
		lambda_overall_vec <- c(overall1 = -6 * dx)
		return(list(
			a_overall_matrix = a_overall_matrix,
			lambda_overall_vec = lambda_overall_vec
		))
	}

	# Axially symmetric
	if (abs(dx - dy) <= tol * max(1, abs(dx), abs(dy))) {
		a_overall_matrix <- cbind(
			rowMeans(a1 * b1),
			rowMeans(a2 * b2 + a5 * b5),
			rowMeans(a3 * b3 + a4 * b4)
		)
		colnames(a_overall_matrix) <- paste0("overall", 1:3)
		lambda_overall_vec <- c(
			overall1 = -3 * (dx + dy),
			overall2 = -(2 * dx + 4 * dz),
			overall3 = -(5 * dx + dz)
		)
		return(list(
			a_overall_matrix = a_overall_matrix,
			lambda_overall_vec = lambda_overall_vec
		))
	}

	# Fully anisotropic
	block_a <- 3 * (dx + dy)
	block_b <- sr3 * (dx - dy)
	block_c <- dx + dy + 4 * dz
	block_trace_half <- 0.5 * (block_a + block_c)
	block_diff_half <- 0.5 * (block_a - block_c)
	block_radius <- sqrt(block_diff_half^2 + block_b^2)

	lambda_overall_vec <- -c(
		block_trace_half + block_radius,
		block_trace_half - block_radius,
		dx + 4 * dy + dz,
		4 * dx + dy + dz,
		dx + dy + 4 * dz
	)

	theta <- 0.5 * atan2(2 * block_b, block_a - block_c)
	cos_theta <- cos(theta)
	sin_theta <- sin(theta)

	mode1_a <- cos_theta * a1 + sin_theta * a2
	mode2_a <- -sin_theta * a1 + cos_theta * a2
	mode1_b <- cos_theta * b1 + sin_theta * b2
	mode2_b <- -sin_theta * b1 + cos_theta * b2
	a_overall_matrix <- cbind(
		rowMeans(mode1_a * mode1_b),
		rowMeans(mode2_a * mode2_b),
		rowMeans(a3 * b3),
		rowMeans(a4 * b4),
		rowMeans(a5 * b5)
	)
	colnames(a_overall_matrix) <- paste0("overall", 1:5)
	names(lambda_overall_vec) <- colnames(a_overall_matrix)

	list(
		a_overall_matrix = a_overall_matrix,
		lambda_overall_vec = lambda_overall_vec
	)
}

#' Calculate relaxation rates from spectral density functions with internal and overall motion
#'
#' Calculate per-pair relaxation rates by combining internal and overall motion
#' into a spectral density function and summing the resulting spectral density
#' terms with user-supplied coefficients.
#'
#' For pair \eqn{p}, this function evaluates
#' \deqn{R_p = \sum_t c_{pt} J_p(\omega_{pt})}
#' where \eqn{t} indexes the supplied spectral density terms, with coefficients
#' \eqn{c_{pt}} and frequencies \eqn{\omega_{pt}} from `spec_den_term_array`.
#' The spectral density is calculated as
#' \deqn{J_p(\omega) = -\sum_i \sum_j
#' a^{\mathrm{int}}_{pj} a^{\mathrm{overall}}_{pi}
#' \frac{\lambda^{\prime}_{ij}}{(\lambda^{\prime}_{ij})^2 + \omega^2}}
#' with combined decay rates
#' \deqn{\lambda^{\prime}_{ij} = \lambda^{\mathrm{int}}_j + \lambda^{\mathrm{overall}}_i.}
#' In the spectral density summation, the minus sign reflects the use of
#' negative decay eigenvalues.
#'
#' @param a_int_matrix (pairs, eigenvalues) matrix of internal motion amplitudes
#' @param lambda_int_vec internal motion eigenvalues associated with a_int_matrix columns
#' @param a_overall_matrix (pairs, eigenvalues) matrix of overall rotational diffusion amplitudes
#' @param lambda_overall_vec overall rotational diffusion eigenvalues associated with a_overall_matrix columns
#' @param spec_den_term_array array `(pairs, terms, components)` with
#'    spectral-density coefficients in component 1 and frequencies in
#'    component 2
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return Numeric vector of relaxation rates, one value per atom pair.
#'
#' @seealso [a_matrix_to_relax_backprop()]
#'
#' @export
a_matrix_to_relax <- function(a_int_matrix, lambda_int_vec, a_overall_matrix, lambda_overall_vec, spec_den_term_array, gradient = FALSE) {
	stopifnot(
		is.matrix(a_int_matrix),
		is.matrix(a_overall_matrix),
		is.array(spec_den_term_array),
		length(dim(spec_den_term_array)) == 3,
		ncol(a_int_matrix) == length(lambda_int_vec),
		ncol(a_overall_matrix) == length(lambda_overall_vec),
		nrow(a_int_matrix) == nrow(a_overall_matrix),
		dim(spec_den_term_array)[1] == nrow(a_int_matrix),
		dim(spec_den_term_array)[3] == 2,
		lambda_int_vec <= 0,
		lambda_overall_vec <= 0
	)

	spec_den_coef_matrix <- spec_den_term_array[, , 1, drop = FALSE]
	dim(spec_den_coef_matrix) <- dim(spec_den_term_array)[1:2]
	spec_den_freq_matrix <- spec_den_term_array[, , 2, drop = FALSE]
	dim(spec_den_freq_matrix) <- dim(spec_den_term_array)[1:2]
	value <- numeric(nrow(a_int_matrix))
	if (gradient) {
		d_value_d_a_int_matrix <- matrix(0, nrow = nrow(a_int_matrix), ncol = ncol(a_int_matrix))
	}

	for (i in seq_along(lambda_overall_vec)) {
		lambda_prime_vec <- lambda_int_vec + lambda_overall_vec[i]
		for (j in seq_along(lambda_prime_vec)) {
			term_vec <- rowSums(
				spec_den_coef_matrix * (
					-a_overall_matrix[, i] * lambda_prime_vec[j] /
						(lambda_prime_vec[j]^2 + spec_den_freq_matrix^2)
				)
			)
			value <- value + a_int_matrix[, j] * term_vec
			if (gradient) {
				d_value_d_a_int_matrix[, j] <- d_value_d_a_int_matrix[, j] + term_vec
			}
		}
	}

	if (gradient) {
		attr(value, "gradient") <- d_value_d_a_int_matrix
	}

	value
}

#' Back-propagate energy derivative from relaxation rates to a matrix
#'
#' For pair \eqn{p} and internal correlation function component \eqn{j}, 
#' the chain rule gives
#' \deqn{\frac{\partial E}{\partial a_{pj}} =
#' \frac{\partial E}{\partial R_p}
#' \frac{\partial R_p}{\partial a_{pj}}.}
#'
#' @param d_relax_d_a_matrix matrix (pairs, internal eigenvalues) from `gradient` attribute of
#'    `a_matrix_to_relax()`
#' @param d_energy_d_relax vector (pairs)
#'
#' @return Matrix `(pairs, internal eigenvalues)` with `d_energy_d_a_matrix`.
#'
#' @seealso [a_matrix_to_relax()]
a_matrix_to_relax_backprop <- function(d_relax_d_a_matrix, d_energy_d_relax) {

	d_relax_d_a_matrix * d_energy_d_relax
}

#' Calculate restraint energy from group norm squared values
#'
#' @param g current group norm squared values
#' @param g0 target group norm squared values
#' @param k force constant
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return restraint energy calculated using \eqn{k*(g-g0)^2}
#'
#' @export
g_to_energy <- function(g, g0, k=1, gradient=FALSE) {

	attr(g, "gradient") <- NULL
    expr1 <- g - g0
    value <- k * expr1^2
    if (gradient) {
		grad <- k * (2 * expr1)
		attr(value, "gradient") <- grad
    } else {
    	attr(value, "gradient") <- NULL
    }
    value
}

#' Quadratic loss function
#'
#' @param x current values
#' @param x0 target values
#' @param k force constant
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return \deqn{k(x - x_0)^2}
#'
#' @export
quadratic_loss <- function(x, x0, k = 1, gradient = FALSE) {

	expr1 <- x - x0
	value <- k * expr1^2
	if (gradient) {
		attr(value, "gradient") <- k * (2 * expr1)
	}
	value
}

#' Scale real or back calculated values towards zero by taking a power
#'
#' @param x numerical values to be rescaled
#' @param p numeric power to raise the values to, usually 1 or less
#'
#' @return \deqn{ ( ( |x| + 1 )^p - 1 ) \operatorname{sgn}(x) }
power_scale <- function(x, p) {
	( ( abs(x) + 1 )^p - 1 ) * sign(x)
}

#' Traditional squared loss function scaled by a power
#'
#' @param x current values
#' @param x0 target values
#' @param p numeric power to raise the values to, usually 1 or less
#' @param k force constant
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return \deqn{ k \left ( ( ( |x| + 1 )^p - 1 ) \operatorname{sgn}(x) - ( ( |x_0| + 1 )^p - 1 ) \operatorname{sgn}(x_0) \right )^2 }
#'
#' @examples
#' x <- seq(-10, 10, by=0.1)
#' loss_x <- sapply(x, power_scaled_loss, 2, 0.25)
#' loss_x_grad <- power_scaled_loss(x, 2, 0.25, gradient=TRUE)
#' oldpar <- par(mfrow = c(2, 1))
#' plot(x, loss_x, type="l", ylab="loss")
#' plot(x, attr(loss_x_grad, "gradient"), type="l", ylab="dloss/dx")
#' points(x[-1]-mean(diff(x))/2, diff(loss_x)/mean(diff(x)), type="l", col="blue")
#' abline(h = 0, col="gray")
#' legend("bottomright",
#'   legend=c("Analytical", "Finite Difference"),
#'   bty="n", lwd=1, col=c("black", "blue"))
#' par(oldpar)
#'
#' @export
power_scaled_loss <- function(x, x0, p=1, k=1, gradient=FALSE) {
	
	expr1 <- ( ( abs(x) + 1 )^p - 1 ) * sign(x) - ( ( abs(x0) + 1 )^p - 1 ) * sign(x0)
	value <- sum(k * expr1^2)
	if (gradient) {
		# D[(((Abs[x] + 1)^p - 1)*Sign[x] - ((Abs[x0] + 1)^p - 1)*Sign[x0])^2, x]
		grad <- k * (2 * expr1) * p * ( abs(x) + 1 )^(p-1) #* sign(x) * sign(x)
		attr(value, "gradient") <- grad
	}
	value
}

#' Calculate group norm squared values from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`)
#' @param grouping_list list of lists of integer vectors giving groupings of models to average interaction tensors
#'
#' @return matrix (pairs, groupings) with group norm squared values
#'
#' @export
coord_array_to_g_matrix <- function(coord_array, atom_pairs, grouping_list) {

	r_array <- coord_array_to_r_array(coord_array, atom_pairs)
	
	d_array <- r_array_to_d_array(r_array, gradient=FALSE)
	
	g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=FALSE))
	
	simplify2array(g_list)
}

#' Calculate g value restraint energy from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`)
#' @param grouping_list list of lists of integer vectors giving groupings of models to average interaction tensors
#' @param g0 target group norm squared values
#' @param k force constant
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return total restraint energy calculated using \eqn{k*(g-g0)^2}
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array (atoms, xyz, models).
#'
#' @export
coord_array_to_g_energy <- function(coord_array, atom_pairs, grouping_list, g0, k, gradient=FALSE) {

	# calculate internuclear vectors
	r_array <- coord_array_to_r_array(coord_array, atom_pairs)
	
	# calculate dipole-dipole interaction tensors
	d_array <- r_array_to_d_array(r_array, gradient=gradient)
	
	# calculate norm squared for different groupings of dipole-dipole interaction tensors
	g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=gradient))
	
	# convert the list above into a matrix (pairs, groupings)
	g_matrix <- simplify2array(g_list)
	
	# calculate energies from the norm squared values
	energy_matrix <- g_to_energy(g_matrix, g0, k, gradient=gradient)
	
	# return the sum of all the individual restraint energies
	value <- sum(energy_matrix)
	
	if (gradient) {
	
		# calculate de/dd = dg/dd * de/dg for all individual interaction tensor components
		d_energy_d_d_array <- 0
		# sum the contributions from the different norm squared values
		for (i in seq_along(g_list)) {
			d_energy_d_d_array <- d_energy_d_d_array + attr(g_list[[i]], "gradient")*attr(energy_matrix, "gradient")[,i]
		}
		
		# calculate de/dr = dd/dr * de/dd for each xyz component of the internuclear vectors
		# coerce d_energy_d_d_array to a vector so it gets repeated over dd/dr
		d_energy_d_r_array_all <- attr(d_array, "gradient")*as.vector(d_energy_d_d_array)
		# sum the individual interaction tensor component derivatives associated with x, y, and z
		d_energy_d_r_array <- apply(d_energy_d_r_array_all, c(1,2,4), sum)
		
		# initialize empty gradient with dimensions equal to input coordinates
		grad <- coord_array
		grad[] <- 0
		
		# propagate the internuclear vector derivatives back onto the atomic coordinates
		for (i in seq_len(dim(d_energy_d_r_array)[1])) {
			
			pair_grad <- t(d_energy_d_r_array[i,,])
			# first component is subtracted so derivative is too
			grad[atom_pairs[i,1],,] <- grad[atom_pairs[i,1],,] - pair_grad
			# second component is added so derivative is too
			grad[atom_pairs[i,2],,] <- grad[atom_pairs[i,2],,] + pair_grad
		}
		
		attr(value, "gradient") <- grad
	}
	
	value
}

#' Calculate g value restraint energy from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array`)
#' @param grouping_mat integer matrix (groups, models) giving groupings of models to average interaction tensors
#' @param g0 target group norm squared values
#' @param loss_func loss function to use to calculate energy
#' @param ... additional parameters passed to `loss_func`
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return total restraint energy calculated using `loss_func`
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array (atoms, xyz, models).
#'
#' @export
coord_array_to_g_energy_refactored <- function(coord_array, atom_pairs, grouping_mat, g0, loss_func = power_scaled_loss, ..., gradient=FALSE) {

	# calculate internuclear vectors
	r_array <- coord_array_to_r_array(coord_array, atom_pairs)
	
	# calculate dipole-dipole interaction tensors
	d_array <- r_array_to_d_array(r_array, gradient=gradient)
	
	# calculate norm squared for different groupings of dipole-dipole interaction tensors
	g_matrix <- d_array_to_g_matrix(d_array, grouping_mat, gradient=gradient)
	
	# calculate energies from the norm squared values
	energy_matrix <- loss_func(g_matrix, g0, ..., gradient=gradient)
	
	# return the sum of all the individual restraint energies
	value <- sum(energy_matrix)
	
	if (gradient) {
	
		# initialize empty gradient with dimensions equal to input coordinates
		d_energy_d_coord_array <- coord_array
		d_energy_d_coord_array[] <- 0

		d_energy_d_g_matrix <- attr(energy_matrix, "gradient")

		# back-propagate derivatives from g matrix to d array
		d_energy_d_d_array <- d_array_to_g_matrix_backprop(attr(g_matrix, "gradient"), d_energy_d_g_matrix)
		
		# back-propagate derivatives from d array to r array
		d_energy_d_r_array <- r_array_to_d_array_backprop(attr(d_array, "gradient"), d_energy_d_d_array)
		
		# accumulate back-propagated derivatives from r array to coord array gradient
		d_energy_d_coord_array <- coord_array_to_r_array_backprop(d_energy_d_coord_array, atom_pairs, d_energy_d_r_array)
		
		# set coordinate gradient to accumulated values
		attr(value, "gradient") <- d_energy_d_coord_array
	}
	
	value
}

#' Convert one rate block of atom-relax columns into a spectral-density term array
#'
#' Reconstruct a per-rate spectral-density term array from columns in an
#' `*_atom_relax.csv`-style data frame.
#'
#' The expected column naming convention is
#' `"<rate_name>_<freq_name>_<component>"`, where `<component>` is either
#' `"coef"` or `"freq"`. If present, the special columns
#' `"<rate_name>_value"` and `"<rate_name>_k"` store the target relaxation
#' values and optional loss weights associated with that rate block.
#'
#' Parsing is anchored to the supplied `rate_name`, so the rate name itself may
#' contain underscores. The `freq_name` portion must not contain underscores.
#' For each `freq_name`, both `coef` and `freq` columns must be present.
#'
#' @param atom_relax_df data frame containing one or more rate blocks encoded as
#'   atom-relax columns
#' @param rate_name character scalar giving the rate block to extract
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`value`}{numeric vector from the `"<rate_name>_value"` column, or
#'       `NULL` if that column is absent}
#'     \item{`k`}{numeric vector from the `"<rate_name>_k"` column, or `NULL`
#'       if that column is absent}
#'     \item{`spec_den_term_array`}{numeric array `(pairs, terms, 2)` with
#'       second-dimension names taken from `<freq_name>` and third-dimension
#'       names `c("coef", "freq")`}
#'   }
#'
#' @examples
#' atom_relax_df <- data.frame(
#'   r1_value = c(1.0, 2.0),
#'   r1_k = c(4.0, 9.0),
#'   r1_0_coef = c(0.1, 0.2),
#'   r1_0_freq = c(0, 0),
#'   r1_2wH_coef = c(0.3, 0.4),
#'   r1_2wH_freq = c(100, 100),
#'   check.names = FALSE
#' )
#'
#' out <- atom_relax_columns_to_spec_den_term_array(atom_relax_df, "r1")
#' out$value
#' dim(out$spec_den_term_array)
#' dimnames(out$spec_den_term_array)[[2]]
#'
#' @export
atom_relax_columns_to_spec_den_term_array <- function(atom_relax_df, rate_name) {
	stopifnot(
		is.data.frame(atom_relax_df),
		is.character(rate_name),
		length(rate_name) == 1,
		nzchar(rate_name)
	)

	value_col <- paste0(rate_name, "_value")
	k_col <- paste0(rate_name, "_k")

	rate_prefix <- paste0(rate_name, "_")
	rate_cols <- names(atom_relax_df)[startsWith(names(atom_relax_df), rate_prefix)]
	term_cols <- setdiff(rate_cols, c(value_col, k_col))

	if (!length(term_cols)) {
		stop("No spectral-density term columns found for rate `", rate_name, "`")
	}

	col_suffix <- sub(paste0("^", rate_prefix), "", term_cols)
	col_parts <- regexec("^([^_]+)_(coef|freq)$", col_suffix)
	col_matches <- regmatches(col_suffix, col_parts)

	if (any(lengths(col_matches) != 3)) {
		stop(
			"All term columns for rate `", rate_name,
			"` must match `", rate_name, "_<freq_name>_<coef|freq>` with no underscores in <freq_name>"
		)
	}

	freq_names <- vapply(col_matches, `[`, character(1), 2)
	components <- vapply(col_matches, `[`, character(1), 3)

	if (anyDuplicated(term_cols)) {
		stop("Duplicate term columns found for rate `", rate_name, "`")
	}

	freq_levels <- unique(freq_names)
	for (freq_name in freq_levels) {
		freq_components <- components[freq_names == freq_name]
		if (!setequal(freq_components, c("coef", "freq"))) {
			stop(
				"Rate `", rate_name, "` term `", freq_name,
				"` must have both `coef` and `freq` columns"
			)
		}
	}

	spec_den_term_array <- array(
		NA_real_,
		dim = c(nrow(atom_relax_df), length(freq_levels), 2),
		dimnames = list(NULL, freq_levels, c("coef", "freq"))
	)

	for (j in seq_along(term_cols)) {
		spec_den_term_array[, freq_names[j], components[j]] <- atom_relax_df[[term_cols[j]]]
	}

	list(
		value = if (value_col %in% names(atom_relax_df)) atom_relax_df[[value_col]] else NULL,
		k = if (k_col %in% names(atom_relax_df)) atom_relax_df[[k_col]] else NULL,
		spec_den_term_array = spec_den_term_array
	)
}

#' Convert a spectral-density term array into one rate block of atom-relax columns
#'
#' Flatten a per-rate spectral-density term array into columns suitable for an
#' `*_atom_relax.csv`-style data frame.
#'
#' Output columns follow the naming convention
#' `"<rate_name>_<freq_name>_<component>"`, where `<component>` is either
#' `"coef"` or `"freq"`. If `value` or `k` is supplied, the additional
#' `"<rate_name>_value"` and `"<rate_name>_k"` columns store the target
#' relaxation values and optional loss weights.
#'
#' The second dimension of `spec_den_term_array` supplies `<freq_name>` through
#' its dimnames. The third dimension must have length two and corresponds to the
#' `coef` and `freq` components. If third-dimension names are missing, they are
#' assumed to be ordered as `c("coef", "freq")`.
#'
#' @param spec_den_term_array numeric array `(pairs, terms, 2)`
#' @param rate_name character scalar giving the rate block name
#' @param value optional numeric vector of length `dim(spec_den_term_array)[1]`
#'   to store in the `"<rate_name>_value"` column. If `NULL`, the value column
#'   is omitted.
#' @param k optional numeric vector of length `1` or
#'   `dim(spec_den_term_array)[1]` to store in the `"<rate_name>_k"` column.
#'   If `NULL`, the k column is omitted.
#'
#' @return A data frame containing optional `"<rate_name>_value"` and
#'   `"<rate_name>_k"` columns, plus one `coef` and one `freq` column for each
#'   term in `spec_den_term_array`.
#'
#' @examples
#' spec_den_term_array <- array(
#'   c(0.1, 0.2, 0, 0, 0.3, 0.4, 100, 100),
#'   dim = c(2, 2, 2),
#'   dimnames = list(NULL, c("0", "2wH"), c("coef", "freq"))
#' )
#'
#' spec_den_term_array_to_atom_relax_columns(
#'   spec_den_term_array = spec_den_term_array,
#'   rate_name = "r1",
#'   value = c(1.0, 2.0),
#'   k = c(4.0, 9.0)
#' )
#'
#' @export
spec_den_term_array_to_atom_relax_columns <- function(spec_den_term_array, rate_name, value = NULL, k = NULL) {
	stopifnot(
		is.array(spec_den_term_array),
		length(dim(spec_den_term_array)) == 3,
		dim(spec_den_term_array)[3] == 2,
		is.character(rate_name),
		length(rate_name) == 1,
		nzchar(rate_name)
	)

	n_pairs <- dim(spec_den_term_array)[1]
	term_names <- dimnames(spec_den_term_array)[[2]]
	component_names <- dimnames(spec_den_term_array)[[3]]

	if (is.null(term_names)) {
		stop("`spec_den_term_array` must have second-dimension names for term frequencies")
	}
	if (any(grepl("_", term_names, fixed = TRUE))) {
		stop("Term names in `spec_den_term_array` must not contain underscores")
	}

	if (is.null(component_names)) {
		component_names <- c("coef", "freq")
	} else if (!identical(component_names, c("coef", "freq"))) {
		stop("Third-dimension names of `spec_den_term_array` must be `c(\"coef\", \"freq\")`")
	}

	if (!is.null(value)) {
		stopifnot(is.numeric(value), length(value) == n_pairs)
	}
	if (!is.null(k)) {
		stopifnot(is.numeric(k), length(k) %in% c(1, n_pairs))
	}

	out <- data.frame(row.names = seq_len(n_pairs), check.names = FALSE)
	if (!is.null(value)) {
		out[[paste0(rate_name, "_value")]] <- value
	}
	if (!is.null(k)) {
		out[[paste0(rate_name, "_k")]] <- rep(k, length.out = n_pairs)
	}

	for (term_name in term_names) {
		for (component_name in component_names) {
			out[[paste(rate_name, term_name, component_name, sep = "_")]] <-
				spec_den_term_array[, term_name, component_name]
		}
	}

	out
}

#' Convert atom-relax columns into a named list of spectral-density term arrays
#'
#' Reconstruct all rate blocks in an `*_atom_relax.csv`-style data frame as a
#' named list. Each list element contains both the optional target relaxation
#' values and the corresponding spectral-density term array for one rate.
#'
#' Rate names are discovered automatically from columns ending in `_value`,
#' `_k`, `_coef`, or `_freq`. Parsing proceeds from the right, so rate names
#' may contain underscores. Term-frequency names must not contain underscores
#' and must be paired as `<freq_name>_coef` and `<freq_name>_freq`.
#'
#' @param atom_relax_df data frame containing one or more atom-relax rate blocks
#'
#' @return Named list whose elements are the return values of
#'   [atom_relax_columns_to_spec_den_term_array()]. Each element is itself a
#'   list with entries `value`, optional `k`, and `spec_den_term_array`.
#'
#' @examples
#' atom_relax_df <- data.frame(
#'   r1_value = c(1.0, 2.0),
#'   r1_k = c(4.0, 9.0),
#'   r1_0_coef = c(0.1, 0.2),
#'   r1_0_freq = c(0, 0),
#'   r2_0_coef = c(0.3, 0.4),
#'   r2_0_freq = c(10, 10),
#'   check.names = FALSE
#' )
#'
#' atom_relax_df_to_spec_den_term_array_list(atom_relax_df)
#'
#' @export
atom_relax_df_to_spec_den_term_array_list <- function(atom_relax_df) {

	stopifnot(is.data.frame(atom_relax_df))

	if (!ncol(atom_relax_df)) {
		return(list())
	}

	rate_names <- character(0)
	for (col_name in names(atom_relax_df)) {
		if (grepl("_(value|k)$", col_name)) {
			rate_names <- c(rate_names, sub("_(value|k)$", "", col_name))
			next
		}

		col_parts <- regexec("^(.+)_([^_]+)_(coef|freq)$", col_name)
		col_match <- regmatches(col_name, col_parts)[[1]]
		if (length(col_match) == 4) {
			rate_names <- c(rate_names, col_match[2])
		}
	}

	rate_names <- unique(rate_names)
	if (!length(rate_names)) {
		# Re-enable this check once all supported inputs are expected to use the
		# `*_atom_relax.csv` rate-column format.
		# stop("No atom-relax rate columns found")
		return(list())
	}

	out <- lapply(rate_names, function(rate_name) {
		atom_relax_columns_to_spec_den_term_array(atom_relax_df, rate_name)
	})
	names(out) <- rate_names

	out
}

#' Convert a named list of spectral-density term arrays into atom-relax columns
#'
#' Flatten a named list of per-rate spectral-density term arrays into a single
#' data frame suitable for an `*_atom_relax.csv`-style file.
#'
#' Each element of `spec_den_term_array_list` must be named by the desired
#' rate name and must either be:
#' \itemize{
#'   \item a spectral-density term array `(pairs, terms, 2)`, or
#'   \item a list with elements `spec_den_term_array` and optional `value` and
#'   `k`
#' }
#' In the second form, `value` and `k` are written to the corresponding
#' `"<rate_name>_value"` and `"<rate_name>_k"` columns.
#'
#' @param spec_den_term_array_list named list of per-rate spectral-density term
#'   arrays, or named list of lists containing `spec_den_term_array` and
#'   optional `value` and `k`
#'
#' @return Data frame combining all rate blocks side by side using the atom-
#'   relax column naming convention
#'
#' @examples
#' spec_den_term_array_list <- list(
#'   r1 = list(
#'     value = c(1.0, 2.0),
#'     k = c(4.0, 9.0),
#'     spec_den_term_array = array(
#'       c(0.1, 0.2, 0, 0),
#'       dim = c(2, 1, 2),
#'       dimnames = list(NULL, "0", c("coef", "freq"))
#'     )
#'   ),
#'   r2 = list(
#'     spec_den_term_array = array(
#'       c(0.3, 0.4, 10, 10),
#'       dim = c(2, 1, 2),
#'       dimnames = list(NULL, "0", c("coef", "freq"))
#'     )
#'   )
#' )
#'
#' spec_den_term_array_list_to_atom_relax_df(spec_den_term_array_list)
#'
#' @export
spec_den_term_array_list_to_atom_relax_df <- function(spec_den_term_array_list) {

	stopifnot(is.list(spec_den_term_array_list))

	if (!length(spec_den_term_array_list)) {
		return(data.frame(check.names = FALSE))
	}
	if (is.null(names(spec_den_term_array_list)) ||
		any(!nzchar(names(spec_den_term_array_list)))) {
		stop("`spec_den_term_array_list` must be a named list")
	}

	rate_dfs <- lapply(names(spec_den_term_array_list), function(rate_name) {
		rate_obj <- spec_den_term_array_list[[rate_name]]

		if (is.array(rate_obj)) {
			spec_den_term_array_to_atom_relax_columns(
				spec_den_term_array = rate_obj,
				rate_name = rate_name,
				value = NULL,
				k = NULL
			)
		} else if (is.list(rate_obj) &&
			"spec_den_term_array" %in% names(rate_obj)) {
			spec_den_term_array_to_atom_relax_columns(
				spec_den_term_array = rate_obj[["spec_den_term_array"]],
				rate_name = rate_name,
				value = rate_obj[["value"]],
				k = rate_obj[["k"]]
			)
		} else {
			stop(
				"Each element of `spec_den_term_array_list` must be an array or a list ",
				"containing `spec_den_term_array` and optional `value` and `k`"
			)
		}
	})

	n_rows <- vapply(rate_dfs, nrow, integer(1))
	if (length(unique(n_rows)) != 1) {
		stop("All rate blocks must have the same number of rows")
	}

	out <- do.call(cbind, rate_dfs)
	value_cols <- grepl("_value$", names(out))
	out[, c(which(value_cols), which(!value_cols)), drop = FALSE]
}

#' Read data for calculating spectral density functions
#'
#' @param prefix_path to prefix of four CSV files
#'
#' @return a list with elements: `atom_pairs`, `groupings`, `a_coef`, and `lambda_coef`
#'
#' @export
read_spec_den_data <- function(prefix_path) {

	atom_pairs <- utils::read.csv(paste0(prefix_path, "_atom_pairs.csv"))
	groupings <- unname(as.matrix(utils::read.csv(paste0(prefix_path, "_groupings.csv"), header=FALSE, row.names=NULL)))
	a_coef <- as.matrix(utils::read.csv(paste0(prefix_path, "_a_coef.csv"), check.names=FALSE))
	lambda_coef <- as.matrix(utils::read.csv(paste0(prefix_path, "_lambda_coef.csv"), check.names=FALSE, row.names=1))
	
	list(
		atom_pairs = atom_pairs,
		groupings = groupings,
		a_coef = a_coef,
		lambda_coef = lambda_coef
	)
}

#' Read data for calculating spectral-density-based relaxation rates
#'
#' @param prefix_path to prefix of four CSV files
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`atom_pairs`}{data frame of atom-pair metadata and target value
#'       columns, excluding flattened spectral-density term columns ending in
#'       `_coef` or `_freq`}
#'     \item{`unit`}{logical flag indicating whether the first interaction is
#'       encoded by unit-vector identifiers rather than atom identifiers. The
#'       current implementation returns a length-one logical. It is intended to
#'       generalize to one flag per interaction when cross-correlation support
#'       is added.}
#'     \item{`relax_data_list`}{named list of per-rate entries returned by
#'       [atom_relax_df_to_spec_den_term_array_list()]}
#'     \item{`groupings`}{grouping matrix}
#'     \item{`a_int_coef`}{matrix of amplitudes used to construct internal-motion
#'       coefficients}
#'     \item{`lambda_int_coef`}{matrix of rate coefficients used to construct
#'       internal-motion eigenvalues}
#'   }
#'
#' @export
read_spec_den_relax_data <- function(prefix_path) {

	atom_relax_df <- utils::read.csv(paste0(prefix_path, "_atom_relax.csv"), check.names = FALSE)
	groupings <- unname(as.matrix(utils::read.csv(paste0(prefix_path, "_groupings.csv"), header=FALSE, row.names=NULL)))
	a_int_coef <- as.matrix(utils::read.csv(paste0(prefix_path, "_a_coef.csv"), check.names=FALSE))
	lambda_int_coef <- as.matrix(utils::read.csv(paste0(prefix_path, "_lambda_coef.csv"), check.names=FALSE, row.names=1))
	relax_data_list <- atom_relax_df_to_spec_den_term_array_list(atom_relax_df)
	first_two_unit <- grepl("_unit$", names(atom_relax_df)[1:2])
	if (any(first_two_unit) && !all(first_two_unit)) {
		stop("The first two columns of `*_atom_relax.csv` must either both end in `_unit` or neither may")
	}
	unit <- all(first_two_unit)
	atom_pair_cols <- !grepl("_(coef|freq)$", names(atom_relax_df))
	atom_pairs <- atom_relax_df[, atom_pair_cols, drop = FALSE]
	
	list(
		atom_pairs = atom_pairs,
		unit = unit,
		relax_data_list = relax_data_list,
		groupings = groupings,
		a_int_coef = a_int_coef,
		lambda_int_coef = lambda_int_coef
	)
}

#' Write data for calculating spectral-density-based relaxation rates
#'
#' Write a `spec_den_relax_data`-style object to the four CSV files read by
#' [read_spec_den_relax_data()].
#'
#' The `*_atom_relax.csv` file is assembled by combining `atom_pairs` with the
#' flattened columns generated from `relax_data_list`. To preserve the `unit`
#' flag on reread, the first two atom-pair columns are written with names ending
#' in `"_unit"` when `unit` is `TRUE`, and without that suffix when `unit` is
#' `FALSE`.
#'
#' @param spec_den_relax_data list containing at least `atom_pairs`, `unit`,
#'   `relax_data_list`, `groupings`, `a_int_coef`, and `lambda_int_coef`
#' @param prefix_path character scalar giving the output file prefix
#'
#' @return Invisibly returns `prefix_path`
#'
#' @export
write_spec_den_relax_data <- function(spec_den_relax_data, prefix_path) {

	stopifnot(
		is.list(spec_den_relax_data),
		all(c("atom_pairs", "unit", "relax_data_list", "groupings", "a_int_coef", "lambda_int_coef") %in%
			names(spec_den_relax_data)),
		is.character(prefix_path),
		length(prefix_path) == 1,
		nzchar(prefix_path)
	)

	atom_pairs <- spec_den_relax_data[["atom_pairs"]]
	unit <- spec_den_relax_data[["unit"]]
	relax_data_list <- spec_den_relax_data[["relax_data_list"]]
	groupings <- spec_den_relax_data[["groupings"]]
	a_int_coef <- spec_den_relax_data[["a_int_coef"]]
	lambda_int_coef <- spec_den_relax_data[["lambda_int_coef"]]

	stopifnot(
		is.data.frame(atom_pairs),
		ncol(atom_pairs) >= 2,
		is.logical(unit),
		length(unit) == 1,
		is.list(relax_data_list),
		is.matrix(groupings),
		is.matrix(a_int_coef),
		is.matrix(lambda_int_coef)
	)

	atom_pairs_out <- atom_pairs
	first_two_names <- names(atom_pairs_out)[1:2]
	if (unit) {
		names(atom_pairs_out)[1:2] <- ifelse(
			grepl("_unit$", first_two_names),
			first_two_names,
			paste0(first_two_names, "_unit")
		)
	} else {
		names(atom_pairs_out)[1:2] <- sub("_unit$", "", first_two_names)
	}

	relax_df <- spec_den_term_array_list_to_atom_relax_df(relax_data_list)
	if (nrow(relax_df) && nrow(relax_df) != nrow(atom_pairs_out)) {
		stop("`atom_pairs` and `relax_data_list` must have the same number of rows")
	}
	atom_relax_df <- cbind(atom_pairs_out, relax_df)

	utils::write.csv(atom_relax_df, paste0(prefix_path, "_atom_relax.csv"), row.names = FALSE, na = "")
	utils::write.table(groupings, paste0(prefix_path, "_groupings.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
	utils::write.csv(a_int_coef, paste0(prefix_path, "_a_coef.csv"), row.names = FALSE)
	utils::write.csv(lambda_int_coef, paste0(prefix_path, "_lambda_coef.csv"))

	invisible(prefix_path)
}

#' Shift from one array dimension to another
#'
#' @param a array whose dimensions should be shifted
#' @param n integer factor to divide `dfrom` dimension and multiply `dto` dimension
#' @param dnames list with dimension names as alternative to specifying `n`
#' @param dfrom integer with dimension to move from
#' @param dto integer with dimension to move to (should be one greater than `dfrom`)
#'
#' @return Array with the same entries as `a`, reshaped so size is moved from
#'   dimension `dfrom` into dimension `dto`.
array_shift <- function(a, n=NULL, dnames=NULL, dfrom=1, dto=2) {

	if (!is.null(n)) {
	
		if (n == 1) {
			return(a)
		}
	
		new_a <- a
	
		new_dim <- dim(a)
		new_dim[dfrom] <- dim(a)[dfrom]/n
		new_dim[dto] <- dim(a)[dto]*n
		dim(new_a) <- new_dim
		
		new_dimnames <- dimnames(a)
		if (!is.null(dimnames(a)[[dfrom]])) {
			new_dimnames[[dfrom]] <- dimnames(a)[[dfrom]][seq_len(length(dimnames(a)[[dfrom]])/n)]
		}
		if (!is.null(dimnames(a)[[dto]])) {
			new_dimnames[[dto]] <- paste(rep(dimnames(a)[[dto]], each=n), seq_len(n), sep="-")
		}
		dimnames(new_a) <- new_dimnames
	
	} else if (!is.null(dnames)) {
	
		new_a <- a
		dim(new_a) <- sapply(dnames, length)
		dimnames(new_a) <- dnames
	
	} else {
	
		stop()
	}
	
	new_a
}

#' Calculate cross relaxation rates from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param rates named numeric vector with rates
#' @param spec_den_data_list list of data for calculating spectral density functions
#' @param proton_mhz spectrometer proton field strength in MHz
#'
#' @return Named list of numeric sigma vectors, one entry for each block in
#'   `spec_den_data_list`.
#'
#' @export
coord_array_to_sigma <- function(coord_array, rates, spec_den_data_list, proton_mhz) {

	lapply(spec_den_data_list, function(spec_den_data) {
	
		# calculate internuclear vectors (convert from Å^-3 to m^-3)
		r_array <- coord_array_to_r_array(coord_array*1e-10, spec_den_data[["atom_pairs"]][,1:2,drop=FALSE])
	
		# calculate dipole-dipole interaction tensors
		d_array <- r_array_to_d_array(r_array)
	
		# calculate the factor by which the number of models should be expanded
		n_shift <- ncol(spec_den_data[["groupings"]])/dim(coord_array)[3]
	
		# shift tensor components from atom pairs into virtual models
		d_array_shifted <- array_shift(d_array, n_shift)
		
		# calculate matrix of g values
		g_mat <- apply(spec_den_data[["groupings"]], 1, d_array_to_g, d_array=d_array_shifted)
		
		# calculate matrix of a values
		a_mat <- g_matrix_to_a_matrix(g_mat, spec_den_data[["a_coef"]])
		
		# calculate lambda eigenvalues
		lambda_vec <- -colSums(rates[rownames(spec_den_data[["lambda_coef"]])]*spec_den_data[["lambda_coef"]])
		
		# update eigenvalues to account for molecular tumbling
		lambda_prime_vec <- lambda_vec - rates["kc"]
		
		a_matrix_to_sigma(a_mat, lambda_prime_vec, proton_mhz)
	})
}

#' Calculate sigma restraint energy from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param rates named numeric vector with rates
#' @param spec_den_data_list list of data for calculating spectral density functions
#' @param proton_mhz spectrometer proton field strength in MHz
#' @param loss_func loss function to use
#' @param ... additional parameters passed to `loss_func`
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return total restraint energy calculated using `loss_func`
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array 
#' (atoms, xyz, models).
#'
#' @export
coord_array_to_sigma_energy <- function(coord_array, rates, spec_den_data_list, proton_mhz, loss_func = power_scaled_loss, ..., gradient=FALSE) {

	# intermediate derivatives that need to be captured for back-propagation
	d_array_list <- g_matrix_list <- sigma_vector_list <- vector("list", length(spec_den_data_list))
	
	for (i in seq_along(spec_den_data_list)) {

		# calculate internuclear vectors
		r_array <- coord_array_to_r_array(coord_array*1e-10, spec_den_data_list[[i]][["atom_pairs"]][,1:2,drop=FALSE])
	
		# calculate dipole-dipole interaction tensors
		d_array_list[[i]] <- r_array_to_d_array(r_array, gradient=gradient)
	
		# calculate the factor by which the number of models should be expanded
		n_shift <- ncol(spec_den_data_list[[i]][["groupings"]])/dim(coord_array)[3]
	
		# shift tensor components from atom pairs into virtual models
		d_array_shifted <- array_shift(d_array_list[[i]], n_shift)
	
		# calculate norm squared for different groupings of dipole-dipole interaction tensors
		g_matrix_list[[i]] <- d_array_to_g_matrix(d_array_shifted, spec_den_data_list[[i]][["groupings"]], gradient=gradient)
		
		# calculate matrix of a values
		a_matrix <- g_matrix_to_a_matrix(g_matrix_list[[i]], spec_den_data_list[[i]][["a_coef"]])
		
		# calculate lambda eigenvalues
		lambda_vector <- -colSums(rates[rownames(spec_den_data_list[[i]][["lambda_coef"]])]*spec_den_data_list[[i]][["lambda_coef"]])
		
		# update eigenvalues to account for molecular tumbling
		lambda_prime_vector <- lambda_vector - rates["kc"]
		
		# calculate cross relaxation rates
		sigma_vector_list[[i]] <- a_matrix_to_sigma(a_matrix, lambda_prime_vector, proton_mhz, gradient=TRUE)
	}
	
	# combine sigmas into a single vector
	sigma_vector <- unlist(sigma_vector_list)

	# combine target sigmas into a single vector
	sigma0_vector <- unlist(lapply(spec_den_data_list, function(x) {
		idx <- seq_len(nrow(x$atom_pairs)/(ncol(x$groupings)/dim(coord_array)[[3]]))
		x$atom_pairs[idx,"sigma"]
	}))
	
	# calculate energies from the sigma vectors
	energy_vector <- loss_func(sigma_vector, sigma0_vector, ..., gradient=gradient)
	
	# return the sum of all the individual restraint energies
	value <- sum(energy_vector)
	
	if (gradient) {
		
		# calculate starting offsets within sigma gradients
		sigma_vector_list_lengths <- sapply(sigma_vector_list, length)
		start_idx <- cumsum(c(1, sigma_vector_list_lengths[-length(sigma_vector_list)]))
	
		# initialize empty gradient with dimensions equal to input coordinates
		d_energy_d_coord_array <- coord_array
		d_energy_d_coord_array[] <- 0

		for (i in seq_along(spec_den_data_list)) {

			# extract particular segment of energy gradient
			d_energy_d_sigma <- attr(energy_vector, "gradient")[seq(start_idx[i], length.out=length(sigma_vector_list[[i]]))]
			
			# back-propagate derivatives from sigma to a matrix
			d_energy_d_a_matrix <- a_matrix_to_sigma_backprop(attr(sigma_vector_list[[i]], "gradient"), d_energy_d_sigma)
			
			# back-propagate derivatives from a array to g matrix
			d_energy_d_g_matrix <- g_matrix_to_a_matrix_backprop(spec_den_data_list[[i]][["a_coef"]], d_energy_d_a_matrix)
			
			# back-propagate derivatives from g matrix to d array
			d_energy_d_d_array <- d_array_to_g_matrix_backprop(attr(g_matrix_list[[i]], "gradient"), d_energy_d_g_matrix)
			
			# because of the way the arrays are structured shifting back isn't necessary
			# calculate the factor by which the number of models should be expanded
			#n_shift <- ncol(spec_den_data_list[[i]][["groupings"]])/dim(coord_array)[3]
			
			# shift tensor components from virtual models back to atom pairs
			#if (n_shift > 1) {
			#	d_energy_d_d_array <- array_shift(d_energy_d_d_array, dnames=dimnames(d_array_list[[i]]))
			#}
		
			# back-propagate derivatives from d array to r array
			d_energy_d_r_array <- r_array_to_d_array_backprop(attr(d_array_list[[i]], "gradient"), d_energy_d_d_array)
		
			# accumulate back-propagated derivatives from r array into coord array gradient
			d_energy_d_coord_array <- coord_array_to_r_array_backprop(d_energy_d_coord_array, spec_den_data_list[[i]][["atom_pairs"]][,1:2,drop=FALSE], d_energy_d_r_array)
		}
		
		# set coordinate gradient to accumulated values
		attr(value, "gradient") <- d_energy_d_coord_array*1e-10
	}
	
	value
}

#' Calculate relaxation rates from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param rates named numeric vector with rates
#' @param spec_den_relax_data_list list of data for calculating relaxation rates
#'
#' @return Named list of numeric matrices of relaxation rates, one per block in
#'   `spec_den_relax_data_list`, with columns corresponding to relaxation-rate
#'   types.
#'
#' @export
coord_array_to_relax <- function(coord_array, rates, spec_den_relax_data_list) {
	lapply(spec_den_relax_data_list, function(spec_den_relax_data) {
		unit <- spec_den_relax_data[["unit"]]
		if (is.null(unit)) {
			unit <- FALSE
		}
		stopifnot(is.logical(unit), length(unit) == 1)

		relax_data_list <- spec_den_relax_data[["relax_data_list"]]
		if (is.null(relax_data_list)) {
			stop("`spec_den_relax_data` must contain `relax_data_list`")
		}

		# calculate internuclear vectors (convert from Å^-3 to m^-3)
		r_array <- coord_array_to_r_array(coord_array * 1e-10, spec_den_relax_data[["atom_pairs"]][, 1:2, drop = FALSE])

		# calculate dipole-dipole interaction tensors
		d_array <- r_array_to_d_array(r_array, dist = !unit, unit = unit)

		# calculate the factor by which the number of models should be expanded
		n_shift <- ncol(spec_den_relax_data[["groupings"]]) / dim(coord_array)[3]

		# shift tensor components from atom pairs into virtual models
		d_array_shifted <- array_shift(d_array, n_shift)

		# calculate matrix of g values
		g_mat <- apply(spec_den_relax_data[["groupings"]], 1, d_array_to_g, d_array = d_array_shifted)

		# calculate matrix of a values
		a_int_mat <- g_matrix_to_a_matrix(g_mat, spec_den_relax_data[["a_int_coef"]])

		# calculate lambda eigenvalues
		lambda_int_vec <- -colSums(rates[rownames(spec_den_relax_data[["lambda_int_coef"]])] * spec_den_relax_data[["lambda_int_coef"]])

		# calculate unit dipole-dipole interaction tensors if not already done
		d_unit_array <- if (unit) d_array else r_array_to_d_array(r_array, dist = FALSE, unit = TRUE)

		# shift unit tensor components from atom pairs into virtual models
		d_unit_array_shifted <- array_shift(d_unit_array, n_shift)
		
		# calculate overall modes from diffusion constants
		overall_modes <- dxyz_dunit_to_overall_modes(rates[c("Dx", "Dy", "Dz")], d_unit_array_shifted)

		relax_mat <- vapply(
			relax_data_list,
			function(relax_entry) {
				a_matrix_to_relax(
					a_int_mat,
					lambda_int_vec,
					overall_modes[["a_overall_matrix"]],
					overall_modes[["lambda_overall_vec"]],
					relax_entry[["spec_den_term_array"]]
				)
			},
			FUN.VALUE = numeric(nrow(a_int_mat))
		)
		if (is.null(dim(relax_mat))) {
			relax_mat <- matrix(relax_mat, ncol = 1)
		}
		colnames(relax_mat) <- names(relax_data_list)
		relax_mat
	})
}

#' Calculate relaxation rate restraint energy from atomic coordinates
#'
#' @param coord_array 3D array (atoms, xyz, models) with atomic coordinates
#' @param rates named numeric vector with ensemble rates
#' @param spec_den_relax_data_list list of data for calculating relaxation rates
#' @param loss_func loss function to use
#' @param ... additional parameters passed to `loss_func`. If a `k` argument is
#'   not supplied here, `coord_array_to_relax_energy()` looks for an optional
#'   numeric `k` field in each `relax_data_list` entry and uses those values as
#'   per-relaxation-rate force constants. Each such `k` may have length 1 or
#'   the number of relaxation rates in that entry.
#' @param gradient a logical value indicating whether to calculate the derivative
#'
#' @return total restraint energy calculated using `loss_func`
#'
#' The optional derivative is contained in the `"gradient"` attribute. It is a 3D array
#' (atoms, xyz, models).
#'
#' @export
coord_array_to_relax_energy <- function(coord_array, rates, spec_den_relax_data_list, loss_func = power_scaled_loss, ..., gradient = FALSE) {
	# intermediate derivatives that need to be captured for back-propagation
	d_array_list <- g_matrix_list <- relax_vector_list <- vector("list", length(spec_den_relax_data_list))

	for (i in seq_along(spec_den_relax_data_list)) {
		spec_den_relax_data <- spec_den_relax_data_list[[i]]
		unit <- spec_den_relax_data[["unit"]]
		if (is.null(unit)) {
			unit <- FALSE
		}
		stopifnot(is.logical(unit), length(unit) == 1)

		relax_data_list <- spec_den_relax_data[["relax_data_list"]]
		if (is.null(relax_data_list)) {
			stop("`spec_den_relax_data` must contain `relax_data_list`")
		}

		# calculate internuclear vectors (convert from Å to m)
		r_array <- coord_array_to_r_array(coord_array * 1e-10, spec_den_relax_data[["atom_pairs"]][, 1:2, drop = FALSE])

		# calculate dipole-dipole interaction tensors
		d_array_list[[i]] <- r_array_to_d_array(r_array, dist = !unit, unit = unit, gradient = gradient)

		# calculate the factor by which the number of models should be expanded
		n_shift <- ncol(spec_den_relax_data[["groupings"]]) / dim(coord_array)[3]

		# shift tensor components from atom pairs into virtual models
		d_array_shifted <- array_shift(d_array_list[[i]], n_shift)

		# calculate matrix of g values
		g_matrix_list[[i]] <- d_array_to_g_matrix(d_array_shifted, spec_den_relax_data[["groupings"]], gradient = gradient)

		# calculate matrix of internal amplitudes
		a_int_matrix <- g_matrix_to_a_matrix(g_matrix_list[[i]], spec_den_relax_data[["a_int_coef"]])

		# calculate internal lambda eigenvalues
		lambda_int_vec <- -colSums(
			rates[rownames(spec_den_relax_data[["lambda_int_coef"]])] * spec_den_relax_data[["lambda_int_coef"]]
		)

		# calculate unit dipole-dipole interaction tensors if not already done
		d_unit_array <- if (unit) d_array_list[[i]] else r_array_to_d_array(r_array, dist = FALSE, unit = TRUE)

		# shift unit tensor components from atom pairs into virtual models
		d_unit_array_shifted <- array_shift(d_unit_array, n_shift)

		# calculate overall modes from diffusion constants
		overall_modes <- dxyz_dunit_to_overall_modes(rates[c("Dx", "Dy", "Dz")], d_unit_array_shifted)

		relax_vector_list[[i]] <- lapply(
			relax_data_list,
			function(relax_entry) {
				a_matrix_to_relax(
					a_int_matrix,
					lambda_int_vec,
					overall_modes[["a_overall_matrix"]],
					overall_modes[["lambda_overall_vec"]],
					relax_entry[["spec_den_term_array"]],
					gradient = gradient
				)
			}
		)
	}

	# combine calculated relaxation rates into a single vector
	relax_vector <- unlist(
		lapply(relax_vector_list, function(x) unlist(x, use.names = FALSE)),
		use.names = FALSE
	)

	# combine target relaxation rates into a single vector
	relax0_vector <- unlist(
		lapply(spec_den_relax_data_list, function(spec_den_relax_data) {
			unlist(lapply(spec_den_relax_data[["relax_data_list"]], `[[`, "value"), use.names = FALSE)
		}),
		use.names = FALSE
	)

	# combine optional force constants into a single vector
	k_vector <- unlist(
		lapply(spec_den_relax_data_list, function(spec_den_relax_data) {
			unlist(
				lapply(spec_den_relax_data[["relax_data_list"]], function(relax_entry) {
					n_value <- length(relax_entry[["value"]])
					if (!"k" %in% names(relax_entry)) {
						rep(1, n_value)
					} else if (length(relax_entry[["k"]]) == 1) {
						rep(relax_entry[["k"]], n_value)
					} else {
						relax_entry[["k"]]
					}
				}),
				use.names = FALSE
			)
		}),
		use.names = FALSE
	)

	loss_args <- list(...)
	if (!"k" %in% names(loss_args)) {
		loss_args$k <- k_vector
	}
	loss_args$x <- relax_vector
	loss_args$x0 <- relax0_vector
	loss_args$gradient <- gradient

	# calculate energies from the relaxation-rate vectors
	energy_vector <- do.call(loss_func, loss_args)

	# return the sum of all the individual restraint energies
	value <- sum(energy_vector)

	if (gradient) {
		# initialize empty gradient with dimensions equal to input coordinates
		d_energy_d_coord_array <- coord_array
		d_energy_d_coord_array[] <- 0

		energy_grad <- attr(energy_vector, "gradient")
		offset <- 0

		for (i in seq_along(spec_den_relax_data_list)) {
			spec_den_relax_data <- spec_den_relax_data_list[[i]]
			relax_value_list <- relax_vector_list[[i]]
			d_energy_d_a_int_matrix <- matrix(
				0,
				nrow = nrow(g_matrix_list[[i]]),
				ncol = ncol(spec_den_relax_data[["a_int_coef"]])
			)

			for (j in seq_along(relax_value_list)) {
				relax_value <- relax_value_list[[j]]
				d_energy_d_relax <- energy_grad[seq.int(offset + 1, length.out = length(relax_value))]
				offset <- offset + length(relax_value)

				d_energy_d_a_int_matrix <- d_energy_d_a_int_matrix +
					a_matrix_to_relax_backprop(attr(relax_value, "gradient"), d_energy_d_relax)
			}

			# back-propagate derivatives from a matrix to g matrix
			d_energy_d_g_matrix <- g_matrix_to_a_matrix_backprop(spec_den_relax_data[["a_int_coef"]], d_energy_d_a_int_matrix)

			# back-propagate derivatives from g matrix to d array
			d_energy_d_d_array <- d_array_to_g_matrix_backprop(attr(g_matrix_list[[i]], "gradient"), d_energy_d_g_matrix)

			# back-propagate derivatives from d array to r array
			d_energy_d_r_array <- r_array_to_d_array_backprop(attr(d_array_list[[i]], "gradient"), d_energy_d_d_array)

			# accumulate back-propagated derivatives from r array into coord array gradient
			d_energy_d_coord_array <- coord_array_to_r_array_backprop(
				d_energy_d_coord_array,
				spec_den_relax_data[["atom_pairs"]][, 1:2, drop = FALSE],
				d_energy_d_r_array
			)
		}

		# set coordinate gradient to accumulated values
		attr(value, "gradient") <- d_energy_d_coord_array * 1e-10
	}

	value
}
