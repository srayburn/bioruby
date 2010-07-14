# = util/SDI/SDI_rerootable.rb - Binary SDI rerootable trees
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#
# == Description
#
# This file contains an implementation of an algorithm to infer gene duplication
# and speciation events in a rooted binary gene tree where the tree is rerootable. 
# 
# == References
#
# * http://wiki.github.com/srayburn/bioruby
# * http://bioinformatics.oxfordjournals.org/cgi/reprint/17/9/821
# * http://www.phyloxml.org
# * https://www.nescent.org/wg_phyloinformatics/PhyloSoC:PhyloXML_support_in_BioRuby
# 

#  load libraries required
require 'bio/util/SDI/SDI'

module Bio
  # == Description
  #
  # Bio::SDI_rerootable implements Speciation/Duplication inferencing for PhyloXML trees where the root can be changed.
  #
  # == Usage
  # 
  # Inputs to the constructor are a rooted phyloXML formatted binary gene tree and a rooted phyloXML formatted binary species tree.
  # The external nodes of the species tree must contain all of the species in the external nodes
  # of the gene tree. 
  #
  # Example:
  # 
  #   require 'bio'
  #   # Create an instance of sdi algorithm
  #   # gene_tree and species_tree are Bio::PhyloXML::Tree objects (see phyloxml documentation)
  #   sdi = Bio::SDI.new(gene_tree, species_tree)  
  #
  #   # compute the speciation and duplication events:
  #   updated_tree = sdi.compute_speciation_duplications 
  #
  #   # Print number of duplications
  #   # puts sdi.duplications_sum
  #

  class SDIR 
    attr_accessor :_min_dup
	attr_accessor :gene_tree
	attr_accessor :species_tree
	
	def initialize(gene_tree, species_tree)
	  @gene_tree = gene_tree
	  @species_tree = species_tree
	end #initialize
	
	def root_and_infer
	  branch_list = []
	  rooted_trees = []
	  prev_root = nil
	  prev_root_c1 = nil
	  prev_root_c2 = nil
	  duplications = 0
	  @_min_dup = 1000000000000
	  prev_root_was_dup = false
	  # The following is a messy deep copy. Would need to implement clone() for PhyloXML tree object.
	  g = Marshal::load(Marshal.dump(gene_tree))
	  
	  if(g.leaves.length <= 1)
	    set_rooted(g, true)
		@_min_dup = 0
		tree_array = [g]
		return tree_array
	  end #if
	  
	  g.nodes.each { |n|
	    if !g.leaves.include?(n) && (g.children(n).length != 2)
		  raise "Gene tree must be completely binary."
		end #if
	  }
	  species_tree.nodes.each { |n|
	    if !species_tree.leaves.include?(n) && (species_tree.children(n).length != 2)
		  raise "Species tree must be completely binary."
		end #if
	  }
	 # puts g.leaves[0]
	  set_rooted(g, true)
	  reroot(g, g.leaves[0])
	  branches = get_branches_preorder(g)
	  
	  sdi  = Bio::SDI_rerootable.new(g, species_tree)
	  g = sdi.compute_speciation_duplications 
      duplications = sdi.duplications_sum
	  used_root_placements = []
	  
	  branches.each { |b|
	    prev_root = g.root
		prev_root_c1 = g.children(prev_root)[0]
		prev_root_c2 = g.children(prev_root)[1]
		prev_root_was_dup = duplication?(prev_root)
		reroot(g,b)
		duplications = sdi.update_mapping(prev_root_was_dup, prev_root_c1, prev_root_c2)
		puts duplications
		puts @_min_dup
		if !used_root_placements.include?(b)
		  if duplications < @_min_dup
			rooted_trees.clear()
		    rooted_trees.push(Marshal::load(Marshal.dump(g)))
		    @_min_dup = duplications
		  elsif duplications == @_min_dup
		    rooted_trees.push(Marshal::load(Marshal.dump(g)))
		  end #if
		end #if
		used_root_placements.push(b)
	  }
	
	return rooted_trees
	
	end
	
	def set_rooted(tree,value)
	  tree.rooted = value
	end #set_rooted
	
	def reroot(tree, new_root)
	    before = Bio::PhyloXML::Writer.new('before.xml')
		after = Bio::PhyloXML::Writer.new('after.xml')
		before.write(tree)
		node1 = tree.children(tree.root)[0]
		node2 = tree.children(tree.root)[1]
		node3 = tree.parent(new_root)
	    prev_root = tree.root
		
		tree.remove_node(tree.root)
		tree.add_edge(node2, node1)
		new_root_node = Bio::PhyloXML::Node.new
		tree.add_node(new_root_node)
		tree.root = new_root_node
		#puts tree.root
    		tree.add_edge(tree.root, node3)
		tree.add_edge(tree.root, new_root)
		tree.remove_edge(node3, new_root)
		
	after.write(tree)
	end #reroot(new_root)
	
	def get_branches_preorder(tree)
	  #puts tree.nodes
	  nodes = [tree.root]
	  b_list = []
	  while !nodes.empty?
	    current = nodes.pop()
		#puts "nodes length"
		#puts nodes.length
		b_list.push(current)
		#puts "current "
		#puts current
		#puts "leaves"
		#puts tree.leaves
		
		if !tree.leaves.include?(current)
		 # puts "not leaf"
		  nodes.push(tree.children(current)[0])
		  nodes.push(tree.children(current)[1])
		end #if
		#puts "nodes"
		#puts nodes
		#puts "end nodes"
	  end #while

	  b_list.shift
	  b_list.shift
	  b_list.shift
	  return b_list
	end #get_branches_preorder
	
  
    def duplication?(node)
	  if node.events != nil
		return (node.events.duplications > 0)
	  else return false
	  end #if
	end #duplication?
  end

end #module