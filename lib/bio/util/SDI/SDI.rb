# = util/SDI/SDI.rb - Binary SDI
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#
# == Description
#
# This file contains an implementation of an algorithm to infer gene duplication
# and speciation events in a rooted binary gene tree. 
#
#
#
# == References
#
# * http://wiki.github.com/srayburn/bioruby
# * http://bioinformatics.oxfordjournals.org/cgi/reprint/17/9/821
# * http://www.phyloxml.org
# * https://www.nescent.org/wg_phyloinformatics/PhyloSoC:PhyloXML_support_in_BioRuby
# 

#  load libraries required

#require 'bio/db/phyloxml/phyloxml_parser'

module Bio
  # == Description
  #
  # Bio::SDI implements Speciation/Duplication inferencing for PhyloXML trees.
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
  class SDI

    # Bio::PhyloXML::Tree object containing rooted binary gene tree
    attr_accessor :gene_tree
    # Bio::PhyloXML::Tree object containing rooted binary species tree of all species in gene tree
    attr_accessor :species_tree
    # hash mapping from node in gene_tree to integer
    attr_accessor :gene_mapping 
    # hash mapping from node in species_tree to integer 
    attr_accessor :species_numbering
    # hash mapping from integer to node in species_tree
    attr_accessor :spec_node_map
	# integer count of duplications
	attr_accessor :duplications_sum
    
	# Create a new Bio::SDI object
	#   
	#    sdi = Bio::SDI.new(gene_tree, species_tree)
	#    #gene_tree and species_tree are Bio::PhyloXML::Tree objects (see phyloxml documentation)
    #   
	# The initialization verifies that trees are non-empty and rooted, and that the species tree leaves
	# include all of the species represented in the gene tree and raises an exception if any of the 
	# above conditions fail.
	# ---
    # *Arguments*:
	# * (required) _gene_tree_: Bio::PhyloXML::Tree object
	# * (required) _species_tree_: Bio::PhyloXML::Tree object
	#
    def initialize(gene_tree, species_tree)
      # Bio::PhyloXML::Tree objects
      @gene_tree = gene_tree
      @species_tree = species_tree

      # Raise exception for empty tree
      if (@gene_tree == nil) || (@species_tree == nil)
        raise "Gene and Species trees must be non-empty."
      end #if

      # raise exception for unrooted tree
      if !is_rooted?(@gene_tree) || !is_rooted?(@species_tree)
        raise "Gene and Species trees must be rooted."
      end # if
    
      @gene_mapping = {}
      @species_numbering = {}
      @spec_node_map = {}
	  @duplications_sum = 0
      @taxonomy_comparison_basis = _determine_taxonomy_comparison_basis
	  
    end #initialize
    
    # Wrapper method for the algorithm. Calls initialization and mapping computation.
	#    
	#   updated_tree = sdi.compute_speciation_duplications
	#
	# ---
	# *Returns*:: Bio::PhyloXML::Tree object
	#
    def compute_speciation_duplications
      initialize_mapping!
      return compute_mapping
    end #compute_speciation_duplications

    # Computes the mapping between gene and species trees and updates each clade in 
    # gene tree to have either a speciation or duplication event. Traverses
	# the gene tree in a post order fashion. Calls recursive method
    # _compute_mapping(node). Called by compute_speciation_duplications()
	# 
	# Note: Should ONLY be called by wrapper compute_speciation_duplications()
	# Note to self: Should I go ahead and make this private?
	# ---
	# *Returns*:: Bio::PhyloXML::Tree object
    #
    def compute_mapping 
 
      root = @gene_tree.root
      next_nodes = @gene_tree.children(root)
      next_nodes.each { |n|
        index = _compute_mapping(n) 
      } #end each

      a = @gene_mapping[next_nodes[0]]
      b = @gene_mapping[next_nodes[1]]
    
	  while a != b
	    if a > b
		  a = @species_numbering[@species_tree.parent(@spec_node_map[a])]
		else
		  b = @species_numbering[@species_tree.parent(@spec_node_map[b])]
		end #if
	  end #while
      @gene_mapping[root] = a
      if (@gene_mapping[root] == @gene_mapping[next_nodes[0]]) || (@gene_mapping[root] == @gene_mapping[next_nodes[1]])
        if root.events == nil
		  root.events = Bio::PhyloXML::Events.new
		end #if
        root.events.duplications = 1
      else
	    if root.events == nil
          root.events = Bio::PhyloXML::Events.new
	    end #if
        root.events.speciations = 1
      end #if
      
      return @gene_tree

    end #compute_mapping
  
    # Recursive, private helper to compute_mapping(). Post order traversal of gene tree.
	# Only called by compute_mapping.
	# ---
	# *Arguments*:
	# * (required) _node_: Bio::PhyloXML::Node object
	#
    def _compute_mapping(node)  		 #:doc:
     
      next_nodes = @gene_tree.children(node)
      next_nodes.each { |n|
        index = _compute_mapping(n) 
      } #end each
      if next_nodes[0] != nil:
        a = @gene_mapping[next_nodes[0]]
        if next_nodes[1] != nil:
          b = @gene_mapping[next_nodes[1]]
        else 
          b = @gene_mapping[next_nodes[0]]
        end #if  
   	   
	    while a != b
	      if a > b
		    a = @species_numbering[@species_tree.parent(@spec_node_map[a])]
		  else
		    b = @species_numbering[@species_tree.parent(@spec_node_map[b])]
		  end #if
		end #while
		@gene_mapping[node] = a
        if (@gene_mapping[node] == @gene_mapping[next_nodes[0]]) || (@gene_mapping[node] == @gene_mapping[next_nodes[1]])
          if node.events == nil
			node.events = Bio::PhyloXML::Events.new
		  end #if
          node.events.duplications = 1
		  @duplications_sum += 1
        else
		  if node.events == nil
            node.events = Bio::PhyloXML::Events.new
		  end #if
          node.events.speciations = 1
        end #if

      end #if

    end #_compute_mapping
  
    # Numbers nodes of species tree in a preorder fashion.
    # Calls private helper _initialize_species_map()
    # modifies instance variables @species_numbering and @spec_node_map
	# Only called by initialize_mapping!
	#  
    def initialize_species_map!
      root = @species_tree.root
      index = 1
    
      @species_numbering[root] = index
      @spec_node_map[index] = root
      index += 1
    
      next_nodes = species_tree.children(root)
      next_nodes.each { |n|
        index = _initialize_species_map!( n, index)
      }
   
    end #initialize_species_map!
 
    # recursive, private helper of initialize_species_map()
	# modifies instance variables @species_numbering and @spec_node_map
	# Only called by initialize_species_map!
	#
    def _initialize_species_map!( node, index)         # :doc:
   
      @species_numbering[node] = index
      @spec_node_map[index] = node
      index = index + 1
   
      next_nodes = species_tree.children(node)
      next_nodes.each{ |n|
        index = _initialize_species_map!( n, index)
      }
      return index

    end #_initialize_species_map!
  
    # Implements initialization of algorithm.
    # Numbers species tree nodes in preorder traversal
    # and sets mapping of leaf nodes in gene tree to
    # the mapping of a matching leaf node in the species tree.
	# Calls initialize_species_map.
	#
    def initialize_mapping!
     
      initialize_species_map!
	  
      index = 0
	  spec_leaves = {}
	  @species_tree.leaves.each{ |n|
	    spec_leaves[taxonomy_to_string(n)] = n
	  }
	  
	
      @gene_tree.leaves.each{ |n|
	    key = spec_leaves[taxonomy_to_string(n)]
		if key == nil
		  raise "All taxonomies in gene tree must be represented in the species tree."
		end #if
	    @gene_mapping[n] = @species_numbering[key]
      }
      
    end #initialize_mapping!

    # Tests to see if nodes are equivalent. Checks taxonomy id first, then code, then scientific name, then common name.
    # Raises fatal exception if nodes do not have enough information to match.
	# ---
	# *Arguments*:
	# * (required) _node1_: Bio::PhyloXML::Node object
	# * (required) _node2_: Bio::PhyloXML::Node object
	# *Returns*:: Boolean
	#
    def node_equal?(node1, node2)

	  
	  if (!node1.taxonomies.empty?) && (!node2.taxonomies.empty?)
	  # compare taxonomy id if exists and provider is the same
	  if (node1.taxonomies[0].taxonomy_id != nil) && (node2.taxonomies[0].taxonomy_id != nil)
		if (node1.taxonomies[0].taxonomy_id.provider != nil) && (node2.taxonomies[0].taxonomy_id.provider != nil)
		  if node1.taxonomies[0].taxonomy_id.value == node2.taxonomies[0].taxonomy_id.value
		   return true
		  else 
			return false
		  end # if
		end #if
	
	  # otherwise compare code
	  elsif (node1.taxonomies[0].code != nil) && (node2.taxonomies[0].code != nil) 
		if node1.taxonomies[0].code == node2.taxonomies[0].code
		  return true
		else 
		  return false
		end #if
  
	  # otherwise compare scientific name
	  elsif (node1.taxonomies[0].scientific_name != nil) && (node2.taxonomies[0].scientific_name != nil)
		if node1.taxonomies[0].scientific_name == node2.taxonomies[0].scientific_name
		  return true
		else 
		  return false
		end #if
  
	  # otherwise compare common names
	  elsif (node1.taxonomies[0].common_names[0] != nil) && (node2.taxonomies[0].common_names[0] != nil)
		if node1.taxonomies[0].common_names[0] == node2.taxonomies[0].common_names[0]
		  return true
		else 
		  return false
		end # if
	  end #if  
	  end #if
	  
	  # otherwise, not enough in common to compare
      raise "Nodes must share common attributes in order to be compared"
    end #node_equal?

    # Tests tree for root.
	# ---
	# *Arguments*:
	# * (required) _tree_: Bio::PhyloXML::Tree object
	# *Returns*:: Boolean
	#
    def is_rooted?(tree)
      if tree.root == nil:
        return false
      else
        return true
      end #if
    
    end #is_rooted?

    # Determines what taxonomy element is shared by all external nodes in gene & species trees. 
	# Helper for taxonomy_to_string
	# ---
	# *Returns*:: Symbol
	#
	def _determine_taxonomy_comparison_basis
	  all_have_id = true
	  all_have_code = true
	  all_have_sn = true
	  all_have_cn = true
	  
	  @species_tree.leaves.each{ |n|
	    if !n.taxonomies.empty?
          if n.taxonomies[0].taxonomy_id == nil
		    all_have_id = false
		  end
		  if n.taxonomies[0].code == nil
		    all_have_code = false
		  end
		  if n.taxonomies[0].scientific_name == nil
		    all_have_sn = false
		  end
		  if n.taxonomies[0].common_names.empty?
		    all_have_cn = false
		  end
		else 
		  raise "Species tree node has no taxonomic information"
		end #if
	  } #end species_tree.leaves.each
	  @gene_tree.leaves.each{ |n|
	    if !n.taxonomies.empty?
          if n.taxonomies[0].taxonomy_id == nil
		    all_have_id = false
		  end
		  if n.taxonomies[0].code == nil
		    all_have_code = false
		  end
		  if n.taxonomies[0].scientific_name == nil
		    all_have_sn = false
		  end
		  if n.taxonomies[0].common_names.empty?
		    all_have_cn = false
		  end
		else 
		  raise "Gene tree node has no taxonomic information"
		end #if
	  } #end gene_tree.leaves.each
	  if all_have_id
	    return :id
	  elsif all_have_code
	    return :code
	  elsif all_have_sn
	    return :sn
	  elsif all_have_cn
	    return :cn
	  else
	    raise "Gene and Species trees have incomparable taxonomies"
	  end #if
	end #_determine_taxonomy_comparison_basis
	
	# Converts taxonomic information to string for use as key in performance enhancement.
	# ---
	# *Arguments*:
	# *(required) _node_: Bio::PhyloXML::Node object
	# *Returns*:: String
	#
	def taxonomy_to_string(node)
	  case @taxonomy_comparison_basis
	    when :id
		  return node.taxonomies[0].taxonomy_id.to_s
		when :code
		  return node.taxonomies[0].code.to_s
		when :sn
		  return node.taxonomies[0].scientific_name.to_s
		when :cn
		  return node.taxonomies[0].common_names[0].to_s
		else 
		  raise "Invalid taxonomy comparison code"
	  end #case
	end #taxonomy_to_string
	
	
	
    private :_initialize_species_map!
    private :_compute_mapping
	private :_determine_taxonomy_comparison_basis

  end
  # == Description
  #
  # Bio::SDI_rerootable implements and extension to Speciation/Duplication inferencing for PhyloXML trees 
  # where the root can be changed. It inherits from Bio::SDI. The extension to the algorithm takes advantage
  # of the observation that only a few nodes need to change if the root changes.
  #
  # == Usage
  # 
  # Inputs to the constructor are a rooted phyloXML formatted binary gene tree and a rooted phyloXML formatted binary species tree.
  # The external nodes of the species tree must contain all of the species in the external nodes
  # of the gene tree. 
  #
  # This object is meant to be used by the SDIR object and not directly by the user. 
  # Note to self: are there private classes in ruby??
  #
  class SDI_rerootable < SDI
    # Updates Events for a rerooted gene tree
	# ---
	# *Arguments*:
	# * (required) _prev_root_was_dup_: Boolean
	# * (required) _prev_root_c1: Bio::PhyloXML::Node object
	# * (required) _prev_root_c2: Bio::PhyloXML::Node object
	# *Returns*:: Integer
	#
    def update_mapping(prev_root_was_dup, prev_root_c1, prev_root_c2)
	  if (@gene_tree.children(@gene_tree.root)[0] == prev_root_c1) || (@gene_tree.children(@gene_tree.root)[1] == prev_root_c1)
	    calculate_mapping_for_node( prev_root_c1)
	  else
	    calculate_mapping_for_node(prev_root_c2)
	  end #if
	  
	  if prev_root_was_dup
	    if @gene_tree.root.events == nil
		  @gene_tree.root.events = Bio::PhyloXML::Events.new
		end #if
	    @gene_tree.root.events.duplications = 1
	  else
	    if @gene_tree.root.events == nil
		  @gene_tree.root.events = Bio::PhyloXML::Events.new
		end #if
	    @gene_tree.root.events.speciations = 1
	  end #if
	  
	  calculate_mapping_for_node(@gene_tree.root)
	  return @duplications_sum
	end #update_mapping
	# Helper method for update_mapping
	# ---
	# *Arguments*:
	# * (required) _node_: Bio::PhyloXML::Node object
	#
	def calculate_mapping_for_node(node)
	  if !@gene_tree.leaves.include?(node)
	    if node.events != nil
		  if node.events.duplications != nil && node.events.duplications > 0
		    was_duplication = true
		  end #if
		else 
		  was_duplication = false
		end #if
	    a = @gene_mapping[@gene_tree.children(node)[0]]
        b = @gene_mapping[@gene_tree.children(node)[1]]
    
        @gene_mapping[node] = @species_numbering[@species_tree.lowest_common_ancestor(@spec_node_map[a], @spec_node_map[b])]
 
        if (@gene_mapping[node] == @gene_mapping[@gene_tree.children(node)[0]]) || (@gene_mapping[node] == @gene_mapping[@gene_tree.children(node)[1]])
          if node.events == nil
			node.events = Bio::PhyloXML::Events.new
		  end #if
          node.events.duplications = 1
		  if !was_duplication
		    @duplications_sum += 1
		  end #if
        else
		  if node.events == nil
            node.events = Bio::PhyloXML::Events.new
          end #if
		  node.events.speciations = 1
        end #if
      end #if
	end  #calculate_mapping_for_node(node)
  end #class
  
	
end #module
